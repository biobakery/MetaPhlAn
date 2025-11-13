#! /usr/bin/env python
__author__ = 'Michal Puncochar (michal.puncochar@unitn.it)'
__version__ = '4.1.1'
__date__ = '4 Nov 2024'


import subprocess as sp
import argparse as ap
import multiprocessing.pool as mpp
import pathlib
import traceback
from typing import Sequence

from tqdm.auto import tqdm

from ..utils import global_flags

try:
    import tomllib  # Python 3.11+
except ModuleNotFoundError:
    import tomli as tomllib  # https://github.com/hukkin/tomli

from . import info
from .multistrain.pipeline import run
from .multistrain.utils import ArgumentType, error
from .database_controller import StrainphlanDatabaseController


class ArgTypes:
    input: Sequence[pathlib.Path] | None
    input_list: pathlib.Path | None
    output_dir: pathlib.Path
    database: pathlib.Path
    threads: int
    config: pathlib.Path
    target: str
    reuse: str
    save_bam_file: bool
    debug: bool
    output_suffix: str



def read_params():
    p = ap.ArgumentParser(description="StrainPhlAn sample2markers multi-strain version.",
                          formatter_class=ap.ArgumentDefaultsHelpFormatter)

    # TODO: consistency --input vs. --samples across strainphlan, sample2markers, sample2markers_multistrain
    input_g = p.add_mutually_exclusive_group(required=True)
    input_g.add_argument('-i', '--input', type=ArgumentType.existing_file, nargs='+',
                         help="Paths to the SAM files, separate multiple samples by space")
    input_g.add_argument('--input_list', type=ArgumentType.existing_file,
                         help='Path to a txt file, on each line there\'s a path to a SAM file as in --input')
    p.add_argument('-o', '--output_dir', type=ArgumentType.creatable_dir, required=True,
                   help="Path to the output directory")
    p.add_argument('-d', '--database', type=ArgumentType.existing_file, required=True,
                   help="Path to the MetaPhlAn database pkl file")
    p.add_argument('-t', '--threads', type=ArgumentType.positive_int, default=1, help="Number of threads")
    # TODO: config or arguments?
    p.add_argument('--config', type=ArgumentType.existing_file, default=None, help="Path to a config file")
    p.add_argument('--reuse', type=str, default='bam', choices=['none', 'bam', 'pileup', 'all'],
                   help="Which intermediate files to reuse if present. None will force to re-run everything.")
    p.add_argument('--save_bam_file', action='store_true', default=False,
                   help="Whether to keep the preprocessed BAM file")
    p.add_argument('--debug', action='store_true', default=False,
                   help="Store intermediate files for debugging including pileup and BAM file")
    p.add_argument('--output_suffix', type=str, default='', help="Suffix to append to json file names.")

    return p


def check_params(argp: ap.ArgumentParser):
    args = argp.parse_args(namespace=ArgTypes())

    if args.input_list is not None:
        with open(args.input_list) as f:
            args.input = [pathlib.Path(line.strip()) for line in f if line.strip()]
        info(f'Found {len(args.input)} samples in the specified input file')

    if args.config is None:
        args.config = pathlib.Path(__file__).parent / 'multistrain' / 'config-default.toml'

    return args


def check_samtools():
    samtools_version_output = sp.run("samtools", capture_output=True).stderr.decode().split('\n')
    for line in samtools_version_output:
        if line.startswith('Version:'):
            samtools_version_info = line
            break
    else:
        error("Couldn't get information about the samtools version")
        exit(1)

    samtools_version = samtools_version_info.rstrip('\n').split(' ')[1]
    samtools_version_parts = [int(x) for x in samtools_version.split('.')]
    if samtools_version_parts < [1, 17]:
        error(f'The samtools version required is >= 1.17, installed is {samtools_version}')
        exit(1)

    info(f'Using samtools version {samtools_version}')


def try_run(args):
    try:
        return run(*args, marker_to_clade=try_run.marker_to_clade, marker_to_ext=try_run.marker_to_ext,
                   clade_to_markers=try_run.clade_to_markers)
    except Exception as e:
        error(f'Error running sample {args[0]}')
        error(str(e))
        error(traceback.format_exc())
        return False


def main():
    info('Starting sample2markers multi-strain version')

    check_samtools()


    argp = read_params()
    args = check_params(argp)

    with open(args.config, 'rb') as f:
        config = tomllib.load(f)

    if args.debug:
        global_flags.debug = True


    mp_db_controller = StrainphlanDatabaseController(args.database)
    mp_db_controller.load_database()
    info('Getting information from the MetaPhlAn database')
    db_name = mp_db_controller.get_database_name()
    marker_to_clade = mp_db_controller.get_markers2clade()
    marker_to_ext = mp_db_controller.get_markers2ext()
    clade_to_markers = mp_db_controller.get_clade2markers()


    ss_args = [(sample_path, args.output_dir, config, args.save_bam_file, args.reuse, db_name, args.output_suffix)
               for sample_path in args.input]

    n_threads = min(args.threads, len(ss_args))

    info(f'Running on {len(ss_args)} samples using {n_threads} processes')

    # bind constant data to the function so that they are shared on fork (copy-on-write)
    try_run.marker_to_clade = marker_to_clade
    try_run.clade_to_markers = clade_to_markers
    try_run.marker_to_ext = marker_to_ext

    if n_threads == 1:
        successes = [try_run(ss_arg) for ss_arg in ss_args]
    else:
        with mpp.Pool(n_threads) as pool:
            successes = list(tqdm(pool.imap_unordered(try_run, ss_args), total=len(ss_args)))

    if all(successes):
        info(f'Successfully finished sample2markers multi-strain, results are stored at {args.output_dir}')
    else:
        error('Finished sample2markers with some errors')


if __name__ == '__main__':
    main()
