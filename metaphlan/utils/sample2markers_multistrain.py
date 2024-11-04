#! /usr/bin/env python
__author__ = 'Michal Puncochar (michal.puncochar@unitn.it)'
__version__ = '4.1.1'
__date__ = '4 Nov 2024'


import subprocess as sp
import argparse as ap
import multiprocessing.pool as mpp
import pathlib
import traceback

from tqdm.auto import tqdm

try:
    import tomllib  # Python 3.11+
except ModuleNotFoundError:
    import tomli as tomllib  # https://github.com/hukkin/tomli

from metaphlan.utils.multistrain.pipeline import run
from metaphlan.utils.multistrain.utils import ArgumentType, error
from metaphlan.utils.database_controller import StrainphlanDatabaseController
from metaphlan.utils import info


def read_params():
    p = ap.ArgumentParser(description="Strain sharing pipeline: After-StrainPhlAn calculations.",
                          formatter_class=ap.ArgumentDefaultsHelpFormatter)

    # TODO: consistency --input vs. --samples across strainphlan, sample2markers, sample2markers_multistrain
    p.add_argument('-i', '--input', type=ArgumentType.existing_file, nargs='+', required=True,
                   help="Paths to the SAM files, separate multiple samples by space")
    p.add_argument('-o', '--output_dir', type=ArgumentType.creatable_dir, required=True,
                   help="Path to the output directory")
    p.add_argument('-d', '--database', type=ArgumentType.existing_file, default='latest',
                   help="Path to the MetaPhlAn database pkl file")
    p.add_argument('-t', '--threads', type=ArgumentType.positive_int, default=1, help="Number of threads")
    # TODO: quasi markers behavior
    # TODO: target (what to calculate) + save_files (what not to remove) + reuse (what to potentially reuse)
    p.add_argument('--config', type=ArgumentType.existing_file, default=None, help="Path to a config file")
    p.add_argument('--reuse', type=str, default='all', choices=['none', 'bam', 'pileup', 'all'],
                   help="Which intermediate files to reuse. None re-runs everything. "
                        "Bam will reuse the sorted bam file.")
    p.add_argument('--save_intermediate_files', action='store_true', default=False,
                   help="Whether to store the intermediate files (potentially big)")

    return p


def check_params(argp):
    args = argp.parse_args()

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


def try_run(*args):
    try:
        run(*args, marker_to_clade=try_run.marker_to_clade)
        return True
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


    mp_db_controller = StrainphlanDatabaseController(args.database)
    marker_to_clade = mp_db_controller.get_markers2clade()
    mp_version = mp_db_controller.get_database_name()


    ss_args = []
    for sample_path in args.input:
        if sample_path.suffix == '.bz2':
            sample_name = sample_path.with_suffix('').stem
        else:
            sample_name = sample_path.stem
        output_dir = args.output_dir / sample_name
        ss_args.append((sample_path, output_dir, config, args.save_intermediate_files, args.reuse, mp_version))

    info(f'Running on {len(ss_args)} samples')
    if args.threads == 1 or len(ss_args) == 1:
        try_run.marker_to_clade = marker_to_clade

        successes = [try_run(*ss_arg) for ss_arg in ss_args]
    else:
        # bind marker_to_clade to the function, so it's accessible in all threads
        # (this is surprisingly without much overhead)
        def init(fu, mtc):
            fu.marker_to_clade = mtc

        with mpp.Pool(args.threads, initializer=init, initargs=(try_run, marker_to_clade)) as pool:
            successes = list(tqdm(pool.starmap(try_run, ss_args), total=len(ss_args)))

    if all(successes):
        info(f'Successfully finished sample2markers multi-strain, results are stored at {args.output_dir}')
    else:
        error('Finished sample2markers with some errors')


if __name__ == '__main__':
    main()
