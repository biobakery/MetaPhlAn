#! /usr/bin/env python
import argparse as ap
import base64
import io
import json
import ctypes
import multiprocessing as mp
import multiprocessing.pool as mpp
import os
import pathlib
import zipfile
from collections import defaultdict

import numpy as np
import pandas as pd
from tqdm import tqdm


try:
    import tomllib  # Python 3.11+

except ModuleNotFoundError:
    import tomli as tomllib  # https://github.com/hukkin/tomli


from . import info, error
from .database_controller import StrainphlanDatabaseController
from .multistrain.utils import ArgumentType, MetaphlanDBInfo
from .multistrain.strain_tracking_commons import load_sample, load_sample_parallel, process_samples_argument


class Arguments:
    samples: list[pathlib.Path]
    sample_list: pathlib.Path
    database: pathlib.Path
    config: pathlib.Path
    clean: bool
    threads: int


def read_params():
    p = ap.ArgumentParser(description="StrainPhlAn multi-strain: define rare alleles",
                          formatter_class=ap.ArgumentDefaultsHelpFormatter)

    group_input = p.add_mutually_exclusive_group(required=True)
    group_input.add_argument('--samples', type=ArgumentType.existing_dir, nargs='+',
                             help="Paths to the allele counts files")
    group_input.add_argument('--sample_list', type=ArgumentType.file_list_of_dirs,
                             help="File with one sample per line, each line is a allele counts file")
    p.add_argument('--database', type=ArgumentType.existing_file, required=True,
                   help="Path to the MetaPhlAn database .pkl")
    p.add_argument('--config', type=ArgumentType.existing_file, help="A path to config file", required=True)
    p.add_argument('--clean', action='store_true', help="Don't reuse anything")
    p.add_argument('--threads', type=ArgumentType.positive_int, default=1, help="Number of threads")

    return p


def check_params(argp: ap.ArgumentParser):
    args = argp.parse_args(namespace=Arguments())

    if args.samples is None:
        args.samples = args.sample_list

    return args


def count_alleles_shared(arg):
    sample_name, sample_ac_path = arg

    config = count_alleles_shared.config
    metaphlan_db = count_alleles_shared.metaphlan_db
    sgb_to_samples = count_alleles_shared.sgb_to_samples
    allele_counts_shared_raw = count_alleles_shared.allele_counts_shared_raw
    marker_to_offset = count_alleles_shared.marker_to_offset
    sgb_to_lock = count_alleles_shared.sgb_to_lock

    _, sample_allele_counts = load_sample(config, sample_ac_path, None, metaphlan_db)

    for m, ac in sample_allele_counts.items():
        sgb_id = metaphlan_db.marker_to_clade[m]
        if sgb_id not in sgb_to_samples.keys() or sample_name not in sgb_to_samples[sgb_id]:
            continue
        with sgb_to_lock[sgb_id]:
            a_np = np.frombuffer(allele_counts_shared_raw, offset=marker_to_offset[m],
                                 count=metaphlan_db.marker_to_len[m]*4, dtype=np.int32).reshape((4, -1))
            a_np += (ac >= config['min_allele_cov'])


def main():
    argp = read_params()
    args = check_params(argp)

    with open(args.config, 'rb') as f:
        config = tomllib.load(f)

    db_dir = pathlib.Path(config["db_dir"])
    db_dir.mkdir(parents=True, exist_ok=True)



    sample_name_to_ac_path = process_samples_argument(args.samples)


    mp_db_controller = StrainphlanDatabaseController(args.database)
    mp_db_controller.load_database()

    info('Processing MetaPhlAn database')
    mp_db_info = MetaphlanDBInfo.from_mp_controller(mp_db_controller)



    sgb_to_samples_file = db_dir / 'sgb_to_samples.tsv'

    if sgb_to_samples_file.exists() and not args.clean:
        info('Reusing SGB to samples file')
        sgb_to_samples = pd.read_csv(sgb_to_samples_file, sep='\t', index_col=0).squeeze('columns') \
            .str.replace("'", '"').map(json.loads).to_dict()
    else:
        info('Counting the prevalence of SGBs in samples')

        load_sample_parallel.config = config
        load_sample_parallel.mp_db_info = mp_db_info

        sgb_to_samples = defaultdict(list)
        with mpp.Pool(args.threads) as pool:
            for sample_name, res in tqdm(pool.imap_unordered(load_sample_parallel, sample_name_to_ac_path.items()),
                                         total=len(sample_name_to_ac_path)):
                if res is None:
                    continue
                clades_present, _ = res
                for clade in clades_present:
                    sgb_to_samples[clade].append(sample_name)

        pd.Series(sgb_to_samples).to_csv(sgb_to_samples_file, sep='\t')


    sgb_to_samples = {sgb_id: samples for sgb_id, samples in sgb_to_samples.items()
                      if len(samples) >= config['min_sgb_prevalence_training']}

    info(f'  Remaining {len(sgb_to_samples)} prevalent SGBs')


    info('Counting alleles')
    info('  Allocating array')
    tot_size = 0
    marker_to_offset = {}
    sgb_to_lock = {}
    for sgb_id in tqdm(sgb_to_samples.keys()):
        sgb_to_lock[sgb_id] = mp.Lock()
        for m in mp_db_info.clade_to_markers[sgb_id]:
            marker_to_offset[m] = tot_size
            tot_size += mp_db_info.marker_to_len[m] * 4

    allele_counts_shared_raw = mp.Array(ctypes.c_int32, tot_size, lock=False)



    count_alleles_shared.config = config
    count_alleles_shared.metaphlan_db = mp_db_info
    count_alleles_shared.sgb_to_samples = sgb_to_samples
    count_alleles_shared.allele_counts_shared_raw = allele_counts_shared_raw
    count_alleles_shared.marker_to_offset = marker_to_offset
    count_alleles_shared.sgb_to_lock = sgb_to_lock

    info('  Running counting')
    with mpp.Pool(args.threads) as pool:
        list(tqdm(pool.imap_unordered(count_alleles_shared, sample_name_to_ac_path.items()),
                  total=len(sample_name_to_ac_path)))


    info('  Converting')
    allele_counts = defaultdict(dict)
    for sgb_id in tqdm(sgb_to_samples.keys()):
        for m in mp_db_info.clade_to_markers[sgb_id]:
            ac = np.frombuffer(allele_counts_shared_raw, offset=marker_to_offset[m],
                               count=mp_db_info.marker_to_len[m]*4, dtype=np.int32).reshape((4, -1))
            sgb_id = mp_db_info.marker_to_clade[m]
            allele_counts[sgb_id][m] = ac




    info('Saving allele counts')
    db_metadata = {
        'database_name': 'mpa_vJun23_CHOCOPhlAnSGB_202307'
    }

    for sgb_id, sgb_allele_counts in tqdm(allele_counts.items()):
        allele_counts_file = db_dir / f'{sgb_id}.zip'
        allele_counts_file.parent.mkdir(exist_ok=True)
        allele_counts_file_tmp = allele_counts_file.with_name(allele_counts_file.name + '.tmp')
        with zipfile.ZipFile(allele_counts_file_tmp, 'w', compression=zipfile.ZIP_DEFLATED) as fo:
            with fo.open('metadata.json', 'w') as fm:
                json.dump(db_metadata, io.TextIOWrapper(fm))

            for m, c in sgb_allele_counts.items():
                with fo.open(m + '.npy', 'w') as fm:
                    np.save(fm, c)

        allele_counts_file_tmp.rename(allele_counts_file)



    info('Removing temporary files')
    os.unlink(sgb_to_samples_file)

    info('Done.')



if __name__ == '__main__':
    main()
