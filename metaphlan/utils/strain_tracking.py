#! /usr/bin/env python
import argparse as ap
import json
import multiprocessing.pool as mpp
import pathlib
import zipfile
from collections import defaultdict, Counter

import numpy as np
import pandas as pd
from tqdm import tqdm

try:
    import tomllib  # Python 3.11+
except ModuleNotFoundError:
    import tomli as tomllib  # https://github.com/hukkin/tomli


from . import info, error, warning
from .database_controller import StrainphlanDatabaseController
from .multistrain.utils import ArgumentType, MetaphlanDBInfo
from .multistrain.strain_tracking_commons import load_sample, load_sample_parallel, load_sample_parallel_only_clades, \
    process_samples_argument

OUTPUT_HEADER = ['sample_1', 'sample_2', 'sgb_id', 'pop_matches', 'pop_positions',
                 'rare_matches', 'rare_positions_1', 'rare_positions_2',
                 'rare_mathces_afs_mean_1', 'rare_matches_afs_std_1',
                 'rare_mathces_afs_mean_2', 'rare_matches_afs_std_2']


class Arguments:
    samples: list[pathlib.Path]
    sample_list: pathlib.Path
    in_memory: bool
    database_rare_alleles: pathlib.Path
    database: pathlib.Path
    clade: str
    min_allele_cov: int|None
    output_dir: pathlib.Path
    skip_input_check: bool
    config: pathlib.Path
    threads_cpu: int
    threads_io: int


def read_params():
    p = ap.ArgumentParser(description="Strain tracking using rare alleles population ANI",
                          formatter_class=ap.ArgumentDefaultsHelpFormatter)

    group_input = p.add_mutually_exclusive_group(required=True)
    group_input.add_argument('--samples', type=ArgumentType.existing_dir, nargs='+',
                             help="Paths to the multi-strain output directory (containing .rare_alleles.zip file) "
                                  "separated by space")
    group_input.add_argument('--sample_list', type=ArgumentType.file_list_of_paths,
                             help="File with one sample per line, each line is a multi-strain directory output "
                                  "directory (containing pileup.tsv.gz file)")

    p.add_argument('--in_memory', action='store_true', help="Whether to pre-load the data in-memory (pay attention!)")
    p.add_argument('--database_rare_alleles', type=ArgumentType.existing_dir, required=True,
                   help="Path to the database of rare alleles (directory)")
    p.add_argument('--database', type=ArgumentType.existing_file, required=True,
                   help="Path to the MetaPhlAn database .pkl")
    p.add_argument('--clade', type=str, help="Restrict to only this clade")
    p.add_argument('--min_allele_cov', type=int, default=None, help="Minimum allele coverage to consider")
    p.add_argument('--output_dir', type=ArgumentType.creatable_dir, required=True,
                   help="Path to the output directory")
    p.add_argument('--skip_input_check', action='store_true', help="Skip checking input files")
    p.add_argument('--config', type=ArgumentType.existing_file, default=None, help="A path to config file")
    p.add_argument('--threads_cpu', type=ArgumentType.positive_int, default=1, help="Number of threads")
    p.add_argument('--threads_io', type=ArgumentType.positive_int, default=1, help="Number of threads")

    return p


def check_params_config(argp: ap.ArgumentParser):
    args = argp.parse_args(namespace=Arguments())

    if args.config is None:
        args.config = pathlib.Path(__file__).parent / 'multistrain' / 'config-default-strain-tracking.toml'

    if not args.output_dir.exists():
        args.output_dir.mkdir()

    if args.samples is None:
        args.samples = args.sample_list

    with open(args.config, 'rb') as f:
        config = tomllib.load(f)

    if args.min_allele_cov is not None:
        config['min_allele_cov'] = args.min_allele_cov

    return args, config


def compute(sample_name_1, sample_name_2, clades_present_1, clades_present_2, allele_counts_1, allele_counts_2,
            rare_alleles_db, prevalent_clades, config, metaphlan_db):
    markers_1 = set(allele_counts_1.keys())
    markers_2 = set(allele_counts_2.keys())

    pop_ani_n = Counter()
    pop_ani_d = Counter()
    rare_ani_n = Counter()
    rare_ani_d_1 = Counter()
    rare_ani_d_2 = Counter()
    rare_match_afs_1 = defaultdict(list)
    rare_match_afs_2 = defaultdict(list)

    clades_common = set(clades_present_1).intersection(clades_present_2).intersection(prevalent_clades)
    markers_common = markers_1.intersection(markers_2)

    for m in markers_common:
        sgb_id = metaphlan_db.marker_to_clade[m]
        if sgb_id not in clades_common:
            continue

        ac1 = allele_counts_1[m]
        ac2 = allele_counts_2[m]

        ac1p = ac1 >= config['min_allele_cov']
        ac2p = ac2 >= config['min_allele_cov']

        # Pop
        d1_pop = np.logical_or.reduce(ac1p, axis=0)
        d2_pop = np.logical_or.reduce(ac2p, axis=0)

        pop_d = d1_pop & d2_pop
        pop_n = np.logical_or.reduce(ac1p & ac2p, axis=0)

        pop_ani_n[sgb_id] += np.count_nonzero(pop_n)
        pop_ani_d[sgb_id] += np.count_nonzero(pop_d)

        # Rare
        if m in rare_alleles_db:
            ac1p_rare = ac1p & rare_alleles_db[m]
            ac2p_rare = ac2p & rare_alleles_db[m]

            rare_d_1 = np.logical_or.reduce(ac1p_rare, axis=0) & d2_pop  # rare in 1 and covered in 2
            rare_d_2 = np.logical_or.reduce(ac2p_rare, axis=0) & d1_pop  # rare in 2 and covered in 1

            rare_ani_d_1[sgb_id] += np.count_nonzero(rare_d_1)
            rare_ani_d_2[sgb_id] += np.count_nonzero(rare_d_2)

            mask_rare_matches = np.logical_or.reduce(ac1p_rare & ac2p_rare, axis=0)
            n_rare_matches = np.count_nonzero(mask_rare_matches)

            if n_rare_matches > 0:
                ac1_rare = ac1[:, mask_rare_matches]
                ac2_rare = ac2[:, mask_rare_matches]
                bc1_rare = ac1_rare.sum(axis=0)
                bc2_rare = ac2_rare.sum(axis=0)
                afs1 = (rare_alleles_db[m][:, mask_rare_matches] * ac1_rare).sum(axis=0) / bc1_rare
                afs2 = (rare_alleles_db[m][:, mask_rare_matches] * ac2_rare).sum(axis=0) / bc2_rare

                rare_ani_n[sgb_id] += n_rare_matches
                rare_match_afs_1[sgb_id].extend(afs1)
                rare_match_afs_2[sgb_id].extend(afs2)

    result_rows = []
    for sgb_id in pop_ani_n.keys():
        if ((rare_ani_d_1[sgb_id] < config['min_rare_total'] and rare_ani_d_2[sgb_id] < config['min_rare_total'])
           or rare_ani_n[sgb_id] < config['min_rare_matches']) and (pop_ani_d[sgb_id] < config['min_pop_aln_len']):
            continue

        match_afs_1 = np.array(rare_match_afs_1[sgb_id])
        match_afs_mean_1 = match_afs_1.mean() if len(match_afs_1) > 0 else None
        match_afs_std_1 = match_afs_1.std() if len(match_afs_1) > 1 else None
        match_afs_2 = np.array(rare_match_afs_2[sgb_id])
        match_afs_mean_2 = match_afs_2.mean() if len(match_afs_2) > 0 else None
        match_afs_std_2 = match_afs_2.std() if len(match_afs_2) > 1 else None

        result_rows.append([
            sample_name_1,
            sample_name_2,
            sgb_id,
            pop_ani_n[sgb_id],
            pop_ani_d[sgb_id],
            rare_ani_n[sgb_id],
            rare_ani_d_1[sgb_id],
            rare_ani_d_2[sgb_id],
            match_afs_mean_1,
            match_afs_std_1,
            match_afs_mean_2,
            match_afs_std_2,
        ])

        result_rows.append([
            sample_name_2,
            sample_name_1,
            sgb_id,
            pop_ani_n[sgb_id],
            pop_ani_d[sgb_id],
            rare_ani_n[sgb_id],
            rare_ani_d_2[sgb_id],
            rare_ani_d_1[sgb_id],
            match_afs_mean_2,
            match_afs_std_2,
            match_afs_mean_1,
            match_afs_std_1,
        ])

    return result_rows


def run(sample_name_1):
    metaphlan_db = run.metaphlan_db
    prevalent_clades = run.prevalent_clades
    config = run.config
    sample_name_to_ac_path = run.sample_name_to_ac_path
    rare_alleles_db = run.rare_alleles_db

    ac_path_1 = sample_name_to_ac_path[sample_name_1]
    r1 = load_sample(config, ac_path_1, None, metaphlan_db)
    if r1 is None:
        return
    clades_present_1, allele_counts_1 = r1

    result_rows = []
    for sample_name_2, ac_path_2 in sample_name_to_ac_path.items():
        if sample_name_2 == sample_name_1:
            break  # run only the triangle, avoid self-comparisons

        r2 = load_sample(config, ac_path_2, None, metaphlan_db)
        if r2 is None:
            continue
        clades_present_2, allele_counts_2 = r2

        result_rows_batch = compute(sample_name_1, sample_name_2, clades_present_1, clades_present_2,
                                    allele_counts_1, allele_counts_2, rare_alleles_db, prevalent_clades,
                                    config, metaphlan_db)
        result_rows.extend(result_rows_batch)

    return result_rows


def run_in_memory(sample_name_1):
    metaphlan_db = run.metaphlan_db
    prevalent_clades = run.prevalent_clades
    config = run.config
    sample_name_to_ac_path = run.sample_name_to_ac_path
    rare_alleles_db = run.rare_alleles_db

    sample_to_clades_present = run_in_memory.sample_to_clades_present
    sample_to_allele_counts = run_in_memory.sample_to_allele_counts

    clades_present_1 = sample_to_clades_present[sample_name_1]
    allele_counts_1 = sample_to_allele_counts[sample_name_1]


    result_rows = []
    for sample_name_2, ac_path_2 in sample_name_to_ac_path.items():
        if sample_name_2 == sample_name_1:
            break  # run only triangle, avoid self-comparisons and duplicates

        clades_present_2 = sample_to_clades_present[sample_name_2]
        allele_counts_2 = sample_to_allele_counts[sample_name_2]

        result_rows_batch = compute(sample_name_1, sample_name_2, clades_present_1, clades_present_2,
                                    allele_counts_1, allele_counts_2, rare_alleles_db, prevalent_clades,
                                    config, metaphlan_db)
        result_rows.extend(result_rows_batch)

    return result_rows


def main():
    info('Starting strain tracking')

    argp = read_params()
    args, config = check_params_config(argp)


    info('Checking input samples')
    sample_name_to_ac_path = process_samples_argument(args.samples, check=not args.skip_input_check)
    info(f'Will run strain tracking on {len(sample_name_to_ac_path)} samples')


    mp_db_controller = StrainphlanDatabaseController(args.database)
    mp_db_controller.load_database()

    info('Processing MetaPhlAn database')
    mp_db_info = MetaphlanDBInfo.from_mp_controller(mp_db_controller)


    # Bind constant data to the function so child processes can access them without copying the data
    #   Works on unix as fork uses copy-on-write so the data won't be copied
    load_sample_parallel.config = config
    load_sample_parallel.mp_db_info = mp_db_info


    if args.in_memory:
        info(f'Loading samples into memory using {args.threads_io} thread(s)')
        sample_to_clades_present = {}
        sample_to_allele_counts = {}
        with mpp.Pool(args.threads_io) as pool:
            for sample_name, res in tqdm(pool.imap_unordered(load_sample_parallel, sample_name_to_ac_path.items()),
                                         total=len(sample_name_to_ac_path)):
                if res is None:
                    continue
                clades_present, sample_allele_counts = res
                sample_to_clades_present[sample_name] = clades_present
                sample_to_allele_counts[sample_name] = sample_allele_counts
    else:
        sample_to_allele_counts = None
        info(f'Counting SGB prevalence using {args.threads_io} thread(s)')
        sample_to_clades_present = {}
        with mpp.Pool(args.threads_io) as pool:
            for res in tqdm(pool.imap_unordered(load_sample_parallel_only_clades, sample_name_to_ac_path.items()),
                            total=len(sample_name_to_ac_path)):
                if res is None:
                    continue
                sample_name, clades_present = res
                sample_to_clades_present[sample_name] = clades_present


    clade_to_samples = defaultdict(list)
    for sample_name, sample_clades in sample_to_clades_present.items():
        for clade in sample_clades:
            clade_to_samples[clade].append(sample_name)


    prevalent_clades = set((clade for clade, clade_samples in clade_to_samples.items()
                            if len(clade_samples) >= config['min_cohort_sgb_prevalence']))

    if args.clade is not None:
        if args.clade not in prevalent_clades:
            warning(f'The clade {args.clade} is not prevalent in the samples, nothing to do')
            return
        prevalent_clades = {args.clade}
    else:
        info(f'There are {len(prevalent_clades)} prevalent SGBs')


    info('Loading rare alleles database')
    rare_alleles_db_dir = args.database_rare_alleles

    if args.clade is not None:
        f = rare_alleles_db_dir / f'{args.clade}.zip'
        if not f.exists():
            error(f'The clade {args.clade} is not in the rare alleles database')
            return
        fs = [f]
    else:
        fs = list(rare_alleles_db_dir.iterdir())

    rare_alleles_db = {}
    prevalent_clades_rare = set()
    for ac_path in tqdm(fs):
        suffix = '.zip'
        assert ac_path.name.endswith(suffix)
        sgb_id = ac_path.name[:-len(suffix)]
        if sgb_id not in prevalent_clades:
            continue

        prevalent_clades_rare.add(sgb_id)
        with zipfile.ZipFile(ac_path, 'r') as f:
            for fi in f.infolist():
                if fi.filename == 'metadata.json':
                    with f.open(fi, 'r') as fm:
                        metadata = json.load(fm)
                    if 'database_name' in metadata and metadata['database_name'] != mp_db_info.db_name:
                        raise Exception(f'Sample {sample_name} has been generated with a different MetaPhlAn database'
                                        f' {metadata["database_name"]} and not the one provided {mp_db_info.db_name}')
                elif fi.filename.endswith('.npy'):
                    m = fi.filename[:-len('.npy')]
                    with f.open(fi, 'r') as fm:
                        allele_prevalences_abs = np.load(fm)
                        tot_counts = allele_prevalences_abs.sum(axis=0)
                        with np.errstate(invalid='ignore'):  # 0/0 case will make NaN and be filtered as False below
                            allele_prevalences_rel = allele_prevalences_abs / tot_counts
                        rare_alleles_db[m] = (allele_prevalences_rel < config['max_rare_allele_freq']) & \
                                             (tot_counts >= config['min_db_pos_coverage'])
                else:
                    warning(f'Ignoring unrecognized file in the allele counts zip file {fi.filename}')


    info(f'There are {len(prevalent_clades)} prevalent SGBs out of which {len(prevalent_clades_rare)} are also '
         f'in the rare database.')



    info(f'Calculating all-vs-all matches using {args.threads_cpu} thread(s)')
    run.metaphlan_db = mp_db_info
    run.prevalent_clades = prevalent_clades
    run.config = config
    run.sample_name_to_ac_path = sample_name_to_ac_path
    run.rare_alleles_db = rare_alleles_db

    result_rows = []
    if args.in_memory:
        run_in_memory.sample_to_clades_present = sample_to_clades_present
        run_in_memory.sample_to_allele_counts = sample_to_allele_counts
        with mpp.Pool(args.threads_cpu) as pool:
            for result_rows_batch in tqdm(pool.imap_unordered(run_in_memory, sample_name_to_ac_path.keys()),
                                          total=len(sample_name_to_ac_path)):
                result_rows.extend(result_rows_batch)
    else:
        with mpp.Pool(args.threads_cpu) as pool:
            for result_rows_batch in tqdm(pool.imap_unordered(run, sample_name_to_ac_path.keys()),
                                          total=len(sample_name_to_ac_path)):
                result_rows.extend(result_rows_batch)


    info(f'Creating a dataframe with {len(result_rows):,} rows')
    df_results = pd.DataFrame(result_rows, columns=OUTPUT_HEADER)


    info('Saving results file as parquet')
    df_results.to_parquet(args.output_dir / 'results.parquet', index=False)

    # info('Saving results file as tsv')
    # df_results.to_csv(args.output_dir / 'results.tsv.gz', sep='\t', index=False)


    info(f'Finished strain tracking. The results are stored at {args.output_dir}')



if __name__ == '__main__':
    main()
