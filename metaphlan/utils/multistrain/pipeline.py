import base64
import bz2
import gzip
import json
import os
import pathlib
import pickle
from collections import Counter
from typing import Sequence

import numpy as np
import numpy.typing as npt
import pandas as pd
import pysam

from ...utils import info, error, warning, info_debug, global_flags
from . import utils
from .genotype_resolving import compute_genotypes
from .get_strainphlan_markers import get_strainphlan_markers
from .linkage import calculate_linkage, linkage_merging
from .model_fitting import fit_model
from .prepare_sam import prepare_sam
from .snp_calling import run_pileup, filter_loci_snp_call, PileupResult


def try_load_bam(bam_path):
    """

    :param pathlib.Path bam_path:
    :return:
    """

    try:
        return pysam.AlignmentFile(bam_path)
    except ValueError as e:
        if 'file has no sequences defined' in str(e):
            error(f'There are no reads in the BAM file {bam_path}, skipping the sample')
            return None
        else:
            raise e


def step_bam(sample_path, bam_path, bai_path, output_dir, db_name, config):
    """

    :param pathlib.Path sample_path:
    :param pathlib.Path bam_path:
    :param pathlib.Path bai_path:
    :param pathlib.Path output_dir:
    :param dict config:
    :return:
    """
    read_lens = prepare_sam(sample_path, bam_path, output_dir, db_name, config)

    bai_path.unlink(missing_ok=True)
    pysam.index(str(bam_path))

    sam_file = try_load_bam(bam_path)

    return sam_file, read_lens


def save_bam(read_lens_path, read_lens):
    """

    :param pathlib.Path read_lens_path:
    :param Sequence[int] read_lens:
    :return:
    """
    read_lens_path_tmp = read_lens_path.with_name(read_lens_path.name + '.tmp')
    with gzip.open(read_lens_path_tmp, 'wt') as f:
        f.write('\n'.join(map(str, read_lens)))
    read_lens_path_tmp.rename(read_lens_path)


def load_bam(bam_path, read_lens_path):
    """

    :param pathlib.Path bam_path:
    :param pathlib.Path read_lens_path:
    :return:
    """
    with gzip.open(read_lens_path, 'rt') as f:
        read_lens = [int(line.strip()) for line in f if line.strip()]

    sam_file = try_load_bam(bam_path)

    return sam_file, read_lens


def save_pileup(pr, pileup_path, seq_error_path, null_err_rate_path):
    """

    :param PileupResult pr:
    :param pathlib.Path pileup_path:
    :param pathlib.Path seq_error_path:
    :param pathlib.Path null_err_rate_path:
    :return:
    """
    pileup_path_tmp = pileup_path.with_name(pileup_path.name + '.tmp')
    with gzip.open(pileup_path_tmp, 'wb') as f:
        for m in pr.base_frequencies.keys():
            bfs = {b: base64.b64encode(a) for b, a in pr.base_frequencies[m].items()}
            er = base64.b64encode(pr.avg_error_rates[m])
            line = b'\t'.join([m.encode()] + [bfs[b] for b in 'ACTG'] + [er]) + b'\n'
            f.write(line)

    pileup_path_tmp.rename(pileup_path)

    seq_error_path_tmp = seq_error_path.with_name(seq_error_path.name + '.tmp')
    pr.df_seq_errors.reset_index().to_csv(seq_error_path_tmp, index=False)
    seq_error_path_tmp.rename(seq_error_path)

    null_err_rate_path_tmp = null_err_rate_path.with_name(null_err_rate_path.name + '.tmp')
    with open(null_err_rate_path_tmp, 'w') as f:
        f.write('null_err_rate\t' + str(pr.null_err_rate) + '\n')
        f.write('substitutions\t' + json.dumps(pr.substitutions) + '\n')
        f.write('substitutions_per_position\t' + json.dumps(pr.substitutions_per_position) + '\n')
        f.write('contexts\t' + json.dumps(pr.contexts) + '\n')
    null_err_rate_path_tmp.rename(null_err_rate_path)


def load_pileup(pileup_path, null_err_rate_path):
    """

    :param pathlib.Path pileup_path:
    :param pathlib.Path null_err_rate_path:
    :return:
    """
    with gzip.open(pileup_path, 'rb') as f:
        base_frequencies = {}
        avg_error_rates = {}
        marker_to_length = {}
        for line in f:
            bf = {}
            m, bf['A'], bf['C'], bf['T'], bf['G'], er = line.strip().split(b'\t')
            m = m.decode()
            base_frequencies[m] = {b: np.frombuffer(base64.decodebytes(x), dtype=np.uint16) for b, x in bf.items()}
            avg_error_rates[m] = np.frombuffer(base64.decodebytes(er), dtype=np.float64)
            marker_to_length[m] = len(avg_error_rates[m])

    with open(null_err_rate_path, 'r') as f:
        desc, null_err_rate = next(f).rstrip('\n').split('\t')
        assert desc == "null_err_rate"
        if null_err_rate == "None":
            null_err_rate = None
        else:
            null_err_rate = float(null_err_rate)

        desc, subs = next(f).rstrip('\n').split('\t')
        assert desc == "substitutions"
        if subs == "None":
            subs = None
        else:
            subs = {k: Counter(v) for k, v in json.loads(subs).items()}

        desc, subs_pos = next(f).rstrip('\n').split('\t')
        assert desc == "substitutions_per_position"
        if subs_pos == "None":
            subs_pos = None
        else:
            subs_pos = {k: {k2: Counter({int(k3): int(v3) for k3, v3 in v2.items()}) for k2, v2 in v.items()}
                        for k, v in json.loads(subs_pos).items()}

        desc, contexts = next(f).rstrip('\n').split('\t')
        assert desc == "contexts"
        if contexts == "None":
            contexts = None
        else:
            contexts = {k: {k2: Counter(v2) for k2, v2 in v.items()} for k, v in json.loads(contexts).items()}

    return PileupResult(base_frequencies, avg_error_rates, marker_to_length, None, null_err_rate, subs, subs_pos,
                        contexts)


def markers_and_species_filtering(base_frequencies, marker_to_clade_db, marker_to_ext, clade_to_n_markers_db, config):
    """

    :param dict[str, dict[str, npt.NDArray]] base_frequencies:
    :param dict[str, str] marker_to_clade_db:
    :param dict[str, Sequence[str]] marker_to_ext:
    :param dict[str, int] clade_to_n_markers_db:
    :param dict config:
    :return:
    """
    markers_without_clade = [m for m in base_frequencies.keys() if m not in marker_to_clade_db]
    if len(markers_without_clade) > 0:
        raise Exception(f'There are {len(markers_without_clade)} markers without a clade:'
                        f' {utils.shorten_text(",".join(markers_without_clade))}')

    # Filter quasi markers with external hits (using >= 1 hit threshold)
    all_markers = sorted(base_frequencies.keys())
    clade_to_n_markers = Counter([marker_to_clade_db[m] for m in all_markers])
    markers_non_quasi = [m for m in all_markers if len(marker_to_ext[m]) == 0]
    markers_quasi = [m for m in all_markers if len(marker_to_ext[m]) > 0]
    markers_quasi_filtered = [m for m in markers_quasi
                              if not any(clade_to_n_markers[ext_sgb] / clade_to_n_markers_db[ext_sgb] >=
                                         config['quasi_marker_frac'] for ext_sgb in marker_to_ext[m])]
    markers_filtered = markers_non_quasi + markers_quasi_filtered

    # Filter markers for breadth
    marker_to_breadth = {m: (np.array([base_frequencies[m][b] for b in 'ACTG']).sum(axis=0) > 0).mean()
                         for m in markers_filtered}
    markers_present = set([m for m, breadth in marker_to_breadth.items() if breadth >= config['min_marker_breadth']])

    # Filter clades / SGBs with enough markers
    clades_to_n_markers_present = Counter([marker_to_clade_db[m] for m in markers_present])
    clades_present = [clade for clade, c in clades_to_n_markers_present.items()
                      if c >= config['min_markers_abs']
                      and c / clade_to_n_markers_db[clade] >= config['min_markers_rel']]

    return markers_present, clades_present


def step_reconstructed_markers(output_dir, config, sam_file, pr, read_lens, marker_to_clade_db, marker_to_ext):
    """

    :param pathlib.Path output_dir:
    :param dict config:
    :param sam_file:
    :param PileupResult pr:
    :param Sequence read_lens:
    :param dict[str, str] marker_to_clade_db:
    :param dict[str, Sequence[str]] marker_to_ext:
    :return:
    """


    s_marker_to_clade_db = pd.Series(marker_to_clade_db)
    clade_to_markers_db = s_marker_to_clade_db.groupby(s_marker_to_clade_db).groups
    clade_to_n_markers_db = {k: len(v) for k, v in clade_to_markers_db.items()}

    markers_present, clades_present = markers_and_species_filtering(pr.base_frequencies, marker_to_clade_db,
                                                                    marker_to_ext, clade_to_n_markers_db, config)


    result_rows = {}
    consensuses_maj = {}
    consensuses_min = {}
    for sgb_id in clades_present:
        info_debug(sgb_id, 'Preparing loci')
        loci_rows = []
        for m in clade_to_markers_db[sgb_id]:
            if m not in markers_present:
                continue

            for pos in range(pr.marker_to_length[m]):
                bfs = {b: pr.base_frequencies[m][b][pos] for b in 'ACTG'}
                bfs = Counter({k: v for k, v in bfs.items() if v > 0})
                if bfs.total() == 0:
                    continue

                if pr.null_err_rate is not None:
                    pos_err_rate = pr.null_err_rate
                else:
                    pos_err_rate = pr.avg_error_rates[m][pos]

                loci_rows.append({
                    'marker': m,
                    'pos': pos + 1,
                    'error_rate': pos_err_rate,
                    'base_frequencies': bfs,
                })
        if len(loci_rows) == 0:
            warning('No loci rows despite all the filtering (possibly a bug in the code)')
            continue

        df_loci_sgb = pd.DataFrame(loci_rows)
        marker_to_length_sgb = {m: pr.marker_to_length[m] for m in
                                df_loci_sgb['marker'].unique()}  # downsize the possibly huge dict
        result_row, data, df_loci_sgb, df_loci_sgb_polyallelic = filter_loci_snp_call(df_loci_sgb, marker_to_length_sgb,
                                                                                      read_lens, config)

        df_loci_sgb_filtered = df_loci_sgb.query('filtered').copy()

        if global_flags.debug:
            result_row['est_err_rate'] = 1 - df_loci_sgb['max_frequency'].sum() / df_loci_sgb['base_coverage'].sum()
            result_row['est_err_rate_filtered'] = 1 - df_loci_sgb_filtered['max_frequency'].sum() / \
                                                  df_loci_sgb_filtered['base_coverage'].sum()
        if 'n_markers_polyallelic_significant' in result_row \
                and result_row['n_markers_polyallelic_significant'] >= config['multistrain_min_markers_abs'] \
                and result_row['n_markers_polyallelic_significant'] / clade_to_n_markers_db[sgb_id] \
                >= config['multistrain_min_markers_rel']:
            result_row['multi_strain'] = True
            info_debug(sgb_id, 'is multistrain')
            df_loci_biallelic_significant = df_loci_sgb_polyallelic.query('biallelic_significant')
            info_debug(sgb_id, 'Creating and merging linkage')
            sgb_linkage = calculate_linkage(df_loci_biallelic_significant, sam_file,
                                            config)  # (m1, p1, m2, p2) => (b1 + b2) => count
            merging_results_before, merging_results = linkage_merging(df_loci_biallelic_significant, sgb_linkage,
                                                                      config)
            sgb_nps = merging_results[0]
            info_debug(sgb_id, 'Fitting model')
            fit_model(sgb_id, df_loci_sgb_filtered, result_row, config)
        else:
            result_row['multi_strain'] = False
            info_debug(sgb_id, 'is not multistrain')
            sgb_nps = []
            merging_results_before = ([], None)


        info_debug(sgb_id, 'Generating genotype by maximizing per-position probabilities')
        consensuses_major, consensuses_minor, qualities_major, qualities_minor =\
            compute_genotypes(df_loci_sgb_filtered, sgb_nps, result_row, marker_to_length_sgb)

        strain_resolved_markers_sgb = [{
            'sgb_id': sgb_id,
            'marker': m,
            'sequence_maj': consensuses_major[m].decode(),
            'sequence_min': consensuses_minor[m].decode(),
            'qualities_maj': utils.qualities_to_phred(qualities_major[m]),
            'qualities_min': utils.qualities_to_phred(qualities_minor[m]),
        } for m in consensuses_major.keys()]


        if global_flags.debug:
            info_debug(sgb_id, 'Saving debug files')
            output_dir_per_sgb = output_dir / 'per_sgb_data'
            output_dir_per_sgb.mkdir(parents=True, exist_ok=True)
            df_loci_sgb.to_csv(output_dir_per_sgb / f'df_loci_{sgb_id}.tsv', sep='\t', index=False)
            with open(output_dir_per_sgb / f'data_{sgb_id}.pic', 'wb') as f:
                pickle.dump(data, f)

        consensuses_maj_sgb, consensuses_min_sgb = get_strainphlan_markers(strain_resolved_markers_sgb, result_row,
                                                                           merging_results_before, config)
        consensuses_maj.update(consensuses_maj_sgb)
        consensuses_min.update(consensuses_min_sgb)

        result_rows[sgb_id] = result_row


    df_results = pd.DataFrame.from_dict(result_rows, orient='index')

    return df_results, consensuses_maj, consensuses_min


def save_reconstructed_markers(output_results, output_major, output_minor, mp_version, df_results, consensuses_maj,
                               consensuses_min):
    """

    :param pathlib.Path output_results:
    :param pathlib.Path output_major:
    :param pathlib.Path output_minor:
    :param str mp_version:
    :param pd.DataFrame df_results:
    :param dict consensuses_maj:
    :param dict consensuses_min:
    :return:
    """

    consensuses_maj = list(consensuses_maj.values())
    consensuses_min = list(consensuses_min.values())

    with bz2.open(output_minor, 'wt') as f:
        json.dump({'database_name': mp_version, 'consensus_markers': consensuses_min}, f, indent=2)

    with bz2.open(output_major, 'wt') as f:
        json.dump({'database_name': mp_version, 'consensus_markers': consensuses_maj}, f, indent=2)

    df_results.to_csv(output_results, sep='\t')



def run(sample_path, output_dir, config, target, save_bam_file, reuse, db_name, marker_to_clade,
        marker_to_ext):
    """

    :param pathlib.Path sample_path:
    :param pathlib.Path output_dir:
    :param dict config:
    :param str target:
    :param bool save_bam_file:
    :param str reuse:
    :param str db_name:
    :param dict[str, str] marker_to_clade:
    :param dict[str, Sequence[str]] marker_to_ext:
    :return:
    """


    if sample_path.suffix == '.bz2':
        sample_name = sample_path.with_suffix('').stem
    else:
        sample_name = sample_path.stem

    output_dir = output_dir / sample_name

    pileup_path_before = output_dir / 'pileup_nonfilt.tsv.gz'
    pileup_path_after = output_dir / 'pileup.tsv.gz'
    read_lens_path = output_dir / 'read_lengths.txt.gz'
    null_err_rate_path_before = output_dir / 'null_err_rate_nonfilt.txt'
    null_err_rate_path_after = output_dir / 'null_err_rate.txt'
    seq_error_path_before = output_dir / 'seq_error_profile_nonfilt.csv'
    seq_error_path_after = output_dir / 'seq_error_profile.csv'
    bam_path = output_dir / f'{sample_name}.sorted.bam'
    bai_path = bam_path.with_suffix(bam_path.suffix + '.bai')
    output_results = output_dir / 'results.tsv'
    output_major = output_dir / f'{sample_name}_major.json.bz2'
    output_minor = output_dir / f'{sample_name}_minor.json.bz2'

    info(f'Starting sample {sample_name}')
    output_dir.mkdir(parents=True, exist_ok=True)

    if global_flags.debug or save_bam_file:
        files_to_remove = []
    else:
        files_to_remove = [bam_path, bai_path]
    dirs_to_remove = []


    pileup_target_files = [pileup_path_before, pileup_path_after, null_err_rate_path_before,
                           null_err_rate_path_after]  # files to consider the step to be satisfied
    bam_target_files = [bam_path, bai_path, read_lens_path]  # files to consider the step to be satisfied
    pileup_output_files = [pileup_path_before, pileup_path_after, null_err_rate_path_before, null_err_rate_path_after,
                           read_lens_path, bam_path, bai_path]  # files to be able to reuse this step
    rm_target_files = [output_results, output_major, output_minor]  # files to consider the step to be satisfied

    def meta_step_pileup():
        if reuse in ['all', 'bam'] and all(map(pathlib.Path.is_file, bam_target_files)):
            info(f'Reusing existing preprocessed BAM file for sample {sample_name}')
            sam_file_, read_lens_ = load_bam(bam_path, read_lens_path)
        else:
            info(f'Preprocessing the BAM file for sample {sample_name}')
            sam_file_, read_lens_ = step_bam(sample_path, bam_path, bai_path, output_dir, db_name, config)
            save_bam(read_lens_path, read_lens_)

        if sam_file_ is None:
            return None, None, None, None, None, None


        info(f'Calculating pileup for sample {sample_name}')
        pr_before_, pr_after_ = run_pileup(sam_file_, config)
        info_debug('Saving pileup')
        save_pileup(pr_before_, pileup_path_before, seq_error_path_before, null_err_rate_path_before)
        save_pileup(pr_after_, pileup_path_after, seq_error_path_after, null_err_rate_path_after)

        return sam_file_, read_lens_, pr_after_


    if target == 'pileup':
        if reuse in ['all', 'pileup'] and all(map(pathlib.Path.is_file, pileup_target_files)):
            info(f'Pileup for sample {sample_name} is already present, nothing to do')
            return True
        else:
            sam_file, *_ = meta_step_pileup()
            if sam_file is None:
                return False
    elif target == 'reconstructed_markers':
        if reuse == 'all' and all(map(pathlib.Path.is_file, rm_target_files)):
            info(f'Reconstructed markers for sample {sample_name} are already present, nothing to do.')
            return True
        else:
            if reuse in ['all', 'pileup'] and all(map(pathlib.Path.is_file, pileup_output_files)):
                info(f'Loading existing pileup for sample {sample_name}')
                sam_file, read_lens = load_bam(bam_path, read_lens_path)
                if sam_file is None:
                    return False
                pr = load_pileup(pileup_path_after, null_err_rate_path_after)
            else:
                sam_file, read_lens, pr = meta_step_pileup()

                if sam_file is None:
                    return False


            info(f'Running marker reconstruction for sample {sample_name}')
            df_results, consensuses_maj, consensuses_min = \
                step_reconstructed_markers(output_dir, config, sam_file, pr, read_lens, marker_to_clade, marker_to_ext)

            save_reconstructed_markers(output_results, output_major, output_minor, db_name, df_results,
                                       consensuses_maj, consensuses_min)


    info(f'Cleaning intermediate files for sample {sample_name}')
    for f in files_to_remove:
        os.unlink(f)

    for d in dirs_to_remove:
        os.rmdir(d)



    info(f'Finished sample {sample_name}')
    return True
