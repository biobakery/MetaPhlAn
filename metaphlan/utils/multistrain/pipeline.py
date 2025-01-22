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

from metaphlan.utils import info, error, warning
from . import utils
from .genotype_resolving import calculate_strain_base_probabilities, resolve_strains
from .get_strainphlan_markers import get_strainphlan_markers
from .linkage import calculate_linkage, linkage_merging
from .model_fitting import fit_model
from .prepare_sam import prepare_sam
from .snp_calling import run_pileup, filter_loci_snp_call
from ..util_fun import info_debug, global_flags


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


def step_pileup(sample_path, bam_path, bai_path, output_dir, config):
    """

    :param pathlib.Path sample_path:
    :param pathlib.Path bam_path:
    :param pathlib.Path bai_path:
    :param pathlib.Path output_dir:
    :param dict config:
    :return:
    """
    info('Preparing the BAM file')
    read_lens = prepare_sam(sample_path, bam_path, output_dir, config)

    bai_path.unlink(missing_ok=True)
    pysam.index(str(bam_path))

    sam_file = try_load_bam(bam_path)

    if sam_file is None:
        return None, None, None, None

    info('Running pileup')
    base_frequencies, error_rates, marker_to_length = run_pileup(sam_file, config)

    return sam_file, base_frequencies, error_rates, marker_to_length, read_lens


def save_pileup(pileup_path, base_frequencies, error_rates, read_lens_path, read_lens):
    """

    :param pathlib.Path pileup_path:
    :param dict[str, dict[str, npt.NDArray]] base_frequencies:
    :param dict[str, npt.NDArray] error_rates:
    :param pathlib.Path read_lens_path:
    :param Sequence[int] read_lens:
    :return:
    """
    pileup_path_tmp = pileup_path.with_name(pileup_path.name + '.tmp')
    with gzip.open(pileup_path_tmp, 'wb') as f:
        for m in base_frequencies.keys():
            bfs = {b: base64.b64encode(a) for b, a in base_frequencies[m].items()}
            er = base64.b64encode(error_rates[m])
            line = b'\t'.join([m.encode()] + [bfs[b] for b in 'ACTG'] + [er]) + b'\n'
            f.write(line)

    pileup_path_tmp.rename(pileup_path)

    read_lens_path_tmp = read_lens_path.with_name(read_lens_path.name + '.tmp')
    with gzip.open(read_lens_path_tmp, 'wt') as f:
        f.write('\n'.join(map(str, read_lens)))
    read_lens_path_tmp.rename(read_lens_path)


def load_pileup(pileup_path, read_lens_path):
    """

    :param pathlib.Path pileup_path:
    :param pathlib.Path read_lens_path:
    :return:
    """
    with gzip.open(pileup_path, 'rb') as f:
        base_frequencies = {}
        error_rates = {}
        marker_to_length = {}
        for line in f:
            bf = {}
            m, bf['A'], bf['C'], bf['T'], bf['G'], er = line.strip().split(b'\t')
            m = m.decode()
            base_frequencies[m] = {b: np.frombuffer(base64.decodebytes(x), dtype=np.uint16) for b, x in bf.items()}
            error_rates[m] = np.frombuffer(base64.decodebytes(er), dtype=np.float64)
            marker_to_length[m] = len(error_rates[m])

    with gzip.open(read_lens_path, 'rt') as f:
        read_lens = [int(line.strip()) for line in f if line.strip()]


    return base_frequencies, error_rates, marker_to_length, read_lens


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


def step_reconstructed_markers(output_dir_linkages, config, sam_file, base_frequencies, error_rates,
                               read_lens, marker_to_length, marker_to_clade_db, marker_to_ext):
    """

    :param pathlib.Path output_dir_linkages:
    :param dict config:
    :param sam_file:
    :param dict[str, dict[str, npt.NDArray]] base_frequencies:
    :param dict[str, npt.NDArray] error_rates:
    :param Sequence read_lens:
    :param dict[str, int] marker_to_length:
    :param dict[str, str] marker_to_clade_db:
    :param dict[str, Sequence[str]] marker_to_ext:
    :return:
    """


    output_dir_linkages.mkdir(parents=True, exist_ok=True)

    s_marker_to_clade_db = pd.Series(marker_to_clade_db)
    clade_to_markers_db = s_marker_to_clade_db.groupby(s_marker_to_clade_db).groups
    clade_to_n_markers_db = {k: len(v) for k, v in clade_to_markers_db.items()}

    markers_present, clades_present = markers_and_species_filtering(base_frequencies, marker_to_clade_db, marker_to_ext,
                                                                    clade_to_n_markers_db, config)


    result_rows = {}
    consensuses_maj = {}
    consensuses_min = {}
    for sgb_id in clades_present:
        info_debug(sgb_id, 'Preparing loci')
        loci_rows = []
        for m in clade_to_markers_db[sgb_id]:
            if m not in markers_present:
                continue

            for pos in range(marker_to_length[m]):
                bfs = {b: base_frequencies[m][b][pos] for b in 'ACTG'}
                bfs = Counter({k: v for k, v in bfs.items() if v > 0})
                if bfs.total() == 0:
                    continue
                loci_rows.append({
                    'marker': m,
                    'pos': pos + 1,
                    'error_rate': error_rates[m][pos],
                    'base_frequencies': bfs,
                })
        if len(loci_rows) == 0:
            warning('No loci rows despite all the filtering (possibly a bug in the code)')
            continue

        df_loci_sgb = pd.DataFrame(loci_rows)
        marker_to_length_sgb = {m: marker_to_length[m] for m in
                                df_loci_sgb['marker'].unique()}  # downsize the possibly huge dict
        result_row, data, df_loci_sgb, df_loci_sgb_polyallelic = filter_loci_snp_call(df_loci_sgb, marker_to_length_sgb,
                                                                                      read_lens, config)

        if global_flags.debug:
            df_loci_sgb_filtered = df_loci_sgb.query('filtered')
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
            fit_model(sgb_id, df_loci_sgb, result_row, config)
        else:
            result_row['multi_strain'] = False
            info_debug(sgb_id, 'is not multistrain')
            sgb_nps = []
            merging_results_before = ([], None)

        
        info_debug(sgb_id, 'Calculating probas')
        df_loci_sgb = calculate_strain_base_probabilities(df_loci_sgb, sgb_nps, result_row, config)
        info_debug(sgb_id, 'Generating genotypes')
        consensuses_major, consensuses_minor, qualities_major, qualities_minor = resolve_strains(df_loci_sgb,
                                                                                                 marker_to_length_sgb)
        strain_resolved_markers_sgb = [{
            'sgb_id': sgb_id,
            'marker': m,
            'sequence_maj': consensuses_major[m].decode(),
            'sequence_min': consensuses_minor[m].decode(),
            'qualities_maj': utils.qualities_to_phred(qualities_major[m]),
            'qualities_min': utils.qualities_to_phred(qualities_minor[m]),
        } for m in consensuses_major.keys()]


        if global_flags.debug:
            info_debug(sgb_id, 'Saving files')
            df_loci_sgb.to_csv(output_dir_linkages / f'df_loci_{sgb_id}.tsv', sep='\t', index=False)
            with open(output_dir_linkages / f'data_{sgb_id}.pic', 'wb') as f:
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



def run(sample_path, output_dir, config, target, save_bam_file, reuse, mp_version, marker_to_clade,
        marker_to_ext):
    """

    :param pathlib.Path sample_path:
    :param pathlib.Path output_dir:
    :param dict config:
    :param str target:
    :param bool save_bam_file:
    :param str reuse:
    :param str mp_version:
    :param dict[str, str] marker_to_clade:
    :param dict[str, Sequence[str]] marker_to_ext:
    :return:
    """


    if sample_path.suffix == '.bz2':
        sample_name = sample_path.with_suffix('').stem
    else:
        sample_name = sample_path.stem

    output_dir = output_dir / sample_name

    pileup_path = output_dir / 'pileup.tsv.gz'
    read_lens_path = output_dir / 'read_lengths.txt.gz'
    bam_path = output_dir / f'{sample_name}.sorted.bam'
    bai_path = bam_path.with_suffix(bam_path.suffix + '.bai')
    output_dir_linkages = output_dir / 'per_sgb_linkages'
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


    pileup_target_files = [pileup_path]  # files to consider the step to be satisfied
    pileup_output_files = [pileup_path, read_lens_path, bam_path, bai_path]  # files to be able to reuse this step
    rm_target_files = [output_results, output_major, output_minor]  # files to consider the step to be satisfied

    if target == 'pileup':
        if reuse in ['all', 'pileup'] and all(map(pathlib.Path.is_file, pileup_target_files)):
            info(f'Pileup for sample {sample_name} is already present, nothing to do')
            return True

        sam_file, base_frequencies, error_rates, marker_to_length, read_lens = step_pileup(sample_path, bam_path,
                                                                                           bai_path, output_dir, config)
        if sam_file is None:
            return False
        save_pileup(pileup_path, base_frequencies, error_rates, read_lens_path, read_lens)

    elif target == 'reconstructed_markers':
        if reuse == 'all' and all(map(pathlib.Path.is_file, rm_target_files)):
            info(f'Reconstructed markers for sample {sample_name} are already present, nothing to do.')
            return True

        if reuse in ['all', 'pileup'] and all(map(pathlib.Path.is_file, pileup_output_files)):
            info(f'Loading existing pileup for sample {sample_name}')
            base_frequencies, error_rates, marker_to_length, read_lens = load_pileup(pileup_path, read_lens_path)
            sam_file = try_load_bam(bam_path)
            if sam_file is None:
                return False
        else:
            info(f'Calculating pileup for sample {sample_name}')
            sam_file, base_frequencies, error_rates, marker_to_length, read_lens = step_pileup(sample_path, bam_path,
                                                                                               bai_path, output_dir,
                                                                                               config)
            if sam_file is None:
                return False

            save_pileup(pileup_path, base_frequencies, error_rates, read_lens_path, read_lens)



        info(f'Running marker reconstruction for sample {sample_name}')
        df_results, consensuses_maj, consensuses_min = \
            step_reconstructed_markers(output_dir_linkages, config, sam_file, base_frequencies, error_rates,
                                       read_lens, marker_to_length, marker_to_clade, marker_to_ext)

        save_reconstructed_markers(output_results, output_major, output_minor, mp_version, df_results, consensuses_maj,
                                   consensuses_min)

        if not global_flags.debug:
            for sgb in df_results.index:
                files_to_remove.append(output_dir_linkages / f'mrb_{sgb}.pic')
            dirs_to_remove.append(output_dir_linkages)


    info(f'Cleaning intermediate files for sample {sample_name}')
    for f in files_to_remove:
        os.unlink(f)

    for d in dirs_to_remove:
        os.rmdir(d)



    info(f'Finished sample {sample_name}')
    return True
