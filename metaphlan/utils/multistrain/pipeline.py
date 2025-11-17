import base64
import bz2
import gzip
import io
import json
import pathlib
import pickle
import zipfile
from collections import Counter
from typing import Sequence

import numpy as np
import numpy.typing as npt
import pandas as pd
import pysam

from .utils import MetaphlanDBInfo
from ...utils import info, error, warning, info_debug, global_flags
from . import utils
from .genotype_resolving import compute_genotypes
from .get_strainphlan_markers import get_strainphlan_markers
from .linkage import calculate_linkage, linkage_merging
from .model_fitting import fit_model
from .prepare_sam import prepare_sam, get_read_lens_from_bam
from .snp_calling import run_pileup, PileupResult, filter_loci_snp_call


ACTG = 'ACTG'


def try_load_bam(bam_path):
    """

    :param pathlib.Path bam_path:
    :return:
    """

    try:
        return pysam.AlignmentFile(bam_path)
    except ValueError as e:
        se = str(e)
        if 'file has no sequences defined' in se or 'file does not contain alignment data' in se:
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
    :param str db_name:
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



def load_read_lens(read_lens_path):
    """

    :param pathlib.Path read_lens_path:
    :return:
    """
    with gzip.open(read_lens_path, 'rt') as f:
        read_lens = [int(line.strip()) for line in f if line.strip()]

    return read_lens


def load_bam(bam_path):
    """

    :param pathlib.Path bam_path:
    :return:
    """
    sam_file = try_load_bam(bam_path)

    return sam_file


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
            bfs = {b: base64.b64encode(pr.base_frequencies[m][i]) for i, b in enumerate(ACTG)}
            er = base64.b64encode(pr.avg_error_rates[m])
            line = b'\t'.join([m.encode()] + [bfs[b] for b in ACTG] + [er]) + b'\n'
            f.write(line)

    pileup_path_tmp.rename(pileup_path)

    if pr.df_seq_errors is not None:
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
            base_frequencies[m] = np.array([np.frombuffer(base64.decodebytes(bf[b]), dtype=np.uint16) for b in ACTG])
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


def markers_and_species_filtering(base_frequencies, mp_db_info, config):
    """

    :param dict[str, npt.NDArray] base_frequencies:
    :param MetaphlanDBInfo mp_db_info:
    :param dict config:
    :return:
    """
    markers_without_clade = [m for m in base_frequencies.keys() if m not in mp_db_info.marker_to_clade]
    if len(markers_without_clade) > 0:
        raise Exception(f'There are {len(markers_without_clade)} markers without a clade:'
                        f' {utils.shorten_text(",".join(markers_without_clade))}')

    # Filter quasi markers with external hits (using >= 1 hit threshold)
    all_markers = sorted(base_frequencies.keys())
    clade_to_n_markers = Counter([mp_db_info.marker_to_clade[m] for m in all_markers])
    markers_non_quasi = [m for m in all_markers if len(mp_db_info.marker_to_ext[m]) == 0]
    markers_quasi = [m for m in all_markers if len(mp_db_info.marker_to_ext[m]) > 0]
    markers_quasi_filtered = [m for m in markers_quasi
                              if not any(clade_to_n_markers[ext_sgb] / mp_db_info.clade_to_n_markers[ext_sgb] >=
                                         config['quasi_marker_frac'] for ext_sgb in mp_db_info.marker_to_ext[m])]
    markers_filtered = markers_non_quasi + markers_quasi_filtered

    # Filter markers for breadth
    marker_to_breadth = {m: (base_frequencies[m].sum(axis=0) > 0).mean()
                         for m in markers_filtered}
    markers_present = [m for m, breadth in marker_to_breadth.items() if breadth >= config['min_marker_breadth']]

    # Filter clades / SGBs with enough markers
    clades_to_n_markers_present = Counter([mp_db_info.marker_to_clade[m] for m in markers_present])
    clades_present = [clade for clade, c in clades_to_n_markers_present.items()
                      if c >= config['min_markers_abs']
                      and c / mp_db_info.clade_to_n_markers[clade] >= config['min_markers_rel']]

    # Keep only markers from kept clades
    markers_present = set([m for m in markers_present if mp_db_info.marker_to_clade[m] in clades_present])

    return markers_present, clades_present


def step_reconstructed_markers(output_dir, config, sam_file, pr, read_lens, mp_db_info, allele_counts_only):
    """
    :param pathlib.Path output_dir:
    :param dict config:
    :param sam_file:
    :param PileupResult pr:
    :param Sequence[int] read_lens:
    :param MetaphlanDBInfo mp_db_info:
    :param bool allele_counts_only:

    :return:
    """

    assert sam_file is not None


    avg_read_len: float = np.mean(read_lens)

    markers_present, clades_present = markers_and_species_filtering(pr.base_frequencies, mp_db_info, config)

    consensuses_maj = {}
    consensuses_min = {}
    result_rows = {}
    bfs_filtered = {}
    for sgb_id in clades_present:
        info_debug(sgb_id, 'Preparing loci')

        clade_markers_present = [m for m in mp_db_info.clade_to_markers[sgb_id] if m in markers_present]
        assert len(clade_markers_present) > 0

        bfs = {}
        err_rates = {}
        marker_to_length_sgb = {}
        for m in clade_markers_present:
            bfs[m] = pr.base_frequencies[m]

            le = bfs[m].shape[1]
            marker_to_length_sgb[m] = le

            if config['error_rate'] == 'estimate' and pr.null_err_rate is not None:
                err_rates[m] = np.array([pr.null_err_rate] * le)
            elif config['error_rate'] == 'phred' or (config['error_rate'] == 'estimate' and pr.null_err_rate is None):
                err_rates[m] = pr.avg_error_rates[m]
            else:
                err_rates[m] = np.array([float(config['error_rate'])] * le)

        result_row, position_mask, polyallelic_significant_masks = filter_loci_snp_call(bfs, err_rates, avg_read_len,
                                                                                        config)

        result_rows[sgb_id] = result_row

        if 'n_markers_polyallelic_significant' in result_row \
                and result_row['n_markers_polyallelic_significant'] >= config['multistrain_min_markers_abs'] \
                and result_row['n_markers_polyallelic_significant'] / mp_db_info.clade_to_n_markers[sgb_id] \
                >= config['multistrain_min_markers_rel']:
            result_row['multi_strain'] = True
            info_debug(sgb_id, 'is multistrain')
        else:
            result_row['multi_strain'] = False
            info_debug(sgb_id, 'is not multistrain')

        for m in clade_markers_present:
            bfs_m = bfs[m]
            # replace with zeros minor alleles that are not in polyallelic significant positions
            with np.errstate(invalid='ignore'):  # 0/0 will be NaN which will produce False down the line
                afs = bfs_m / bfs_m.sum(axis=0)
            mask = (afs > 0.5) | polyallelic_significant_masks[m]  # major base or significant polyallelic position
            bfm_filtered = bfs_m.copy()
            bfm_filtered[~mask] = 0
            bfs_filtered[m] = bfm_filtered


        if not allele_counts_only:
            loci_rows = []
            for m in clade_markers_present:
                for pos in range(pr.marker_to_length[m]):
                    bfs_m = pr.base_frequencies[m][:, pos]
                    base_coverage = bfs_m.sum()
                    max_frequency = bfs_m.max()

                    if base_coverage < config['min_output_base_coverage']:
                        continue

                    pos_err_rate = err_rates[m][pos]
                    allelism = np.count_nonzero(bfs_m > 0)
                    polyallelic_significant = polyallelic_significant_masks[m][pos]
                    biallelic_significant = (allelism == 2) and polyallelic_significant
                    filtered = position_mask[m][pos]

                    loci_rows.append({
                        'marker': m,
                        'pos': pos + 1,
                        'error_rate': pos_err_rate,
                        'base_frequencies': bfs_m,
                        'base_coverage': base_coverage,
                        'max_frequency': max_frequency,
                        'allelism': allelism,
                        'filtered': filtered,
                        'polyallelic_significant': polyallelic_significant,
                        'biallelic_significant': biallelic_significant,
                    })

            df_loci_sgb = pd.DataFrame(loci_rows)
            output_dir_per_sgb = output_dir / 'per_sgb_data'
            if global_flags.debug:
                info_debug(sgb_id, 'Saving debug files')
                output_dir_per_sgb.mkdir(parents=True, exist_ok=True)
                df_loci_sgb.to_csv(output_dir_per_sgb / f'df_loci_{sgb_id}.tsv', sep='\t', index=False)

            df_loci_sgb_filtered = df_loci_sgb.query('filtered').copy()
            del df_loci_sgb

            if global_flags.debug:
                result_row['est_err_rate'] = 1 - df_loci_sgb_filtered['max_frequency'].sum() / \
                                             df_loci_sgb_filtered['base_coverage'].sum()


            if len(df_loci_sgb_filtered) == 0:
                warning(f'No loci left after filtering. Skipping SGB {sgb_id}, {output_dir.name}')
                continue

            if result_row['multi_strain']:
                info_debug(sgb_id, 'Creating and merging linkage')
                df_loci_biallelic_significant = df_loci_sgb_filtered.query('biallelic_significant')
                sgb_linkage = calculate_linkage(df_loci_biallelic_significant, sam_file,
                                                config)  # (m1, p1, m2, p2) => (b1 + b2) => count
                merging_results_before, merging_results = linkage_merging(df_loci_biallelic_significant, sgb_linkage,
                                                                          config)
                sgb_nps = merging_results[0]
                info_debug(sgb_id, 'Fitting model')
                fit_model(sgb_id, df_loci_sgb_filtered, result_row, config)
            else:
                sgb_nps = []
                merging_results_before = ([], None)


            info_debug(sgb_id, 'Generating genotype by maximizing per-position probabilities')
            consensuses_major, consensuses_minor, qualities_major, qualities_minor, log_probas_switch = \
                compute_genotypes(df_loci_sgb_filtered, sgb_nps, result_row, marker_to_length_sgb)

            strain_resolved_markers_sgb = [{
                'marker': m,
                'sequence_maj': consensuses_major[m].decode(),
                'sequence_min': consensuses_minor[m].decode(),
                'log_p_maj': qualities_major[m],
                'log_p_min': qualities_minor[m],
                'polyallelic_significant': polyallelic_significant_masks[m],
                'log_probas_switch': log_probas_switch[m],
            } for m in consensuses_major.keys()]

            if global_flags.debug:
                info_debug('Saving strain resolved markers')
                with open(output_dir_per_sgb / f'strain_resolved_markers_{sgb_id}.pic', 'wb') as f:
                    pickle.dump(strain_resolved_markers_sgb, f)

            consensuses_maj_sgb, consensuses_min_sgb = get_strainphlan_markers(strain_resolved_markers_sgb, result_row,
                                                                               merging_results_before, avg_read_len,
                                                                               config)
            consensuses_maj.update(consensuses_maj_sgb)
            consensuses_min.update(consensuses_min_sgb)

    df_results = pd.DataFrame.from_dict(result_rows, orient='index')

    return df_results, bfs_filtered, consensuses_maj, consensuses_min


def save_reconstructed_markers(output_results, output_major, output_minor, db_name, df_results, consensuses_maj,
                               consensuses_min, bfs_filtered, allele_counts_file, allele_counts_only):
    """

    :param pathlib.Path output_results:
    :param pathlib.Path output_major:
    :param pathlib.Path output_minor:
    :param str db_name:
    :param pd.DataFrame df_results:
    :param dict consensuses_maj:
    :param dict consensuses_min:
    :param dict[str, np.ndarray] bfs_filtered:
    :param pathlib.Path allele_counts_file:
    :param bool allele_counts_only:
    :return:
    """


    # Reconstructed genotypes
    if not allele_counts_only:
        consensuses_maj = list(consensuses_maj.values())
        consensuses_min = list(consensuses_min.values())

        with bz2.open(output_minor, 'wt') as f:
            json.dump({'database_name': db_name, 'consensus_markers': consensuses_min}, f, indent=2)

        with bz2.open(output_major, 'wt') as f:
            json.dump({'database_name': db_name, 'consensus_markers': consensuses_maj}, f, indent=2)


    # Results dataframe
    df_results.to_csv(output_results, sep='\t')


    # Allele counts
    allele_counts_metadata = {
        'database_name': db_name
    }
    pileup_path_tmp = allele_counts_file.with_name(allele_counts_file.name + '.tmp')
    with zipfile.ZipFile(pileup_path_tmp, 'w', compression=zipfile.ZIP_DEFLATED) as f:
        with f.open('metadata.json', 'w') as fm:
            json.dump(allele_counts_metadata, io.TextIOWrapper(fm))

        for m, bfs in bfs_filtered.items():
            with f.open(m + '.npy', 'w') as fm:
                np.save(fm, bfs)

    pileup_path_tmp.rename(allele_counts_file)



def run(sample_path, output_dir, config, save_bam_file, reuse, output_suffix, allele_counts_only, mp_db_info):
    """

    :param pathlib.Path sample_path:
    :param pathlib.Path output_dir:
    :param dict config:
    :param bool save_bam_file:
    :param str reuse:
    :param str output_suffix:
    :param bool allele_counts_only:
    :param MetaphlanDBInfo mp_db_info:
    :return:
    """

    if sample_path.suffix == '.bz2':
        sample_name = sample_path.with_suffix('').stem
    else:
        sample_name = sample_path.stem

    output_dir = output_dir / sample_name

    # Main outputs
    bam_path = output_dir / f'{sample_name}.sorted.bam'
    bai_path = bam_path.with_suffix(bam_path.suffix + '.bai')
    read_lens_path = output_dir / f'{sample_name}.read_lengths.txt.gz'

    results_path = output_dir / f'{sample_name}.results.tsv'
    output_major_path = output_dir / f'{sample_name}.major{output_suffix}.json.bz2'
    output_minor_path = output_dir / f'{sample_name}.minor{output_suffix}.json.bz2'
    output_allele_counts_path = output_dir / f'{sample_name}.allele_counts.zip'

    # Debug outputs
    pileup_path_before = output_dir / 'full_pileup_nonfilt.tsv.gz'
    pileup_path_after = output_dir / 'full_pileup_filt.tsv.gz'
    null_err_rate_path_before = output_dir / 'null_err_rate_nonfilt.txt'
    null_err_rate_path_after = output_dir / 'null_err_rate.txt'
    seq_error_path_before = output_dir / 'seq_error_profile_nonfilt.csv'
    seq_error_path_after = output_dir / 'seq_error_profile.csv'


    info(f'Starting sample {sample_name}')
    output_dir.mkdir(parents=True, exist_ok=True)

    if global_flags.debug or save_bam_file:
        files_to_remove = []
    else:
        files_to_remove = [bam_path, bai_path, read_lens_path]
    dirs_to_remove = []


    bam_target_files = [bam_path, bai_path]
    pileup_target_files = [pileup_path_before, pileup_path_after, null_err_rate_path_before, null_err_rate_path_after] \
                          + bam_target_files

    if reuse in ['bam', 'pileup'] and all(map(pathlib.Path.is_file, bam_target_files)):
        info(f'Loading existing preprocessed BAM file for sample {sample_name}')
        sam_file = load_bam(bam_path)
        if read_lens_path.is_file():
            read_lens = load_read_lens(read_lens_path)
        else:
            info('Getting read lengths from BAM file')
            read_lens = get_read_lens_from_bam(sam_file)
    else:
        info(f'Preprocessing the intput SAM file for sample {sample_name}')
        sam_file, read_lens = step_bam(sample_path, bam_path, bai_path, output_dir, mp_db_info.db_name, config)
        save_bam(read_lens_path, read_lens)

    if sam_file is None:
        return False

    if reuse in ['pileup'] and all(map(pathlib.Path.is_file, pileup_target_files)):
        info(f'Loading existing pileup for sample {sample_name}')
        pr_after = load_pileup(pileup_path_after, null_err_rate_path_after)
    else:
        info(f'Calculating pileup for sample {sample_name}')
        pr_before, pr_after = run_pileup(sam_file, config)
        if global_flags.debug:
            info_debug('Saving pileup')
            save_pileup(pr_before, pileup_path_before, seq_error_path_before, null_err_rate_path_before)
            save_pileup(pr_after, pileup_path_after, seq_error_path_after, null_err_rate_path_after)


    info(f'Running marker reconstruction for sample {sample_name}')
    df_results, bfs_filtered, consensuses_maj, consensuses_min = \
        step_reconstructed_markers(output_dir, config, sam_file, pr_after, read_lens, mp_db_info, allele_counts_only)

    save_reconstructed_markers(results_path, output_major_path, output_minor_path, mp_db_info.db_name, df_results,
                               consensuses_maj, consensuses_min, bfs_filtered, output_allele_counts_path,
                               allele_counts_only)


    info(f'Cleaning intermediate files for sample {sample_name}')
    for f in files_to_remove:
        f.unlink(missing_ok=True)

    for d in dirs_to_remove:
        if d.exists():
            d.rmdir()



    info(f'Finished sample {sample_name}')
    return True
