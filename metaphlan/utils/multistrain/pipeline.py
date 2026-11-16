import bz2
import gzip
import io
import json
import pathlib
import zipfile
from collections import Counter
from typing import Sequence

import numpy as np
import numpy.typing as npt
import pandas as pd
import pysam

from .utils import MetaphlanDBInfo
from ...utils import info, error, warning, info_debug
from . import utils
from .genotype_resolving import compute_genotypes
from .get_strainphlan_markers import get_strainphlan_markers
from .linkage import calculate_linkage, linkage_merging
from .model_fitting import fit_model
from .prepare_sam import prepare_sam, get_read_lens_from_bam, NoAlignments
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
            return None  # no reads; the caller writes empty outputs
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

        if 'n_markers_significantly_polyallelic' in result_row \
                and result_row['n_markers_significantly_polyallelic'] >= config['multistrain_min_markers_abs'] \
                and result_row['n_markers_significantly_polyallelic'] / mp_db_info.clade_to_n_markers[sgb_id] \
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
            _reconstruct_sgb(sgb_id, clade_markers_present, pr, err_rates, position_mask, polyallelic_significant_masks,
                             marker_to_length_sgb, result_row, sam_file, avg_read_len, config, output_dir,
                             consensuses_maj, consensuses_min)

    df_results = pd.DataFrame.from_dict(result_rows, orient='index')

    return df_results, bfs_filtered, consensuses_maj, consensuses_min


# ==============================================================================================
# Per-position table of an SGB's markers, as numpy columns
# ==============================================================================================
class Loci:
    """Column store of the per-position rows of an SGB (marker, position, counts, coverage, error rate, flags)."""

    def __init__(self, markers, mk, pos, bf, cov, maxf, err, allelism, poly_sig, filtered):
        self.markers = markers  # marker names (row 'marker' = markers[mk])
        self.mk = mk
        self.pos = pos  # 1-based
        self.bf = bf  # (N, 4) uint16, ACTG
        self.cov = cov  # uint64 like bfs_m.sum()
        self.maxf = maxf  # uint16 like bfs_m.max()
        self.err = err
        self.allelism = allelism
        self.poly_sig = poly_sig
        self.bi_sig = (allelism == 2) & poly_sig
        self.filtered = filtered

    def subset(self, mask):
        return Loci(self.markers, self.mk[mask], self.pos[mask], self.bf[mask], self.cov[mask], self.maxf[mask],
                    self.err[mask], self.allelism[mask], self.poly_sig[mask], self.filtered[mask])

    def __len__(self):
        return len(self.pos)

    def marker_col(self):
        return np.array(self.markers, dtype=object)[self.mk]

    def fit_frame(self):
        """The columns fit_model reads."""
        return pd.DataFrame({
            'max_frequency': self.maxf,
            'base_coverage': self.cov,
            'error_rate': self.err,
            'polyallelic_significant': self.poly_sig,
        })

    def frame(self, with_bf=True):
        """Full DataFrame (used for the small biallelic-significant subset)."""
        d = {
            'marker': self.marker_col(),
            'pos': self.pos,
            'error_rate': self.err,
        }
        if with_bf:
            col = np.empty(len(self.pos), dtype=object)
            for i in range(len(self.pos)):
                col[i] = self.bf[i]
            d['base_frequencies'] = col
        d.update({
            'base_coverage': self.cov,
            'max_frequency': self.maxf,
            'allelism': self.allelism,
            'filtered': self.filtered,
            'polyallelic_significant': self.poly_sig,
            'biallelic_significant': self.bi_sig,
        })
        return pd.DataFrame(d)


def build_loci(clade_markers_present, base_frequencies, err_rates, position_mask, polyallelic_significant_masks,
               marker_to_length, min_output_base_coverage):
    mks, poss, bfs, errs, pms, fms = [], [], [], [], [], []
    for i, m in enumerate(clade_markers_present):
        bfm = base_frequencies[m]  # (4, L) uint16
        L = marker_to_length[m]
        cov = bfm.sum(axis=0, dtype=np.uint64)
        keep = ~(cov < min_output_base_coverage)  # NaN keeps all
        idx = np.nonzero(keep)[0]
        mks.append(np.full(len(idx), i, dtype=np.int64))
        poss.append(idx + 1)
        bfs.append(bfm[:, idx].T)
        errs.append(np.asarray(err_rates[m])[idx])
        pms.append(polyallelic_significant_masks[m][idx])
        fms.append(position_mask[m][idx])
    if not mks:
        z = np.zeros(0, dtype=np.int64)
        return Loci(list(clade_markers_present), z, z, np.zeros((0, 4), np.uint16), z.astype(np.uint64),
                    z.astype(np.uint16), z.astype(float), z, z.astype(bool), z.astype(bool))
    bf = np.ascontiguousarray(np.concatenate(bfs))
    cov = bf.sum(axis=1, dtype=np.uint64)
    maxf = bf.max(axis=1)
    allelism = (bf > 0).sum(axis=1).astype(np.int64)
    return Loci(list(clade_markers_present), np.concatenate(mks), np.concatenate(poss).astype(np.int64), bf, cov, maxf,
                np.concatenate(errs).astype(np.float64), allelism, np.concatenate(pms).astype(bool),
                np.concatenate(fms).astype(bool))

def _reconstruct_sgb(sgb_id, clade_markers_present, pr, err_rates, position_mask, polyallelic_significant_masks,
                     marker_to_length_sgb, result_row, sam_file, avg_read_len, config, output_dir,
                     consensuses_maj, consensuses_min):
    """Linkage, model fit, genotypes and output markers of one SGB."""
    loci_all = build_loci(clade_markers_present, pr.base_frequencies, err_rates, position_mask,
                          polyallelic_significant_masks, pr.marker_to_length, config['min_output_base_coverage'])
    loci = loci_all.subset(loci_all.filtered)
    del loci_all

    if len(loci) == 0:
        warning(f'No loci left after filtering. Skipping SGB {sgb_id}, {output_dir.name}')
        return

    if result_row['multi_strain']:
        info_debug(sgb_id, 'Creating and merging linkage')
        df_loci_biallelic_significant = loci.subset(loci.bi_sig).frame()
        sgb_linkage = calculate_linkage(df_loci_biallelic_significant, sam_file, config)
        merging_results_before, merging_results = linkage_merging(df_loci_biallelic_significant, sgb_linkage, config)
        sgb_nps = merging_results[0]
        info_debug(sgb_id, 'Fitting model')
        fit_model(sgb_id, loci.fit_frame(), result_row, config)
    else:
        sgb_nps = []
        merging_results_before = ([], None)

    info_debug(sgb_id, 'Generating genotype by maximizing per-position probabilities')
    consensuses_major, consensuses_minor, qualities_major, qualities_minor, log_probas_switch = \
        compute_genotypes(loci, sgb_nps, result_row, marker_to_length_sgb)

    strain_resolved_markers_sgb = [{
        'marker': m,
        'sequence_maj': consensuses_major[m].decode(),
        'sequence_min': consensuses_minor[m].decode(),
        'log_p_maj': qualities_major[m],
        'log_p_min': qualities_minor[m],
        'polyallelic_significant': polyallelic_significant_masks[m],
        'log_probas_switch': log_probas_switch[m],
    } for m in consensuses_major.keys()]

    consensuses_maj_sgb, consensuses_min_sgb = get_strainphlan_markers(strain_resolved_markers_sgb, result_row,
                                                                       merging_results_before, avg_read_len, config)
    result_row['n_markers_reconstructed_major'] = len(consensuses_maj_sgb)
    result_row['n_markers_reconstructed_minor'] = len(consensuses_min_sgb)
    consensuses_maj.update(consensuses_maj_sgb)
    consensuses_min.update(consensuses_min_sgb)


def save_reconstructed_markers(output_results, output_major, output_minor, db_name, df_results, consensuses_maj,
                               consensuses_min, bfs_filtered, allele_counts_file, allele_counts_only, config):
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
    :param dict config:
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
        'database_name': db_name,
        'base_order': ACTG,
        'min_mapping_quality': config['min_mapping_quality'],
        'min_read_pident': config['min_read_pident'],
        'min_base_coverage': config['min_base_coverage'],
        'min_base_quality': config['min_base_quality']
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
    :param str reuse: 'bam' to reuse the preprocessed BAM of an earlier run if present, 'none' to redo everything
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

    info(f'Starting sample {sample_name}')
    output_dir.mkdir(parents=True, exist_ok=True)

    if save_bam_file:
        files_to_remove = []
    else:
        files_to_remove = [bam_path, bai_path, read_lens_path]

    bam_target_files = [bam_path, bai_path]

    def empty_outputs(reason):
        warning(f'Sample {sample_name}: {reason}. Writing empty outputs (no markers in the major and minor JSONs)')
        save_reconstructed_markers(results_path, output_major_path, output_minor_path, mp_db_info.db_name,
                                   pd.DataFrame(), {}, {}, {}, output_allele_counts_path, allele_counts_only, config)
        for f in files_to_remove:
            f.unlink(missing_ok=True)
        info(f'Finished sample {sample_name}')
        return True

    if reuse == 'bam' and all(map(pathlib.Path.is_file, bam_target_files)):
        info(f'Loading existing preprocessed BAM file for sample {sample_name}')
        sam_file = load_bam(bam_path)
        if sam_file is None:
            return empty_outputs('the preprocessed BAM has no reads')
        if read_lens_path.is_file():
            read_lens = load_read_lens(read_lens_path)
        else:
            info('Getting read lengths from BAM file')
            read_lens = get_read_lens_from_bam(sam_file)
    else:
        info(f'Preprocessing the input file of sample {sample_name}')
        try:
            sam_file, read_lens = step_bam(sample_path, bam_path, bai_path, output_dir, mp_db_info.db_name, config)
        except NoAlignments as e:
            return empty_outputs(str(e))
        if sam_file is None:
            return empty_outputs('the preprocessed BAM has no reads')
        save_bam(read_lens_path, read_lens)

    info(f'Calculating pileup for sample {sample_name}')
    pr_before, pr_after = run_pileup(sam_file, config)


    info(f'Running marker reconstruction for sample {sample_name}')
    df_results, bfs_filtered, consensuses_maj, consensuses_min = \
        step_reconstructed_markers(output_dir, config, sam_file, pr_after, read_lens, mp_db_info, allele_counts_only)

    save_reconstructed_markers(results_path, output_major_path, output_minor_path, mp_db_info.db_name, df_results,
                               consensuses_maj, consensuses_min, bfs_filtered, output_allele_counts_path,
                               allele_counts_only, config)
    if not allele_counts_only and not consensuses_maj:
        warning(f'Sample {sample_name}: no marker was reconstructed ({len(df_results)} species processed); the major '
                f'and minor JSONs are empty')


    info(f'Cleaning intermediate files for sample {sample_name}')
    for f in files_to_remove:
        f.unlink(missing_ok=True)


    info(f'Finished sample {sample_name}')
    return True
