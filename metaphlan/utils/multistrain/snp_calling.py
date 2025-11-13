from collections import Counter, defaultdict, OrderedDict
from dataclasses import dataclass

import numpy as np

DTYPE_ERROR_RATES = np.float64

DTPYE_BASE_FREQUENCY = np.uint16
import numpy.typing as npt
import pandas as pd
import scipy.stats as sps
import statsmodels.stats.multitest as smsm
from typing import Sequence, Any

from ...utils import info, info_debug, warning

ACTGactg = 'ACTGactg'
ACTG = 'ACTG'
ACTG_to_i = dict(zip(ACTG, range(len(ACTG))))

BAM_FUNMAP = 4
BAM_FSECONDARY = 256
BAM_FQCFAIL = 512
BAM_FDUP = 1024


@dataclass
class PileupResult:
    base_frequencies: dict
    avg_error_rates: dict
    marker_to_length: dict
    df_seq_errors: pd.DataFrame | None
    null_err_rate: float | None
    substitutions: dict[Counter]
    substitutions_per_position: dict[dict[Counter]]
    contexts:  dict[dict[Counter]]


read_types = ['R1', 'R2', 'UN']

rev_c = {}
for b1, b2 in ['CG', 'AT']:
    rev_c[b1] = b2
    rev_c[b2] = b1

for b1 in 'ACTG':
    for b2 in 'ACTG':
        rev_c[b1 + b2] = rev_c[b1] + rev_c[b2]



def pu_to_r(pu):
    if pu.alignment.is_paired:
        assert (not pu.alignment.is_read1) or (not pu.alignment.is_read2)
        # from metaphlan we can get pairs but we don't know which is R1 and which R2
        return 'R1' if pu.alignment.is_read1 else 'R2' if pu.alignment.is_read2 else 'UN'
    else:
        return 'UN'


def run_pileup_inner(sam_file, config, drop_read_positions=None):
    all_markers = sam_file.references
    marker_lengths = sam_file.lengths
    marker_to_length = dict(zip(all_markers, marker_lengths))


    base_frequencies = {}
    avg_error_rates = {}

    all_positions = {r: Counter() for r in read_types}
    minor_positions = {r: Counter() for r in read_types}
    per_position_eror_rates_sum = {r: Counter() for r in read_types}
    subs = {r: Counter() for r in read_types}
    subs_pos = {r: defaultdict(Counter) for r in read_types}
    contexts = {r: defaultdict(Counter) for r in read_types}
    err_rates_positions = 0

    flag_filter = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL
    if config['filter_duplicates']:
        flag_filter |= BAM_FDUP

    for base_pileup in sam_file.pileup(contig=None, stepper='samtools', min_base_quality=config['min_base_quality'],
                                       min_mapping_quality=config['min_mapping_quality'], ignore_orphans=False,
                                       flag_filter=flag_filter):
        marker = base_pileup.reference_name
        pos = base_pileup.reference_pos  # 0 indexed

        query_sequences = base_pileup.get_query_sequences()
        query_qualities = base_pileup.get_query_qualities()
        query_positions = base_pileup.get_query_positions()
        rs = [pu_to_r(pu) for pu in base_pileup.pileups]

        # positions from the start of the read (if base is lower-case, the read had been reverse-complemented and the
        #   position has to be counted from the end)
        bqppur = [(b.upper(), q, p if b.isupper() else pu.alignment.qlen - p - 1, pu, r)
                  for b, q, p, pu, r in zip(query_sequences, query_qualities, query_positions, base_pileup.pileups, rs)
                  if b in ACTGactg]

        if drop_read_positions is not None:
            bqppur = [(b, q, p, pu, r) for b, q, p, pu, r in bqppur if (r, p) not in drop_read_positions]

        if not bqppur:  # not covered position
            continue

        bases, qualities, positions, pups, rs = zip(*bqppur)
        base_coverage = len(bases)

        if base_coverage < config['min_base_coverage']:
            continue

        # base is uppercase for fwd, lowercase for reverse mapping
        base_frequency = Counter(bases)

        qualities = np.array(qualities)
        err_rates = 10 ** (-qualities / 10)
        avg_error_rate = np.mean(err_rates)

        if marker not in base_frequencies:
            le = marker_to_length[marker]
            base_frequencies[marker] = np.zeros((4, le), DTPYE_BASE_FREQUENCY)
            avg_error_rates[marker] = np.zeros(le, DTYPE_ERROR_RATES)

        for b in base_frequency.keys():
            base_frequencies[marker][ACTG_to_i[b], pos] = base_frequency[b]
        avg_error_rates[marker][pos] = avg_error_rate

        if base_coverage >= config['seq_error_stats_min_cov']:
            major_base = base_frequency.most_common(1)[0][0]
            err_rates_positions += 1
            for b, p, e, pu, r in zip(bases, positions, err_rates, pups, rs):
                all_positions[r][p] += 1
                per_position_eror_rates_sum[r][p] += e

                if b != major_base:
                    minor_positions[r][p] += 1

                    bb = major_base + b
                    if not pu.alignment.is_forward:
                        bb = rev_c[bb]

                    subs[r][bb] += 1
                    subs_pos[r][bb][p] += 1

                    qp = pu.query_position
                    if qp > 0 and qp + 1 < pu.alignment.query_length:
                        # ctx = pu.alignment.get_reference_sequence()[qp-1:qp+2].upper()
                        ctx = pu.alignment.query_alignment_sequence[qp - 1:qp + 2]
                        if not all(b_c in ACTG for b_c in ctx):
                            continue

                        if not pu.alignment.is_forward:
                            ctx = ''.join(rev_c[x] for x in reversed(ctx))

                        assert ctx[1] == bb[1]

                        contexts[r][bb][ctx] += 1


    minor_positions_tot = sum(minor_positions.values(), Counter())
    all_positions_tot = sum(all_positions.values(), Counter())


    if all_positions_tot.total() > 0:
        a = config['null_err_rate_pseudocount']
        null_err_rate = (minor_positions_tot.total() + a) / (all_positions_tot.total() + a)
        info_debug(f'Null err. rate: {null_err_rate}')

        if null_err_rate > config['max_null_err_rate']:
            if config['max_null_err_rate'] > 0:
                warning(f'Null error rate is higher than expected ({null_err_rate}), will not use it')
            null_err_rate = None

        dfs = []
        for r in read_types:
            s_minor_positions = pd.Series(minor_positions[r], name='minor_bases')
            s_all_positions = pd.Series(all_positions[r], name='all_bases')
            s_per_position_eror_rates_sum = pd.Series(per_position_eror_rates_sum[r], name='sum_error_rates')
            df_seq_errors_r = pd.concat([s_minor_positions, s_all_positions, s_per_position_eror_rates_sum], axis=1) \
                .rename_axis('read_position').sort_index()
            df_seq_errors_r['r'] = r
            if len(df_seq_errors_r) > 0:
                dfs.append(df_seq_errors_r)

        if len(dfs) > 0:
            df_seq_errors = pd.concat(dfs).reset_index().set_index(['r', 'read_position'])

            a = config['null_err_rate_excess_pseudocount']
            df_seq_errors['null_err_rate_excess'] = (df_seq_errors['minor_bases'] + a) / (df_seq_errors['sum_error_rates'] + a)
        else:
            warning(f'No positions to get seq. errors statistics')
            df_seq_errors = None
    else:
        info('No positions to establish null error rate, will use default')
        df_seq_errors = None
        null_err_rate = None



    return PileupResult(base_frequencies, avg_error_rates, marker_to_length, df_seq_errors, null_err_rate, subs,
                        subs_pos, contexts)


def run_pileup(sam_file, config):
    """
    Returns a dataframe with all loci (filtered for SGB markers and min base coverage), with calculated p-values
    :param sam_file: pysam opened file
    :param dict config:
    :return:
    """

    pr_before = run_pileup_inner(sam_file, config)

    pr_after = pr_before
    if config['discard_read_ends'] and pr_before.df_seq_errors is not None:
        df_se = pr_before.df_seq_errors.query(f'all_bases>{config["discard_read_ends_support"]}')
        if len(df_se) > 0:
            c = np.log(df_se['null_err_rate_excess'].dropna())
            q25 = c.quantile(0.25)
            q75 = c.quantile(0.75)
            iqr = q75 - q25
            mx = q75 + 1.5 * iqr
            m = c > mx
            drop_read_positions = c.index[m]

            if len(drop_read_positions) > 0:
                warning(f'Will not use the following positions of the reads, '
                        f'there might be adapter contamination '
                        f'or other issue with sequencing: {list(drop_read_positions)}')
                pr_after = run_pileup_inner(sam_file, config, drop_read_positions)
        else:
            info_debug('Not dropping any read positions, not enough support')

    return pr_before, pr_after


def mask_iqr(x, q25, q75):
    if q25 is not None and q75 is not None:
        iqr = q75 - q25
        return (x > q75 + 1.5 * iqr) | (x < q25 - 1.5 * iqr)
    else:
        return np.zeros(x.shape, dtype=bool)


def fit_nb(q50, avg_read_len_i, z, p_t):
    if q50 is None:
        return None, None

    q50s = np.arange(avg_read_len_i + 1) / avg_read_len_i * q50  # taper the median linearly at the ends of the markers

    # convert median to mean assuming poisson dist.
    a = 50
    b = 50 / 3 - 50 * q50s
    c = -1
    ks = (-b + np.sqrt(b ** 2 - 4 * a * c)) / 2 / a

    # get Negative binomial parameters, z is the variance:mean ratio
    rs = ks / (z - 1)
    ps = rs / (rs + ks)
    d_nb = sps.nbinom(rs, ps)

    min_cs, max_cs = d_nb.isf(1 - p_t), d_nb.isf(p_t) - 1

    return min_cs, max_cs


def mask_nb(bc, min_cs, max_cs):
    """

    :param npt.NDArray[int] bc:
    :param npt.NDArray[int] min_cs:
    :param npt.NDArray[int] max_cs:
    :return:
    """
    mask = np.ones(len(bc), dtype=bool)
    if min_cs is not None:
        crop = min(len(min_cs) - 1, len(bc))
        mask[:crop] = (bc[:crop] >= min_cs[:crop]) & (bc[:crop] <= max_cs[:crop])
        mask[-crop:] = (bc[-crop:] >= min_cs[:crop:][::-1]) & (bc[-crop:] <= max_cs[:crop][::-1])
        mask[crop:-crop] = (bc[crop:-crop] >= min_cs[-1]) & (bc[crop:-crop] <= max_cs[-1])

    return mask



def filter_loci_snp_call(bfs, err_rates, avg_read_len, config):
    """

    :param dict[str, npt.NDArray[int]] bfs:
    :param dict[str, npt.NDArray[float]] err_rates:
    :param float avg_read_len:
    :param dict config:
    :return:
    """
    result_row = {}

    base_coverages = {m: bf.sum(axis=0) for m, bf in bfs.items()}
    avg_read_len_i = int(avg_read_len)
    assert avg_read_len_i > 0

    base_coverages_fit = {m: bc[avg_read_len_i:-avg_read_len_i] for m, bc in base_coverages.items()
                          if len(bc) > 2 * avg_read_len_i}

    q25, q50, q75 = None, None, None
    if len(base_coverages_fit) > 0:
        base_coverage_fit = np.concatenate(list(base_coverages_fit.values()))
        base_coverage_fit = base_coverage_fit[base_coverage_fit > 0]
        if len(base_coverage_fit) > 0:
            q25 = np.quantile(base_coverage_fit, .25)
            q50 = np.median(base_coverage_fit)
            q75 = np.quantile(base_coverage_fit, .75)


    min_cs, max_cs = fit_nb(q50, avg_read_len_i, config['outlier_coverage_overdispersion'], config['outlier_coverage_p'])
    position_mask = {m: mask_nb(bc, min_cs, max_cs) for m, bc in base_coverages.items()}
    outlier_coverage = sum(np.count_nonzero(~x) for x in position_mask.values()) /\
                       sum(len(x) for x in position_mask.values())
    marker_to_outlier_coverage = {m: (~x).mean() for m, x in position_mask.items()}
    markers_outlier_covered = set((m for m, oc in marker_to_outlier_coverage.items()
                                   if oc > config['outlier_coverage_marker']))


    result_row['coverage_25'] = q25
    result_row['coverage_50'] = q50
    result_row['coverage_75'] = q75
    result_row['outlier_coverage'] = outlier_coverage
    result_row['outlier_covered_markers'] = len(markers_outlier_covered)


    if config['trim_marker_ends'] > 0:
        t = config['trim_marker_ends']
        for m, pm in position_mask.items():
            pm[:t] = False
            pm[-t:] = False

    marker_covered_positions = [np.count_nonzero(bc[position_mask[m]]) for m, bc in base_coverages.items()]

    result_row['n_markers_total'] = len(bfs)
    result_row['n_positions_total'] = sum((np.count_nonzero(bc) for bc in base_coverages.values()))
    result_row['n_positions_filtered'] = sum(marker_covered_positions)
    result_row['n_markers_filtered'] = np.count_nonzero(marker_covered_positions)

    allelism = {m: (bf > 0).sum(axis=0) for m, bf in bfs.items()}
    allelism_filtered = {m: a * position_mask[m] for m, a in allelism.items()}
    result_row['n_biallelic_total'] = sum((np.count_nonzero(a == 2) for a in allelism_filtered.values()))
    result_row['n_triallelic_total'] = sum((np.count_nonzero(a == 3) for a in allelism_filtered.values()))
    result_row['n_quadallelic_total'] = sum((np.count_nonzero(a == 4) for a in allelism_filtered.values()))


    # Polyallelism test with FDR
    m_to_pvals = OrderedDict()
    m_to_poly_mask = {}
    for m in base_coverages.keys():
        polyallelic_mask = allelism_filtered[m] > 1
        if polyallelic_mask.sum() > 0:
            bf_poly = bfs[m][:, polyallelic_mask]
            bc_poly = base_coverages[m][polyallelic_mask]
            er_poly = err_rates[m][polyallelic_mask]
            poly_max_f = bf_poly.max(axis=0)
            p_values = sps.binom.cdf(poly_max_f, bc_poly, 1 - er_poly)
            m_to_pvals[m] = p_values
        m_to_poly_mask[m] = polyallelic_mask

    m_to_qvals = {}
    if len(m_to_pvals) > 0:
        all_pvalues = np.concatenate(list(m_to_pvals.values()))
        _, all_qvalues, _, _ = smsm.multipletests(all_pvalues, method="fdr_bh")

        # split concatenated array to per-marker arrays
        spacer = 0
        for m, ps in m_to_pvals.items():
            m_to_qvals[m] = all_qvalues[spacer: spacer+len(ps)]
            spacer += len(ps)
        assert spacer == len(all_qvalues)

    polyallelic_significant_masks = {}
    for m in base_coverages.keys():
        polyallelic_mask = m_to_poly_mask[m]
        mask_sig_polyallelic_full = np.zeros(len(base_coverages[m]), dtype=bool)
        if m in m_to_qvals:
            q_vals = m_to_qvals[m]
            mask_sig_polyallelic = q_vals < config['polymorphism_q_alpha']
            mask_sig_polyallelic_full[polyallelic_mask] = mask_sig_polyallelic
        polyallelic_significant_masks[m] = mask_sig_polyallelic_full


    # Polyallelism test without FDR (just p-values)
    # polyallelic_significant_masks = {}
    # for m in base_coverages.keys():
    #     polyallelic_mask = allelism_filtered[m] > 1
    #     mask_sig_polyallelic_full = np.zeros(len(base_coverages[m]), dtype=bool)
    #     if polyallelic_mask.sum() > 0:
    #         bf_poly = bfs[m][:, polyallelic_mask]
    #         bc_poly = base_coverages[m][polyallelic_mask]
    #         er_poly = err_rates[m][polyallelic_mask]
    #         poly_max_f = bf_poly.max(axis=0)
    #         p_values = sps.binom.cdf(poly_max_f, bc_poly, 1 - er_poly)
    #         mask_sig_polyallelic = p_values < config['polymorphism_p_alpha']
    #         mask_sig_polyallelic_full[polyallelic_mask] = mask_sig_polyallelic
    #     polyallelic_significant_masks[m] = mask_sig_polyallelic_full

    result_row['n_polyallelic_significant'] = sum((psm.sum() for psm in polyallelic_significant_masks.values()))
    result_row['n_markers_polyallelic_significant'] = sum((psm.sum() > 0 for psm in polyallelic_significant_masks.values()))

    return result_row, position_mask, polyallelic_significant_masks
