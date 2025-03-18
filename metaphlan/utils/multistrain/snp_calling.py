from collections import Counter, defaultdict
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import pandas as pd
import scipy.stats as sps
import statsmodels.stats.multitest as smsm
from typing import Sequence

from ...utils import info_debug, warning

ACTGactg = 'ACTGactg'
ACTG = 'ACTG'

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
            base_frequencies[marker] = {b: np.zeros(le, np.uint16) for b in ACTG}
            avg_error_rates[marker] = np.zeros(le, np.float64)

        for b in base_frequency.keys():
            base_frequencies[marker][b][pos] = base_frequency[b]
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
    a = config['null_err_rate_pseudocount']
    null_err_rate = (minor_positions_tot.total() + a) / (all_positions_tot.total() + a)
    info_debug(f'Null err. rate: {null_err_rate}')

    if null_err_rate > config['max_null_err_rate']:
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
        dfs.append(df_seq_errors_r)

    df_seq_errors = pd.concat(dfs).reset_index().set_index(['r', 'read_position'])

    a = config['null_err_rate_excess_pseudocount']
    df_seq_errors['null_err_rate_excess'] = (df_seq_errors['minor_bases'] + a) / (df_seq_errors['sum_error_rates'] + a)


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
    if config['discard_read_ends']:
        df_se = pr_before.df_seq_errors.query(f'all_bases>{config["discard_read_ends_support"]}')
        if len(df_se) > 0:
            c = df_se['null_err_rate_excess'].dropna()
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


def filter_loci_snp_call(df_loci_sgb, marker_to_length, read_lens, config):
    """

    :param pd.DataFrame df_loci_sgb:
    :param dict[str, int] marker_to_length:
    :param Sequence[int] read_lens:
    :param dict config:
    :return:
    """
    df_loci_sgb['filtered'] = False
    df_loci_sgb['polyallelic_significant'] = False
    df_loci_sgb['biallelic_significant'] = False

    df_loci_sgb['end_pos'] = np.minimum(df_loci_sgb['pos'],
                                        df_loci_sgb['marker'].map(marker_to_length) - df_loci_sgb['pos'] + 1)
    df_loci_sgb['base_coverage'] = df_loci_sgb['base_frequencies'].map(lambda bf: bf.total())
    df_loci_sgb['max_frequency'] = df_loci_sgb['base_frequencies'].map(lambda bf: max(bf.values()))
    df_loci_sgb['allelism'] = df_loci_sgb['base_frequencies'].map(lambda bf: sum(x > 0 for x in bf.values()))

    result_row = {}
    data = {}


    # filter sites with outlier coverage (outside sgb-wise 1.5 IQR)
    # to determine the IQR consider only positions at least avg-read-length away from the ends
    avg_read_len = np.mean(read_lens)
    base_coverage_fit = df_loci_sgb.query(f'end_pos >= {avg_read_len}')['base_coverage']
    base_coverage_all = df_loci_sgb['base_coverage']
    q25 = base_coverage_fit.quantile(.25)
    q75 = base_coverage_fit.quantile(.75)
    iqr = q75 - q25
    mask_fit = (base_coverage_fit > q75 + 1.5 * iqr) | (base_coverage_fit < q25 - 1.5 * iqr)
    mask_all = (base_coverage_all > q75 + 1.5 * iqr) | (base_coverage_all < q25 - 1.5 * iqr)
    result_row['coverage_25'] = q25
    result_row['coverage_50'] = base_coverage_fit.median()
    result_row['coverage_75'] = q75
    result_row['outlier_coverage'] = mask_fit.mean()
    marker_to_outlier_coverage = mask_fit.groupby(df_loci_sgb['marker']).mean()
    markers_outlier_covered = marker_to_outlier_coverage.index[marker_to_outlier_coverage > config['outlier_coverage_marker']]
    result_row['outlier_covered_markers'] = len(markers_outlier_covered)
    data['outlier_covered_markers'] = markers_outlier_covered

    df_loci_sgb_filtered = df_loci_sgb[(~mask_all) & (~df_loci_sgb['marker'].isin(markers_outlier_covered))]

    # trim the ends of markers
    df_loci_sgb_filtered = df_loci_sgb_filtered[df_loci_sgb_filtered['end_pos'] >= config['trim_marker_ends']]


    result_row['n_markers_total'] = len(df_loci_sgb['marker'].unique())
    result_row['n_positions_total'] = len(df_loci_sgb)
    result_row['n_positions_filtered'] = len(df_loci_sgb_filtered)
    result_row['n_markers_filtered'] = len(df_loci_sgb_filtered['marker'].unique())
    df_loci_sgb.loc[df_loci_sgb_filtered.index, 'filtered'] = True


    # poly-allelic
    df_loci_sgb_polyallelic = df_loci_sgb.query('filtered and allelism > 1').copy()
    result_row['n_biallelic_total'] = (df_loci_sgb_polyallelic['allelism'] == 2).sum()
    result_row['n_triallelic_total'] = (df_loci_sgb_polyallelic['allelism'] == 3).sum()
    result_row['n_quadallelic_total'] = (df_loci_sgb_polyallelic['allelism'] == 4).sum()

    if len(df_loci_sgb_polyallelic) > 0:
        # calculate p and q values of being polymorphic
        df_loci_sgb_polyallelic['p_value'] = df_loci_sgb_polyallelic.apply(
            lambda r: sps.binom.cdf(r['max_frequency'], r['base_coverage'], 1 - r['error_rate']), axis=1)
        _, q_values, _, _ = smsm.multipletests(df_loci_sgb_polyallelic['p_value'], method='fdr_bh')
        df_loci_sgb_polyallelic_significant = df_loci_sgb_polyallelic[q_values < config['polymorphism_q_alpha']]
        df_loci_sgb_biallelic_significant = df_loci_sgb_polyallelic_significant.query('allelism == 2')
        result_row['n_polyallelic_significant'] = len(df_loci_sgb_polyallelic_significant)
        result_row['n_biallelic_significant'] = len(df_loci_sgb_biallelic_significant)
        result_row['n_markers_polyallelic_significant'] = len(df_loci_sgb_polyallelic_significant['marker'].unique())


        if len(df_loci_sgb_polyallelic_significant) > 0:
            # Filter markers by the number of sites
            gb_marker = df_loci_sgb_polyallelic_significant.groupby('marker')
            # Test whether the SNPs are not uniformly distributed across the marker
            ms = []
            ps = []
            for m, idx_m in gb_marker.indices.items():
                df_loci_marker = df_loci_sgb_polyallelic_significant.iloc[idx_m]
                if len(df_loci_marker) < 2:
                    p_ks = 1.0
                else:
                    pos = df_loci_marker['pos']
                    d_uniform = sps.uniform(config['trim_marker_ends'], marker_to_length[m] - config['trim_marker_ends'])
                    _, p_ks = sps.kstest(pos, d_uniform.cdf)
                ms.append(m)
                ps.append(p_ks)

            _, qs, _, _ = smsm.multipletests(ps, method='fdr_bh')
            qs = pd.Series(index=ms, data=qs)
            non_uniformly_covered_markers = qs.index[qs < config['q_uniformity_alpha']]
            result_row['n_markers_non_uniform_coverage'] = len(non_uniformly_covered_markers)
            data['non_uniformly_covered_markers'] = non_uniformly_covered_markers

            df_loci_sgb.loc[df_loci_sgb_polyallelic_significant.index, 'polyallelic_significant'] = True
            df_loci_sgb_polyallelic.loc[df_loci_sgb_polyallelic_significant.index, 'polyallelic_significant'] = True

            if len(df_loci_sgb_biallelic_significant) > 0:
                df_loci_sgb.loc[df_loci_sgb_biallelic_significant.index, 'biallelic_significant'] = True
                df_loci_sgb_polyallelic.loc[df_loci_sgb_biallelic_significant.index, 'biallelic_significant'] = True

    return result_row, data, df_loci_sgb, df_loci_sgb_polyallelic

