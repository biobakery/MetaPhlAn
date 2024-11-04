from collections import Counter

import numpy as np
import pandas as pd
import scipy.stats as sps
import statsmodels.stats.multitest as smsm


def run_pileup(sam_file, config):
    """
    Returns a dataframe with all loci (filtered for SGB markers and min base coverage), with calculated p-values
    :param sam_file: pysam opened file
    :param config:
    :return:
    """

    ACTGactg = 'ACTGactg'
    ACTG = 'ACTG'

    all_markers = sam_file.references
    marker_lengths = sam_file.lengths
    marker_to_length = dict(zip(all_markers, marker_lengths))


    base_frequencies = {}
    error_rates = {}

    for base_pileup in sam_file.pileup(contig=None, stepper='samtools', min_base_quality=config['min_base_quality'],
                                       min_mapping_quality=config['min_mapping_quality'], ignore_orphans=False):
        marker = base_pileup.reference_name

        pos = base_pileup.pos  # 0 indexed

        query_sequences = base_pileup.get_query_sequences()
        query_qualities = base_pileup.get_query_qualities()

        bq = [(b.upper(), q) for b, q in zip(query_sequences, query_qualities) if b in ACTGactg]

        if not bq:  # not covered position (?)
            continue

        bases, qualities = zip(*bq)
        base_coverage = len(bases)

        if base_coverage < config['min_base_coverage']:
            continue

        base_frequency = Counter(bases)

        qualities = np.array(qualities)
        error_rate = np.mean(10 ** (-qualities / 10))

        if marker not in base_frequencies:
            le = marker_to_length[marker]
            base_frequencies[marker] = {b: np.zeros(le, np.uint16) for b in ACTG}
            error_rates[marker] = np.zeros(le, np.float64)

        for b in base_frequency.keys():
            base_frequencies[marker][b][pos] = base_frequency[b]
        error_rates[marker][pos] = error_rate


    return base_frequencies, error_rates, marker_to_length


def filter_loci(df_loci_sgb, marker_to_length, config):
    df_loci_sgb['filtered'] = False
    df_loci_sgb['biallelic_significant'] = False

    df_loci_sgb['end_pos'] = np.minimum(df_loci_sgb['pos'], df_loci_sgb['marker'].map(marker_to_length) - df_loci_sgb['pos'] + 1)
    # df_loci_sgb['base_frequencies'] = df_loci_sgb.apply(lambda r: {b: r['base_frequency_' + b] for b in ACTG}, axis=1)
    df_loci_sgb['base_coverage'] = df_loci_sgb['base_frequencies'].map(lambda bf: bf.total())
    df_loci_sgb['max_frequency'] = df_loci_sgb['base_frequencies'].map(lambda bf: max(bf.values()))
    df_loci_sgb['allelism'] = df_loci_sgb['base_frequencies'].map(lambda bf: sum(x > 0 for x in bf.values()))

    result_row = {}
    data = {}

    # trim the ends of markers
    df_loci_sgb_filtered = df_loci_sgb[df_loci_sgb['end_pos'] >= config['trim_marker_ends']]

    # filter markers with low breath of coverage
    marker_num_positions = df_loci_sgb_filtered.groupby('marker').size()
    marker_lens = pd.Series(index=marker_num_positions.index,
                            data=marker_num_positions.index.map(marker_to_length))
    marker_to_breath_rel = marker_num_positions / (marker_lens - 2 * config['trim_marker_ends'])
    low_breadth_markers = marker_to_breath_rel.index[marker_to_breath_rel < config['min_marker_breadth']]
    result_row['n_markers_total'] = len(marker_num_positions)
    result_row['n_markers_low_breadth'] = len(low_breadth_markers)
    data['low_breadth_markers'] = low_breadth_markers

    df_loci_sgb_filtered = df_loci_sgb_filtered[~df_loci_sgb_filtered['marker'].isin(low_breadth_markers)]

    # filter sites with unusual coverage (outside sgb-wise 1.5 IQR)
    base_coverage = df_loci_sgb_filtered['base_coverage']
    q25 = base_coverage.quantile(.25)
    q75 = base_coverage.quantile(.75)
    iqr = q75 - q25
    mask = (base_coverage > q75 + 1.5 * iqr) | (base_coverage < q25 - 1.5 * iqr)
    result_row['coverage_25'] = q25
    result_row['coverage_50'] = base_coverage.median()
    result_row['coverage_75'] = q75
    result_row['unusual_coverage'] = mask.mean()
    marker_to_unusual_coverage = mask.groupby(df_loci_sgb_filtered['marker']).mean()
    markers_unusually_covered = marker_to_unusual_coverage.index[marker_to_unusual_coverage > config['unusual_coverage_marker']]
    result_row['unusually_covered_markers'] = len(markers_unusually_covered)
    data['unusually_covered_markers'] = markers_unusually_covered

    result_row['n_positions_total'] = len(df_loci_sgb)
    df_loci_sgb_filtered = df_loci_sgb_filtered[(~mask) & (~df_loci_sgb_filtered['marker'].isin(markers_unusually_covered))]
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
        result_row['n_biallelic_significant'] = len(df_loci_sgb_biallelic_significant)

        if len(df_loci_sgb_biallelic_significant) > 0:
            # Filter markers by the number of sites
            gb_marker = df_loci_sgb_biallelic_significant.groupby('marker')
            # Test whether the SNPs are not uniformly distributed across the marker
            ms = []
            ps = []
            for m, idx_m in gb_marker.indices.items():
                df_loci_marker = df_loci_sgb_biallelic_significant.iloc[idx_m]
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
            df_loci_sgb.loc[df_loci_sgb_biallelic_significant.index, 'biallelic_significant'] = True
            df_loci_sgb_polyallelic.loc[df_loci_sgb_biallelic_significant.index, 'biallelic_significant'] = True

    return result_row, data, df_loci_sgb, df_loci_sgb_polyallelic
