from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import scipy.stats as sps
import statsmodels.stats.multitest as smsm


def reuse_error_rate(error_rate_info_file):
    with open(error_rate_info_file) as f:
        try:
            line = next(f)
        except StopIteration:
            line = ''
        if line.startswith('PCR_duplicates'):
            next(f)
            mismatches_1 = int(next(f).strip().split('\t')[1])
            matches_1 = int(next(f).strip().split('\t')[1])
            next(f)
            next(f)
            next(f)
            next(f)
            try:
                line = next(f)
            except StopIteration:
                line = ''
        else:
            mismatches_1 = 0
            matches_1 = 0

        if line.startswith('Pair-end mates overlap'):
            next(f)
            mismatches_2 = int(next(f).strip().split('\t')[1])
            matches_2 = int(next(f).strip().split('\t')[1])
            next(f)
            next(f)
            next(f)
        else:
            mismatches_2 = 0
            matches_2 = 0

    return matches_1, matches_2, mismatches_1, mismatches_2


def save_error_rate(matches_1, matches_2, mismatches_1, mismatches_2, qual_c, config, error_rate_info_file):
    n_positions_1 = matches_1 + mismatches_1
    with open(error_rate_info_file, 'w') as f:
        if n_positions_1 > 0:
            err_rate_1 = mismatches_1 / n_positions_1
            if err_rate_1 > 0:
                e_phred_1 = -10 * np.log10(err_rate_1)
            else:
                e_phred_1 = 'MAX'
            qual_s = '{' + ', '.join(f'{q}: {v / qual_c.total():.1%}' for q, v in qual_c.items()) + '}'
            f.write('PCR_duplicates:\n')
            f.write(f'evaluated_positions\t{n_positions_1}\n')
            f.write(f'mismatches\t{mismatches_1}\n')
            f.write(f'matches\t{matches_1}\n')
            f.write(f'error_rate\t{err_rate_1}\n')
            f.write(f'empirical_phred\t{e_phred_1}\n')
            f.write(f'base_qualities\t{qual_s}\n')
            f.write(f'base_quality_threshold\t{config["min_base_quality"]}\n')

        n_positions_2 = matches_2 + mismatches_2
        if n_positions_2 > 0:
            err_rate_2 = mismatches_2 / n_positions_2
            if err_rate_2 > 0:
                e_phred_2 = -10 * np.log10(err_rate_2)
            else:
                e_phred_2 = 'MAX'
            f.write('Pair-end mates overlap:\n')
            f.write(f'evaluated_positions\t{n_positions_2}\n')
            f.write(f'mismatches\t{mismatches_2}\n')
            f.write(f'matches\t{matches_2}\n')
            f.write(f'error_rate\t{err_rate_2}\n')
            f.write(f'empirical_phred\t{e_phred_2}\n')
            f.write(f'base_quality_threshold\t{config["min_base_quality"]}\n')


def estimate_error_rate_pcr_duplicates(sam_file, config):
    matches = 0
    mismatches = 0
    qual_c = Counter()

    last_pos = -1
    group_reads = []
    for read in sam_file.fetch():
        if not read.is_paired or not read.mate_is_mapped: # only pair-end with mapped mates
            continue

        if read.pos == last_pos:
            group_reads.append(read)
        else:
            if len(group_reads) > 1:
                gb = defaultdict(list)
                for r in group_reads:
                    gb[(r.rnext, r.pnext, r.rlen)].append(r)

                for k, reads in gb.items():
                    if len(reads) < 2:
                        continue

                    seqs = [[b if q >= config['min_base_quality'] else '-'
                             for b, q in zip(r.query_sequence, r.query_qualities)] for r in reads]
                    for r in reads:
                        for q in r.query_qualities:
                            qual_c[q] += 1

                    for b in zip(*seqs):
                        b = [x for x in b if x != '-']
                        if not len(b):
                            continue
                        if all(x == b[0] for x in b):
                            matches += len(b)
                        else:
                            m = Counter(b).most_common(1)[0][1]
                            mm = len(b) - m
                            matches += m
                            mismatches += mm

            group_reads = [read]
            last_pos = read.pos

    return matches, mismatches, qual_c


def estimate_error_rate_overlaps(sam_file, config):
    ACTG = 'ACTG'

    matches = 0
    mismatches = 0

    waiting_for_mate = {}
    for read in sam_file.fetch():
        if read.qname in waiting_for_mate:
            r1 = waiting_for_mate[read.qname]
            r2 = read

            del waiting_for_mate[read.qname]

            ref_pos_1 = r1.get_reference_positions(full_length=True)
            read_seq_1 = r1.query_alignment_sequence
            read_quals_1 = r1.query_qualities

            s1 = {}
            for p, b, q in zip(ref_pos_1, read_seq_1, read_quals_1):
                if q < config['min_base_quality']:
                    continue
                if b not in ACTG:
                    continue
                s1[p] = b

            ref_pos_2 = r2.get_reference_positions(full_length=True)
            read_seq_2 = r2.query_alignment_sequence
            read_quals_2 = r2.query_qualities

            for p, b2, q in zip(ref_pos_2, read_seq_2, read_quals_2):
                if q < config['min_base_quality']:
                    continue
                if p not in s1:
                    continue
                if b2 not in ACTG:
                    continue

                b1 = s1[p]

                if b1 == b2:
                    matches += 2
                else:
                    matches += 1
                    mismatches += 1

            continue

        if not read.is_paired or not read.mate_is_mapped or (read.mate_is_forward == read.is_forward):  # proper pair
            continue
        if read.reference_id != read.next_reference_id:  # same reference
            continue
        if read.pos > read.pnext:  # mate is after this one (process the pair only once)
            continue
        if read.pos + read.qlen < read.pnext:  # mate is not too far (overlap impossible)
            continue

        waiting_for_mate[read.qname] = read

    return matches, mismatches


def run_pileup(sam_file, config, error_rate):
    """
    Returns a dataframe with all loci (filtered for SGB markers and min base coverage), with calculated p-values
    :param sam_file: pysam opened file
    :param config:
    :param error_rate: Error rate, either float or "phred" to be calculated
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

        if error_rate == 'phred':
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
