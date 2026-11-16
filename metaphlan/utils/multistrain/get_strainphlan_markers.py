import numpy as np

from . import utils
from .linkage import eval_pair

ACTG = 'ACTG'


# ==============================================================================================
# StrainPhlAn markers of the major and minor strain: quality masking, statistics and the release decision
# ==============================================================================================
ACTG_B = np.frombuffer(b'ACTG', dtype=np.uint8)
GAP = ord('-')


def _phred_ints(log_p):
    """ord(c) - 33 for every character produced by utils.qualities_to_phred(log_p)."""
    q = np.array(log_p)
    q = np.minimum(q, 0)
    q = -10 * utils.log1mexp(q) / np.log(10)
    out = np.where(np.isnan(q), 0, np.minimum(q, 93))
    return np.trunc(out).astype(np.int64)


def get_strainphlan_markers(strain_resolved_markers_sgb, result_row, merging_results, avg_read_len, config):
    """Major and minor consensus markers to output for one SGB (vectorised over positions); fills result_row."""
    ACTGset = 'ACTG'
    t = config['trim_marker_ends']
    t = max(t, int(avg_read_len / 2))
    min_q = config['min_output_quality']

    def calc_breadth_arr(s):
        L = len(s)
        if t > 0:
            lo, hi = t, L - t
            if hi <= lo:
                return 0
            n_t = hi - lo
            gaps = int(np.count_nonzero(s[lo:hi] == GAP))
        else:
            n_t = L
            if n_t == 0:
                return 0
            gaps = int(np.count_nonzero(s == GAP))
        return 100 * (1 - gaps / n_t)

    snp_rate_n = 0
    snp_rate_d = 0
    polyallelic_masked = 0
    p_switch = 0
    avg_quality_min_n = 0
    avg_quality_maj_n = 0
    avg_quality_min_d = 0
    avg_quality_maj_d = 0
    p_switch_d = 0
    consensuses_sgb_maj = {}
    consensuses_sgb_min = {}
    breadths_min = []
    masked_is_np = False  # the count becomes numpy.int64 once it adds a numpy.False_ ("ps and ..."), as in MetaPhlAn's code
    for srm in strain_resolved_markers_sgb:
        q_maj = _phred_ints(srm['log_p_maj'])
        q_min = _phred_ints(srm['log_p_min'])
        smaj = np.frombuffer(srm['sequence_maj'].encode('ascii'), dtype=np.uint8).copy()
        smin = np.frombuffer(srm['sequence_min'].encode('ascii'), dtype=np.uint8).copy()
        smaj[~(q_maj >= min_q)] = GAP  # a base is kept if q >= min_output_quality (NaN gaps all)
        smin[~(q_min >= min_q)] = GAP
        maj_called = smaj != GAP
        min_called = smin != GAP
        ps = np.asarray(srm['polyallelic_significant'], dtype=bool)
        polyallelic_masked += int(np.count_nonzero(ps & ~maj_called))
        masked_is_np = masked_is_np or not ps.all()
        lmin = np.asarray(srm['log_p_min'])
        lmaj = np.asarray(srm['log_p_maj'])
        lsw = np.asarray(srm['log_probas_switch'])
        avg_quality_min_n += np.exp(np.logaddexp.reduce(lmin[min_called]))
        avg_quality_min_d += int(np.count_nonzero(min_called))
        avg_quality_maj_n += np.exp(np.logaddexp.reduce(lmaj[maj_called]))
        avg_quality_maj_d += int(np.count_nonzero(maj_called))
        p_switch += np.exp(np.logaddexp.reduce(lsw[maj_called]))
        p_switch_d += int(np.count_nonzero(maj_called))

        breadth_maj = calc_breadth_arr(smaj)
        breadth_min = calc_breadth_arr(smin)
        breadths_min.append(breadth_min)

        sequence_maj = smaj.tobytes().decode('ascii')
        sequence_min = smin.tobytes().decode('ascii')
        if breadth_maj >= config['min_output_breadth']:
            consensuses_sgb_maj[srm['marker']] = {'marker': srm['marker'], 'sequence': sequence_maj,
                                                   'breadth': breadth_maj}
        if breadth_min >= config['min_output_breadth']:
            consensuses_sgb_min[srm['marker']] = {'marker': srm['marker'], 'sequence': sequence_min,
                                                   'breadth': breadth_min}

        isb_maj = np.isin(smaj, ACTG_B)
        isb_min = np.isin(smin, ACTG_B)
        both = isb_maj & isb_min
        snp_rate_n += int(np.count_nonzero(smaj[both] != smin[both]))
        snp_rate_d += int(np.count_nonzero(both))

    snp_rate_reconstructed = snp_rate_n / snp_rate_d if snp_rate_d > 0 else 0
    avg_quality_min = avg_quality_min_n / avg_quality_min_d if avg_quality_min_d > 0 else None
    avg_quality_maj = avg_quality_maj_n / avg_quality_maj_d if avg_quality_maj_d > 0 else None

    if masked_is_np:
        polyallelic_masked = np.int64(polyallelic_masked)

    result_row.update({
        'snp_rate_reconstructed': snp_rate_reconstructed,
        'n_positions_polyallelic_masked': polyallelic_masked,
        'p_switch': p_switch,
        'p_switch_d': p_switch_d,
        'avg_quality_maj': avg_quality_maj,
        'avg_quality_min': avg_quality_min,
    })

    if result_row['multi_strain'] and result_row['r_fit_var'] is not None and result_row['snp_rate_var'] is not None:
        snp_rate_predicted = result_row['snp_rate']
        r_var = result_row['r_fit_var']
        snp_rate_var = result_row['snp_rate_var']

        pred_completeness = 1 - np.maximum(0, snp_rate_predicted - snp_rate_reconstructed) / snp_rate_predicted
        pred_recall = pred_completeness * np.median(breadths_min) / 100

        node_pairs, g = merging_results
        node_pairs_reconstructed = []
        for np_ in node_pairs:
            assert len(np_.marker_pos) == 1
            m, p = np_.marker_pos[0]
            if m not in consensuses_sgb_min.keys():
                continue
            b_major = consensuses_sgb_maj[m]['sequence'][p - 1]
            b_minor = consensuses_sgb_min[m]['sequence'][p - 1]
            if b_major not in ACTGset or b_minor not in ACTGset:
                continue
            node_pairs_reconstructed.append(np_)

        xs = []
        for np1, np2 in zip(node_pairs_reconstructed[:-1], node_pairs_reconstructed[1:]):
            r = eval_pair(g, np1, np2)
            w_along = sum(x > 0 for x in r['w_along'])
            w_across = sum(x > 0 for x in r['w_across'])
            xs.append((w_along, w_across))

        good = 0
        bad = 0
        for x in xs:
            if x == (1, 0) or x == (0, 0):
                continue
            if x == (2, 0):
                good += 1
            elif x[1] == 1:
                bad += 1
            elif x[1] == 2:
                bad += 2
                good -= 1
            else:
                assert False

        pred_precision_support = good + bad
        pred_precision = max(0, good / pred_precision_support if pred_precision_support > 0 else 0)

        if snp_rate_reconstructed >= config['min_output_snp_rate'] / 100 \
                and pred_recall >= config['min_recall'] / 100 \
                and pred_precision >= config['min_precision'] / 100 \
                and pred_precision_support >= config['min_precision_support'] \
                and (np.isnan(config['max_snp_rate_var']) or 0 <= snp_rate_var <= config['max_snp_rate_var']) \
                and (np.isnan(config['max_r_var']) or 0 <= r_var <= config['max_r_var']):
            output_minor_strain = True
            consensuses_min = consensuses_sgb_min
        else:
            output_minor_strain = False
            consensuses_min = {}

        result_row.update({
            'pred_recall': pred_recall,
            'pred_precision': pred_precision,
            'pred_precision_support': pred_precision_support,
            'output_minor_strain': output_minor_strain,
        })
    else:
        consensuses_min = {}

    return consensuses_sgb_maj, consensuses_min
