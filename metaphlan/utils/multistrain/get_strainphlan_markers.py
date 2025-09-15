import numpy as np

from . import utils
from .linkage import eval_pair

ACTG = 'ACTG'


def get_strainphlan_markers(strain_resolved_markers_sgb, result_row, merging_results, avg_read_len, config):
    """
    Assesses the quality of reconstructed markers and produces the final genotypes to output

    :param strain_resolved_markers_sgb:
    :param dict result_row:
    :param merging_results:
    :param float avg_read_len:
    :param config:
    :return:
    """

    def calc_breadth(s):
        t = config['trim_marker_ends']
        t = max(t, int(avg_read_len/2))
        s_trimmed = s[t:-t] if t > 0 else s
        return 100 * (1 - s_trimmed.count('-') / len(s_trimmed)) if len(s_trimmed) > 0 else 0

    snp_rate_n = 0
    snp_rate_d = 0
    consensuses_sgb_maj = {}
    consensuses_sgb_min = {}
    breadths_min = []
    for srm in strain_resolved_markers_sgb:
        phred_maj = utils.qualities_to_phred(srm['log_p_maj'])
        phred_min = utils.qualities_to_phred(srm['log_p_min'])

        q_maj = [ord(q) - 33 for q in phred_maj]
        q_min = [ord(q) - 33 for q in phred_min]

        sequence_maj = ''.join(b if q >= config['min_output_quality'] else '-' for b, q in zip(srm['sequence_maj'], q_maj))
        sequence_min = ''.join(b if q >= config['min_output_quality'] else '-' for b, q in zip(srm['sequence_min'], q_min))
        breadth_maj = calc_breadth(sequence_maj)
        breadth_min = calc_breadth(sequence_min)
        breadths_min.append(breadth_min)

        if breadth_maj >= config['min_output_breadth']:
            consensuses_sgb_maj[srm['marker']] = {
                'marker': srm['marker'],
                'sequence': sequence_maj,
                'breadth': breadth_maj,
            }

        if breadth_min >= config['min_output_breadth']:
            consensuses_sgb_min[srm['marker']] = {
                'marker': srm['marker'],
                'sequence': sequence_min,
                'breadth': breadth_min,
            }

        for a, b in zip(sequence_maj, sequence_min):
            if a not in ACTG or b not in ACTG:
                continue
            snp_rate_n += a != b
            snp_rate_d += 1

    snp_rate_reconstructed = snp_rate_n / snp_rate_d if snp_rate_d > 0 else 0

    result_row.update({
        'snp_rate_reconstructed': snp_rate_reconstructed,
    })

    if result_row['multi_strain'] and result_row['r_fit_var'] is not None and result_row['snp_rate_var'] is not None:
        snp_rate_predicted = result_row['snp_rate']
        r_var = result_row['r_fit_var']
        snp_rate_var = result_row['snp_rate_var']

        # recall
        pred_completeness = 1 - np.maximum(0, snp_rate_predicted - snp_rate_reconstructed) / snp_rate_predicted
        pred_recall = pred_completeness * np.median(breadths_min) / 100

        # precision
        node_pairs, g = merging_results
        node_pairs_reconstructed = []
        for np_ in node_pairs:
            assert len(np_.marker_pos) == 1  # this should be before merging
            m, p = np_.marker_pos[0]
            if m not in consensuses_sgb_min.keys():
                continue

            b_major = consensuses_sgb_maj[m]['sequence'][p - 1]
            b_minor = consensuses_sgb_min[m]['sequence'][p - 1]
            if b_major not in ACTG or b_minor not in ACTG:  # this node pair wasn't actually reconstructed
                continue

            node_pairs_reconstructed.append(np_)

        xs = []
        for i, (np1, np2) in enumerate(zip(node_pairs_reconstructed[:-1], node_pairs_reconstructed[1:])):
            r = eval_pair(g, np1, np2)

            w_along = sum(x > 0 for x in r['w_along'])
            w_across = sum(x > 0 for x in r['w_across'])
            x = (w_along, w_across)
            xs.append(x)

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

        # filter
        if snp_rate_reconstructed >= config['min_output_snp_rate'] / 100 \
                and pred_recall >= config['min_recall'] / 100 \
                and pred_precision >= config['min_precision'] / 100 \
                and pred_precision_support >= config['min_precision_support'] \
                and (np.isnan(config['max_snp_rate_var']) or 0 <= snp_rate_var <= config['max_snp_rate_var']) \
                and (np.isnan(config['max_r_var']) or 0 <= r_var <= config['max_r_var']):
            consensuses_min = consensuses_sgb_min
        else:
            consensuses_min = {}

        result_row.update({
            'pred_recall': pred_recall,
            'pred_precision': pred_precision,
            'pred_precision_support': pred_precision_support,
        })
    else:
        consensuses_min = {}


    return consensuses_sgb_maj, consensuses_min
