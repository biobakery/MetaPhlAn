import itertools as it

import numpy as np
import pandas as pd
import scipy.stats as sps

from typing import Iterable
from .linkage import NodePair
from .utils import log1mexp

ACTG = 'ACTG'



def get_probas_1_mono_vect(freq, cov, eps):
    log_facts_mono = {}
    for b in ACTG:
        log_p_evidence = sps.multinomial.logpmf(freq, cov, p=np.array([1 - eps if x == b else eps / 3 for x in ACTG]).T)
        log_fact = log_p_evidence + np.log(1 / 4)  # use equal prior
        log_facts_mono[b] = log_fact

    log_margin_mono = np.logaddexp.reduce(list(log_facts_mono.values()))

    return log_margin_mono, log_facts_mono


def get_probas_1_bi_vect(freq, cov, r, eps):
    # model: bi-allelic with bases b_maj, b_min, ratio maj/min = r and noise eps
    # p(B_maj = b_maj, B_min = b_min | data, r, position is biallelic) =
    #   = p(data | b_maj, b_min, r) * p(b_maj, b_min | r) / margin
    log_facts_bi = {a: {} for a in ACTG}
    for b_maj, b_min in it.permutations(ACTG, 2):
        # without noise p[maj, min, other, other] = [r, 1-r, 0, 0]
        # apply one step of noise:
        #   p = [(1-e) * r, e/3 * r, e/3 * r, e/3 * r] + [e/3 * (1-r), (1-e) * (1-r), e/3 * (1-r), e/3 * (1-r)] =
        #     = [(1-e) * r + e/3 * (1-r), e/3 * r + (1-e) * (1-r), e/3, e/3]
        p_multinom = np.array([
            r * (1 - eps) + (1 - r) * eps / 3 if x == b_maj else
            (1 - r) * (1 - eps) + r * eps / 3 if x == b_min else
            eps / 3 for x in ACTG
        ]).T

        log_p_evidence = sps.multinomial.logpmf(freq, cov, p_multinom)
        log_fact = log_p_evidence + np.log(1 / 12)  # use equal prior
        log_facts_bi[b_maj][b_min] = log_fact

    # sum the facts in log-space
    log_margin_bi = np.logaddexp.reduce([y for x in log_facts_bi.values() for y in x.values()])

    return log_facts_bi, log_margin_bi


def get_probas_2_vect(log_margin_mono, log_facts_mono, log_facts_bi, log_margin_bi, freq, cov, eps, snp_rate, r):
    # H0 ... mono-allelic with eps noise
    # H1 ... bi-allelic with r
    log_p_biallelic_prior = np.log(snp_rate)
    log_p_not_biallelic_prior = np.log(1 - snp_rate)

    r_adj = r * (1 - eps) + (1 - r) * eps / 3
    log_p_not_biallelic_ev = np.logaddexp.reduce([sps.binom.logpmf(cov_i, cov, 1 - eps) for cov_i in freq.T])
    log_p_biallelic_ev = np.logaddexp.reduce([sps.binom.logpmf(cov_i, cov, r_adj) for cov_i in freq.T])

    log_p_not_biallelic_fact = log_p_not_biallelic_ev + log_p_not_biallelic_prior
    log_p_biallelic_fact = log_p_biallelic_ev + log_p_biallelic_prior

    # The p_biallelic_fact/p_not_biallelic_fact would give us the posterior odds.
    # Assuming H0 and H1 are exclusive and their union complete, we can convert to posterior probability and use
    #   them as weights
    log_margin = np.logaddexp(log_p_not_biallelic_fact, log_p_biallelic_fact)
    log_p_not_biallelic_post = log_p_not_biallelic_fact - log_margin
    log_p_biallelic_post = log_p_biallelic_fact - log_margin

    # Then the probability of a base coming from major/minor strain is a weighted combination of the mono/bi-allelic
    #   posterior probabilities
    log_p_major_minor = {}
    for bmaj in ACTG:
        log_p_major_minor[bmaj] = {}
        for bmin in ACTG:
            if bmaj == bmin:
                logp = log_p_not_biallelic_post + log_facts_mono[bmaj] - log_margin_mono
            else:
                logp = log_p_biallelic_post + log_facts_bi[bmaj][bmin] - log_margin_bi
            log_p_major_minor[bmaj][bmin] = logp

    log_p_major = {b: np.logaddexp.reduce(list(log_p_major_minor[b].values())) for b in ACTG}
    log_p_minor = {b: np.logaddexp.reduce([v[b] for v in log_p_major_minor.values()]) for b in ACTG}

    return log_p_major_minor, log_p_major, log_p_minor


def get_probas_switch_vect(max_f, cov):
    p_bi = 0.01  # prior that a position is bi-allelic (prior SNP rate)

    # p(r) = snp_rate * uniform for r in [0.5..1), 1-snp_rate for r=1
    # r is logarithmically sampled for 100 points
    rs = np.log10(np.logspace(0, 0.5, 100)) + 0.5
    rs_prior = np.ones(len(rs)) * p_bi / (len(rs) - 1)
    rs_prior[-1] = 1 - p_bi

    max_f = np.expand_dims(max_f, 0)
    cov = np.expand_dims(cov, 0)
    rs = np.expand_dims(rs, 1)
    rs_prior = np.expand_dims(rs_prior, 1)

    log_p_rs = np.log(rs_prior)

    min_f = cov - max_f

    # probability of switch ... true minor base becomes observed major base
    # p(switch | data) = int_r [ p(switch | r, data) * p(r | data) ]
    #                  = int_r [ p(switch | r, data) * p(data | r) * p(r) / p(data) ]
    #                  = int_r [ Binom(min_f, cov, r) * Binom(max_f, cov, r) * p(r) / p(data) ]

    log_p1s = sps.binom.logpmf(max_f, cov, rs) + log_p_rs
    log_p2s = sps.binom.logpmf(min_f, cov, rs)

    z = np.logaddexp.reduce(log_p1s)

    log_p_switch = np.logaddexp.reduce(log_p1s - z + log_p2s)

    log_q = log1mexp(log_p_switch)

    return log_q


def get_probas_switch_multistrain_vect(max_f, cov, eps, r):
    r_adj = r * (1 - eps) + (1 - r) * eps / 3
    min_f = cov - max_f

    log_p1s = sps.binom.logpmf(max_f, cov, r_adj)
    log_p2s = sps.binom.logpmf(min_f, cov, r_adj)
    log_p_switch = log_p2s - np.logaddexp(log_p1s, log_p2s)

    return log_p_switch


# ==============================================================================================
# Genotypes of the major and minor strain at every position (vectorised over positions)
# ==============================================================================================
ACTG_B = np.frombuffer(b'ACTG', dtype=np.uint8)


def compute_genotypes(loci, node_pairs, result_row, marker_to_length):
    """Major/minor consensus, their per-position log probabilities and the switch probabilities, per marker.
    loci: pipeline.Loci of the SGB's filtered positions."""
    consensuses_major, consensuses_minor, qualities_major, qualities_minor, log_probas_switch = {}, {}, {}, {}, {}

    marker_col = loci.marker_col()
    # set built from the unique markers in row order (its iteration order is the marker order of the output)
    markers = set(pd.unique(marker_col))
    for consensuses in [consensuses_major, consensuses_minor]:
        consensuses.update({m: bytearray(b'-' * marker_to_length[m]) for m in markers})
    for qualities in [qualities_major, qualities_minor]:
        qualities.update({m: np.zeros(marker_to_length[m]) for m in markers})
    log_probas_switch.update({m: np.full(marker_to_length[m], -np.inf) for m in markers})

    n = len(loci)
    if n == 0:
        return consensuses_major, consensuses_minor, qualities_major, qualities_minor, log_probas_switch

    base_frequencies = loci.bf  # uint16, like np.array(list(df['base_frequencies'].map(list).values))
    covs = loci.cov
    epss = loci.err

    maj_code = np.zeros(n, dtype=np.int64)
    min_code = np.zeros(n, dtype=np.int64)
    maj_lp = np.zeros(n)
    min_lp = np.zeros(n)
    has_minor = False

    if result_row['multi_strain']:
        has_minor = True
        r = result_row['r_fit']
        snp_rate = result_row['snp_rate']
        log_margin_mono, log_facts_mono = get_probas_1_mono_vect(base_frequencies, covs, epss)
        log_facts_bi, log_margin_bi = get_probas_1_bi_vect(base_frequencies, covs, r, epss)
        log_p_maj_min, log_p_maj, log_p_min = get_probas_2_vect(log_margin_mono, log_facts_mono, log_facts_bi,
                                                                log_margin_bi, base_frequencies, covs, epss,
                                                                snp_rate, r)
        lsw = get_probas_switch_multistrain_vect(loci.maxf, covs, epss, r)

        # 16-state choice with max() semantics over the states in order (bmaj in ACTG, then bmin in ACTG): the first
        # entry is kept unless a later one is strictly greater, so ties go to the first maximum and NaN never wins
        # unless it is the first entry
        best = np.zeros(n, dtype=np.int64)
        best_val = np.array(np.broadcast_to(log_p_maj_min[ACTG[0]][ACTG[0]], (n,)), dtype=np.float64)
        for k in range(1, 16):
            v = np.broadcast_to(log_p_maj_min[ACTG[k // 4]][ACTG[k % 4]], (n,))
            gt = v > best_val
            best[gt] = k
            best_val[gt] = v[gt]
        maj_code = best // 4
        min_code = best % 4
        lpm = np.stack([np.broadcast_to(log_p_maj[b], (n,)) for b in ACTG], axis=1)
        lpn = np.stack([np.broadcast_to(log_p_min[b], (n,)) for b in ACTG], axis=1)
        maj_lp = lpm[np.arange(n), maj_code]
        min_lp = lpn[np.arange(n), min_code]

        # linked node pairs (singletons with phaser 'none', blocks with graph merging)
        if node_pairs:
            key_to_row = {(marker_col[i], int(loci.pos[i])): i for i in range(n)}
            for node_pair in node_pairs:
                marker_pos = node_pair.marker_pos
                gt_a = node_pair.genotype_maj
                gt_b = node_pair.genotype_min
                rows = [key_to_row[(m, p)] for m, p in marker_pos]
                log_p_ab = 0
                log_p_ba = 0
                for row, ba, bb in zip(rows, gt_a, gt_b):
                    f_ab = log_facts_bi[ba][bb][row]
                    f_ba = log_facts_bi[bb][ba][row]
                    log_margin = np.logaddexp(f_ab, f_ba)
                    log_p_ab += f_ab - log_margin
                    log_p_ba += f_ba - log_margin
                log_p_ab /= len(marker_pos)
                log_p_ba /= len(marker_pos)
                for row, ba, bb in zip(rows, gt_a, gt_b):
                    # max({ba: log_p_ab, bb: log_p_ba}.items()) for the major and
                    # max(reversed({ba: log_p_ba, bb: log_p_ab}.items())) for the minor; in both the first item is
                    # kept unless the second is strictly greater (ties and NaN keep major = ba, minor = bb)
                    if log_p_ba > log_p_ab:
                        maj_code[row], maj_lp[row] = ACTG.index(bb), log_p_ba
                        min_code[row], min_lp[row] = ACTG.index(ba), log_p_ba
                    else:
                        maj_code[row], maj_lp[row] = ACTG.index(ba), log_p_ab
                        min_code[row], min_lp[row] = ACTG.index(bb), log_p_ab
    else:
        # single strain: quality from unique (max_frequency, base_coverage) pairs
        pairs = np.stack([loci.maxf.astype(np.uint64), covs.astype(np.uint64)], axis=1)
        uniq, inv = np.unique(pairs, axis=0, return_inverse=True)
        lp_u = get_probas_switch_vect(uniq[:, 0], uniq[:, 1])
        maj_lp = np.asarray(lp_u)[inv.ravel()]
        maj_code = base_frequencies.argmax(axis=1)
        lsw = None

    # scatter into per-marker arrays
    order = np.argsort(loci.mk, kind='stable')
    mk_sorted = loci.mk[order]
    bounds = np.flatnonzero(np.diff(mk_sorted)) + 1
    for seg in np.split(order, bounds):
        if len(seg) == 0:
            continue
        m = loci.markers[loci.mk[seg[0]]]
        p = loci.pos[seg] - 1
        cm = np.frombuffer(consensuses_major[m], dtype=np.uint8).copy()
        cm[p] = ACTG_B[maj_code[seg]]
        consensuses_major[m] = bytearray(cm.tobytes())
        qualities_major[m][p] = maj_lp[seg]
        if has_minor:
            cn = np.frombuffer(consensuses_minor[m], dtype=np.uint8).copy()
            cn[p] = ACTG_B[min_code[seg]]
            consensuses_minor[m] = bytearray(cn.tobytes())
            qualities_minor[m][p] = min_lp[seg]
            log_probas_switch[m][p] = lsw[seg]

    return consensuses_major, consensuses_minor, qualities_major, qualities_minor, log_probas_switch
