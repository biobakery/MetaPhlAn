import itertools as it

import numpy as np
import pandas as pd
import scipy.stats as sps

from typing import Iterable
from .linkage import NodePair
from .utils import FrozenCounter, log1mexp

ACTG = 'ACTG'


def get_probas_1_mono(get_probas_1_mono_cache, freq, cov, eps):
    """
    Get the mono-allelic probability

    :param dict get_probas_1_mono_cache:
    :param FrozenCounter freq:
    :param int cov:
    :param float eps:
    :return:
    """

    def get_probas_1_mono_inner():
        # model: mono-allelic position with base b and noise eps
        # p(B = b | data, position is monoallelic) = p(data | b) * p(b) / margin
        log_facts_mono = {}
        for b in ACTG:
            d = sps.multinomial(cov, p=[1 - eps if x == b else eps / 3 for x in ACTG])
            log_p_evidence = d.logpmf([freq[x] for x in ACTG])
            log_fact = log_p_evidence + np.log(1 / 4)  # use equal prior
            log_facts_mono[b] = log_fact

        log_margin_mono = np.logaddexp.reduce(list(log_facts_mono.values()))

        return log_margin_mono, log_facts_mono

    key = (freq, cov, eps)
    if key not in get_probas_1_mono_cache:
        get_probas_1_mono_cache[key] = get_probas_1_mono_inner()
    return get_probas_1_mono_cache[key]


def get_probas_1_bi(get_probas_1_bi_cache, freq, cov, r, eps):
    """
    Get the bi-allelic probability

    :param dict get_probas_1_bi_cache:
    :param FrozenCoutner freq:
    :param int cov:
    :param float r:
    :param float eps:
    :return:
    """

    def get_probas_1_bi_inner():

        # model: bi-allelic with bases b_maj, b_min, ratio maj/min = r and noise eps
        # p(B_maj = b_maj, B_min = b_min | data, r, position is biallelic) =
        #   = p(data | b_maj, b_min, r) * p(b_maj, b_min | r) / margin
        log_facts_bi = {a: {b: -np.inf for b in ACTG} for a in ACTG}
        for b_maj, b_min in it.permutations(ACTG, 2):
            # without noise p[maj, min, other, other] = [r, 1-r, 0, 0]
            # apply one step of noise:
            #   p = [(1-e) * r, e/3 * r, e/3 * r, e/3 * r] + [e/3 * (1-r), (1-e) * (1-r), e/3 * (1-r), e/3 * (1-r)] =
            #     = [(1-e) * r + e/3 * (1-r), e/3 * r + (1-e) * (1-r), e/3, e/3]
            p_multinom = [
                r * (1 - eps) + (1 - r) * eps / 3 if x == b_maj else
                (1 - r) * (1 - eps) + r * eps / 3 if x == b_min else
                eps / 3 for x in ACTG
            ]

            d = sps.multinomial(cov, p=p_multinom)
            log_p_evidence = d.logpmf([freq[x] for x in ACTG])
            log_fact = log_p_evidence + np.log(1 / 12)  # use equal prior
            log_facts_bi[b_maj][b_min] = log_fact

        # sum the facts in log-space
        log_margin_bi = np.logaddexp.reduce([y for x in log_facts_bi.values() for y in x.values()])

        return log_facts_bi, log_margin_bi

    key = (freq, cov, r, eps)
    if key not in get_probas_1_bi_cache:
        get_probas_1_bi_cache[key] = get_probas_1_bi_inner()
    return get_probas_1_bi_cache[key]



def get_probas_2(get_probas_2_cache, log_margin_mono, log_facts_mono, log_facts_bi, log_margin_bi, freq, cov, eps,
                 snp_rate, r):
    """

    Combines the mono and bi allelic probabilities into a final posterior probability

    :param get_probas_2_cache:
    :param float log_margin_mono:
    :param dict[str, float] log_facts_mono:
    :param dict[str, dict[str, float] log_facts_bi:
    :param float log_margin_bi:
    :param FrozenCounter freq:
    :param int cov:
    :param float eps:
    :param float snp_rate:
    :param float r:
    :return:
    """

    def get_probas_2_inner():
        # H0 ... mono-allelic with eps noise
        # H1 ... bi-allelic with r
        log_p_biallelic_prior = np.log(snp_rate)
        log_p_not_biallelic_prior = np.log(1 - snp_rate)

        r_adj = r * (1 - eps) + (1 - r) * eps / 3
        log_p_not_biallelic_ev = np.logaddexp.reduce([sps.binom.logpmf(cov_i, cov, 1 - eps) for cov_i in freq.values()])
        log_p_biallelic_ev = np.logaddexp.reduce([sps.binom.logpmf(cov_i, cov, r_adj) for cov_i in freq.values()])

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

    key = (freq, cov, r, eps, snp_rate)
    if key not in get_probas_2_cache:
        get_probas_2_cache[key] = get_probas_2_inner()
    return get_probas_2_cache[key]


def get_probas_switch(get_probas_switch_cache, max_f, cov):
    def get_probas_switch_inner():
        p_bi = 0.01  # prior that a position is bi-allelic (prior SNP rate)

        rs = np.log10(np.logspace(0, 0.5, 100)) + 0.5
        rs_prior = np.ones(len(rs)) * p_bi / (len(rs) - 1)
        rs_prior[-1] = 1 - p_bi

        log_p_rs = np.log(rs_prior)

        min_f = cov - max_f

        log_p1s = sps.binom.logpmf(max_f, cov, rs) + log_p_rs
        log_p2s = sps.binom.logpmf(min_f, cov, rs)

        z = np.logaddexp.reduce(log_p1s)

        log_p_switch = np.logaddexp.reduce(log_p1s - z + log_p2s)

        log_q = log1mexp(log_p_switch)

        return log_q

    key = (max_f, cov)
    if key not in get_probas_switch_cache:
        get_probas_switch_cache[key] = get_probas_switch_inner()

    return get_probas_switch_cache[key]


def compute_genotypes(df_loci_sgb_filtered, node_pairs, result_row, marker_to_length):
    """
    Calculates per-position probabilities of each base for both major and minor in the multi-strain case or
        single-strain case. Then maximize those probabilities to generate genotypes

    :param pd.DataFrame df_loci_sgb_filtered:
    :param Iterable[NodePair] node_pairs:
    :param dict result_row:
    :param dict marker_to_length:
    :return:
    """


    consensuses_major = {}
    consensuses_minor = {}
    qualities_major = {}
    qualities_minor = {}

    markers = set(df_loci_sgb_filtered['marker'].unique())

    for consensuses in [consensuses_major, consensuses_minor]:
        consensuses.update({m: bytearray(b'-' * marker_to_length[m]) for m in markers})
    for qualities in [qualities_major, qualities_minor]:
        qualities.update({m: np.zeros(marker_to_length[m]) for m in markers})

    df_loci_sgb_filtered['base_frequencies'] = df_loci_sgb_filtered['base_frequencies'].map(FrozenCounter)

    get_probas_1_mono_precomputed = {}
    res = df_loci_sgb_filtered.apply(lambda row_: get_probas_1_mono(
        get_probas_1_mono_precomputed, row_['base_frequencies'], row_['base_coverage'], row_['error_rate']), axis=1)
    res_a = np.array(res.values.tolist())
    log_margin_mono, log_facts_mono = res_a[:, 0], res_a[:, 1]
    df_loci_sgb_filtered['log_facts_mono'] = log_facts_mono
    df_loci_sgb_filtered['log_margin_mono'] = log_margin_mono


    if result_row['multi_strain']:
        r = result_row['r_fit']
        snp_rate = result_row['snp_rate']

        get_probas_1_bi_precomputed = {}
        res = df_loci_sgb_filtered.apply(lambda row_: get_probas_1_bi(
            get_probas_1_bi_precomputed, row_['base_frequencies'], row_['base_coverage'], r, row_['error_rate']),
                                         axis=1)
        res_a = np.array(res.values.tolist())
        log_facts_bi, log_margin_bi = res_a[:, 0], res_a[:, 1]
        df_loci_sgb_filtered['log_facts_bi'] = log_facts_bi
        df_loci_sgb_filtered['log_margin_bi'] = log_margin_bi

        get_probas_2_precomputed = {}
        args = ['log_margin_mono', 'log_facts_mono', 'log_facts_bi', 'log_margin_bi', 'base_frequencies',
                'base_coverage', 'error_rate']
        res = df_loci_sgb_filtered.apply(lambda row_: get_probas_2(
            get_probas_2_precomputed, *(row_[arg] for arg in args), snp_rate, r), axis=1)
        res_a = np.array(res.values.tolist())
        log_p_maj_min, log_p_maj, log_p_min = res_a[:, 0], res_a[:, 1], res_a[:, 2]
        df_loci_sgb_filtered['log_p_maj_min'] = log_p_maj_min
        df_loci_sgb_filtered['log_p_maj'] = log_p_maj
        df_loci_sgb_filtered['log_p_min'] = log_p_min

        df_loci_idx = df_loci_sgb_filtered.set_index(['marker', 'pos']).sort_index()
        df_loci_idx['log_p_maj_linked'] = pd.NA
        df_loci_idx['log_p_min_linked'] = pd.NA

        # link the probabilities across genotypes

        for node_pair in node_pairs:
            marker_pos = node_pair.marker_pos
            gt_a = node_pair.genotype_maj
            gt_b = node_pair.genotype_min

            log_p_ab = 0
            log_p_ba = 0
            for mp, ba, bb in zip(marker_pos, gt_a, gt_b):
                log_facts_bi = df_loci_idx.loc[mp, 'log_facts_bi']
                # use odds (coditioned on those two options) => renormalize to 1
                log_margin_bi = np.logaddexp(log_facts_bi[ba][bb], log_facts_bi[bb][ba])
                log_p_ab += log_facts_bi[ba][bb] - log_margin_bi
                log_p_ba += log_facts_bi[bb][ba] - log_margin_bi

            # geometric mean of probabilities to make sure all the positions are assigned together
            log_p_ab /= len(marker_pos)
            log_p_ba /= len(marker_pos)

            log_p_majs = []
            log_p_mins = []
            for mp, ba, bb in zip(marker_pos, gt_a, gt_b):
                log_p_maj = {ba: log_p_ab, bb: log_p_ba}
                log_p_min = {ba: log_p_ba, bb: log_p_ab}
                log_p_majs.append(log_p_maj)
                log_p_mins.append(log_p_min)

            df_loci_idx.loc[marker_pos, 'log_p_maj_linked'] = log_p_majs
            df_loci_idx.loc[marker_pos, 'log_p_min_linked'] = log_p_mins

        df_loci_sgb_filtered = df_loci_idx.reset_index()

        for _, row in df_loci_sgb_filtered.iterrows():
            pos = row['pos'] - 1
            marker = row['marker']

            if not pd.isna(row['log_p_maj_linked']):
                log_p_maj = row['log_p_maj_linked']
                log_p_min = row['log_p_min_linked']
                maj_b, maj_log_p = max(log_p_maj.items(), key=lambda x: x[1])
                # if 50:50, make sure to pick the other bas
                min_b, min_log_p = max(reversed(log_p_min.items()), key=lambda x: x[1])
            else:
                # Take the most probable combination of b_maj, b_min, then the quality will be marginalized probability
                log_p_maj_min = row['log_p_maj_min']
                maj_b, min_b, _ = max(((bmaj, bmin, logp) for bmaj, v in log_p_maj_min.items()
                                       for bmin, logp in v.items()), key=lambda x: x[2])
                maj_log_p = row['log_p_maj'][maj_b]
                min_log_p = row['log_p_min'][min_b]

            consensuses_major[marker][pos] = ord(maj_b)
            qualities_major[marker][pos] = maj_log_p
            consensuses_minor[marker][pos] = ord(min_b)
            qualities_minor[marker][pos] = min_log_p

    else:
        get_probas_switch_cache = {}
        df_loci_sgb_filtered['log_p'] = df_loci_sgb_filtered.apply(
            lambda row_: get_probas_switch(get_probas_switch_cache, row_['max_frequency'], row_['base_coverage']), axis=1)

        # df_loci_sgb_filtered['log_p'] = df_loci_sgb_filtered.apply(
        #     lambda row_: {b: row_['log_facts_mono'][b] - row_['log_margin_mono'] for b in ACTG}, axis=1)

        for _, row in df_loci_sgb_filtered.iterrows():
            pos = row['pos'] - 1
            marker = row['marker']
            # maj_b, maj_log_p = max(row['log_p'].items(), key=lambda x: x[1])
            maj_b = row['base_frequencies'].most_common()[0][0]
            maj_log_p = row['log_p']
            consensuses_major[marker][pos] = ord(maj_b)
            qualities_major[marker][pos] = maj_log_p

    return consensuses_major, consensuses_minor, qualities_major, qualities_minor
