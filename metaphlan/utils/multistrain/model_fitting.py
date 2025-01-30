import numpy as np
import numpy.typing as npt
import scipy.optimize as spo
import scipy.special as spsp
import scipy.stats as sps
from typing import Sequence

from .utils import defaultdict_with_args
from metaphlan.utils import warning
from ..util_fun import info_debug


def expit_d(x):
    """d/dx expit(x)"""
    sx = spsp.expit(x)
    return sx * (1 - sx)


def expit_d2(x):
    """d2/dx2 expit(x)"""
    sx = spsp.expit(x)
    return sx * (1 - sx) * (1 - 2*sx)


def binom_pmf_f(args):
    """pmf = c(n,k) * p**k * (1-p)**(n-k)"""
    max_frequency, base_coverage, r = args
    return sps.binom(base_coverage, r).pmf(max_frequency)


def binom_logpmf_f(args):
    """pmf = c(n,k) * p**k * (1-p)**(n-k)"""
    max_frequency, base_coverage, r = args
    return sps.binom(base_coverage, r).logpmf(max_frequency)


def get_calc_fs(df_vc):
    """
    a: probability of max-frequency k/n given the base is biallelic with strain ratio r and error rate eps
    a1 = Binom(n, k,   p=r*(1-eps) + (1-r)*eps/3)
    a2 = Binom(n, n-k, p=r*(1-eps) + (1-r)*eps/3)

    b: probability of max-frequency k/n given the base is mono allelic with error rate eps
    b1 = Binom(n, k,   p=1-eps)
    b2 = Binom(n, n-k, p=1-eps)

    a = a1 + a2, b = b1 + b2 ... folded Binomial distributions (folded because the true major base can be minor)

    da1, da2, dda1, dda2 are the first and second derivatives wrt. r of a1, a2
    """

    binom_pmf_cache = defaultdict_with_args(binom_pmf_f)
    binom_logpmf_cache = defaultdict_with_args(binom_logpmf_f)

    def binom_pmf_d_f(args):
        """
        d/dp pmf = c(n,k) * p**(k-1) * (1-p)**(n-k-1) * (k-np)
                 = pmf * (k - np) / p / (1-p)
        """
        k, n, p = args
        pmf = binom_pmf_cache[args]
        return pmf * (k - n * p) / p / (1 - p)


    def binom_pmf_d2_f(args):
        """d2/dp2 pmf = c(n,k) * ((k-1)*k * p**(k-2) * (1-p)**(n-k) - 2*k*(n-k) * p**(k-1) * (1-p)**(n-k-1) + (n-k-1)*(n-k) * p**k * (1-p)**(n-k-2))"""
        k, n, p = args
        pmf = binom_pmf_cache[args]
        # return math.comb(n,k) * ((k-1)*k * p**(k-2) * (1-p)**(n-k) - 2*k*(n-k) * p**(k-1) * (1-p)**(n-k-1) + (n-k-1)*(n-k) * p**k * (1-p)**(n-k-2))
        return pmf * ((k - 1) * k / p ** 2 - 2 * k * (n - k) / p / (1 - p) + (n - k - 1) * (n - k) / (1 - p) ** 2)


    def calc_ab(r):
        a1 = df_vc.apply(lambda row: binom_pmf_cache[(row['max_frequency'], row['base_coverage'], r * (1 - row['error_rate']) + (1 - r) * row['error_rate'] / 3)], axis=1)
        a2 = df_vc.apply(lambda row: binom_pmf_cache[(row['base_coverage'] - row['max_frequency'], row['base_coverage'], r * (1 - row['error_rate']) + (1 - r) * row['error_rate'] / 3)], axis=1)
        b1 = df_vc.apply(lambda row: binom_pmf_cache[(row['max_frequency'], row['base_coverage'], 1 - row['error_rate'])], axis=1)
        b2 = df_vc.apply(lambda row: binom_pmf_cache[(row['base_coverage'] - row['max_frequency'], row['base_coverage'], 1 - row['error_rate'])], axis=1)
        return a1 + a2, b1 + b2


    def calc_logab(r):
        a1 = df_vc.apply(lambda row: binom_logpmf_cache[(row['max_frequency'], row['base_coverage'], r * (1 - row['error_rate']) + (1 - r) * row['error_rate'] / 3)], axis=1)
        a2 = df_vc.apply(lambda row: binom_logpmf_cache[(row['base_coverage'] - row['max_frequency'], row['base_coverage'], r * (1 - row['error_rate']) + (1 - r) * row['error_rate'] / 3)], axis=1)
        b1 = df_vc.apply(lambda row: binom_logpmf_cache[(row['max_frequency'], row['base_coverage'], 1 - row['error_rate'])], axis=1)
        b2 = df_vc.apply(lambda row: binom_logpmf_cache[(row['base_coverage'] - row['max_frequency'], row['base_coverage'], 1 - row['error_rate'])], axis=1)
        return np.logaddexp(a1, a2), np.logaddexp(b1, b2)


    def calc_dab(r):
        da1 = df_vc.apply(lambda row: (1-row['error_rate'] - row['error_rate']/3)*binom_pmf_d_f((row['max_frequency'], row['base_coverage'], r * (1-row['error_rate']) + (1-r) * row['error_rate']/3)), axis=1)
        da2 = df_vc.apply(lambda row: (1-row['error_rate'] - row['error_rate']/3)*binom_pmf_d_f((row['base_coverage'] - row['max_frequency'], row['base_coverage'], r * (1-row['error_rate']) + (1-r) * row['error_rate']/3)), axis=1)
        return da1 + da2


    def calc_d2ab(r):
        dda1 = df_vc.apply(lambda row: (1-4*row['error_rate']/3)**2 * binom_pmf_d2_f((row['max_frequency'], row['base_coverage'], r * (1-row['error_rate']) + (1-r) * row['error_rate']/3)), axis=1)
        dda2 = df_vc.apply(lambda row: (1-4*row['error_rate']/3)**2 * binom_pmf_d2_f((row['base_coverage'] - row['max_frequency'], row['base_coverage'], r * (1-row['error_rate']) + (1-r) * row['error_rate']/3)), axis=1)
        return dda1 + dda2

    return calc_ab, calc_logab, calc_dab, calc_d2ab


def neg_log_lik(x, df_vc, ab_cache, logab_cache, dab_cache, d2ab_cache):
    logit_snp_rate, logit_r = x
    snp_rate = spsp.expit(logit_snp_rate)
    r = spsp.expit(logit_r)

    if r < 0.5:
        r = 1 - r

    a, b = ab_cache[r]
    lik = a * snp_rate + b * (1 - snp_rate)
    nll = -np.log(lik)

    return (nll * df_vc['count']).mean()


def neg_log_lik_logspace(x, df_vc, ab_cache, logab_cache, dab_cache, d2ab_cache):
    logit_snp_rate, logit_r = x
    snp_rate = spsp.expit(logit_snp_rate)
    r = spsp.expit(logit_r)

    if r < 0.5:
        r = 1 - r

    loga, logb = logab_cache[r]
    nll = -np.logaddexp(loga + np.log(snp_rate), logb + np.log(1 - snp_rate))

    return (nll * df_vc['count']).mean()


def neg_log_lik_jac(x, df_vc, ab_cache, logab_cache, dab_cache, d2ab_cache):
    logit_snp_rate, logit_r = x
    snp_rate = spsp.expit(logit_snp_rate)
    r = spsp.expit(logit_r)

    if r < 0.5:
        r = 1 - r

    a, b = ab_cache[r]
    da = dab_cache[r]

    lik = a * snp_rate + b * (1 - snp_rate)

    lik_ds = (a - b) * expit_d(logit_snp_rate)
    lik_dr = da * snp_rate * expit_d(logit_r)

    nllik_ds = - lik_ds / lik
    nllik_dr = - lik_dr / lik

    nll_ds = (nllik_ds * df_vc['count']).mean()
    nll_dr = (nllik_dr * df_vc['count']).mean()

    return np.array([nll_ds, nll_dr])


def neg_log_lik_d2(x, ab_cache, logab_cache, dab_cache, d2ab_cache):
    logit_snp_rate, logit_r = x
    snp_rate = spsp.expit(logit_snp_rate)
    r = spsp.expit(logit_r)

    if r < 0.5:
        r = 1 - r

    a, b = ab_cache[r]
    da = dab_cache[r]
    dda = d2ab_cache[r]

    lik = a * snp_rate + b * (1 - snp_rate)

    lik_ds = (a - b) * expit_d(logit_snp_rate)
    lik_dr = da * snp_rate * expit_d(logit_r)

    lik_dss = (a - b) * expit_d2(logit_snp_rate)
    lik_drr = dda * snp_rate * expit_d(logit_r) ** 2 + da * snp_rate * expit_d2(logit_r)
    lik_dsr = da * expit_d(logit_snp_rate) * expit_d(logit_r)

    # d2 ll = (lik * lik'' - (lik')**2) / lik**2
    # (f * d/dxdy f - d/dx f * d/dy f) / f**2
    # nllik_dss = - (lik * lik_dss - lik_ds ** 2) / lik ** 2
    # nllik_drr = - (lik * lik_drr - lik_dr ** 2) / lik ** 2
    # nllik_drs = - (lik * lik_dsr - lik_dr * lik_ds) / lik ** 2

    nllik_dss = - (lik_dss / lik - (lik_ds / lik) ** 2)
    nllik_drr = - (lik_drr / lik - (lik_dr / lik) ** 2)
    nllik_drs = - (lik_dsr - lik_dr * lik_ds / lik) / lik

    return nllik_dss, nllik_drr, nllik_drs


def div_in_logspace(a, log_b):
    return np.sign(a) * np.exp(np.log(np.abs(a)) - log_b)


def neg_log_lik_d2_logspace(x, ab_cache, logab_cache, dab_cache, d2ab_cache):
    logit_snp_rate, logit_r = x
    snp_rate = spsp.expit(logit_snp_rate)
    r = spsp.expit(logit_r)

    if r < 0.5:
        r = 1 - r

    loga, logb = logab_cache[r]
    ll = np.logaddexp(loga + np.log(snp_rate), logb + np.log(1 - snp_rate))
    a, b = np.exp(loga), np.exp(logb)
    da = dab_cache[r]
    dda = d2ab_cache[r]

    lik = np.exp(ll)

    lik_ds = (a - b) * expit_d(logit_snp_rate)
    lik_dr = da * snp_rate * expit_d(logit_r)

    lik_dss = (a - b) * expit_d2(logit_snp_rate)
    lik_drr = dda * snp_rate * expit_d(logit_r) ** 2 + da * snp_rate * expit_d2(logit_r)
    lik_dsr = da * expit_d(logit_snp_rate) * expit_d(logit_r)

    # d2 ll = (lik * lik'' - (lik')**2) / lik**2
    # (f * d/dxdy f - d/dx f * d/dy f) / f**2
    # nllik_dss = - (lik * lik_dss - lik_ds ** 2) / lik ** 2
    # nllik_drr = - (lik * lik_drr - lik_dr ** 2) / lik ** 2
    # nllik_drs = - (lik * lik_dsr - lik_dr * lik_ds) / lik ** 2

    # nllik_dss = - (lik_dss / lik - (lik_ds / lik) ** 2)
    # nllik_drr = - (lik_drr / lik - (lik_dr / lik) ** 2)
    # nllik_drs = - (lik_dsr - lik_dr * lik_ds / lik) / lik

    nllik_dss = - (div_in_logspace(lik_dss, ll) - div_in_logspace(lik_ds, ll) ** 2)
    nllik_drr = - (div_in_logspace(lik_drr, ll) - div_in_logspace(lik_dr, ll) ** 2)
    nllik_drs = - div_in_logspace(lik_dsr - div_in_logspace(lik_dr * lik_ds, ll), ll)

    return nllik_dss, nllik_drr, nllik_drs


def neg_log_lik_hess(x, df_vc, ab_cache, logab_cache, dab_cache, d2ab_cache):
    nllik_dss, nllik_drr, nllik_drs = neg_log_lik_d2(x, ab_cache, logab_cache, dab_cache, d2ab_cache)

    nllik_dss = (nllik_dss * df_vc['count']).mean()
    nllik_drr = (nllik_drr * df_vc['count']).mean()
    nllik_drs = (nllik_drs * df_vc['count']).mean()

    return np.array([[nllik_dss, nllik_drs], [nllik_drs, nllik_drr]])


def fisher_information(x, df_vc, ab_cache, logab_cache, dab_cache, d2ab_cache):
    """
    Observed fisher information matrix => asymptotically corresponds to the inverse of variance of the MLE estimate

    :param npt.NDArray x:
    :return:
    """
    nllik_dss, nllik_drr, nllik_drs = neg_log_lik_d2_logspace(x, ab_cache, logab_cache, dab_cache, d2ab_cache)

    counts = df_vc['count'].values

    i_s = np.dot(nllik_dss, counts)
    i_r = np.dot(nllik_drr, counts)
    i_sr = np.dot(nllik_drs, counts)

    return np.array([[i_s, i_sr], [i_sr, i_r]])


def run_optimization(x0, args, method, config):
    """

    :param npt.NDArray x0:
    :param Sequence args:
    :param str method:
    :param dict config:
    :return:
    """
    if method == 'Newton-CG':
        f = neg_log_lik
        jac = neg_log_lik_jac
        hess = neg_log_lik_hess
    elif method == 'BFGS':
        f = neg_log_lik
        jac = neg_log_lik_jac
        hess = None
    elif method == 'Nelder-Mead':
        f = neg_log_lik_logspace
        jac = None
        hess = None
    else:
        assert False

    res = spo.minimize(f, x0, jac=jac, hess=hess, args=args, method=method, options=dict(maxiter=config['max_iter']))
    res.method = method

    if not res.success or np.isnan(res.x).any():
        warning(f'SNP rate and strain ratio fitting failed ({res.message}) after {res.nit} iterations')
        return None

    info_debug(f'Optimization took', res.nit, 'iterations')
    fisher_mat = fisher_information(res.x, *args)
    try:
        fisher_mat_inv = np.linalg.inv(fisher_mat)
    except np.linalg.LinAlgError as e:
        if str(e) == 'Singular matrix':
            warning(f'Fisher information matrix is singular')
            return None
        else:
            raise e

    res.snp_rate_var, res.r_var = np.array([fisher_mat_inv[0, 0], fisher_mat_inv[1, 1]])
    if res.snp_rate_var < 0 or res.r_var < 0:
        warning(f'Fisher information matrix gave negative variance')
        return None

    res.snp_rate, res.r = spsp.expit(res.x)
    if res.r < 0.5:
        info_debug('Optimization of r converged below 0.5, flipping it')
        res.r = 1 - res.r

    return res


def fit_model(sgb_id, df_loci_sgb, result_row, config):
    """

    :param str sgb_id:
    :param pd.DataFrame df_loci_sgb:
    :param dict result_row:
    :param dict config:
    :return:
    """
    df_loci_sgb_filtered = df_loci_sgb.query('filtered')

    if len(df_loci_sgb_filtered) == 0:
        result_row['snp_rate'] = None
        result_row['r_fit'] = None
        return

    df_vc = df_loci_sgb_filtered[['max_frequency', 'base_coverage', 'error_rate']].value_counts()\
        .rename('count').to_frame().reset_index()

    calc_ab, calc_logab, calc_dab, calc_d2ab = get_calc_fs(df_vc)

    ab_cache = defaultdict_with_args(calc_ab)
    logab_cache = defaultdict_with_args(calc_logab)
    dab_cache = defaultdict_with_args(calc_dab)
    d2ab_cache = defaultdict_with_args(calc_d2ab)

    # determine starting guess
    df_loci_polyallelic = df_loci_sgb_filtered.query('polyallelic_significant')
    assert len(df_loci_polyallelic) > 0
    r0 = min(config['init_r_max'], (df_loci_polyallelic['max_frequency'] / df_loci_polyallelic['base_coverage']).mean())
    snp_rate0 = max(config['init_snp_rate_min'], len(df_loci_polyallelic) / len(df_loci_sgb_filtered))
    info_debug(f'[{sgb_id}] Picking initial guess', r0, snp_rate0)

    x0 = (spsp.logit(snp_rate0), spsp.logit(r0))
    args = (df_vc, ab_cache, logab_cache, dab_cache, d2ab_cache)


    res = run_optimization(x0, args, 'Newton-CG', config)
    if res is None:
        warning(f'[{sgb_id}] SNP rate and strain ratio fitting failed, '
                f'will restart with a 1st order optimization method')
        res = run_optimization(x0, args, 'BFGS', config)
    if res is None:
        warning(f'[{sgb_id}] SNP rate and strain ratio fitting failed, '
                f'will restart with a 0th order optimization method')
        res = run_optimization(x0, args, 'Nelder-Mead', config)
    if res is not None:
        snp_rate = res.snp_rate
        r = res.r
        snp_rate_var = res.snp_rate_var
        r_var = res.r_var
    else:
        warning(f'[{sgb_id}] SNP rate and strain ratio fitting failed, will use fallback values')
        snp_rate, r = None, None
        snp_rate_var, r_var = None, None

    result_row['snp_rate'] = snp_rate
    result_row['r_fit'] = r
    result_row['snp_rate_var'] = snp_rate_var
    result_row['r_fit_var'] = r_var
    result_row['fit_method'] = res.method
