import numpy as np
import numpy.typing as npt
import scipy.optimize as spo
import scipy.special as spsp
import scipy.stats as sps
from typing import Sequence

from ...utils import info_debug, warning


def transform_r(logit_r):
    r = 0.5 + spsp.expit(logit_r) / 2
    assert 0.5 <= r <= 1
    return r


def transform_r_rev(r):
    assert 0.5 <= r <= 1
    logit_r = spsp.logit(2 * r - 1)
    return logit_r


def expit_d(x):
    """d/dx expit(x)"""
    sx = spsp.expit(x)
    return sx * (1 - sx)


def expit_d2(x):
    """d2/dx2 expit(x)"""
    sx = spsp.expit(x)
    return sx * (1 - sx) * (1 - 2 * sx)


def binom_pmf_d(k, n, p, pmf):
    """
    d/dp pmf = c(n,k) * p**(k-1) * (1-p)**(n-k-1) * (k-np)
             = pmf * (k - np) / p / (1-p)
    """
    return pmf * (k - n * p) / p / (1 - p)


def binom_pmf_d2(k, n, p, pmf):
    """
    d2/dp2 pmf = c(n,k) * ((k-1)*k * p**(k-2) * (1-p)**(n-k) - 2*k*(n-k) * p**(k-1) * (1-p)**(n-k-1) +
    (n-k-1)*(n-k) * p**k * (1-p)**(n-k-2))
    """
    # return math.comb(n,k) * ((k-1)*k * p**(k-2) * (1-p)**(n-k) - 2*k*(n-k) * p**(k-1) * (1-p)**(n-k-1) +
    # (n-k-1)*(n-k) * p**k * (1-p)**(n-k-2))
    return pmf * ((k - 1) * k / p ** 2 - 2 * k * (n - k) / p / (1 - p) + (n - k - 1) * (n - k) / (1 - p) ** 2)


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

    bcs = df_vc['base_coverage'].values
    max_fs = df_vc['max_frequency'].values
    min_fs = bcs - max_fs
    eps = df_vc['error_rate'].values

    b1 = sps.binom.pmf(max_fs, bcs, 1 - eps)
    b2 = sps.binom.pmf(min_fs, bcs, 1 - eps)
    b = b1 + b2

    log_b1 = sps.binom.logpmf(max_fs, bcs, 1 - eps)
    log_b2 = sps.binom.logpmf(min_fs, bcs, 1 - eps)
    log_b = np.logaddexp(log_b1, log_b2)

    def get_r_adj(r):
        return r * (1 - eps) + (1 - r) * eps / 3

    def calc_logab(r):
        r_adj = get_r_adj(r)
        log_a1 = sps.binom.logpmf(max_fs, bcs, r_adj)
        log_a2 = sps.binom.logpmf(min_fs, bcs, r_adj)
        log_a = np.logaddexp(log_a1, log_a2)

        return log_a, log_b

    def calc_all(r):
        r_adj = get_r_adj(r)
        a1 = sps.binom.pmf(max_fs, bcs, r_adj)
        a2 = sps.binom.pmf(min_fs, bcs, r_adj)
        a = a1 + a2


        da1 = (1 - eps - eps / 3) * binom_pmf_d(max_fs, bcs, r_adj, a1)
        da2 = (1 - eps - eps / 3) * binom_pmf_d(min_fs, bcs, r_adj, a2)
        da = da1 + da2

        dda1 = (1 - 4 * eps / 3) ** 2 * binom_pmf_d2(max_fs, bcs, r_adj, a1)
        dda2 = (1 - 4 * eps / 3) ** 2 * binom_pmf_d2(min_fs, bcs, r_adj, a2)
        dda = dda1 + dda2

        return a, b, da, dda

    return calc_all, calc_logab


def neg_log_lik(x, counts, calc_all, calc_logab):
    logit_snp_rate, logit_r = x
    snp_rate = spsp.expit(logit_snp_rate)
    r = transform_r(logit_r)

    # L(x) = p(bi) * p(x | bi) + p(~bi) * p(x | ~bi)
    #    p(bi) = snp_rate
    #    p(~bi) = 1 - snp_rate
    #    p(x | bi) = a(r) ... bi-allelic with ratio r
    #    p(x | ~bi) = b   ... mono-allelic independent of r

    a, b, da, dda = calc_all(r)
    lik = a * snp_rate + b * (1 - snp_rate)
    nll = -np.log(lik)

    return np.dot(nll, counts)


def neg_log_lik_logspace(x, counts, calc_all, calc_logab):
    logit_snp_rate, logit_r = x
    snp_rate = spsp.expit(logit_snp_rate)
    r = transform_r(logit_r)

    loga, logb = calc_logab(r)
    nll = -np.logaddexp(loga + np.log(snp_rate), logb + np.log(1 - snp_rate))

    return np.dot(nll, counts)


def neg_log_lik_jac(x, counts, calc_all, calc_logab):
    logit_snp_rate, logit_r = x
    snp_rate = spsp.expit(logit_snp_rate)
    r = transform_r(logit_r)

    a, b, da, dda = calc_all(r)

    lik = a * snp_rate + b * (1 - snp_rate)

    # ls = logit(snp_rate)
    # lr = transform_inv(r) ... logit after mapping from [0.5, 1] => [0, 1]
    # L(x) = p(bi) * p(x | bi) + p(~bi) * p(x | ~bi)
    #      = snp_rate * a(r) + (1 - snp_rate) * b
    #      = expit(ls) * a(r) + (1 - expit(ls) * b     ... used for differentiate w.r.t ls
    #      = snp_rate * a(t(lr)) + (1 - snp_rate) * b  ... used for differentiate w.r.t lr
    # dL / d ls = expit_d(ls) * a(r) - expit_d(ls) * b
    # dL / d lr = snp_rate * d_a(r) * d_t(lr)


    lik_ds = (a - b) * expit_d(logit_snp_rate)
    lik_dr = da * snp_rate * expit_d(logit_r) / 2

    # d logL / dx = 1 / L * dL/dx
    nllik_ds = - lik_ds / lik
    nllik_dr = - lik_dr / lik

    nll_ds = np.dot(nllik_ds, counts)
    nll_dr = np.dot(nllik_dr, counts)

    return np.array([nll_ds, nll_dr])


def neg_log_lik_d2(x, calc_all, calc_logab):
    logit_snp_rate, logit_r = x
    snp_rate = spsp.expit(logit_snp_rate)
    r = transform_r(logit_r)

    a, b, da, dda = calc_all(r)

    lik = a * snp_rate + b * (1 - snp_rate)

    lik_ds = (a - b) * expit_d(logit_snp_rate)
    lik_dr = da * snp_rate * expit_d(logit_r) / 2

    lik_dss = (a - b) * expit_d2(logit_snp_rate)
    lik_drr = dda * snp_rate * (expit_d(logit_r) / 2) ** 2 + da * snp_rate * expit_d2(logit_r) / 2
    lik_dsr = da * expit_d(logit_snp_rate) * expit_d(logit_r) / 2

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


def neg_log_lik_d2_logspace(x, calc_all, calc_logab):
    logit_snp_rate, logit_r = x
    snp_rate = spsp.expit(logit_snp_rate)
    r = transform_r(logit_r)

    a, b, da, dda = calc_all(r)
    loga, logb = calc_logab(r)
    ll = np.logaddexp(loga + np.log(snp_rate), logb + np.log(1 - snp_rate))
    lik = np.exp(ll)
    a, b = np.exp(loga), np.exp(logb)

    lik_ds = (a - b) * expit_d(logit_snp_rate)
    lik_dr = da * snp_rate * expit_d(logit_r) / 2

    lik_dss = (a - b) * expit_d2(logit_snp_rate)
    lik_drr = dda * snp_rate * (expit_d(logit_r) / 2) ** 2 + da * snp_rate * expit_d2(logit_r) / 2
    lik_dsr = da * expit_d(logit_snp_rate) * expit_d(logit_r) / 2

    # d ll = lik' / lik ... see above (derivative of log)
    # d2 ll = (lik * lik'' - (lik')**2) / lik**2
    # (f * d/dxdy f - d/dx f * d/dy f) / f**2
    # nllik_dss = - (lik * lik_dss - lik_ds ** 2) / lik ** 2
    # nllik_drr = - (lik * lik_drr - lik_dr ** 2) / lik ** 2
    # nllik_drs = - (lik * lik_dsr - lik_dr * lik_ds) / lik ** 2


    nllik_dss = - (lik_dss / lik - (lik_ds / lik) ** 2)
    nllik_drr = - (lik_drr / lik - (lik_dr / lik) ** 2)
    nllik_drs = - (lik_dsr - lik_dr * lik_ds / lik) / lik

    # nllik_dss = - (div_in_logspace(lik_dss, ll) - div_in_logspace(lik_ds, ll) ** 2)
    # nllik_drr = - (div_in_logspace(lik_drr, ll) - div_in_logspace(lik_dr, ll) ** 2)
    # nllik_drs = - div_in_logspace(lik_dsr - div_in_logspace(lik_dr * lik_ds, ll), ll)


    return nllik_dss, nllik_drr, nllik_drs


def neg_log_lik_d2_logspace_untransformed(x, calc_all, calc_logab):
    logit_snp_rate, logit_r = x
    snp_rate = spsp.expit(logit_snp_rate)
    r = transform_r(logit_r)

    a, b, da, dda = calc_all(r)
    loga, logb = calc_logab(r)
    ll = np.logaddexp(loga + np.log(snp_rate), logb + np.log(1 - snp_rate))
    a, b = np.exp(loga), np.exp(logb)
    lik = np.exp(ll)

    lik_ds = a - b
    lik_dr = da * snp_rate

    lik_dss = 0
    lik_drr = dda * snp_rate
    lik_dsr = da

    # d ll = lik' / lik ... see above (derivative of log)
    # d2 ll = (lik * lik'' - (lik')**2) / lik**2
    # (f * d/dxdy f - d/dx f * d/dy f) / f**2
    # nllik_dss = - (lik * lik_dss - lik_ds ** 2) / lik ** 2
    # nllik_drr = - (lik * lik_drr - lik_dr ** 2) / lik ** 2
    # nllik_drs = - (lik * lik_dsr - lik_dr * lik_ds) / lik ** 2


    nllik_dss = - (lik_dss / lik - (lik_ds / lik) ** 2)
    nllik_drr = - (lik_drr / lik - (lik_dr / lik) ** 2)
    nllik_drs = - (lik_dsr - lik_dr * lik_ds / lik) / lik

    # nllik_dss = - (div_in_logspace(lik_dss, ll) - div_in_logspace(lik_ds, ll) ** 2)
    # nllik_drr = - (div_in_logspace(lik_drr, ll) - div_in_logspace(lik_dr, ll) ** 2)
    # nllik_drs = - div_in_logspace(lik_dsr - div_in_logspace(lik_dr * lik_ds, ll), ll)


    return nllik_dss, nllik_drr, nllik_drs


def neg_log_lik_hess(x, counts, calc_all, calc_logab):
    nllik_dss, nllik_drr, nllik_drs = neg_log_lik_d2(x, calc_all, calc_logab)

    i_s = np.dot(nllik_dss, counts)
    i_r = np.dot(nllik_drr, counts)
    i_sr = np.dot(nllik_drs, counts)

    return np.array([[i_s, i_sr], [i_sr, i_r]])


def neg_log_lik_hess_untransformed(x, counts, calc_all, calc_logab):
    nllik_dss, nllik_drr, nllik_drs = neg_log_lik_d2_logspace_untransformed(x, calc_all, calc_logab)

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


    if not res.success and method == 'BFGS' and \
            'Desired error not necessarily achieved due to precision loss' in res.message:
        if np.linalg.norm(res.jac) < config['acceptable_tol']:
            warning('SNP rate and strain ratio fitting failed was not precise, but acceptable')
            res.success = True

    if not res.success or np.isnan(res.x).any():
        warning(f'SNP rate and strain ratio fitting failed ({res.message}) after {res.nit} iterations')
        return None

    info_debug(f'Optimization took', res.nit, 'iterations')
    fisher_mat = neg_log_lik_hess(res.x, *args)
    try:
        fisher_mat_inv = np.linalg.inv(fisher_mat)
    except np.linalg.LinAlgError as e:
        if str(e) == 'Singular matrix':
            warning(f'Fisher information matrix is singular ({fisher_mat})')
            # fisher_mat_inv = [[np.nan, np.nan], [np.nan, np.nan]]
            return None
        else:
            raise e

    fisher_mat_ut = neg_log_lik_hess_untransformed(res.x, *args)
    try:
        fisher_mat_ut_inv = np.linalg.inv(fisher_mat_ut)
    except np.linalg.LinAlgError as e:
        if str(e) == 'Singular matrix':
            warning(f'Fisher information matrix is singular ({fisher_mat_ut})')
            fisher_mat_ut_inv = [[np.nan, np.nan], [np.nan, np.nan]]
        else:
            raise e

    res.snp_rate_var, res.r_var = fisher_mat_inv[0, 0], fisher_mat_inv[1, 1]
    res.snp_rate_var_ut, res.r_var_ut = fisher_mat_ut_inv[0, 0], fisher_mat_ut_inv[1, 1]
    if res.snp_rate_var < -config['negative_variance_tol'] or res.r_var < -config['negative_variance_tol']:
        warning(f'Fisher information matrix gave negative variance ({res.snp_rate_var}, {res.r_var})')
        return None

    res.snp_rate = spsp.expit(res.x[0])
    res.r = transform_r(res.x[1])

    return res


def fit_model(sgb_id, df_loci_sgb_filtered, result_row, config):
    """

    :param str sgb_id:
    :param pd.DataFrame df_loci_sgb_filtered:
    :param dict result_row:
    :param dict config:
    :return:
    """
    if len(df_loci_sgb_filtered) == 0:
        result_row['snp_rate'] = None
        result_row['r_fit'] = None
        return

    df_vc = df_loci_sgb_filtered[['max_frequency', 'base_coverage', 'error_rate']].value_counts()\
        .rename('count').to_frame().reset_index()
    counts = df_vc['count'].values

    calc_all, calc_logab = get_calc_fs(df_vc)

    # determine starting guess
    df_loci_polyallelic = df_loci_sgb_filtered.query('polyallelic_significant')
    assert len(df_loci_polyallelic) > 0
    r0 = min(config['init_r_max'], (df_loci_polyallelic['max_frequency'] / df_loci_polyallelic['base_coverage']).mean())
    r0 = max(r0, config['init_r_min'])
    snp_rate0 = max(config['init_snp_rate_min'], len(df_loci_polyallelic) / len(df_loci_sgb_filtered))
    info_debug(f'[{sgb_id}] Picking initial guess', r0, snp_rate0)

    x0 = (spsp.logit(snp_rate0), transform_r_rev(r0))
    args = (counts, calc_all, calc_logab)


    # Try three optimization methods in succession
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
        result_row['snp_rate'] = res.snp_rate
        result_row['r_fit'] = res.r
        result_row['snp_rate_var'] = res.snp_rate_var
        result_row['r_fit_var'] = res.r_var
        result_row['snp_rate_var_ut'] = res.snp_rate_var_ut
        result_row['r_fit_var_ut'] = res.r_var_ut
        result_row['fit_method'] = res.method
    else:
        warning(f'[{sgb_id}] SNP rate and strain ratio fitting failed, will use fallback values')
        result_row['snp_rate'] = config['snp_rate_fallback']
        result_row['r_fit'] = config['r_fallback']
        result_row['snp_rate_var'] = None
        result_row['r_fit_var'] = None
        result_row['snp_rate_var_ut'] = None
        result_row['r_fit_var_ut'] = None
        result_row['fit_method'] = 'fallback_values'
