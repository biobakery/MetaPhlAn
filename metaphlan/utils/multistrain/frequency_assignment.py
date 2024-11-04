from collections import Counter

import numpy as np
import scipy.optimize as spo
import scipy.stats as sps
from tqdm.auto import tqdm

from .utils import logit, inv_logit
from metaphlan.utils import warning


def binary_search(arr, x, low=None, high=None):
    if low is None:
        low = 0
    if high is None:
        high = len(arr) - 1

    if high >= low:
        mid = (high + low) // 2
        if np.isclose(arr[mid], x):
            return mid
        elif arr[mid] > x:
            return binary_search(arr, x, low, mid - 1)
        else:
            return binary_search(arr, x, mid + 1, high)
    else:
        return high


def cvm_exact(samples, cdf_f):
    """
    Calculates the Cramer-von Mises criterion (w^2 * 10,000)
    :param samples: Samples constituting the empirical distribution
    :param cdf_f: Cumulative distribution function
    :return:
    """
    n_steps = 1000
    domain_min = .5
    domain_max = 1.0
    # dx = (domain_max - domain_min) / n_steps
    domain = np.linspace(domain_min, domain_max, n_steps + 1)[:-1]  # dont' include 1
    cdf_1 = cdf_f(domain)

    pmf_2 = np.zeros(len(domain))
    i = 0
    for x in sorted(samples):
        while x > domain[i]:
            i += 1
        pmf_2[i] += 1
    pmf_2 /= len(samples)
    cdf_2 = np.cumsum(pmf_2)


    dy = np.diff(cdf_1)
    # cdf_1 = (cdf_1[1:] + cdf_1[:-1]) / 2
    # cdf_2 = (cdf_2[1:] + cdf_2[:-1]) / 2
    cdf_1 = cdf_1[1:]
    cdf_2 = cdf_2[1:]
    loss = (cdf_1 - cdf_2)**2

    return np.dot(loss, dy) * 100*100


def cvm_test(w_sq, n):
    T = w_sq * n / 100 / 100
    p = max(0, 1. - sps._hypotests._cdf_cvm(T, n))
    return p


class SumOfBinoms:
    domains_full = {}
    domains_folded = {}

    @staticmethod
    def get_domain_full(domain_bc):
        domain_full = set()
        for bc in domain_bc:
            if bc not in SumOfBinoms.domains_full:
                SumOfBinoms.domains_full[bc] = set()
                for i in range(1, bc):
                    x = i / bc
                    assert 0 < x < 1, "assertion failed for {} = {} / {}".format(x, i, bc)
                    SumOfBinoms.domains_full[bc].add(x)
            domain_full.update(SumOfBinoms.domains_full[bc])
        return np.array(sorted(domain_full))

    @staticmethod
    def get_domain_folded(domain_bc):
        domain_folded = set()
        for bc in domain_bc:
            if bc not in SumOfBinoms.domains_folded:
                SumOfBinoms.domains_folded[bc] = set()
                for i in range((bc + 1) // 2, bc):
                    x = i / bc
                    assert .5 <= x < 1, "assertion failed for {} = {} / {}".format(x, i, bc)
                    SumOfBinoms.domains_folded[bc].add(x)
            domain_folded.update(SumOfBinoms.domains_folded[bc])
        return np.array(sorted(domain_folded))

    def __init__(self):
        self.r = None
        self.domain_to_cdf = None
        self.cdf = None
        self.domain_full = None
        self.domain_folded = None
        self.base_coverage_pdf = None
        self.domain_bc = None


    def fit(self, base_coverage, r_obs=None):
        # generate empirical pdf of the coverage
        domain_bc = np.arange(0, base_coverage.max() + 1)
        self.domain_bc = domain_bc

        base_coverage_counter = Counter(base_coverage)
        base_coverage_pdf = np.array([base_coverage_counter[x] for x in domain_bc])
        base_coverage_pdf = base_coverage_pdf / base_coverage_pdf.sum()
        self.base_coverage_pdf = base_coverage_pdf

        # generate domain for the ratios
        self.domain_folded = SumOfBinoms.get_domain_folded(domain_bc)
        self.domain_full = SumOfBinoms.get_domain_full(domain_bc)

        # assert np.isclose(self.domain_full, 1 - self.domain_full[::-1]).all()

        if r_obs is not None:
            # optimize the log-likelihood to find the latent ratio
            def neg_ll_f(r_logit):
                return -self.ll_f(r_obs, r_logit)

            res = spo.minimize_scalar(neg_ll_f)
            if not res.success:
                warning(f'Unsuccessful fit of SOB distribution {res.message}')
                return False

            self.r = inv_logit(np.square(res.x))

            pmf = self.pmf_folded(self.domain_folded)
            cdf = np.cumsum(pmf)

            self.cdf = cdf
            self.domain_to_cdf = dict(zip(self.domain_folded, self.cdf))

        else:
            self.r = r_obs

        return True

    def _pmf_full_unnorm(self, x, r):
        """
        x ... observed ratio(s)
        r ... true ratio

        P(x) = sum_bc P(x | bc) * P(bc)
        P(bc) ... probability of the base covered bc times ... coming from empirical distribution or fitted neg. binom.
        if x * bc is not integer => result is 0 probability (only descrete values are allowed)
        P(x*bc | bc) ... binomial distribution at x*bc out of bc successes with success probability r
        x = 0 and x = 1 are removed (probability 0)

        this function returns *un-normalized* probability (because of the masking of 0 and 1)
        """
        x = np.atleast_1d(x)
        result_p = np.zeros(len(x))
        for bc in self.domain_bc:  # sum for all possible base coverages
            p_bc = self.base_coverage_pdf[bc]
            if p_bc == 0:
                continue
            xbc = x * bc
            xbc_round = np.round(xbc)
            mask = np.isclose(xbc, xbc_round)  # mask saying wheter x*bc is an integer
            p_bin = sps.binom.pmf(xbc_round.astype(int), n=bc, p=r)
            result_p += p_bc * p_bin * mask

        result_p = result_p * (x < 1) * (x > 0)  # mask x=0 and x=1 cases

        return result_p

    def pmf_full(self, x, r=None):
        """
        Normalize the probability so it sums to one over the domain
        """
        if r is None:
            r = self.r

        p_x = self._pmf_full_unnorm(x, r)
        denom = self._pmf_full_unnorm(self.domain_full, r).sum()
        return p_x / denom

    def pmf_folded(self, x, r=None):
        """
        Fold the distribution around 0.5
        """
        if r is None:
            r = self.r

        x = np.atleast_1d(x)
        p_right = self._pmf_full_unnorm(x, r) * (x >= 0.5)
        p_left = self._pmf_full_unnorm(1 - x, r) * (x > 0.5)
        denom = self._pmf_full_unnorm(self.domain_full, r).sum()
        return (p_right + p_left) / denom

    def cdf_folded(self, x):
        x = np.atleast_1d(x)
        result = []
        for xx in x:
            if xx in self.domain_to_cdf:
                result.append(self.domain_to_cdf[xx])
            else:
                i = binary_search(self.domain_folded, xx)
                result.append(self.cdf[i])
        return np.array(result)

    def cvm(self, sample):
        sample_c = Counter(sample)
        folded_domain_set = set(self.domain_folded)

        for x in list(sample_c.keys()):
            if x not in folded_domain_set:
                i = binary_search(self.domain_folded, x)
                sample_c[self.domain_folded[i]] = sample_c[x]

        # sample_pmf_denom = sum(sample_c.values())
        sample_pmf_denom = len(sample)

        pmf_dist = self.pmf_folded(self.domain_folded)

        r = 0
        cdf_dist = 0
        cdf_data = 0
        for i, d in enumerate(self.domain_folded):
            cdf_dist += pmf_dist[i]
            cdf_data += sample_c[d] / sample_pmf_denom
            z = (cdf_dist - cdf_data) ** 2
            dy = pmf_dist[i]
            r += z * dy

        return r * len(sample)

    def ll_f(self, r_obs, logit_r):
        """
        log-likelihood function of the data to be minimized
        logit_r can range over all real numbers
        """
        r = inv_logit(np.square(logit_r))
        probas = self.pmf_folded(r_obs, r)
        return np.log(probas).sum()


class FoldedT:
    def __init__(self, df, loc=0, scale=1):
        self.df = df
        self.args = [df, loc, scale]
        self.kwds = dict(df=df, loc=loc, scale=scale)

    def pdf(self, x):
        x = np.atleast_1d(x)
        p_right = sps.t.pdf(x, *self.args) * (x >= 0)
        p_left = sps.t.pdf(-x, *self.args) * (x > 0)
        return p_right + p_left

    def cdf(self, x):
        x = np.atleast_1d(x)

        d_r = sps.t(*self.args)
        args_rev = list(self.args)
        args_rev[1] *= -1  # negate the loc param to get the flipped left tail as the right tail
        d_l = sps.t(*args_rev)

        return (d_r.cdf(x) - d_r.cdf(0) + d_l.cdf(x) - d_l.cdf(0)) * (x >= 0)

    def logpdf(self, x):
        x = np.atleast_1d(x)
        return np.log(self.pdf(x))

    @classmethod
    def fit(cls, x, fit_df=True, x0=(1, 1), x0_df=2, tol=1e-3, floc=0):
        assert floc == 0, "Only folding around 0 is supported"

        def transform_params(params):
            params = np.array(list(params))
            params = np.square(params)
            if fit_df:
                params[0] += 1
            return params

        def inv_transform_params(params):
            params = np.array(list(params))
            if fit_df:
                params[0] -= 1
            params = np.sqrt(params)
            return params

        d = cls(df=x0_df)

        def neg_ll_f(params):
            params = transform_params(params)
            if fit_df:
                d.args = params
            else:
                d.args = [x0_df] + list(params)

            return -np.sum(d.logpdf(x))

        if fit_df:
            x0 = [x0_df] + list(x0)
        x0 = inv_transform_params(x0)

        res = spo.minimize(neg_ll_f, x0=x0, tol=tol)
        if not res.success:
            warning(f'Unsuccessful fit of t-distribution ({res.message})')

        res_x = transform_params(res.x)
        if fit_df:
            args = res_x
        else:
            args = [x0_df] + list(res_x)

        return np.array(args)


def fit_frequency_distributions(df_loci, datas, results, config):
    if config['frequency_distribution'] == "normal":
        folded_dist_cls = sps.foldnorm
        unfolded_dist_cls = sps.norm

        def fold_to_unfold(params):
            return [params[0] * params[2], params[2]]
    elif config['frequency_distribution'] == "cauchy":
        folded_dist_cls = sps.foldcauchy
        unfolded_dist_cls = sps.cauchy

        def fold_to_unfold(params):
            return [params[0] * params[2], params[2]]
    elif config['frequency_distribution'] == "t":
        folded_dist_cls = FoldedT
        unfolded_dist_cls = sps.t

        def fold_to_unfold(params):
            return params
    else:
        raise Exception(f'Unrecognized value "{config["frequency_distribution"]}" '
                        f'for config parameter "frequency_distribution"')

    df_loci_biallelic_significant = df_loci.query('biallelic_significant')
    gb_sgb = df_loci_biallelic_significant.groupby('sgb').groups
    for sgb in tqdm(gb_sgb.keys()):
        data = datas[sgb]
        result = results[sgb]

        df_loci_sgb = df_loci_biallelic_significant.loc[gb_sgb[sgb]]

        if len(df_loci_sgb) < config['min_significant_biallelic_loci']:
            continue

        # fit a folded-normal distribution to the logit of the ratio of the major base
        max_ratio = df_loci_sgb['max_ratio'].values
        max_ratio_logit = logit(max_ratio)

        folded_dist_params = folded_dist_cls.fit(max_ratio_logit, floc=0)
        folded_dist = folded_dist_cls(*folded_dist_params)

        unfolded_dist_params = fold_to_unfold(folded_dist_params)
        unfolded_dist = unfolded_dist_cls(*unfolded_dist_params)

        mu_fn = folded_dist.kwds['loc']
        data['foldnorm_dist'] = folded_dist
        data['norm_dist'] = unfolded_dist

        r_foldnorm = inv_logit(mu_fn)
        cvm_foldnorm = cvm_exact(max_ratio, lambda x: folded_dist.cdf(logit(x)))
        cvm_foldnorm_p = cvm_test(cvm_foldnorm, len(max_ratio))
        result['r_foldnorm'] = r_foldnorm
        result['cvm_foldnorm'] = cvm_foldnorm
        result['cvm_foldnorm_p'] = cvm_foldnorm_p

        # fit a "sum-of-binomials" distribution
        r_sob = np.nan
        cvm_sob = np.inf
        cvm_sob_p = np.nan
        base_coverage = df_loci_sgb['base_coverage']
        if base_coverage.mean() < config['sob_mean_base_coverage']:
            sob = SumOfBinoms()
            if sob.fit(base_coverage, max_ratio):
                data['sob_dist'] = sob
                r_sob = sob.r
                cvm_sob = cvm_exact(max_ratio, sob.cdf_folded)
                cvm_sob_p = cvm_test(cvm_sob, len(max_ratio))
        result['r_sob'] = r_sob
        result['cvm_sob'] = cvm_sob
        result['cvm_sob_p'] = cvm_sob_p

        if min(cvm_foldnorm, cvm_sob) > config['max_cvm']:
            continue

        if cvm_foldnorm <= cvm_sob:
            result['major_minor_ratio'] = r_foldnorm
        else:
            result['major_minor_ratio'] = r_sob
