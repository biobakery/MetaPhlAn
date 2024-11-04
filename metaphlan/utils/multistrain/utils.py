import argparse as ap
from collections import Counter, defaultdict
import os
import pathlib
import subprocess as sp
import shlex

import numpy as np
import scipy.special as spsp

from metaphlan.utils import error


def get_mem_usage():
    return int(sp.run(shlex.split('ps u -p %d' % os.getpid()), capture_output=True).stdout.decode()
               .split('\n')[1].split()[5]) // 1024


class FrozenCounter(Counter):
    # https://stackoverflow.com/questions/10045562/hashing-an-immutable-dictionary-in-python
    def __hash__(self):
        """Implements hash(self) -> int"""
        if not hasattr(self, '_hash'):
            self._hash = hash(frozenset(self.items()))
        return self._hash


class defaultdict_with_args(defaultdict):
    def __missing__(self, key):
        if self.default_factory:
            dict.__setitem__(self, key, self.default_factory(key))
            return self[key]
        else:
            defaultdict.__missing__(self, key)


class ArgumentType:
    @staticmethod
    def existing_file(path):
        path = pathlib.Path(path).resolve()
        if not os.path.isfile(path):
            raise ap.ArgumentTypeError('The file does not exist (%s).' % path)
        return path

    @classmethod
    def list_in_file(cls, path):
        if not path:
            return []

        path = cls.existing_file(path)
        with open(path) as f:
            r = [x.rstrip('\r\n') for x in f]
            r = [x for x in r if x]
            return r

    @staticmethod
    def existing_dir(path):
        path = pathlib.Path(path).resolve()
        if not path.is_dir():
            raise ap.ArgumentTypeError('The directory does not exist (%s).' % path)
        return path

    @staticmethod
    def creatable_dir(path):
        path = pathlib.Path(path).resolve()
        if path.is_dir() or path.parent.resolve().is_dir():
            return path
        raise ap.ArgumentTypeError('Neither the directory nor its parent exist (%s).' % path)

    @staticmethod
    def creatable_file(path):
        path = pathlib.Path(path).resolve()
        if path.is_file() or path.parent.resolve().is_dir():
            return path
        raise ap.ArgumentTypeError('Neither the file nor its parent exist (%s).' % path)

    @staticmethod
    def percentage(x):
        x = float(x)
        if x < 0 or x > 100:
            raise ap.ArgumentTypeError('The percentage must be between 0.0 and 100.0')

    @staticmethod
    def fraction(x):
        x = float(x)
        if x < 0 or x > 1:
            raise ap.ArgumentTypeError('The fraction must be between 0.0 and 1.0')

    @staticmethod
    def positive_int(x):
        x = int(x)
        if x <= 0:
            raise ap.ArgumentTypeError('The number must be greater than 0')
        return x


def cluster_get_queue_max_runs(queue):
    command = 'qmgr -c "print queue {}"'.format(queue)
    r = sp.run(command, shell=True, check=True, capture_output=True)
    lines = r.stdout.decode().split('\n')
    queue_attributes = {}
    for line in lines:
        if not line.startswith('set queue {}'.format(queue)):
            continue
        import re

        k, x, v = re.split(r' (\+?=) ', line[len('set queue {} '.format(queue)):])
        if x == '=':
            queue_attributes[k] = [v]
        elif x == '+=':
            queue_attributes[k].append(v)
        else:
            raise Exception('Failed to parse')
    if 'max_run' in queue_attributes:
        for max_run in queue_attributes['max_run']:
            if max_run.startswith('[u:PBS_GENERIC='):
                assert max_run.endswith(']')
                max_jobs = int(max_run[len('[u:PBS_GENERIC='):-1])
                return max_jobs
        raise Exception('Failed to parse')


def cluster_get_running_jobs(queue):
    command = 'qstat -u $USER {}'.format(queue)
    r = sp.run(command, shell=True, check=True, capture_output=True)
    stdout = r.stdout.decode()
    return max(0, len(stdout.split('\n')) - 6)


logit = spsp.logit   # np.log(x / (1-x))
inv_logit = spsp.expit  # exp_x / (exp_x + 1)   ... converts log-odds to probability


def logit_pdf(x, f_pdf):  # f(logit(x)) * d/dx logit(x) = f(logit(x)) * 1 / ( x*(1-x) )
    return f_pdf(logit(x)) / x / (1 - x)


def log1mexp(x):
    # https://stackoverflow.com/a/20888225
    with np.errstate(divide='ignore'):  # probabilities 1 will get division by zero warning and return -inf => supress
        return np.where(x < np.log(.5), np.log1p(-np.exp(x)), np.log(-np.expm1(x)))



def qualities_to_phred(q):
    q = np.array(q)
    q = np.minimum(q, 0)  # prevent numerical issues where log-probas go above zero
    q = -10 * log1mexp(q) / np.log(10)  # -10 * log10( 1 - exp(q) )
    return ''.join(chr(33 + int(min(qq, 93))) if not np.isnan(qq) else chr(33) for qq in q)


def run_bash(command):
    try:
        return sp.run(command, check=True, shell=True, capture_output=True)
    except sp.CalledProcessError as e:
        error(f'Failed to run a bash command {command}')
        error(f'stdout: {e.stdout.decode()}')
        error(f'stderr: {e.stderr.decode()}')
        exit(1)


def shorten_text(text, limit=500, over_limit_suffix='...'):
    if len(text) > limit:
        return text[:limit] + over_limit_suffix
    return text
