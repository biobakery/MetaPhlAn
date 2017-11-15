#!/usr/bin/env python

import sys
import os
import bz2
import gzip

import glob


p2 = float(sys.version_info[0]) < 3.0


def ignore_spaces(l):
    if (l[0] == '@') or (l[0] == '>'):
        return l.replace(' ', '_')

    return l


def fopen(fn):
    fileName, fileExtension = os.path.splitext(fn)

    if fileExtension == '.bz2':
        return (bz2.open(fn, "rt") if not p2 else bz2.BZ2File(fn, "r"))

    if fileExtension == '.gz':
        return gzip.open(fn, "rt")

    return open(fn)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        if sys.argv[1] in ['-h', '--help', '-v', '--version']:
            sys.exit(0)

    if len(sys.argv) < 2:
        for l in sys.stdin:
            sys.stdout.write(ignore_spaces(l))

    if len(sys.argv) > 1:
        for a in sys.argv[1:]:
            for f in a.split(','):
                if os.path.isdir(f):
                    for ff in glob.iglob(os.path.join(f, "*fastq*")):
                        with fopen(ff) as inf:
                            for l in inf:
                                sys.stdout.write(ignore_spaces(l))
                else:
                    with fopen(f) as inf:
                        for l in inf:
                            sys.stdout.write(ignore_spaces(l))
