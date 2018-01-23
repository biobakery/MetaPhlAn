#!/usr/bin/env python


import sys
import os
import bz2
import gzip
import glob
from Bio import SeqIO
try:
    import StringIO as uio
except ImportError:
    import io as uio


p2 = float(sys.version_info[0]) < 3.0


def ignore_spaces(l, forced=False):
    if (l[0] == '@') or (l[0] == '>') or forced:
        return l.replace(' ', '_')

    return l


def fastx(l):
    if l[0] == '@':
        return 'fastq'

    if l[0] == '>':
        return 'fasta'

    sys.stderr.write("\nError, input data has to be in fastq or fasta format\n")
    sys.exit(-1)


def fopen(fn):
    fileName, fileExtension = os.path.splitext(fn)

    if fileExtension == '.bz2':
        return (bz2.open(fn, "rt") if not p2 else bz2.BZ2File(fn, "r"))

    if fileExtension == '.gz':
        return gzip.open(fn, "rt")

    return open(fn)


def read_and_write_raw_int(fd, min_len=None):
    fmt = None

    if min_len:
        r = []

        while True:
            l = fd.readline()

            if not fmt:
                fmt = fastx(l)
                readn = (4 if fmt == 'fastq' else 2)

            r.append(l)

            if len(r) == readn:
                break

        for record in SeqIO.parse(uio.StringIO("".join(r)), fmt):
            if len(record) < min_len:
                continue

            record.id = ignore_spaces(record.description, forced=True)
            record.description = ""
            SeqIO.write(record, sys.stdout, fmt)

        for record in SeqIO.parse(fd, fmt):
            if len(record) < min_len:
                continue

            record.id = ignore_spaces(record.description, forced=True)
            record.description = ""
            SeqIO.write(record, sys.stdout, fmt)
    else:
        for l in fd:
            sys.stdout.write(ignore_spaces(l))


def read_and_write_raw(fd, opened=False, min_len=None):
    if opened:
        read_and_write_raw_int(fd, min_len=min_len)
    else:
        with fopen(fd) as inf:
            read_and_write_raw_int(inf, min_len=min_len)


if __name__ == '__main__':
    min_len = None
    args = []

    if len(sys.argv) > 1:
        for l in sys.argv[1:]:
            if l in ['-h', '--help', '-v', '--version']:
                sys.stderr.write("Help message for " +
                                 os.path.basename(sys.argv[0]) + "\n")
                sys.exit(0)

            if min_len == 'next':
                min_len = int(l)
            elif l in ['-l', '--min_len']:
                min_len = 'next'
            else:
                args.append(l)

    if len(args) == 0:
        read_and_write_raw(sys.stdin, opened=True, min_len=min_len)
    else:
        files = []

        for a in args:
            for f in a.split(','):
                if os.path.isdir(f):
                    files += list(glob.iglob(os.path.join(f, "*fastq*")))
                else:
                    files += [f]

        for f in files:
            read_and_write_raw(f, opened=False, min_len=min_len)
