#!/usr/bin/env python


import sys
import os
import bz2
import gzip
import glob
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser
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
    if l:
        if l[0] == '@':
            return 'fastq'

        if l[0] == '>':
            return 'fasta'

    sys.stderr.write("\nError, input data has to be in fastq or fasta format\n\n")
    sys.exit(-1)

def print_record(description, sequence, qual, fmt):
    if fmt == 'fastq':
        return '@{}\n{}\n+\n{}\n'.format(description, sequence, qual)

    if fmt == 'fasta':
        return '>{}\n{}\n'.format(description, sequence)

def fopen(fn):
    fileName, fileExtension = os.path.splitext(fn)

    if fileExtension == '.bz2':
        return (bz2.open(fn, "rt") if not p2 else bz2.BZ2File(fn, "r"))
    if fileExtension == '.gz':
        return gzip.open(fn, "rt")

    return open(fn)


def read_and_write_raw_int(fd, min_len=None):
    fmt = None
    avg_read_length = 0
    nreads = 0
    discarded = 0
    #if min_len:
    r = []

    while True:
        l = fd.readline()

        if not fmt:
            fmt = fastx(l)
            parser = FastqGeneralIterator if fmt == 'fastq' else SimpleFastaParser
            readn = 4 if fmt == 'fastq' else 2

        r.append(l)

        if len(r) == readn:
            break

    for record in parser(uio.StringIO("".join(r))):
        if readn == 4:
            description, sequence, qual = record
        else:
            qual = None
            description, sequence = record

        if len(sequence) >= min_len:
            description = ignore_spaces(description, forced=True)
            avg_read_length = len(sequence) + avg_read_length
            _ = sys.stdout.write(print_record(description, sequence, qual, fmt))
        else:
            discarded = discarded + 1
    
    for idx, record in enumerate(parser(fd),2):
        if readn == 4:
            description, sequence, qual = record
        else:
            qual = None
            description, sequence = record

        if len(sequence) >= min_len:
            avg_read_length = len(sequence) + avg_read_length
            description = ignore_spaces(description, forced=True)
            _ = sys.stdout.write(print_record(description, sequence, qual, fmt))
        else:
            discarded = discarded + 1
    # else:
    #     for idx, l in enumerate(fd,1):
    #         avg_read_length = len(l) + avg_read_length
    #         _ = sys.stdout.write(ignore_spaces(l))

    nreads = idx - discarded
    avg_read_length /= nreads
    return (nreads, avg_read_length)

def read_and_write_raw(fd, opened=False, min_len=None):
    if opened:  #fd is stdin
        nreads, avg_read_length = read_and_write_raw_int(fd, min_len=min_len)
    else:
        with fopen(fd) as inf:
            nreads, avg_read_length = read_and_write_raw_int(inf, min_len=min_len)
    return  (nreads, avg_read_length)


def main():
    min_len = 0
    args = []
    nreads = None
    avg_read_length = None

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
        nreads, avg_read_length = read_and_write_raw(sys.stdin, opened=True, min_len=min_len)
    else:
        files = []
        avg_read_length, nreads = 0, 0
        for a in args:
            for f in a.split(','):
                if os.path.isdir(f):
                    files += list(glob.iglob(os.path.join(f, "*fastq*")))
                else:
                    files += [f]

        for f in files:
            f_nreads, f_avg_read_length = read_and_write_raw(f, opened=False, min_len=min_len)
            nreads += f_nreads
            avg_read_length += f_avg_read_length
        avg_read_length = avg_read_length / len(files)

    if nreads and avg_read_length:
        sys.stderr.write('{}\t{}'.format(nreads,avg_read_length) )
    else:
        exit(1)

if __name__ == '__main__':
    main()