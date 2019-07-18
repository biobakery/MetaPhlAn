#!/usr/bin/env python
#Author: Duy Tin Truong (duytin.truong@unitn.it)
#        at CIBIO, University of Trento, Italy

__author__ = ('Duy Tin Truong (duytin.truong@unitn.it), '
              'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it)')
__version__ = '0.2'
__date__    = '10 Jul 19'

import sys
import argparse as ap
import bz2
import gzip
import tarfile
import ooSubprocess

PYTHON_VERSION = float(sys.version_info[0])

def read_params():
    p = ap.ArgumentParser()
    p.add_argument('--input_file', required=True, default=None, type=str)

    return vars(p.parse_args())


def dump_file(ifn):
    file_ext = ''
    if ifn.endswith('.tar.bz2'):
        ifile = tarfile.open(ifn, 'r:bz2')
        file_ext = '.tar.bz2'
    elif ifn.endswith('.tar.gz'):
        ifile = tarfile.open(ifn, 'r:gz')
        file_ext = '.tar.gz'
    elif ifn.endswith('.bz2'):
        if(PYTHON_VERSION < 3):
            ifile = bz2.BZ2File(ifn, 'r')
        else:    
            ifile = bz2.open(ifn, 'rt')
        file_ext = '.bz2'
    elif ifn.endswith('.gz'):
        if(PYTHON_VERSION < 3):
            ifile = gzip.BZ2File(ifn, 'r')
        else:    
            ifile = gzip.open(ifn, 'rt')
        file_ext = '.gz'
    elif ifn.endswith('.fastq'):
        ifile = open(ifn, 'r')
        file_ext = '.fastq'
    elif ifn.endswith('.sam'):
        ifile = open(ifn, 'r')
        file_ext = '.sam'
    elif ifn.endswith('.sra'):
        oosp = ooSubprocess.ooSubprocess()
        ifile = oosp.ex(
                        'fastq-dump',
                        args=[
                                '-Z', ifn,
                                '--split-spot'],
                        get_out_pipe=True)
        file_ext = '.sra'
    else:
        raise Exception('Unrecognized format! The format should be .bz2, .gz'\
                '.tar.bz2, .tar.gz, .sra, .sam.bz2, .sam, or .fastq\n')

    try:
        if file_ext in ['.tar.bz2', '.tar.gz']:
            for tar_info in ifile:
                ifile2 = ifile.extractfile(tar_info)
                if ifile2 != None:
                    for line in ifile2:
                        sys.stdout.write(line)
        else:
            for line in ifile:
                sys.stdout.write(line)
    except:
        sys.stderr.write('Error while dumping file %s\n'%ifn)
        raise


if __name__ == "__main__":
    args = read_params()
    dump_file(args['input_file'])
