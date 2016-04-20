#!/usr/bin/env python
#Author: Duy Tin Truong (duytin.truong@unitn.it)
#        at CIBIO, University of Trento, Italy

__author__  = 'Duy Tin Truong (duytin.truong@unitn.it)'
__version__ = '0.1'
__date__    = '17 Sep 2015'

import sys
import os
import argparse 
import numpy
from Bio import SeqIO
import glob


def read_params():
    p = argparse.ArgumentParser()
    p.add_argument(
        '--ifn_alignments', 
        nargs='+',
        required=True, 
        default=None, 
        type=str,
        help='The alignment file.')
    p.add_argument(
        '--log_ofn', 
        required=True, 
        default=None, 
        type=str,
        help='The log file.')
    p.add_argument(
        '--nprocs', 
        required=True, 
        default=None, 
        type=int,
        help='Number of processors.')
    p.add_argument(
        '--bootstrap_raxml', 
        required=False, 
        default=0, 
        type=int,
        help='The number of runs for bootstraping when building the tree. '\
             'Default 0.')
    p.add_argument(
        '--verbose', 
        required=False, 
        dest='quiet',
        action='store_false',
        help='Show all information. Default "not set".')
    p.set_defaults(quiet=True)

    return p.parse_args()


def run(cmd):
    print cmd
    os.system(cmd)


def main(args):
    lfile = open(args.log_ofn, 'w')
    for ifn_alignment in args.ifn_alignments:
        if 'remove_' in ifn_alignment:
            continue
        sample2polrate = {}
        ifn_polymorphic = ifn_alignment.replace('.fasta', '.polymorphic')
        singles = set([])
        with open(ifn_polymorphic, 'r') as ifile:
            for line in ifile:
                if line[0] == '#':
                    continue
                line = line.strip().split()
                val = float(line[1])
                if line[0][:3] in ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__', 't__']:
                    singles.add(line[0])
                    continue
                sample2polrate[line[0]] = val
        median = numpy.median(sample2polrate.values())
        std = numpy.std(sample2polrate.values())
        for s in sample2polrate:
            if sample2polrate[s] <= median + std:
                singles.add(s)

        if len(sample2polrate):
            log_line = '%s\t%d\t%d\t%f\n'%\
                        (os.path.basename(ifn_polymorphic).replace('.polymorphic', ''), 
                        len(singles), 
                        len(sample2polrate), 
                        float(len(singles)) / len(sample2polrate))
        else:
            log_line = '%s\t%d\t%d\t%f\n'%\
                        (os.path.basename(ifn_polymorphic).replace('.polymorphic', ''), 
                        len(singles), 
                        len(sample2polrate), 
                        0)
        lfile.write(log_line)

        ifn_alignment2 = ifn_alignment.replace('.fasta', '') + '.remove_multiple_strains.fasta'
        with open(ifn_alignment2, 'w') as ofile:
            for rec in SeqIO.parse(open(ifn_alignment, 'r'), 'fasta'):
                if rec.name in singles:
                    SeqIO.write(rec, ofile, 'fasta')

        with open(ifn_alignment2 + '.log', 'w') as ofile:
            ofile.write(log_line)

        output_suffix = os.path.basename(ifn_alignment2).replace('.polymorphic', '').replace('.fasta', '')
        output_suffix += '.tree'
        if args.bootstrap_raxml:
            cmd = 'raxmlHPC-PTHREADS-SSE3 '
            cmd += '-f a '
            cmd += '-m GTRGAMMA '
            #cmd += '-b 1234 '
            cmd += '-x 1234 '
            cmd += '-N %d '%(args.bootstrap_raxml)
            cmd += '-s %s '%os.path.abspath(ifn_alignment2)
            cmd += '-w %s '%os.path.abspath(os.path.dirname(ifn_alignment2))
            cmd += '-n %s '%output_suffix 
            cmd += '-p 1234 '
        else:
            cmd = 'raxmlHPC-PTHREADS-SSE3 '
            cmd += '-m GTRCAT '
            cmd += '-s %s '%os.path.abspath(ifn_alignment2)
            cmd += '-w %s '%os.path.abspath(os.path.dirname(ifn_alignment2))
            cmd += '-n %s '%output_suffix
            cmd += '-p 1234 '
        if args.nprocs:
            cmd += '-T %d '%(args.nprocs)
        raxfns = glob.glob('%s/RAxML_*%s*'%(os.path.dirname(ifn_alignment2), output_suffix))
        for fn in raxfns:
            os.remove(fn)
        '''
        if len(raxfns) == 0:
            run(cmd)
        '''
        run(cmd)
    lfile.close()
     




if __name__ == "__main__":
    args = read_params()
    main(args)
