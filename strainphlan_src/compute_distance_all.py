#!/usr/bin/env python
#Author: Duy Tin Truong (duytin.truong@unitn.it)
#        at CIBIO, University of Trento, Italy

__author__  = 'Duy Tin Truong (duytin.truong@unitn.it)'
__version__ = '0.1'
__date__    = '9 Feb 2015'

import sys
import os
ABS_PATH = os.path.abspath(sys.argv[0])
MAIN_DIR = os.path.dirname(ABS_PATH)
os.environ['PATH'] += ':%s'%MAIN_DIR
sys.path.append(MAIN_DIR)
import argparse as ap
from Bio import SeqIO, Seq, SeqRecord
from collections import defaultdict
import numpy
from compute_distance import compute_dist_matrix
from ooSubprocess import parallelize


def read_params():
    p = ap.ArgumentParser()
    p.add_argument('--ifn_alignments', nargs='+', required=True, default=None, type=str)
    p.add_argument('--nprocs', required=True, default=None, type=int)
    p.add_argument('--count_gaps', 
                   required=False,
                   dest='ignore_gaps', 
                   action='store_false')
    p.set_defaults(ignore_gaps=True)


    return vars(p.parse_args())



def compute_dist_matrix_wrapper(args):
    compute_dist_matrix(
                        args['ifn_alignment'], 
                        args['ofn_prefix'],
                        args['ignore_gaps'],
                        overwrite=True) 
    


def main(args):
    args_list = []
    for i in range(len(args['ifn_alignments'])):
        args_list.append({})
        args_list[i]['ifn_alignment'] = args['ifn_alignments'][i]
        args_list[i]['ofn_prefix'] = args['ifn_alignments'][i] 
        if not args['ignore_gaps']:
            args_list[i]['ofn_prefix'] += '.count_gaps'
        args_list[i]['ignore_gaps'] = args['ignore_gaps']

    parallelize(compute_dist_matrix_wrapper, args_list, args['nprocs'])



if __name__ == "__main__":
    args = read_params()
    main(args)
