#!/usr/bin/env python
#Author: Duy Tin Truong (duytin.truong@unitn.it)
#        at CIBIO, University of Trento, Italy

__author__  = 'Duy Tin Truong (duytin.truong@unitn.it)'
__version__ = '0.1'
__date__    = '18 Jul 2015'

import sys
import os
import argparse 


def read_params():
    p = argparse.ArgumentParser()
    p.add_argument(
        '--input_file', 
        required=False, 
        default=None, 
        type=str,
        help='The input sam file.')
    p.add_argument(
        '--min_align_score', 
        required=True, 
        default=None, 
        type=int,
        help='The sam records with alignment score smaller than this value '
             'will be discarded.')
    p.add_argument(
        '--verbose', 
        required=False, 
        dest='quiet',
        action='store_false',
        help='Show all information. Default "not set".')
    p.set_defaults(quiet=True)

    return p.parse_args()


def main(args):
    if args.input_file == None:
        ifile = sys.stdin
    else:
        ifile = open(args.input_file, 'r')
    for line in ifile:
        if line[0] == '@':
            sys.stdout.write(line)
        else:
            spline = line.split()
            read_length = len(spline[9])
            align_score = float(spline[11].replace('AS:i:', ''))
            if align_score >= args.min_align_score * read_length / 100.0:
                sys.stdout.write(line)



if __name__ == "__main__":
    args = read_params()
    main(args)
