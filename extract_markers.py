#!/usr/bin/env python
#Author: Duy Tin Truong (duytin.truong@unitn.it)
#        at CIBIO, University of Trento, Italy

__author__  = 'Duy Tin Truong (duytin.truong@unitn.it)'
__version__ = '0.1'
__date__    = '27th Nov 2014'

import sys
import os
import argparse as ap
import pickle
import bz2


def read_params():
    p = ap.ArgumentParser()
    p.add_argument('--clades', nargs='+', required=True, type=str)
    p.add_argument('--mpa_pkl', required=True, type=str)
    p.add_argument('--ofn_markers', required=True, type=str)

    return vars(p.parse_args())


def main(args):
    with open(args['mpa_pkl'], 'rb') as ifile:
        db = pickle.loads(bz2.decompress(ifile.read()))
    with open(args['ofn_markers'], 'w') as ofile:
        for marker in db['markers']:
            for clade in args['clades']:
                if clade in db['markers'][marker]['taxon']:
                    ofile.write('%s\n'%marker)
                    break



if __name__ == "__main__":
    args = read_params()
    main(args)
