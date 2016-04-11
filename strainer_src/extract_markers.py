#!/usr/bin/env python
#Author: Duy Tin Truong (duytin.truong@unitn.it)
#        at CIBIO, University of Trento, Italy

__author__  = 'Duy Tin Truong (duytin.truong@unitn.it)'
__version__ = '0.1'
__date__    = '1 Sep 2014'

import sys
import os
import argparse as ap
import pickle
import bz2
from Bio import SeqIO, Seq, SeqRecord

def read_params():
    p = ap.ArgumentParser()
    p.add_argument('--mpa_pkl', required=True, default=None, type=str)
    p.add_argument('--ifn_markers', required=True, default=None, type=str)
    p.add_argument('--clade', required=True, default=None, type=str)
    p.add_argument('--ofn_markers', required=True, default=None, type=str)
    return vars(p.parse_args())


def extract_markers(mpa_pkl, ifn_markers, clade, ofn_markers):
    with open(mpa_pkl, 'rb') as ifile:
        db = pickle.loads(bz2.decompress(ifile.read()))
    markers = set([])
    for marker in db['markers']:
        if clade == db['markers'][marker]['taxon'].split('|')[-1]:
            markers.add(marker)
    print 'number of markers', len(markers)
    with open(ofn_markers, 'w') as ofile:
        for rec in SeqIO.parse(open(ifn_markers, 'r'), 'fasta'):
            if rec.name in markers:
                SeqIO.write(rec, ofile, 'fasta')


if __name__ == "__main__":
    args = read_params()
    extract_markers(
                    mpa_pkl=args['mpa_pkl'],
                    ifn_markers=args['ifn_markers'],
                    clade=args['clade'],
                    ofn_markers=args['ofn_markers'])
