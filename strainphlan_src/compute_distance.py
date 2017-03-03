#!/usr/bin/env python
#Author: Duy Tin Truong (duytin.truong@unitn.it)
#        at CIBIO, University of Trento, Italy

__author__  = 'Duy Tin Truong (duytin.truong@unitn.it)'
__version__ = '0.1'
__date__    = '1 Sep 2014'

import sys
import os
ABS_PATH = os.path.abspath(sys.argv[0])
MAIN_DIR = os.path.dirname(ABS_PATH)
os.environ['PATH'] += ':%s'%MAIN_DIR
sys.path.append(MAIN_DIR)

from mixed_utils import dist2file, statistics
import argparse as ap
from Bio import SeqIO, Seq, SeqRecord
from collections import defaultdict
import numpy
from ooSubprocess import ooSubprocess


'''
SUBST = {
        'A':{'A':0.0, 'C':1.0, 'G':1.0, 'T':1.0, '-':1.0},
        'C':{'A':1.0, 'C':0.0, 'G':1.0, 'T':1.0, '-':1.0}, 
        'G':{'A':1.0, 'C':1.0, 'G':0.0, 'T':1.0, '-':1.0}, 
        'T':{'A':1.0, 'C':1.0, 'G':1.0, 'T':0.0, '-':1.0},
        '-':{'A':1.0, 'C':1.0, 'G':1.0, 'T':1.0, '-':0.0}}
'''


def read_params():
    p = ap.ArgumentParser()
    p.add_argument('--ifn_alignment', required=True, default=None, type=str)
    p.add_argument('--ofn_prefix', required=True, default=None, type=str)
    p.add_argument('--count_gaps', 
                   required=False,
                   dest='ignore_gaps', 
                   action='store_false')
    p.set_defaults(ignore_gaps=True)
    p.add_argument('--overwrite', 
                   required=False,
                   dest='overwrite', 
                   action='store_true')
    p.set_defaults(overwrite=False)

    return vars(p.parse_args())


def get_dist(seq1, seq2, ignore_gaps):
    if len(seq1) != len(seq2):
        print >> sys.stderr, 'Error: Two sequences have different lengths!'
        print >> sys.stderr, 'Cannot compute the distance!'
        exit(1)

    abs_dist = 0.0
    abs_snp = 0
    for i in range(len(seq1)):
        if seq1[i] != '-' and seq2[i] != '-':
            abs_snp += 1
        if seq1[i] != seq2[i]:
            if ignore_gaps:
                if seq1[i] != '-' and seq2[i] != '-':
                    abs_dist += 1.0
            else:
                abs_dist += 1.0

    abs_sim = len(seq1) - abs_dist
    rel_dist = abs_dist / float(len(seq1))
    rel_sim = 1.0 - rel_dist
    rel_snp = abs_snp / float(len(seq1))
    return abs_dist, rel_dist, abs_sim, rel_sim, abs_snp, rel_snp
    

def compute_dist_matrix(ifn_alignment, ofn_prefix, ignore_gaps, overwrite):
    ofn_abs_dist = ofn_prefix + '.abs_dist'
    if (not overwrite) and os.path.isfile(ofn_abs_dist.replace('.abs_dist', '.rel_dist')):
        print 'File %s exists, skip!'%ofn_abs_dist
        return
    else:
        print 'Compute dist_matrix for %s'%ofn_abs_dist
    #print 'Compute dist_matrix for %s'%ofn_abs_dist

    recs = [rec for rec in SeqIO.parse(open((ifn_alignment), 'r'), 'fasta')]
    abs_dist = numpy.zeros((len(recs), len(recs)))
    abs_dist_flat = []
    rel_dist = numpy.zeros((len(recs), len(recs)))
    rel_dist_flat = []

    abs_sim = numpy.zeros((len(recs), len(recs)))
    abs_sim_flat = []
    rel_sim = numpy.zeros((len(recs), len(recs)))
    rel_sim_flat = []

    abs_snp = numpy.zeros((len(recs), len(recs)))
    abs_snp_flat = []
    rel_snp = numpy.zeros((len(recs), len(recs)))
    rel_snp_flat = []

    for i in range(len(recs)):
        for j in range(i, len(recs)):
            abs_d, rel_d, abs_s, rel_s, abs_sp, rel_sp = get_dist(recs[i].seq, 
                                                                recs[j].seq,
                                                                ignore_gaps)

            abs_dist[i][j] = abs_d
            abs_dist[j][i] = abs_d
            abs_dist_flat.append(abs_d)
            
            rel_dist[i][j] = rel_d
            rel_dist[j][i] = rel_d
            rel_dist_flat.append(rel_d)

            abs_sim[i][j] = abs_s
            abs_sim[j][i] = abs_s
            abs_sim_flat.append(abs_s)
            
            rel_sim[i][j] = rel_s
            rel_sim[j][i] = rel_s
            rel_sim_flat.append(rel_s)

            abs_snp[i][j] = abs_sp
            abs_snp[j][i] = abs_sp
            abs_snp_flat.append(abs_sp)
            
            rel_snp[i][j] = rel_sp
            rel_snp[j][i] = rel_sp
            rel_snp_flat.append(rel_sp)

    labels = [rec.name for rec in recs]
    oosp = ooSubprocess()

    ofn_abs_dist = ofn_prefix + '.abs_dist'
    dist2file(abs_dist, labels, ofn_abs_dist)
    with open(ofn_abs_dist + '.info', 'w') as ofile:
        ofile.write(statistics(abs_dist_flat)[1])
    '''
    if len(abs_dist_flat) > 0:
        oosp.ex('hclust2.py',
        args=['-i', ofn_abs_dist,
              '-o', ofn_abs_dist + '.png',
              '--f_dist_f', 'euclidean',
              '--s_dist_f', 'euclidean',
              '-l', '--dpi', '300',
              '--flabel_size', '5',
              '--slabel_size', '5',
              '--max_flabel_len', '200'])
    '''

    ofn_rel_dist = ofn_prefix + '.rel_dist'
    dist2file(rel_dist, labels, ofn_rel_dist)
    with open(ofn_rel_dist + '.info', 'w') as ofile:
        ofile.write(statistics(rel_dist_flat)[1])
    '''
    if len(rel_dist_flat) > 0:
        oosp.ex('hclust2.py',
        args=['-i', ofn_rel_dist,
              '-o', ofn_rel_dist + '.png',
              '--f_dist_f', 'euclidean',
              '--s_dist_f', 'euclidean',
              '-l', '--dpi', '300',
              '--flabel_size', '5',
              '--slabel_size', '5',
              '--max_flabel_len', '200'])
    '''

    ofn_abs_sim = ofn_prefix + '.abs_sim'
    dist2file(abs_sim, labels, ofn_abs_sim)
    with open(ofn_abs_sim + '.info', 'w') as ofile:
        ofile.write(statistics(abs_sim_flat)[1])
    '''
    if len(abs_sim_flat) > 0:
        oosp.ex('hclust2.py',
        args=['-i', ofn_abs_sim,
              '-o', ofn_abs_sim + '.png',
              '--f_dist_f', 'euclidean',
              '--s_dist_f', 'euclidean',
              '-l', '--dpi', '300',
              '--flabel_size', '5',
              '--slabel_size', '5',
              '--max_flabel_len', '200'])
    '''

    ofn_rel_sim = ofn_prefix + '.rel_sim'
    dist2file(rel_sim, labels, ofn_rel_sim)
    with open(ofn_rel_sim + '.info', 'w') as ofile:
        ofile.write(statistics(rel_sim_flat)[1])
    '''
    if len(rel_sim_flat) > 0:
        oosp.ex('hclust2.py',
        args=['-i', ofn_rel_sim,
              '-o', ofn_rel_sim + '.png',
              '--f_dist_f', 'euclidean',
              '--s_dist_f', 'euclidean',
              '-l', '--dpi', '300',
              '--flabel_size', '5',
              '--slabel_size', '5',
              '--max_flabel_len', '200'])
    '''       

    ofn_abs_snp = ofn_prefix + '.abs_snp'
    dist2file(abs_snp, labels, ofn_abs_snp)
    with open(ofn_abs_snp + '.info', 'w') as ofile:
        ofile.write(statistics(abs_snp_flat)[1])
    ofn_rel_snp = ofn_prefix + '.rel_snp'
    dist2file(rel_snp, labels, ofn_rel_snp)
    with open(ofn_rel_snp + '.info', 'w') as ofile:
        ofile.write(statistics(rel_snp_flat)[1])
 



def main(args):
    compute_dist_matrix(
                        args['ifn_alignment'], 
                        args['ofn_prefix'],
                        args['ignore_gaps'],
                        args['overwrite']) 
    
if __name__ == "__main__":
    args = read_params()
    main(args)
