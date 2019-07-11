#!/usr/bin/env python
#Author: Duy Tin Truong (duytin.truong@unitn.it)
#        at CIBIO, University of Trento, Italy

__author__ = ('Duy Tin Truong (duytin.truong@unitn.it), '
              'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it)')
__version__ = '0.2'
__date__    = '10 Jul 19'


from Bio import SeqIO
import argparse as ap
import sys

def read_params(args):
	p = ap.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
	p.add_argument('--min_len', required = True, default = None, type = int)
	return vars(p.parse_args())
	
if __name__ == '__main__':
	args = read_params(sys.argv)
	min_len = args['min_len']
	with sys.stdout as outf:
		list_r = []
		for r in SeqIO.parse(sys.stdin, "fastq"):
			if len(r) >= min_len:
				list_r.append(r)
		SeqIO.write(list_r, outf, "fastq")
