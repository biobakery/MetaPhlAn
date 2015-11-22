#!/usr/bin/env python
#Author: Duy Tin Truong (duytin.truong@unitn.it)
#        at CIBIO, University of Trento, Italy


import sys
import argparse

def read_params():
    p = argparse.ArgumentParser()
    p.add_argument('--input_file', required=True, default='-', type=str)

    return vars(p.parse_args())


def fix_AF1(ifn):
    if ifn == '-':
        ifile = sys.stdin
    else:
        ifile = open(ifn, 'r')

    for line in ifile:
        if line[0] != '#':
            if 'AF1=0' in line:
                spline = line.split()
                if spline[3] != spline[4] and spline[4].upper() in ['A', 'T', 'C', 'G']:
                    line = line.replace('AF1=0', 'AF1=1')
        sys.stdout.write(line)

    if ifn != '-':
        ifile.close()


if __name__ == "__main__":
    args = read_params()
    fix_AF1(args['input_file'])
