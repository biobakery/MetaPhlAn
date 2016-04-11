#!/usr/bin/env python
# Author: Duy Tin Truong (duytin.truong@unitn.it)
#		at CIBIO, University of Trento, Italy

__author__ = 'Duy Tin Truong (duytin.truong@unitn.it)'
__version__ = '0.1'
__date__ = '1st Sep 2014'

import numpy
import sys

def dist2file(dist, labels, ofn):
    with open(ofn, 'w') as ofile:
        ofile.write('ID')
        for label in labels:
            ofile.write('\t%s'%label)
        ofile.write('\n')
        for i in range(len(labels)):
            ofile.write('%s\t'%labels[i])
            for j in range(len(labels)):
                if j == len(labels) - 1:
                    ofile.write('%f\n'%dist[i][j])
                else:
                    ofile.write('%f\t'%dist[i][j])



def statistics(vals):
    vals = numpy.array(vals)
    result = {}
    if len(vals.shape) == 1:
        num_elems = len(vals)
        nrows = num_elems
        ncols = 1
    else:
        nrows, ncols = vals.shape
        num_elems = nrows * ncols
    if num_elems > 0:
        result['nrows'] = nrows
        result['ncols'] = ncols
        result['size'] = num_elems
        result['average'] = numpy.average(vals)
        result['min'] = numpy.min(vals)
        result['max'] = numpy.max(vals)
        result['median'] = numpy.percentile(vals, 50)
        result['percentile_25'] = numpy.percentile(vals, 25)
        result['percentile_75'] = numpy.percentile(vals, 75)
    else:
        result['nrows'] = nrows
        result['ncols'] = ncols
        result['size'] = num_elems
        result['average'] = 0
        result['min'] = 0
        result['max'] = 0
        result['median'] = 0
        result['percentile_25'] = 0
        result['percentile_75'] = 0

    str_result = ''
    for key in ['nrows',
                'ncols',
                'size',
                'average',
                'min',
                'max',
                'median',
                'percentile_25',
                'percentile_75']:
        str_result += '%s: %s\n'%(key, result[key])

    return result, str_result



def dict2str(dict_var):
    result = ''
    for key in dict_var:
        result += '%s: %s\n'%(key, dict_var[key])
    return result


def openr( fn, mode = "r" ):
    if fn is None:
        return sys.stdin
    return bz2.BZ2File(fn) if fn.endswith(".bz2") else open(fn,mode)
    

def openw( fn ):
    if fn is None:
        return sys.stdout
    return bz2.BZ2File(fn,"w") if fn.endswith(".bz2") else open(fn,"w")
            

def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False
