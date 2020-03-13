#!/usr/bin/env python
__author__ = ('Duy Tin Truong (duytin.truong@unitn.it), '
              'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it)')
__version__ = '3.0'
__date__    = '21 Feb 2020'

import argparse as ap
import pandas
import dendropy
import numpy

def read_params():
    p = ap.ArgumentParser()
    p.add_argument('-t', '--ifn_trees', nargs='+', required=True, default=None, type=str)
    p.add_argument('-f', '--ifn_metadatas', nargs='+', required=True, default=None, type=str)
    p.add_argument('--string_to_remove',
                   required=False, default='', type=str,
                   help='string to be removed in the tree node names')
    p.add_argument('-m', '--metadatas',
                    nargs='+',
                    required=False,
                    default=['all'],
                    type=str,
                    help='The metadata fields that you want to add. '\
                          'Default: add all metadata from the first line.')
    return vars(p.parse_args())


def get_index_col(ifn):
    with open(ifn, 'r') as ifile:
        line = ifile.readline()
        line = line.strip().split()
        for i in range(len(line)):
            if line[i].upper() == 'SAMPLEID':
                return i
    return -1


def main():
    args = read_params()
    add_fields = args['metadatas']
    for ifn_tree in args['ifn_trees']:
        print ('Input:', ifn_tree)
        df_list = []
        samples = []
        for ifn in args['ifn_metadatas']:
            index_col = get_index_col(ifn)
            df = pandas.read_csv(
                ifn,
                sep='\t',
                dtype=numpy.unicode_,
                header=0,
                index_col=index_col)
            df = df.transpose()
            df_list.append(df)
            samples += df.columns.values.tolist()
            if add_fields == ['all']:
                with open(ifn, 'r') as ifile:
                    add_fields = [f for f in ifile.readline().strip().split('\t') \
                                  if f.upper() != 'SAMPLEID']
        print ('number of samples in metadata: %d'%len(samples))
        count = 0
        with open(ifn_tree, 'r') as ifile:
            line = ifile.readline()
        line = line.replace(args['string_to_remove'], '')
        tree = dendropy.Tree.get(stream=open(ifn_tree, 'r'), schema='newick')
        for node in tree.leaf_nodes():
            sample = node.__getattribute__("taxon").__str__().strip("'")
            sample = sample.replace(' ', '_')
            sample = sample.replace(args['string_to_remove'], '')
            prefixes = [prefix for prefix in
                            ['k__', 'p__', 'c__', 'o__',
                             'f__', 'g__', 's__'] \
                        if prefix in sample]

            metadata = sample
            if len(prefixes) == 0:
                count += 1
                for meta in add_fields:
                    old_meta = meta
                    for i in range(len(df_list)):
                        if sample in df_list[i].columns.values.tolist():
                            df = df_list[i]
                            if meta.lower() in df[sample]:
                                meta = meta.lower()
                            elif meta.upper() in df[sample]:
                                meta = meta.upper()
                            elif meta.title() in df[sample]:
                                meta = meta.title()
                            if meta in df[sample]:
                                metadata += '|%s-%s'%(
                                                old_meta,
                                                str(df[sample][meta]).replace(':','_'))
                                break # take the first metadata

            line = line.replace(sample + ':', metadata + ':')

        ofn_tree = ifn_tree + '.metadata'
        print ('Number of samples in tree: %d'%count)
        print ('Output:', ofn_tree)
        with open(ofn_tree, 'w') as ofile:
            ofile.write(line)

if __name__ == '__main__':
    main()