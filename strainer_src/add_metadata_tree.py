#!/usr/bin/env python
#Author: Duy Tin Truong (duytin.truong@unitn.it)
#        at CIBIO, University of Trento, Italy


import sys
import os
import argparse as ap
import pandas
import copy
import ConfigParser
import dendropy
import numpy


def read_params():
    p = ap.ArgumentParser()
    p.add_argument('--ifn_trees', nargs='+', required=True, default=None, type=str)
    p.add_argument('--ifn_metadatas', nargs='+', required=True, default=None, type=str)
    p.add_argument('--string_to_remove', 
                   required=False, default='', type=str,
                   help='string to be removed in the tree node names')
    p.add_argument(
                    '--metadatas', 
                    nargs='+', 
                    required=False, 
                    default=[
                                'subjectID', 
                                'visit_number', 
                                'bodysite', 
                                'gender',
                                'disease',
                                'country',
                                'age',
                                'WMSphase',
                                'dsuvspd',
                                'ethnicity'
                                ],
                    type=str,
                    help='The metadata fields that you want to add. '\
                            'Default: subjectID, visit_number, bodysite, '\
                            'gender, disease, country, age, WMSphase, '\
                            'dsuvspd, ethnicity')
    p.add_argument('--header_first_col', dest='header_first_col', required=False, action='store_true')
    p.set_defaults(header_first_col=False)
    p.add_argument('--add_cluster_info', dest='add_cluster_info', required=False, action='store_true')
    p.set_defaults(add_cluster_info=False)
    p.add_argument('--add_dataset_info', dest='add_dataset_info', required=False, action='store_true')
    p.set_defaults(add_dataset_info=False)
    p.add_argument(
                   '--add_multiple_strains_info',
                   dest='add_multiple_strains_info', 
                   required=False, 
                   action='store_true')
    p.set_defaults(add_multiple_strains_info=False)


    return vars(p.parse_args())


def get_index_col(ifn):
    with open(ifn, 'r') as ifile:
        line = ifile.readline()
        line = line.strip().split()
        for i in range(len(line)):
            if line[i] == 'sampleID':
                return i
    return -1


def main(args):
    if args['add_dataset_info']:
        args['metadatas'].append('dataset')

    for ifn_tree in args['ifn_trees']:
        print 'Input:', ifn_tree

        sample2cluster = {}
        # read clusters
        if args['add_cluster_info']:
            if 'remove_' in ifn_tree:
                clade = ifn_tree.split('.')[-3] + '.' + ifn_tree.split('.')[-2] 
            else:
                clade = ifn_tree.split('.')[-2]
            idir = os.path.dirname(ifn_tree)
            cfn = '%s/%s.clusters'%(idir, clade)
            if os.path.isfile(cfn):
                config = ConfigParser.ConfigParser()
                config.readfp(open(cfn, 'r'))
                sample2cluster = eval(config.get('info', 'node2cluster'))
            else:
                print 'File not found %s'%cfn
                continue
        print 'len(sample2cluster)', len(sample2cluster)

        singles = set([])
        if args['add_multiple_strains_info'] and 'remove_' not in ifn_tree:
            sample2polrate = {}
            clade = ifn_tree.split('.')[-2]
            idir = os.path.dirname(ifn_tree)
            pfn = '%s/%s.polymorphic'%(idir, clade)
            if os.path.isfile(pfn):
                with open(pfn, 'r') as ifile:
                    for pline in ifile:
                        if pline[0] == '#':
                            continue
                        pline = pline.strip().split()
                        val = float(pline[1])
                        if pline[0][:3] in ['s__', 'g__']:
                            singles.add(pline[0])
                            continue
                        sample2polrate[pline[0]] = val
                median = numpy.median(sample2polrate.values())
                std = numpy.std(sample2polrate.values())
                for s in sample2polrate:
                    if sample2polrate[s] <= median + std:
                        singles.add(s)
            else:
                print 'File not found %s'%pfn
                exit(1)

        df_list = []
        samples = []
        for ifn in args['ifn_metadatas']:
            index_col = get_index_col(ifn)
            df = pandas.read_csv(
                ifn,
                sep='\t', 
                dtype=unicode,
                header=0, 
                index_col=index_col)
            if not args['header_first_col']:
                df = df.transpose()
            df_list.append(df)
            samples += df.columns.values.tolist()

        print 'number of samples in metadata: %d'%len(samples)
        count = 0
        with open(ifn_tree, 'r') as ifile:
            line = ifile.readline()
        line = line.replace(args['string_to_remove'], '')
        tree = dendropy.Tree(stream=open(ifn_tree, 'r'), schema='newick')
        for node in tree.leaf_nodes():
            sample = node.get_node_str().strip("'")
            sample = sample.replace(args['string_to_remove'], '')
            prefixes = [prefix for prefix in 
                            ['k__', 'p__', 'c__', 'o__', 
                             'f__', 'g__', 's__'] \
                        if prefix in sample]

            metadata = sample
            if len(prefixes) == 0:
                count += 1
                for meta in args['metadatas']:
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
                                                str(df[sample][meta]).replace(':','_').lower())
                                break # take the first metadata

            if sample in sample2cluster:
                metadata += '|cluster-%d'%sample2cluster[sample]
            if args['add_multiple_strains_info'] and sample[:3] not in ['s__', 'g__'] and len(singles):
                if sample in singles:
                    metadata += '|single-strain'
                else:
                    metadata += '|multiple-strains'

            line = line.replace(sample + ':', metadata + ':')

        ofn_tree = ifn_tree + '.metadata'
        print 'Number of samples in tree: %d'%count
        print 'Output:', ofn_tree
        with open(ofn_tree, 'w') as ofile:
            ofile.write(line)

if __name__ == "__main__":
    args = read_params()
    main(args)
