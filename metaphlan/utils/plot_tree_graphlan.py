#!/usr/bin/env python
__author__ = ('Duy Tin Truong (duytin.truong@unitn.it), '
              'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it)')
__version__ = '3.0'
__date__    = '21 Feb 2020'

import argparse as ap
import dendropy
from io import StringIO
import re
from collections import defaultdict
import matplotlib.colors as colors
import subprocess


def read_params():
    p = ap.ArgumentParser()
    p.add_argument('-t', '--ifn_tree',
                   required=True,
                   default=None,
                   type=str,
                   help='The input tree in newick format.')
    p.add_argument('-m', '--colorized_metadata',
                   required=False,
                   default='unset',
                   type=str,
                   help='The metadata field to colorize. Default "unset".')
    p.add_argument('--fig_size',
                   required=False,
                   default=8,
                   type=float,
                   help='The figure size. Default "8".')
    p.add_argument('--legend_marker_size',
                   required=False,
                   default=20,
                   type=int,
                   help='The legend marker size. Default "20".'
                   )
    p.add_argument('--legend_font_size',
                   required=False,
                   default=10,
                   type=int,
                   help='The legend font size. Default "10".'
                   )
    p.add_argument('--legend_marker_edge_width',
                   required=False,
                   default=0.2,
                   type=float,
                   help='The legend marker edge width. Default "0.2".'
                   )
    p.add_argument('--leaf_marker_size',
                   required=False,
                   default=20,
                   type=int,
                   help='The legend marker size. Default "20".'
                   )
    p.add_argument('--leaf_marker_edge_width',
                   required=False,
                   default=0.2,
                   type=float,
                   help='The legend marker edge width. Default "0.2".'
                   )
    p.add_argument('--dpi',
                   required=False,
                   default=300,
                   type=int,
                   help='The figure dpi.')
    p.add_argument('--figure_extension',
                   required=False,
                   default='.png',
                   type=str,
                   help='The figure extension. Default ".png".')
    p.add_argument('--ofn_prefix',
                   required=False,
                   default=None,
                   type=str,
                   help='The prefix of output files.')
    return p.parse_args()

def run(cmd):
    print (cmd)
    subprocess.call(cmd.split())

def main():
    args = read_params()
    tree = dendropy.Tree.get_from_path(args.ifn_tree, schema='newick',
                                       preserve_underscores=True)
    tree.reroot_at_midpoint()
    count = 0
    metadatas = set([])
    node2metadata = {}
    for node in tree.preorder_node_iter():
        nodestr = node.__getattribute__("taxon").__str__().strip("'")
        if node.is_leaf():
            if '.' in nodestr:
                nodestr = nodestr.replace('.',',')
                node.taxon = dendropy.Taxon(label=nodestr)
            substrs = re.findall(
                         '%s-[a-zA-Z0-9.]*'%args.colorized_metadata,
                          nodestr)
            if substrs:
                md = substrs[0].replace(args.colorized_metadata + '-', '')
                metadatas.add(md)
                node2metadata[nodestr] = md
        else:
            count += 1
            node.taxon = dendropy.Taxon(label='node_%d'%count)
    metadatas = sorted(list(metadatas))
    color_names = list(colors.cnames.keys())
    metadata2color = {}
    for i, md in enumerate(metadatas):
        metadata2color[md] = color_names[i % len(color_names)]

    if not args.ofn_prefix:
        args.ofn_prefix = args.ifn_tree
    ofn_tree = args.ofn_prefix + '.graphlantree'
    tree.write_to_path(ofn_tree, 'newick')
    ofn_annot = args.ofn_prefix + '.annot'
    with open(ofn_annot, 'w') as ofile:
        #ofile.write('clade_separation\t0\n')
        ofile.write('branch_bracket_width\t0\n')
        #ofile.write('clade_separation\t0.15\n')
        ofile.write('branch_bracket_depth\t0\n')
        #ofile.write('branch_thickness\t1.25\n')
        ofile.write('annotation_background_width\t0\n')

        # legend
        ofile.write('#legends\n')
        ofile.write('class_legend_font_size\t%d\n'%args.legend_font_size)

        for md in metadata2color:
            ofile.write('%s\tclade_marker_size\t%d\n'%(md, args.legend_marker_size))
            ofile.write('%s\tclade_marker_color\t%s\n'%(md, metadata2color[md]))
            ofile.write('%s\tclade_marker_edge_width\t%f\n'%(md, args.legend_marker_edge_width))

        # remove intermedate nodes
        for node in tree.preorder_node_iter():
            if not node.is_leaf():
                nodestr = node.__getattribute__("taxon").__str__().strip("'")
                ofile.write('%s\tclade_marker_size\t0\n'%(nodestr))

        # colorize leaf nodes
        for node in tree.seed_node.leaf_nodes():
            nodestr = node.__getattribute__("taxon").__str__().strip("'")
            if nodestr in node2metadata:
                leaf_color = metadata2color[node2metadata[nodestr]]
                ofile.write('%s\tclade_marker_size\t%d\n'%(nodestr, args.leaf_marker_size))
                ofile.write('%s\tclade_marker_color\t%s\n'%(nodestr, leaf_color))
                ofile.write('%s\tclade_marker_edge_width\t%f\n'%(nodestr, args.leaf_marker_edge_width))

    ofn_xml = args.ofn_prefix + '.xml'
    cmd = 'graphlan_annotate.py --annot %s %s %s'%(ofn_annot, ofn_tree, ofn_xml)
    run(cmd)

    ofn_fig = args.ofn_prefix + args.figure_extension
    cmd = 'graphlan.py %s %s --dpi %d --size %f'%(ofn_xml, ofn_fig, args.dpi, args.fig_size)
    run(cmd)

    print ('Output file: %s'%ofn_fig)

if __name__ == '__main__':
    main()
