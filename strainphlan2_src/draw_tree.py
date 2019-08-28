#!/usr/bin/env python
#Authors: Duy Tin Truong (duytin.truong@unitn.it)
#        at CIBIO, University of Trento, Italy

__author__ = ('Duy Tin Truong (duytin.truong@unitn.it), '
              'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it)')
__version__ = '0.2'
__date__    = '10 Jul 19'

import argparse as ap
import dendropy
from io import StringIO
import re
from collections import defaultdict
import matplotlib.colors as colors
import subprocess

GRAPHLAN_PATH="/mnt/d/ScientificWork/CIBIO_repos/graphlan/"

"""
Reads and parses the command line arguments of the script.

:returns: the parsed arguments
"""
def read_params():
    p = ap.ArgumentParser()
    p.add_argument('--tree',
                   required=True,
                   default=None,
                   type=str,
                   help='The input tree in newick format.')
    p.add_argument('--leyend',
                    required=True,
                    default=None,
                    type=str,
                    help='The first lines of the metadata of GraPhlAn.')
    p.add_argument('--metadata',
                   required=True,
                   default=None,
                   type=str,
                   help='The last lines of the metadata of GraPhlAn.')
    return p.parse_args()


"""
Draws a tree using GraPhlAn

:param tree: the tree to draw
:param leyend: leyend of the tree
:param metadata: the metadata of the tree
"""
def draw_tree(args_tree, args_metadata, args_leyend):
    tree = dendropy.Tree.get_from_path(args_tree, schema='newick',
                                       preserve_underscores=True)
    # tree.reroot_at_midpoint()
    mrca = tree.mrca(taxon_labels=["Escherichia_fergusonii_ATCC_35469.CU928158.2"])
    # tree.prune_taxa_with_labels(["Salmonella_enterica_subsp__enterica_serovar_Typhimurium_str__LT2.NC_003197.2"])
    tree.reroot_at_node(mrca, update_bipartitions=False)

    count = 0
    for node in tree.preorder_node_iter():
        nodestr = node.__getattribute__("taxon").__str__().strip("'")
        if node.is_leaf():
            if '.' in nodestr:
                nodestr = nodestr.replace('.',',')
                node.taxon = dendropy.Taxon(label=nodestr)
        else:
            count += 1
            node.taxon = dendropy.Taxon(label='node_%d'%count)

    ofn_prefix = args_tree
    ofn_tree = ofn_prefix + '.graphlantree'
    tree.write_to_path(ofn_tree, 'newick')

    ofn_annot = ofn_prefix + '.annot'
    with open(ofn_annot, 'w') as ofile:
        with open(args_leyend) as file_handle:
            for line in file_handle:
                ofile.write(line)
        ofile.write("\n")
        for node in tree.preorder_node_iter():
            if not node.is_leaf():
                nodestr = node.__getattribute__("taxon").__str__().strip("'")
                ofile.write('%s\tclade_marker_size\t0\n'%(nodestr))

        with open(args_metadata) as file_handle:
            for line in file_handle:
                split = line.split("\t")
                if split[0] in str(tree.seed_node.leaf_nodes()):
                    ofile.write(line)
                    
    ofn_xml = ofn_prefix + '.xml'
    cmd = GRAPHLAN_PATH+'graphlan_annotate.py --annot %s %s %s'%(ofn_annot, ofn_tree, ofn_xml)
    subprocess.call(cmd.split())

    ofn_fig = ofn_prefix + '.png'
    cmd = GRAPHLAN_PATH+'graphlan.py %s %s --dpi %d --size %f'%(ofn_xml, ofn_fig, 300, 8.0)
    subprocess.call(cmd.split())

    print ('Output file: %s'%ofn_fig)



"""
Main function

:param tree: the tree to draw
:param leyend: leyend of the tree
:param metadata: the metadata of the tree
"""
if __name__ == "__main__":
    args = read_params()
    draw_tree(args.tree, args.metadata, args.leyend)
