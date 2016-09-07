#!/usr/bin/env python
#Author: Duy Tin Truong (duytin.truong@unitn.it)
#        at CIBIO, University of Trento, Italy

__author__  = 'Duy Tin Truong (duytin.truong@unitn.it)'
__version__ = '0.1'
__date__    = '7 Sep 2016'

import sys
import os
import argparse 
from ete2 import Tree, TreeStyle, NodeStyle

def read_params():
    p = argparse.ArgumentParser()
    p.add_argument(
        '--ifn', 
        required=True, 
        default=None, 
        type=str,
        help='The input tree file.')
    p.add_argument(
        '--ofn', 
        required=False, 
        default=None, 
        type=str,
        help='The input tree file.')


    return p.parse_args()


def main(args):
    if args.ofn == None:
        args.ofn = args.ifn + '.png'
    ts = TreeStyle()
    ts.show_leaf_name = True
    nstyle = NodeStyle()
    #nstyle["shape"] = "sphere"
    nstyle["size"] = 5
    #nstyle["fgcolor"] = "darkred"
    nstyle["vt_line_width"] = 2
    nstyle["hz_line_width"] = 2

    #ts.mode = "c"
    tree = Tree(args.ifn)
    for n in tree.traverse():
        n.set_style(nstyle)
    tree.render(args.ofn, tree_style=ts, dpi=300, units='in') #, h=20, w=20)


if __name__ == "__main__":
    args = read_params()
    main(args)
