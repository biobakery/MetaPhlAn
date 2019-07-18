#!/usr/bin/env python
#Author: Duy Tin Truong (duytin.truong@unitn.it)
#        at CIBIO, University of Trento, Italy

__author__ = ('Duy Tin Truong (duytin.truong@unitn.it), '
              'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it)')
__version__ = '0.2'
__date__    = '10 Jul 19'


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
    nstyle["size"] = 5
    nstyle["vt_line_width"] = 2
    nstyle["hz_line_width"] = 2

    tree = Tree(args.ifn)
    for n in tree.traverse():
        n.set_style(nstyle)
    tree.render(args.ofn, tree_style=ts, dpi=300, units='in')


if __name__ == "__main__":
    args = read_params()
    main(args)
