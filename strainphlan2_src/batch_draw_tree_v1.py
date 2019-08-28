#!/usr/bin/env python

__author__ = ('Duy Tin Truong (duytin.truong@unitn.it), '
              'Francesco Asnicar (f.asnicar@unitn.it), '
              'Moreno Zolfo (moreno.zolfo@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it), '
              'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it)')
__version__ = '2.0.0'
__date__ = '29 Jul 2019'

import os, sys, pickle, time, shutil
from utils import info, error

if sys.version_info[0] >= 3:
    error("StrainPhlAn2 requires Python 2, your current Python version is {}.{}.{}"
                    .format(sys.version_info[0], sys.version_info[1], 
                        sys.version_info[2]), exit=True)

import argparse as ap
from utils import check_clade
from draw_tree import draw_tree


"""
Reads and parses the command line arguments of the script.

:returns: the parsed arguments
"""
def read_params():
    p = ap.ArgumentParser(description="")
    p.add_argument('-c', '--clades', type=str, 
                   nargs='+', default=[],
                   help="The clades to investigate")
    p.add_argument('-o', '--output_dir', type=str, default=None,
                   help="The output directory")
    p.add_argument('-a', '--metadata', type=str, default=None,
                   help="The metadata to draw the tree")
    
    return p.parse_args()


"""
Checks the mandatory command line arguments of the script.

:param args: the arguments of the script
:returns: the checked args
"""
def check_params(args):
    if not args.clades:
        error('-c (or --clades) must be specified', exit=True, 
            init_new_line=True)
    elif not args.output_dir:
        error('-o (or --output_dir) must be specified', exit=True, 
            init_new_line=True)
    elif not args.metadata:
        error('-a (or --metadata) must be specified', exit=True, 
            init_new_line=True)
    if not args.output_dir.endswith('/'):
        args.output_dir += '/'

    return args


"""
Executes draw_tree for many clades 

:param database: the MetaPhlan2 markers database
:param samples: the folder containing the markers generated with script samples_to_markers.py
:param references: the reference genomes
:param clade: the clade to investigate
:param output_dir: the output directory
:param metadata: the metadata to draw the tree
:param nprocs: the threads used for execution
"""
def batch_draw_tree(clades, output_dir, metadata):
    for c in clades:
        if check_clade(c):
            clade_dir = output_dir+c+"/results_base/"
            if os.path.isfile(clade_dir+"RAxML_bestTree."+c+".tree"):
                draw_tree(clade_dir+"RAxML_bestTree."+c+".tree", metadata)
            clade_dir = output_dir+c+"/results_relaxed2/"
            if os.path.isfile(clade_dir+"RAxML_bestTree."+c+".tree"):    
                draw_tree(clade_dir+"RAxML_bestTree."+c+".tree", metadata)
            clade_dir = output_dir+c+"/results_relaxed3/"
            if os.path.isfile(clade_dir+"RAxML_bestTree."+c+".tree"):
                draw_tree(clade_dir+"RAxML_bestTree."+c+".tree", metadata)


"""
Main function

:param clade: the clade to investigate
:param output_dir: the output directory
:param metadata: the metadata to draw the tree
"""
if __name__ == "__main__":
    info("Start draw_tree batch execution")
    t0 = time.time()
    args = read_params()
    args = check_params(args)
    batch_draw_tree(args.clades, args.output_dir, args.metadata)
    exec_time = time.time() - t0
    info("Finish draw_tree batch execution ("+str(round(exec_time, 2))+
        " seconds): Results are stored at \""+os.getcwd()+"/"+args.output_dir+"\"\n",
         init_new_line=True)
