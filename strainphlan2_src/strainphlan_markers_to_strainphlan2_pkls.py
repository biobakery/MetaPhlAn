#!/usr/bin/env python

__author__ = ('Duy Tin Truong (duytin.truong@unitn.it), '
              'Francesco Asnicar (f.asnicar@unitn.it), '
              'Moreno Zolfo (moreno.zolfo@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it), '
              'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it)')
__version__ = '2.0.0'
__date__ = '29 Jul 2019'

import msgpack, os, time
import argparse as ap
from utils import info, error, optimized_dump, get_breath
from samples_to_markers import BREATH_THRESHOLD


"""
Reads and parses the command line arguments of the script.

:returns: the parsed arguments
"""
def read_params():
    p = ap.ArgumentParser(description="")
    p.add_argument('-m', '--markers', type=str, 
                   nargs='+', default=[],
                   help="The the markers for each sample")
    p.add_argument('-o', '--output_dir', type=str, default=None,
                   help="The output directory")
    
    return p.parse_args()


"""
Checks the mandatory command line arguments of the script.

:param args: the arguments of the script
:returns: the checked args
"""
def check_params(args):
    if not args.markers:
        error('-m (or --markers) must be specified', exit=True, 
            init_new_line=True)
    elif not args.output_dir:
        error('-o (or --output_dir) must be specified', exit=True, 
            init_new_line=True)
    if not args.output_dir.endswith('/'):
        args.output_dir += '/'
    return args 


"""
Converts Strainphlan consensus markers to Strainphlan2 markers Pickle files

:param markers: the Strainphlan consensus marker files
:param output_dir: the output directory
"""
def strainphlan_markers_to_strainphlan2_pkls(markers, output_dir):
    for s in markers:        
        n, _ = os.path.splitext(os.path.basename(s))
        with open(s, 'rb') as ifile:
            consensus = []
            marker2seq = msgpack.unpack(ifile, use_list=False, encoding='utf-8')
            for marker in marker2seq:                
                seq = marker2seq[marker]['seq']
                breath = get_breath(seq)
                if(breath > BREATH_THRESHOLD):
                    consensus.append({"marker":marker, "breath":breath, "sequence":seq})
        markers_pkl = open(output_dir+n+'.pkl', 'wb')
        optimized_dump(markers_pkl, consensus)  


"""
Main function

:param samples: the Strainphlan consensus marker files
:param output_dir: the output directory
"""
if __name__ == "__main__":
    info("Start StrainPhlAn markers to StrainPhlAn2 Pickle markers execution")
    t0 = time.time()
    args = read_params()
    args = check_params(args)
    strainphlan_markers_to_strainphlan2_pkls(args.samples, args.output_dir)
    exec_time = time.time() - t0
    info("Finish StrainPhlAn markers to StrainPhlAn2 Pickle markers execution ("+str(round(exec_time, 2))+
        " seconds): Results are stored at \""+os.getcwd()+"/"+args.output_dir+"\"\n",
         init_new_line=True)