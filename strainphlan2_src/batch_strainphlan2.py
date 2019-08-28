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

if sys.version_info[0] < 3:
    error("StrainPhlAn2 requires Python 3, your current Python version is {}.{}.{}"
                    .format(sys.version_info[0], sys.version_info[1], 
                        sys.version_info[2]), exit=True)

import argparse as ap
from strainphlan2 import strainphlan2
from utils import create_folder, check_clade
from extract_markers import extract_markers


"""
Reads and parses the command line arguments of the script.

:returns: the parsed arguments
"""
def read_params():
    p = ap.ArgumentParser(description="")
    p.add_argument('-d', '--database', type=str, default=None,
                   help="The input MetaPhlAn dtabase")
    p.add_argument('-m', '--clade_markers', type=str, default=None,
                   help="The input MetaPhlAn dtabase")
    p.add_argument('-s', '--samples', type=str, 
                   nargs='+', default=[],
                   help="The the markers for each sample")
    p.add_argument('-r', '--references', type=str, 
                   nargs='+', default=[],
                   help="The reference genomes")
    p.add_argument('-c', '--clades', type=str, 
                   nargs='+', default=[],
                   help="The clades to investigate")
    p.add_argument('-o', '--output_dir', type=str, default=None,
                   help="The output directory")
    p.add_argument('--metadata', type=str, default=None,
                   help="The metadata about references")
    p.add_argument('-n', '--nprocs', type=int, default=None,
                   help="The number of threads to use")
    
    return p.parse_args()


"""
Checks the mandatory command line arguments of the script.

:param args: the arguments of the script
:returns: the checked args
"""
def check_params(args):
    if not (args.database or args.clade_markers):
        error('-d (or --database) or -m (or --clade_markers) must be specified', 
            exit=True, init_new_line=True)
    elif args.database and args.clade_markers:
        error('-d (or --database) and -m (or --clade_markers) can '+
            'not be specified at same time', exit=True, 
            init_new_line=True)
    elif not args.samples:
        error('-s (or --samples) must be specified', exit=True, 
            init_new_line=True)
    elif not args.clades:
        error('-c (or --clades) must be specified', exit=True, 
            init_new_line=True)
    elif not args.output_dir:
        error('-o (or --output_dir) must be specified', exit=True, 
            init_new_line=True)
    elif args.database and not os.path.exists(args.database):
        error('The database does not exist', exit=True, 
            init_new_line=True)
    elif args.clade_markers and not os.path.exists(args.clade_markers):
        error('The clade markers file does not exist', exit=True, 
            init_new_line=True) 
    for s in args.samples:
        if not os.path.exists(s):
            error('The input sample file \"'+s+'\" does not exist', exit=True, 
                init_new_line=True)
    if len(args.samples)+len(args.references) < 4:
        error('The inputs samples + references are less than 4', exit=True, 
            init_new_line=True)
    if not args.output_dir.endswith('/'):
        args.output_dir += '/'
    if not args.nprocs:
        args.nprocs = 1

    return args


"""
Executes StrainPhlAn2 for many clades 

:param database: the MetaPhlan2 markers database
:param samples: the folder containing the markers generated with script samples_to_markers.py
:param clade: the clade to investigate
:param output_dir: the output directory
:param metadata: the metadata about references
:param nprocs: the threads used for execution
"""
def batch_strainphlan2(database, samples, clades, output_dir, metadata, nprocs):
    for c in clades:
        if check_clade(c):
            references = getReferences(c, metadata)
            create_folder(output_dir+c)
            clade_dir = output_dir+c+"/"
            clade_markers = extract_markers(database, c, clade_dir)
            create_folder(output_dir+c+"/results_base/")
            clade_dir = output_dir+c+"/results_base/"
            strainphlan2(None, clade_markers, samples, references, c, clade_dir,
                20, 80, nprocs)

            create_folder(output_dir+c+"/results_relaxed2/")
            clade_dir = output_dir+c+"/results_relaxed2/"
            strainphlan2(None, clade_markers, samples, references, c, clade_dir,
                15, 20, nprocs)

            create_folder(output_dir+c+"/results_relaxed3/")
            clade_dir = output_dir+c+"/results_relaxed3/"
            strainphlan2(None, clade_markers, samples, references, c, clade_dir,
                10, 10, nprocs)


"""
Get the reference files for a specific clade from a metadata file

:param clade: the clade to get the references
:param metadata: the file with the relation between clade and references
:returns: the reference files for the specific clade
"""
def getReferences(clade, metadata):
    csv_reader = open(metadata, "r")  
    references = list()  
    for row in csv_reader:
        row = row.replace("\n","").split(";")
        if row[1] == clade:
            if os.path.exists("/shares/CIBIO-Storage/CM/scratch/data/isolates/CM_carmen_caritro/assembly/"+row[0]+"/"+row[0]+".fasta"):
                references.append("/shares/CIBIO-Storage/CM/scratch/data/isolates/CM_carmen_caritro/assembly/"+row[0]+"/"+row[0]+".fasta")
            elif os.path.exists("/shares/CIBIO-Storage/CM/scratch/tmp_projects/smanara_milkisolates/"+row[0]+"/"+row[0]+".fasta"):
                references.append("/shares/CIBIO-Storage/CM/scratch/tmp_projects/smanara_milkisolates/"+row[0]+"/"+row[0]+".fasta")

    csv_reader.close()
    return references


"""
Main function

:param database: the MetaPhlan2 markers database
:param clade_markers: a FASTA containing the markers of a specific clade
:param samples: the folder containing the markers generated with script samples_to_markers.py
:param clade: the clade to investigate
:param output_dir: the output directory
:param metadata: the metadata about references
:param nprocs: the threads used for execution
"""
if __name__ == "__main__":
    info("Start StrainPhlAn2 batch execution")
    t0 = time.time()
    args = read_params()
    args = check_params(args)
    batch_strainphlan2(args.database, args.samples, args.clades, 
        args.output_dir, args.metadata, args.nprocs)
    exec_time = time.time() - t0
    info("Finish StrainPhlAn2 batch execution ("+str(round(exec_time, 2))+
        " seconds): Results are stored at \""+os.getcwd()+"/"+args.output_dir+"\"\n",
         init_new_line=True)