#!/usr/bin/env python

__author__ = ('Duy Tin Truong (duytin.truong@unitn.it), '
              'Francesco Asnicar (f.asnicar@unitn.it), '
              'Moreno Zolfo (moreno.zolfo@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it), '
              'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it)')
__version__ = '2.0.0'
__date__ = '17 Jul 2019'

import pickle, bz2, os, time
import subprocess as sb
import argparse as ap
from Bio import SeqIO, Seq, SeqRecord
from external_exec import execute, compose_command
from utils import error, info


"""
Reads and parses the command line arguments of the script.

:returns: the parsed arguments
"""
def read_params():
    p = ap.ArgumentParser(description="")
    p.add_argument('-d', '--database', type=str, default=None,
                   help="The input MetaPhlAn dtabase")
    p.add_argument('-c', '--clade', type=str, default=None,
                   help="The clades to investigate")
    p.add_argument('-o', '--output_dir', type=str, default=None,
                   help="The output directory")
    
    return p.parse_args()


"""
Checks the mandatory command line arguments of the script.

:returns: the checked args
"""
def check_params(args):
    if not args.database:
        error('-d (or --database) must be specified', exit=True, 
            init_new_line=True)
    elif not args.clade:
        error('-c (or --clade) must be specified', exit=True, 
            init_new_line=True)
    elif not args.output_dir:
        error('-o (or --output_dir) must be specified', exit=True, 
            init_new_line=True)
    elif not os.path.exists(args.database):
        error('The database does not exist', exit=True, 
            init_new_line=True)
    if not args.output_dir.endswith('/'):
        args.output_dir += '/'
    
    return args


"""
Checks the mandatory programs to execute of the script.

"""
def check_dependencies():
    try:
        sb.check_call("bowtie2-inspect", stdout=sb.DEVNULL, stderr=sb.DEVNULL)
    except Exception as e:
        error('Program not installed or not present in the system path\n'+str(e), 
            init_new_line=True, exit=True)


"""
Extract the markers of a specific clade in a MetaPhlAn database

:param database: the MetaPhlan markers database
:param clade: the clade to extract markers
:param output_dir: the output directory
:returns: the output file with the extracted sequences of the marker
"""
def extract_markers(database, clade, output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir, 755)

    info('\tGenerating DB markers FASTA...', init_new_line=True)
    fasta_markers = generate_markers_fasta(database, output_dir)
    info('\tDone.', init_new_line=True)

    info('\tLoading MetaPhlan database...', init_new_line=True)
    db = pickle.load(bz2.BZ2File(database))
    info('\tDone.',init_new_line=True)
    markers = set([])
    for marker in db['markers']:
        if clade == db['markers'][marker]['taxon'].split('|')[-1]:
            markers.add(marker)
    info('\tNumber of markers for the clade \"'+clade+"\": "+str(len(markers)), init_new_line=True)
    output_file = output_dir+clade+".fna"
    info('\tExporting markers...', init_new_line=True)
    with open(output_file, 'w') as ofile:
        for rec in SeqIO.parse(open(fasta_markers, 'r'), 'fasta'):
            if rec.name in markers:
                SeqIO.write(rec, ofile, 'fasta')
    info('\tDone.', init_new_line=True)
    
    os.remove(fasta_markers)
    return output_file


"""
Generates a FASTA file with the markers form a MetaPhlAn database

:param database: the MetaPhlan markers database
:param output_dir: the output directory
"""
def generate_markers_fasta(database, output_dir):
    db_markers_faa = output_dir+"db_markers.fna"
    bowtie_database, _ = os.path.splitext(database)
    params = {
        "program_name" : "bowtie2-inspect",
        "command_line" : "#program_name# #input# > #output#"
    }
    execute(compose_command(params, input_file=bowtie_database, output_file=db_markers_faa))
    return db_markers_faa


"""
Main function

:param database: the MetaPhlan markers database
:param clade: the clade to extract markers
:param output_dir: the output directory
"""
if __name__ == "__main__":
    info("Start execution: "+format(time.ctime(int(time.time()))))
    args = read_params()
    check_dependencies()
    args = check_params(args)
    extract_markers(args.database, args.clade, args.output_dir)
    info("Finish execution: "+format(time.ctime(int(time.time())))+"\n", init_new_line=True)
