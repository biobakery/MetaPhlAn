#!/usr/bin/env python
__author__ = 'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it)'
__version__ = '4.1.0'
__date__ = '23 Aug 2023'


import bz2, os, pickle, time
from Bio import SeqIO
try:
    from .util_fun import info, error, warning
except ImportError:
    from util_fun import info, error, warning
import argparse as ap

def read_params():
    """ Reads and parses the command line arguments of the script

    Returns:
        namespace: The populated namespace with the command line arguments
    """
    p = ap.ArgumentParser(
        description="", formatter_class=ap.ArgumentDefaultsHelpFormatter)    
    p.add_argument('--in_sgbs', type=str, default=None, help="The path to the file containing the SGBs to keep")
    p.add_argument('--in_pkl', type=str, default=None, help="The input MetaPhlAn PKL database")
    p.add_argument('--in_fna', type=str, default=None, help="The input MetaPhlAn FNA database")
    p.add_argument('--out_pkl', type=str, default=None, help="The output MetaPhlAn PKL database")
    p.add_argument('--out_fna', type=str, default=None, help="The output MetaPhlAn FNA database")
    return p.parse_args()


def check_params(args):
    """Checks the mandatory command line arguments of the script

    Args:
        args (namespace): the arguments to check
    """
    if not args.in_sgbs:
        error('--in_sgbs must be specified', exit=True)
    elif not os.path.exists(args.in_sgbs):
        error('The file {} does not exist'.format(
            args.in_sgbs), exit=True)        
    if not args.in_pkl:
        error('--in_pkl must be specified', exit=True)
    elif not os.path.exists(args.in_pkl):
        error('The file {} does not exist'.format(
            args.in_pkl), exit=True)        
    if not args.in_fna:
        error('--in_fna must be specified', exit=True)
    elif not os.path.exists(args.in_fna):
        error('The file {} does not exist'.format(
            args.in_fna), exit=True)        


def get_to_keep_sgbs(input_file, sgbs_in_db):
    to_keep_sgbs = set()
    with open(input_file, 'r') as rf:
        for line in rf:
            sgb = line.strip()
            if sgb.startswith('t__SGB'):
                if sgb in sgbs_in_db:
                    to_keep_sgbs.add(sgb)
                else:
                    warning('The input SGB {} is not in the database'.format(sgb), exit=False)
            else:
                error('The file containing the SGBs to be kept is not properly formatted!!', exit=True)
    return to_keep_sgbs


def get_sgbs_in_db(taxonomy):
    sgbs_in_db = set()
    for taxa in taxonomy:
        sgbs_in_db.add(taxa.split('|')[-1])
    return sgbs_in_db


def create_toy_database(input_sgbs, input_pkl, input_fna, output_pkl, output_fna):
    info('Loading input PKL database {}...'.format(input_pkl.split('/')[-1]))
    mpa_pkl = pickle.load(bz2.open(input_pkl,'rb'))
    info('Done.')
    sgbs_in_db = get_sgbs_in_db(mpa_pkl['taxonomy'])
    info('Checking input SGBs...')
    to_keep_sgbs = get_to_keep_sgbs(input_sgbs, sgbs_in_db)
    info('Done.')
    info('Filtering input PKL...')
    filtered_taxa = dict()
    filtered_markers = dict()
    filtered_merged = dict()
    for taxa in mpa_pkl['taxonomy']:
        if taxa.split('|')[-1] in to_keep_sgbs:
            filtered_taxa[taxa] = mpa_pkl['taxonomy'][taxa]
    for marker in mpa_pkl['markers']:
        if mpa_pkl['markers'][marker]['clade'] in to_keep_sgbs:
            ext = list()
            for x in mpa_pkl['markers'][marker]['ext']:
                if 't__' + x in to_keep_sgbs:
                    ext.append(x)
            mpa_pkl['markers'][marker]['ext'] = ext
            filtered_markers[marker] = mpa_pkl['markers'][marker] 
    for merge in mpa_pkl['merged_taxon']:
        if merge[0].split('|')[-1] in to_keep_sgbs:
            filtered_merged[merge] = mpa_pkl['merged_taxon'][merge]      
    
    new_mpa_pkl = dict()
    new_mpa_pkl['taxonomy'] = filtered_taxa
    new_mpa_pkl['markers'] = filtered_markers
    new_mpa_pkl['merged_taxon'] = filtered_merged
    info('Done.')
    info('Writting output PKL database {}...'.format(output_pkl.split('/')[-1]))
    with bz2.BZ2File(output_pkl, 'w') as write_file:
        pickle.dump(new_mpa_pkl, write_file, protocol=2)
    info('Done.')
    info('Filtering input FNA database {}...'.format(input_fna))
    filtered_sequences = list()
    if input_fna.endswith('.bz2'):
        handle = bz2.open(input_fna, 'rt')
    else:
        handle = open(input_fna, 'r')    
    for record in SeqIO.parse(handle, "fasta"):
        if record.id in filtered_markers:
            filtered_sequences.append(record)            
    handle.close()
    info('Done.')
    info('Writting output FNA database {}...'.format(output_fna))
    if output_fna.endswith('.bz2'):
        handle = bz2.open(output_fna, 'wt')
    else:
        handle = open(output_fna, 'w')
    SeqIO.write(filtered_sequences, handle, "fasta")
    handle.close()
    info('Done.')

def main():
    t0 = time.time()
    args = read_params()
    info("Start generating toy database")
    check_params(args)
    create_toy_database(args.in_sgbs, args.in_pkl, args.in_fna, args.out_pkl, args.out_fna)
    exec_time = time.time() - t0
    info("Finish generating toy database ({} seconds)".format(round(exec_time, 2)))


if __name__ == '__main__':
    main()
