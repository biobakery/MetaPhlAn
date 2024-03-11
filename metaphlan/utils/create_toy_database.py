#!/usr/bin/env python
__author__ = 'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it)'
__version__ = '4.1.1'
__date__ = '11 Mar 2024'


import bz2, os, pickle, time
from Bio import SeqIO
try:
    from .util_fun import info, error, warning
    from .external_exec import build_bowtie2_db, compress_bz2, generate_markers_fasta
except ImportError:
    from util_fun import info, error, warning
    from external_exec import build_bowtie2_db, compress_bz2, generate_markers_fasta
import argparse as ap


metaphlan_script_install_folder = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DB_FOLDER = os.path.join(metaphlan_script_install_folder, '..',  "metaphlan_databases")
DEFAULT_DB_FOLDER= os.environ.get('METAPHLAN_DB_DIR', DEFAULT_DB_FOLDER)
INDEX = 'latest'

def read_params():
    """ Reads and parses the command line arguments of the script

    Returns:
        namespace: The populated namespace with the command line arguments
    """
    p = ap.ArgumentParser(formatter_class=ap.RawTextHelpFormatter, add_help=False)   
    requiredNamed = p.add_argument_group('requiered arguments')
    requiredNamed.add_argument('--in_sgbs', type=str, default=None, help="The path to the file containing the SGBs to keep")
    requiredNamed.add_argument('--out_dir', type=str, default=None, help="The output folder for the new MetaPhlAn database")
    requiredNamed.add_argument('--out_name', type=str, default=None, help="The name for the new MetaPhlAn database")   
    optionalNamed = p.add_argument_group('optional arguments')
    optionalNamed.add_argument('--bowtie2db', metavar="METAPHLAN_BOWTIE2_DB", type=str, default=DEFAULT_DB_FOLDER,
        help=("Folder containing the MetaPhlAn database. You can specify the location by exporting the DEFAULT_DB_FOLDER variable in the shell."))
    optionalNamed.add_argument('--index', type=str, default=INDEX,
        help=("Specify the id of the database version to use. If \"latest\", MetaPhlAn will get the latest version."))
    optionalNamed.add_argument("-h", "--help", action="help", help="Show this help message and exit")
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
    if not args.out_dir:
        error('--out_dir must be specified', exit=True)
    elif not os.path.exists(args.out_dir):
        error('The file {} does not exist'.format(
            args.out_dir), exit=True)   
    if not args.out_name:
        error('--out_name must be specified', exit=True)


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


def resolve_database(bowtie2db, index):
    if index == 'latest':
        if os.path.exists(os.path.join(bowtie2db, 'mpa_latest')):
            with open(os.path.join(bowtie2db, 'mpa_latest'), 'r') as mpa_latest:
                return '{}/{}'.format(bowtie2db, [line.strip() for line in mpa_latest if not line.startswith('#')][0])
        else:
            error('The default MetaPhlAn database cannot be found at: {}'.format(
                os.path.join(bowtie2db, 'mpa_latest')), exit=True)
    else:
        if os.path.exists('{}/{}.pkl'.format(bowtie2db, index)):
            return '{}/{}'.format(bowtie2db, index)
        else:
            error('The default MetaPhlAn database cannot be found at: {}'.format(
                os.path.join(bowtie2db, index + '.pkl')), exit=True)
            
    
def create_toy_database(input_sgbs, index, bowtie2db, output_folder, output_name):
    info('Loading input PKL database {}...'.format(index))
    input_database = resolve_database(bowtie2db, index) 
    input_pkl = '{}.pkl'.format(input_database)
    mpa_pkl = pickle.load(bz2.open(input_pkl,'rb'))
    sgbs_in_db = get_sgbs_in_db(mpa_pkl['taxonomy'])
    info('Done.')
    
    info('Checking input SGBs...')
    to_keep_sgbs = get_to_keep_sgbs(input_sgbs, sgbs_in_db)
    info('Done.')
    
    info('Filtering input PKL...')
    new_mpa_pkl = {}
    new_mpa_pkl['taxonomy'], new_mpa_pkl['markers'], new_mpa_pkl['merged_taxon'] = {}, {}, {}
    for taxa in mpa_pkl['taxonomy']:
        if taxa.split('|')[-1] in to_keep_sgbs:
            new_mpa_pkl['taxonomy'][taxa] = mpa_pkl['taxonomy'][taxa]
    for marker in mpa_pkl['markers']:
        if mpa_pkl['markers'][marker]['clade'] in to_keep_sgbs:
            ext = []
            for x in mpa_pkl['markers'][marker]['ext']:
                if 't__' + x in to_keep_sgbs:
                    ext.append(x)
            mpa_pkl['markers'][marker]['ext'] = ext
            new_mpa_pkl['markers'][marker] = mpa_pkl['markers'][marker] 
    for merge in mpa_pkl['merged_taxon']:
        if merge[0].split('|')[-1] in to_keep_sgbs:
            new_mpa_pkl['merged_taxon'][merge] = mpa_pkl['merged_taxon'][merge]      
    info('Done.')
    
    info('Writting output PKL database {}...'.format(output_name))
    output_pkl = os.path.join(output_folder, output_name + '.pkl')
    with bz2.BZ2File(output_pkl, 'w') as write_file:
        pickle.dump(new_mpa_pkl, write_file, protocol=2)
    info('Done.')
    
    info('Filtering input FNA database {}...'.format(index))
    input_fna = '{}.fna'.format(input_database)
    remove = False
    if os.path.exists(input_fna):
        handle = open(input_fna, 'r')
    elif os.path.exists(input_fna + '.bz2'):
        handle = bz2.open(input_fna + '.bz2', 'rt')
    else:
        info('Extracting database FASTA file from the Bowtie2 database...')
        for bt2file in ['.1','.2','.3','.4','.rev.1','.rev.2']:
            if not os.path.exists('{}{}.bt2l'.format(input_database, bt2file)):
                error('The MetaPhlAn {} Bowtie2 database cannot (totally or partially) be found at: {}'.format(index, bowtie2db), exit=True)
        input_fna = generate_markers_fasta(input_pkl, output_folder)
        remove = True
        info('Done.')
        handle = open(input_fna, 'r') 
    filtered_sequences = []
    for record in SeqIO.parse(handle, "fasta"):
        if record.id in new_mpa_pkl['markers']:
            filtered_sequences.append(record)            
    handle.close()
    if remove:
        os.remove(input_fna)
    info('Done.')
    
    info('Writting output FNA database {}...'.format(output_name))
    output_fna = os.path.join(output_folder, output_name + '.fna')
    with open(output_fna, 'w') as handle:
        SeqIO.write(filtered_sequences, handle, "fasta")
    info('Done.')
    
    info('Building Bowtie2 database')
    output_bowtie = os.path.join(output_folder, output_name)
    build_bowtie2_db(output_fna, output_bowtie)
    compress_bz2(output_fna)
    info('Done')
    
    
def main():
    t0 = time.time()
    args = read_params()
    info("Start generating toy database")
    check_params(args)
    create_toy_database(args.in_sgbs, args.index, args.bowtie2db, args.out_dir, args.out_name)
    exec_time = time.time() - t0
    info("Finish generating toy database ({} seconds)".format(round(exec_time, 2)))


if __name__ == '__main__':
    main()
