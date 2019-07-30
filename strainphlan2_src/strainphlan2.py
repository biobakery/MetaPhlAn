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
from shutil import copyfile, rmtree
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from external_exec import decompress_bz2, create_blastn_db, execute_blastn
from external_exec import generate_phylophlan_config_file, create_phylophlan_db, execute_phylophlan
from extract_markers import extract_markers
from parallelisation import execute_pool
from utils import optimized_dump, create_folder, check_clade, parse_marker_name, get_breath

# Regular expression to remove comments: \n\"\"\"[^"]+\n\"\"\"

"""
Global variables
"""
MARKER_IN_SAMPLES_THRESHOLD = 80
MARKERS_IN_SAMPLE_THRESHOLD = 20
TRIM_SEQUENCES = 50


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
    p.add_argument('-c', '--clade', type=str, default=None,
                   help="The clade to investigate")
    p.add_argument('-o', '--output_dir', type=str, default=None,
                   help="The output directory")
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
    elif not args.references:
        error('-r (or --references) must be specified', exit=True, 
            init_new_line=True)
    elif not args.clade:
        error('-c (or --clade) must be specified', exit=True, 
            init_new_line=True)
    elif not check_clade(args.clade):
        error('The introduced clade is not at species level', exit=True, 
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
    for s in args.references:
        if not os.path.exists(s):
            error('The reference file \"'+s+'\" does not exist', exit=True, 
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
Gets a binary matrix representing the presence/ausence of the clade 
markers in the uploaded samples

:param database: the MetaPhlan2 markers database
:param clade_markers_file: a FASTA containing the markers of a specific clade
:param samples: the folder containing the markers generated with script samples_to_markers.py
:param clade: the clade to investigate
:param tmp_dir: the temporal output directory
:param nprocs: the threads used for execution
:returns: the clade markers FASTA file and the markers matrix
"""
def get_markers_matrix(database, clade_markers_file, samples, clade, tmp_dir, nprocs):
    if not clade_markers_file:
        clade_markers_file = extract_markers(database, clade, tmp_dir)
    else:
        n, _ = os.path.splitext(os.path.basename(clade_markers_file))
        copyfile(clade_markers_file, tmp_dir+n+".fna")
        clade_markers_file = tmp_dir+n+".fna"

    clade_markers = []
    for rec in SeqIO.parse(open(clade_markers_file, 'r'), 'fasta'):
        clade_markers.append(rec.id)
    
    markers_matrix = execute_pool(((get_matrix_for_sample, s, clade_markers) for s in samples), 
        nprocs)

    return markers_matrix, clade_markers_file


"""
Gets a matrix with the presence / ausence of the clade markers in 
a sample
:param s the path to the sample
:param clade_markers: a list with the names of the clade markers
:returns: the markers matrix for the sample
"""
def get_matrix_for_sample(s, clade_markers):
    sample = pickle.load(open(s, "rb"))
    markers = {"sample": s}
    for m in clade_markers:
        markers.update({m : 0})
    for r in sample:
        if r['marker'] in clade_markers:
            markers.update({r['marker'] : 1})
    return markers


"""
Remove bad markers and samples from the markers matrix:
First, checks if the percentage of markers of a sample sample reaches a
threshold, if not, removes the sample.
Then, checks if the percentage of samples that contain a marker reaches 
a threhold, if not, removes the marker.

:param markers_matrix: the markers matrix
:returns: the filtered markers matrix
"""
def clean_markers_matrix(markers_matrix):    
    # Checks if the percentage of markers of a sample sample reachs a threshold, 
    # if not, removes the sample    
    cleaned_markers_matrix = []
    to_remove = []

    for m in markers_matrix:
        if sum(list(m.values())[1:]) < MARKERS_IN_SAMPLE_THRESHOLD:
            to_remove.append(m)
        else:
            cleaned_markers_matrix.append({'sample': m['sample']})

    for r in to_remove:
        markers_matrix.remove(r)

    # Checks how many samples were deleted
    if len(cleaned_markers_matrix) < 4:
        error("Phylogeny can not be inferred. Too many samples were discarded", 
            exit=True, init_new_line=True) 

    # Checks if the percentage of samples that contain a marker reachs a threhold,
    # if not, removes the marker
    markers = list(markers_matrix[0])[1:]
    for m in markers:
        marker_scores = []
        for s in markers_matrix:
            marker_scores.append(s[m])
        if (sum(marker_scores) * 100) / len(marker_scores) >= MARKER_IN_SAMPLES_THRESHOLD:            
            counter = 0
            for s in markers_matrix:
                cleaned_markers_matrix[counter].update({m : s[m]})
                counter += 1

    # Checks how many markers were deleted
    if sum(list(cleaned_markers_matrix[0].values())[1:]) < MARKERS_IN_SAMPLE_THRESHOLD:
        error("Phylogeny can not be inferred. Too many markers were discarded", 
            exit=True, init_new_line=True) 

    return cleaned_markers_matrix


"""
For each sample, writes the FASTA files with the sequences of the filtered markers

:input cleaned_markers_matrix: a list with the filtered markers
:param samples: the folder containing the markers generated with script samples_to_markers.py
:param references: the FASTA reference files
:param tmp_dir: the output temporal directory
:param nprocs: the threads used for execution
:returns: the folder of the FASTA files
"""
def matrix_markers_to_fasta(cleaned_markers_matrix, samples, references, tmp_dir, nprocs):     
    tmp_dir=tmp_dir+"samples_as_markers/"    
    create_folder(tmp_dir)

    filtered_samples = []
    filtered_names = []
    for s in cleaned_markers_matrix:
        filtered_names.append(os.path.splitext(os.path.basename(s['sample']))[0])
        filtered_samples.append(s['sample'])

    for r in references:
        r_name = os.path.splitext(os.path.basename(r))[0]            
        if r_name in filtered_names:
            copyfile(r, tmp_dir+r_name+".fna")

    execute_pool(((sample_markers_to_fasta, s, filtered_samples, tmp_dir, 
        list(cleaned_markers_matrix[0])) for s in samples), 
        nprocs)
    return tmp_dir


"""
Writes a FASTA file with the filtered clade markers of a sample

:param s: the path to the sample
:param filtered_samples: a list with the filtered samples
:param tmp_dir: the temporal output directory
:param filtered_clade_markers: a list with the filtered clade markers
"""
def sample_markers_to_fasta(s, filtered_samples, tmp_dir, filtered_clade_markers):
    if s in filtered_samples:
        sample_name = os.path.splitext(os.path.basename(s))[0]
        with open(tmp_dir+sample_name+'.fna', 'w') as marker_fna:
            sample = pickle.load(open(s, "rb"))
            for r in sample:
                if r['marker'] in filtered_clade_markers:
                    marker_name = parse_marker_name(r['marker'])
                    seq = SeqRecord(Seq(r['sequence'][TRIM_SEQUENCES:-TRIM_SEQUENCES].replace("*","N"), generic_dna), 
                        id=marker_name, description=marker_name)
                    SeqIO.write(seq, marker_fna, 'fasta')


"""
Writes a FASTA file with the sequences of the filtered clade markers

:param markers: a list with the names of filtered markers
:param tmp_dir: the output temporal directory
:param clade: the threads used for execution
:param clade_markers_file: the FASTA with the clade markers
"""
def cleaned_clade_markers_to_fasta(markers, tmp_dir, clade, clade_markers_file):
    tmp_dir=tmp_dir+clade+"/"
    create_folder(tmp_dir)

    clade_markers = {}
    for rec in SeqIO.parse(open(clade_markers_file, 'r'), 'fasta'):
        clade_markers.update({rec.id: rec.seq})

    for m in list(markers[0])[1:]:
        marker_name = parse_marker_name(m)
        with open(tmp_dir+marker_name+'.fna', 'w') as marker_fna:
            seq = SeqRecord(clade_markers.get(m), id=marker_name, description=marker_name)
            SeqIO.write(seq, marker_fna, 'fasta')


"""
Gets markers from reference files and returns the marker matrix with the
reference markers

:param tmp_dir: the temporal output directory
:param clade_markers_file: the clade markers FASTA file
:param markers_matrix: the markers matrix
:param references: the list of the reference files
:param nprocs: the threads using in the BLASTn executions
:returns: the marker matrix with references
"""
def get_markers_from_references(tmp_dir, clade_markers_file, markers_matrix, references, nprocs): 
    blastn_dir = tmp_dir+"blastn/"
    create_folder(blastn_dir)
    clade_markers = list(markers_matrix[0])[1:]

    results = execute_pool(((process_reference, s, blastn_dir, 
        clade_markers_file, clade_markers) for s in references), 
        nprocs)
    for r in results:
        markers_matrix.append(r)
        
    return markers_matrix 


"""
Processes each reference file and get a markers dictionary to add 
to the markers matrix

:param r: the path to the reference file
:param blastn_dir: the temporal blastn output directory
:param clade_markers_file: the clade markers FASTA file
:param markers_matrix: the markers matrix
:returns: the markers of the reference file as a dictionary
"""
def process_reference(r, blastn_dir, clade_markers_file, clade_markers):
    n, _ = os.path.splitext(os.path.basename(r))               
    _, f = os.path.splitext(r)
    if f == ".bz2":
        r = decompress_bz2(r, blastn_dir)
    blastn_db = create_blastn_db(blastn_dir, r)
    blastn_file = execute_blastn(blastn_dir, clade_markers_file, blastn_db)
    return parse_blastn_results(blastn_dir+n+'.pkl', 
        clade_markers, blastn_file, r)


"""
Parses BLASTn results and gets the presence of the clade markers in
the reference file

:param sample: the name of the reference file
:param clade_markers: a list of the clade_markers_name
:param blastn_file: the BLASTn output file
:param reference: the reference FASTA file
:returns: A dictionary with the presence of the clade markers
"""
def parse_blastn_results(sample, clade_markers, blastn_file, reference):
    markers = {"sample": sample}
    for m in clade_markers:
        markers.update({m : 0})

    reference_sequences = {}
    for rec in SeqIO.parse(open(reference, 'r'), 'fasta'):
        reference_sequences.update({rec.id : rec.seq})

    blastn_result = open(blastn_file, "r")    
    processed_markers = []
    for line in blastn_result:
        if line == '':
            break
        query = line.split("\t")[0]
        if not query in processed_markers:
            processed_markers.append(query)
            markers.update({query : 1})   
    blastn_result.close()
    return markers


"""
Executes PhyloPhlAn2 to compute phylogeny

:param samples_markers_dir: the temporal samples markers directory
:param num_samples: the number of filtered samples
:param tmp_dir: the temporal output directory
:param output_dir: the output_directory
:param clade: the clade
:param nproc: the number of threads to run phylophlan
"""
def compute_phylogeny(samples_markers_dir, num_samples, tmp_dir, output_dir, clade, nprocs):    
    info("\tCreating PhyloPhlAn2 database...", init_new_line=True)
    create_phylophlan_db(tmp_dir, clade)
    info("\tDone.", init_new_line=True)
    info("\tGenerating PhyloPhlAn2 configuration file...", init_new_line=True)
    conf_file = generate_phylophlan_config_file(tmp_dir)
    info("\tDone.", init_new_line=True)
    info("\tProcessing samples...", init_new_line=True)
    min_entries = int(MARKER_IN_SAMPLES_THRESHOLD*num_samples/100)
    execute_phylophlan(samples_markers_dir, conf_file, min_entries, tmp_dir, 
        output_dir, clade, nprocs)
    info("\tDone.", init_new_line=True)


"""
Executes StrainPhlAn2 

:param database: the MetaPhlan2 markers database
:param clade_markers: a FASTA containing the markers of a specific clade
:param samples: the folder containing the markers generated with script samples_to_markers.py
:param references: the reference genomes
:param clade: the clade to investigate
:param output_dir: the output directory
:param nprocs: the threads used for execution
"""
def strainphlan2(database, clade_markers, samples, references, clade, output_dir, nprocs):
    info("Creating temporal directory...", init_new_line=True)
    tmp_dir = output_dir+'tmp/'
    create_folder(tmp_dir)
    info("Done.", init_new_line=True)
    info("Getting markers from sample files...", init_new_line=True)
    markers_matrix, clade_markers_file = get_markers_matrix(database, clade_markers, 
        samples, clade, tmp_dir, nprocs)
    info("Done.", init_new_line=True)    
    info("Getting markers from reference files...", init_new_line=True)
    markers_matrix = get_markers_from_references(tmp_dir, clade_markers_file, 
        markers_matrix, references, nprocs)
    info("Done.", init_new_line=True)    
    info("Removing bad markers / samples...", init_new_line=True)
    cleaned_markers_matrix = clean_markers_matrix(markers_matrix)
    info("Done.", init_new_line=True)    
    info("Writting samples as markers' FASTA files...", init_new_line=True)
    samples_as_markers_dir = matrix_markers_to_fasta(cleaned_markers_matrix, 
        samples, references, tmp_dir, nprocs)
    info("Done.", init_new_line=True)
    info("Writting filtered clade markers as FASTA file...", init_new_line=True)
    cleaned_clade_markers_to_fasta(cleaned_markers_matrix, tmp_dir, clade, clade_markers_file)
    info("Done.", init_new_line=True)
    info("Executing PhyloPhlAn2...", init_new_line=True)
    compute_phylogeny(samples_as_markers_dir, len(cleaned_markers_matrix), tmp_dir, 
        output_dir, clade, nprocs)
    info("Done.", init_new_line=True)    
    info("Removing temporal files...", init_new_line=True)
    rmtree(tmp_dir, ignore_errors=False, onerror=None)
    info("Done.", init_new_line=True)


"""
Main function

:param database: the MetaPhlan2 markers database
:param clade_markers: a FASTA containing the markers of a specific clade
:param samples: the folder containing the markers generated with script samples_to_markers.py
:param references: the reference genomes
:param clade: the clade to investigate
:param output_dir: the output directory
:param nprocs: the threads used for execution
"""
if __name__ == "__main__":
    info("Start StrainPhlAn2 execution")
    t0 = time.time()
    args = read_params()
    args = check_params(args)
    strainphlan2(args.database, args.clade_markers, args.samples, args.references, 
        args.clade, args.output_dir, args.nprocs)
    exec_time = time.time() - t0
    info("Finish StrainPhlAn2 execution ("+str(round(exec_time, 2))+
        " seconds): Results are stored at \""+os.getcwd()+"/"+args.output_dir+"\"\n",
         init_new_line=True)