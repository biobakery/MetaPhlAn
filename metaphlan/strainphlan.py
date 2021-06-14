#!/usr/bin/env python
__author__ = ('Aitor Blanco Miguez (aitor.blancomiguez@unitn.it), '
              'Duy Tin Truong (duytin.truong@unitn.it), '
              'Francesco Asnicar (f.asnicar@unitn.it), '
              'Moreno Zolfo (moreno.zolfo@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)')
__version__ = '3.0.10'
__date__ = '14 Jun 2021'


import sys
try:
    from .utils import *
except ImportError:
    from utils import *
if sys.version_info[0] < 3:
    error("StrainPhlAn " + __version__ + " requires Python 3, your current Python version is {}.{}.{}"
                    .format(sys.version_info[0], sys.version_info[1], 
                        sys.version_info[2]), exit=True)

import os, pickle, time, bz2, numpy, collections, tempfile
import argparse as ap
from shutil import copyfile, rmtree, move
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

DEFAULT_DATABASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
    "metaphlan_databases/mpa_v30_CHOCOPhlAn_201901.pkl")
PHYLOPHLAN_MODES = ['accurate', 'fast']

# Regular expression to remove comments: \n\"\"\"[^"]+\n\"\"\"

"""
Reads and parses the command line arguments of the script.

:returns: the parsed arguments
"""
def read_params():
    p = ap.ArgumentParser(description="", formatter_class=ap.ArgumentDefaultsHelpFormatter)
    p.add_argument('-d', '--database', type=str, default=DEFAULT_DATABASE,
                   help="The input MetaPhlAn " + __version__ + " database")
    p.add_argument('-m', '--clade_markers', type=str, default=None,
                   help="The clade markers as FASTA file")
    p.add_argument('-s', '--samples', type=str, 
                   nargs='+', default=[],
                   help="The reconstructed markers for each sample")
    p.add_argument('-r', '--references', type=str, 
                   nargs='+', default=[],
                   help="The reference genomes")
    p.add_argument('-c', '--clade', type=str, default=None,
                   help="The clade to investigate")
    p.add_argument('-o', '--output_dir', type=str, default=None,
                   help="The output directory")
    p.add_argument('-n', '--nprocs', type=int, default=1,
                   help="The number of threads to use")
    p.add_argument('--secondary_samples', type=str, 
                    nargs='+', default=[],
                    help="The reconstructed markers for each secondary sample")
    p.add_argument('--secondary_references', type=str, 
                    nargs='+', default=[],
                    help="The secondary reference genomes")
    p.add_argument('--trim_sequences', type=int, default=50,
                    help="The number of bases to remove from both ends when trimming markers")
    p.add_argument('--marker_in_n_samples', type=int, default=80,
                    help="Theshold defining the minimum percentage of samples to keep a marker")
    p.add_argument('--sample_with_n_markers', type=int, default=20,
                    help="Threshold defining the minimun number of markers to keep a sample")
    p.add_argument('--secondary_sample_with_n_markers', type=int, default=20,
                    help="Threshold defining the minimun number of markers to keep a secondary sample")
    p.add_argument('--phylophlan_mode', choices=PHYLOPHLAN_MODES, default='fast',
                    help="The presets for fast or accurate phylogenetic analysis")
    p.add_argument('--phylophlan_configuration', type=str, default=None,
                    help="The PhyloPhlAn configuration file")
    p.add_argument('--mutation_rates', action='store_true', default=False,
                   help=("If specified, StrainPhlAn will produce a mutation rates table for each of the aligned markers and a summary table "
                         "for the concatenated MSA. This operation can take long time to finish"))
    p.add_argument('--print_clades_only', action='store_true', default=False,
                   help=("If specified, StrainPhlAn will only print the potential clades and stop the execution"))
    p.add_argument('--debug', action='store_true', default=False,
                   help=("If specified, StrainPhlAn will not remove the temporal folders"))
    
    return p.parse_args()


"""
Checks the mandatory command line arguments of the script.

:param args: the arguments of the script
:returns: the checked args
"""
def check_params(args):
    if args.print_clades_only and args.clade_markers:
        error('-m (or --clade_markers) cannot be specified together with --print_clades_only', 
            exit=True, init_new_line=True)
    elif not args.samples:
        error('-s (or --samples) must be specified', exit=True, 
            init_new_line=True)
    elif not args.print_clades_only and not args.clade:
        error('-c (or --clade) must be specified', exit=True, 
            init_new_line=True)
    elif not args.print_clades_only and not args.output_dir:
        error('-o (or --output_dir) must be specified', exit=True, 
            init_new_line=True)
    elif args.database and not os.path.exists(args.database):
        error('The database does not exist', exit=True, 
            init_new_line=True)
    elif args.clade_markers and not os.path.exists(args.clade_markers):
        error('The clade markers file does not exist', exit=True, 
            init_new_line=True)
    elif args.phylophlan_configuration and not os.path.exists(args.phylophlan_configuration):
        error('The phylophlan configuration file does not exist', exit=True, 
            init_new_line=True)
    for s in args.samples:
        if not os.path.exists(s):
            error('The input sample file \"'+s+'\" does not exist', exit=True, 
                init_new_line=True)
    for s in args.references:
        if not os.path.exists(s):
            error('The reference file \"'+s+'\" does not exist', exit=True, 
                init_new_line=True)
    for s in args.secondary_samples:
        if not os.path.exists(s):
            error('The secondary input sample file \"'+s+'\" does not exist', exit=True, 
                init_new_line=True)
    for s in args.secondary_references:
        if not os.path.exists(s):
            error('The secondary reference file \"'+s+'\" does not exist', exit=True, 
                init_new_line=True)
    if len(args.samples)+len(args.references) < 4:
        error('The main inputs samples + references are less than 4', exit=True, 
            init_new_line=True)
    if not args.print_clades_only and not args.output_dir.endswith('/'):
        args.output_dir += '/'

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
Adds secondary samples to the marker matrix

:param secondary_samples: the paths to the secondary samples
:param cleaned_markers_matrix: the markers matrix of the main samples and references
:param secondary_samples_with_n_markers: threshold defining the minimun number of markers to keep 
    a secondary sample
:param nprocs: number of threads to use
:returns: the filtered markers matrix with the secondary samples
"""
def add_secondary_samples(secondary_samples, cleaned_markers_matrix, 
    secondary_samples_with_n_markers, nprocs):
    clade_markers = list(cleaned_markers_matrix[0].keys())[1:]
    markers_matrix = execute_pool(((get_matrix_for_sample, s, clade_markers) for s in secondary_samples), 
        nprocs)
    for m in markers_matrix:
        if sum(list(m.values())[1:]) >= secondary_samples_with_n_markers:
            cleaned_markers_matrix.append(m)  

    return cleaned_markers_matrix

    
"""
Adds secondary references to the marker matrix

:param secondary_references: the paths to the secondary references
:param cleaned_markers_matrix: the markers matrix of the main samples and references
:param secondary_samples_with_n_markers: threshold defining the minimun number of markers to keep 
    a secondary reference
:param clade_markers_file: a FASTA containing the markers of a specific clade
:param tmp_dir: the temporal output directory
:param nprocs: the number of threads to use
:returns: the filtered markers matrix with the secondary references
"""
def add_secondary_references(secondary_references, cleaned_markers_matrix, 
    secondary_samples_with_n_markers, clade_markers_file, tmp_dir, nprocs):
    clade_markers = list(cleaned_markers_matrix[0])[1:]
    markers_matrix = execute_pool(((process_reference, s, tmp_dir+"blastn/", 
        clade_markers_file, clade_markers) for s in secondary_references), 
        nprocs)
    for m in markers_matrix:
        if sum(list(m.values())[1:]) >= secondary_samples_with_n_markers:
            cleaned_markers_matrix.append(m)    

    return cleaned_markers_matrix 


"""
Adds secondary samples and references to the marker matrix

:param secondary_samples: the paths to the secondary samples
:param secondary_references: the paths to the secondary references
:param cleaned_markers_matrix: the markers matrix of the main samples and references
:param secondary_samples_with_n_markers: threshold defining the minimun number 
    of markers to keep a secondary reference
:param clade_markers_file: a FASTA containing the markers of a specific clade
:param tmp_dir: the temporal output directory
:param nprocs: the number of threads to use
:returns: the filtered markers matrix with the secondary samples and references
"""
def add_secondary_samples_and_references(secondary_samples, secondary_references, 
    cleaned_markers_matrix, secondary_samples_with_n_markers, clade_markers_file, 
    tmp_dir, nprocs):    
    if len(secondary_samples) > 0:
        info("Getting markers from secondary sample files...", init_new_line=True)
        cleaned_markers_matrix = add_secondary_samples(secondary_samples, cleaned_markers_matrix, 
        secondary_samples_with_n_markers, nprocs)
        info("Done.", init_new_line=True)     
    if len(secondary_references) > 0:
        info("Getting markers from secondary reference files...", init_new_line=True)
        cleaned_markers_matrix = add_secondary_references(secondary_references, cleaned_markers_matrix, 
        secondary_samples_with_n_markers, clade_markers_file, tmp_dir, nprocs)
        info("Done.", init_new_line=True) 
    return cleaned_markers_matrix


"""
Gets a matrix with the presence / ausence of the clade markers in 
a sample
:param sample_path the path to the sample
:param clade_markers: a list with the names of the clade markers
:returns: the markers matrix for the sample
"""
def get_matrix_for_sample(sample_path, clade_markers):
    if os.path.splitext(sample_path)[1] == ".bz2":
        sample = pickle.load(bz2.BZ2File(sample_path))
    else:
        sample = pickle.load(open(sample_path, "rb"))
    markers = {"sample": sample_path}
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
:param samples_with_n_markers: threshold defining the minimun number of markers 
    to keep a sample
:param marker_in_n_samples: threshold defining the minimum percentage of 
    samples to keep a marker
:returns: the filtered markers matrix
"""
def clean_markers_matrix(markers_matrix, samples_with_n_markers, 
    marker_in_n_samples, messages = True):    
    # Checks if the percentage of markers of a sample sample reachs a threshold, 
    # if not, removes the sample    
    cleaned_markers_matrix = []
    to_remove = []
    for m in markers_matrix:
        if sum(list(m.values())[1:]) < samples_with_n_markers:
            to_remove.append(m)
        else:
            cleaned_markers_matrix.append({'sample': m['sample']})    
    for r in to_remove:
        markers_matrix.remove(r)

    # Checks how many samples were deleted
    if len(cleaned_markers_matrix) < 4:
        if messages:
            error("Phylogeny can not be inferred. Too many samples were discarded", 
                exit=True, init_new_line=True) 
        else:
            return []

    # Checks if the percentage of samples that contain a marker reachs a threhold,
    # if not, removes the marker
    markers = list(markers_matrix[0])[1:]
    for m in markers:
        marker_scores = []
        for s in markers_matrix:
            marker_scores.append(s[m])
        if (sum(marker_scores) * 100) / len(marker_scores) >= marker_in_n_samples:            
            counter = 0
            for s in markers_matrix:
                cleaned_markers_matrix[counter].update({m : s[m]})
                counter += 1

    # Checks how many markers were deleted
    if len(list(cleaned_markers_matrix[0].values())[1:]) < samples_with_n_markers:        
        if messages:
            error("Phylogeny can not be inferred. Too many markers were discarded", 
                exit=True, init_new_line=True) 
        else:
            return []

    # Checks again if the percentage of markers of a sample sample is at least 1, 
    # if not, removes the sample  
    to_remove = []
    for m in cleaned_markers_matrix:
        if sum(list(m.values())[1:]) < 1:
            to_remove.append(m)
    for r in to_remove:
        cleaned_markers_matrix.remove(r)

    # Checks again how many samples were deleted
    if len(cleaned_markers_matrix) < 4:        
        if messages:
            error("Phylogeny can not be inferred. No enough markers were kept for the samples", 
                exit=True, init_new_line=True)
        else:
            return []

    return cleaned_markers_matrix


"""
For each sample, writes the FASTA files with the sequences of the filtered markers

:param cleaned_markers_matrix: a list with the filtered markers
:param clade: the name of the clade to investigate
:param samples: the folder containing the markers generated with script samples_to_markers.py
:param references: the FASTA reference files
:param trim_sequences: the number of bases to remove when trimming markers
:param tmp_dir: the output temporal directory
:param nprocs: the threads used for execution
:returns: the folder of the FASTA files
"""
def matrix_markers_to_fasta(cleaned_markers_matrix, clade, samples, references, 
    trim_sequences, tmp_dir, nprocs):     
    tmp_dir=tmp_dir+clade+".StrainPhlAn3/"    
    create_folder(tmp_dir)

    filtered_samples = []
    filtered_names = []
    for s in cleaned_markers_matrix:
        filtered_names.append(os.path.splitext(os.path.basename(s['sample']))[0])
        filtered_samples.append(s['sample'])

    for r in references:
        _, f = os.path.splitext(r)
        if f == ".bz2":
            r = decompress_bz2(r, tmp_dir)
        else:
            r_name = os.path.splitext(os.path.basename(r))[0]            
            if r_name in filtered_names:
                copyfile(r, tmp_dir+r_name+".fna")

    execute_pool(((sample_markers_to_fasta, s, filtered_samples, tmp_dir, 
        list(cleaned_markers_matrix[0]), trim_sequences) for s in samples), 
        nprocs)
    return tmp_dir


"""
Writes a FASTA file with the filtered clade markers of a sample

:param sample_path: the path to the sample
:param filtered_samples: a list with the filtered samples
:param tmp_dir: the temporal output directory
:param filtered_clade_markers: a list with the filtered clade markers
:param trim_sequences: the number of bases to remove when trimming markers
"""
def sample_markers_to_fasta(sample_path, filtered_samples, tmp_dir, filtered_clade_markers, trim_sequences):
    if sample_path in filtered_samples:
        sample_name = os.path.splitext(os.path.basename(sample_path))[0].replace(".pkl","")
        with open(tmp_dir+sample_name+'.fna', 'w') as marker_fna:
            if os.path.splitext(sample_path)[1] == ".bz2":
                sample = pickle.load(bz2.BZ2File(sample_path))
            else:
                sample = pickle.load(open(sample_path, "rb"))
            for r in sample:
                if r['marker'] in filtered_clade_markers:
                    marker_name = parse_marker_name(r['marker'])
                    seq = SeqRecord(Seq(r['sequence'][trim_sequences:-trim_sequences].replace("*","N").replace("-","N").strip('N')), id=marker_name, description=marker_name)
                    SeqIO.write(seq, marker_fna, 'fasta')


"""
Writes a FASTA file with the sequences of the filtered clade markers

:param markers: a list with the names of filtered markers
:param tmp_dir: the output temporal directory
:param clade: the threads used for execution
:param clade_markers_file: the FASTA with the clade markers
"""
def cleaned_clade_markers_to_fasta(markers, tmp_dir, clade, clade_markers_file):
    tmp_dir=tmp_dir+clade[:30]+"/"
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
def get_markers_from_references(tmp_dir, clade_markers_file, markers_matrix, 
    references, nprocs): 
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
def process_reference(reference, blastn_dir, clade_markers_file, clade_markers):
    n, _ = os.path.splitext(os.path.basename(reference))               
    _, f = os.path.splitext(reference)
    if f == ".bz2":
        reference = decompress_bz2(reference, blastn_dir)
    blastn_db = create_blastn_db(blastn_dir, reference)
    blastn_file = execute_blastn(blastn_dir, clade_markers_file, blastn_db)
    return parse_blastn_results(blastn_dir+n+'.pkl', 
        clade_markers, blastn_file, reference)


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
Gets PhyloPhlAn configuration

:param phylophlan_mode: the precision of the phylogenetic analysis
:returns: the configuration to create a PhyloPhlAn configuration file
"""
def get_phylophlan_configuration():
    configuration = dict()
    # blastn, tblastn, diamond
    configuration.update({'map':'blastn'})
    # muscle, mafft, opal, upp
    configuration.update({'aligner':'mafft'})
    # trimal
    configuration.update({'trim':'trimal'})
    # fasttree, raxml, iqtree, astral, astrid
    configuration.update({'tree1':'raxml'})
    # configuration.update({'tree2':''})
    
    return configuration


"""
Executes PhyloPhlAn2 to compute phylogeny

:param samples_markers_dir: the temporal samples markers directory
:param num_samples: the number of filtered samples
:param tmp_dir: the temporal output directory
:param output_dir: the output_directory
:param clade: the clade
:param marker_in_n_samples: threshold defining the minimum percentage of samples to keep a marker
:param phylophlan_mode: the precision of the phylogenetic analysis
:param phylophlan_configuration: the PhyloPhlAn configuration file
:param mutation_rates: whether get  the mutation rates for the markers
:param nproc: the number of threads to run phylophlan
"""
def compute_phylogeny(samples_markers_dir, num_samples, tmp_dir, output_dir, clade, 
    marker_in_n_samples, phylophlan_mode, phylophlan_configuration, mutation_rates, nprocs):    
    info("\tCreating PhyloPhlAn 3.0 database...", init_new_line=True)
    create_phylophlan_db(tmp_dir, clade[:30])
    info("\tDone.", init_new_line=True)
    if not phylophlan_configuration:     
        info("\tGenerating PhyloPhlAn 3.0 configuration file...", init_new_line=True)
        conf = get_phylophlan_configuration()
        phylophlan_configuration = generate_phylophlan_config_file(tmp_dir, conf)
        info("\tDone.", init_new_line=True)   
    info("\tProcessing samples...", init_new_line=True)
    min_entries = int(marker_in_n_samples*num_samples/100)
    execute_phylophlan(samples_markers_dir, phylophlan_configuration, min_entries,
        tmp_dir, output_dir, clade, phylophlan_mode, mutation_rates, nprocs)
    if mutation_rates:
        move(output_dir+"mutation_rates.tsv",output_dir+clade+".mutation")
        move(output_dir+"mutation_rates",output_dir+clade+"_mutation_rates")
    info("\tDone.", init_new_line=True)


"""
Generates a file with the polimorfic rates of the species for each sample

:param samples: the folder containing the markers to secondary samples generated with 
    script samples_to_markers.py
:param clade_markers_file: a FASTA containing the markers of a specific clade
:param clade: the clade to investigate
:param output_dir: the output directory
"""
def calculate_polimorfic_rates(samples, clade_markers_file, clade, output_dir):
    clade_markers = []
    for rec in SeqIO.parse(open(clade_markers_file, 'r'), 'fasta'):
        clade_markers.append(rec.id)
    with open(output_dir+clade+".polymorphic", 'w') as polimorfic_file:
        polimorfic_file.write("sample\tpercentage_of_polymorphic_sites\tavg_by_marker\tmedian_by_marker" +
            "\tstd_by_marker\tmin_by_marker\tmax_by_marker\tq25_by_marker\tq75_by_marker")
        for sample_path in samples:
            if os.path.splitext(sample_path)[1] == ".bz2":
                sample = pickle.load(bz2.BZ2File(sample_path))
            else:
                sample = pickle.load(open(sample_path, "rb"))
            p_stats, p_count, m_len = list(), 0, 0
            for m in sample:    
                if(m['marker'] in clade_markers):        
                    p_count += m['sequence'].count('*')
                    m_len += len(m['sequence'])
                    p_stats.append(m['sequence'].count('*')*100/len(m['sequence']))
            if not m_len == 0:
                polimorfic_file.write("\n"+os.path.splitext(os.path.basename(sample_path))[0] +
                    "\t"+str(p_count*100/m_len)+"\t"+str(numpy.average(p_stats)) +
                    "\t"+str(numpy.percentile(p_stats,50)) + "\t"+str(numpy.std(p_stats)) +
                    "\t"+str(numpy.min(p_stats)) + "\t"+str(numpy.max(p_stats)) +
                    "\t"+str(numpy.percentile(p_stats,25)) + "\t"+str(numpy.percentile(p_stats,75)))
    return len(clade_markers)


"""
Writes the information file for the execution

:param cleaned_markers_matrix: a list with the filtered markers
:param num_markers_for_clade: the number of markers in available for the clade
:param clade: the clade to investigate
:param output_dir: the output directory
:param n_p_samples: the main samples
:param n_s_samples: the secondary samples
:param n_p_references: the main references
:param n_s_references: the secondary references
:param trim_sequences: the number of bases to remove when trimming markers
:param samples_with_n_markers: threshold defining the minimun number of markers 
    to keep a main sample
:param marker_in_n_samples: threshold defining the minimum percentage of samples
    to keep a marker
:param secondary_samples_with_n_markers: threshold defining the minimun number of markers 
    to keep a secondary sample
:param phylophlan_mode: the precision of the phylogenetic analysis
:param nprocs: the threads used for execution
"""
def write_info(cleaned_markers_matrix, num_markers_for_clade, clade, output_dir, p_samples, 
    s_samples, p_references, s_references, trim_sequences, samples_with_n_markers, 
    marker_in_n_samples, secondary_samples_with_n_markers, phylophlan_mode, nprocs):
    f_p_samples, f_s_samples, f_p_references, f_s_references = 0, 0, 0, 0
    for c in cleaned_markers_matrix:
        if c['sample'] in p_samples:
            f_p_samples += 1
        elif c['sample'] in s_samples:
            f_s_samples += 1
        else:
            for r in p_references:
                if os.path.splitext(os.path.basename(c['sample']))[0] in r:
                    f_p_references += 1
                    break            
            for r in s_references:
                if os.path.splitext(os.path.basename(c['sample']))[0] in r:
                    f_s_references += 1
                    break
    with open(output_dir+clade+".info", 'w') as info_file:
        info_file.write("Clade: "+clade+"\nNumber of main samples: "+str(len(p_samples)) +
            "\nNumber of secondary samples: " + str(len(s_samples)) +
            "\nNumber of main references: " + str(len(p_references)) +            
            "\nNumber of secondary references: " + str(len(s_references)) +
            "\nNumber of available markers for the clade: " + 
            str(num_markers_for_clade)+"\nFiltering parameters: " +
            "\n\tNumber of bases to remove when trimming markers: "+ str(trim_sequences) +
            "\n\tMinimun number of markers to keep a main sample: "+ str(samples_with_n_markers) +
            "\n\tMinimun number of markers to keep a secondary sample: " + 
            str(secondary_samples_with_n_markers) +
            "\n\tMinimum percentage of samples to keep a marker: "+ str(marker_in_n_samples) +
            "\nNumber of markers selected after filtering: "+str(len(cleaned_markers_matrix[0].keys())-1) +
            "\nNumber of main samples after filtering: "+str(f_p_samples) +
            "\nNumber of secondary samples after filtering: " + str(f_s_samples) +
            "\nNumber of main references after filtering: " + str(f_p_references) +            
            "\nNumber of secondary references after filtering: " + str(f_s_references) +
            "\nPhyloPhlan phylogenetic precision mode: "+ phylophlan_mode +
            "\nNumber of processes used: "+ str(nprocs) + "\n" )


"""
Prints the clades detected in the reconstructed markers

:param database: the MetaPhlan2 markers database
:param samples: the folder containing the markers to main samples generated with script 
    samples_to_markers.py
:param samples_with_n_markers: threshold defining the minimun number of markers 
    to keep a main sample
:param marker_in_n_samples: threshold defining the minimum percentage of samples
    to keep a marker
"""
def print_clades(database, samples, samples_with_n_markers, marker_in_n_samples):
    sample_id = 0
    markers2species = dict()
    species2markers = dict()
    species_markers_matrix = dict()
    species2samples = dict()

    info('Loading MetaPhlan '+ __version__ + ' database...', init_new_line=True)
    db = pickle.load(bz2.BZ2File(database))
    for marker in db['markers']:
        species = db['markers'][marker]['clade']
        markers2species[marker] = species
        if not species in species2markers:
            species2markers[species] = list()
        species2markers[species].append(marker)
    info('Done.',init_new_line=True)

    info('Detecting clades...', init_new_line=True)     
    for sample_path in samples:
        if os.path.splitext(sample_path)[1] == ".bz2":
            sample = pickle.load(bz2.BZ2File(sample_path))
        else:
            sample = pickle.load(open(sample_path, "rb"))
        for r in sample:
            if r['marker'] in markers2species:
                species = markers2species[r['marker']]
                if species not in species_markers_matrix:
                    species_markers_matrix[species] = list()
                    for sp in samples:
                        markers = {"sample": sp}
                        for m in species2markers[species]:
                            markers.update({m : 0})
                        species_markers_matrix[species].append(markers)
                species_markers_matrix[species][sample_id][r['marker']] = 1
        sample_id += 1  
    for species in species_markers_matrix:    
        cleaned_markers_matrix = clean_markers_matrix(species_markers_matrix[species], samples_with_n_markers, marker_in_n_samples, False)
        if len(cleaned_markers_matrix) >= 4:
            species2samples[species] = len(cleaned_markers_matrix)
    info('Done.',init_new_line=True)

    info('Detected clades: ', init_new_line=True)
    sorted_species2samples = collections.OrderedDict(sorted(species2samples.items(), key=lambda kv: kv[1], reverse=True))
    for species in sorted_species2samples:
        info('\t' + species + ': in ' + str(sorted_species2samples[species]) + ' samples.', init_new_line=True)
    info('Done.',init_new_line=True) 


"""
Executes StrainPhlAn 

:param database: the MetaPhlan2 markers database
:param clade_markers: a FASTA containing the markers of a specific clade
:param samples: the folder containing the markers to main samples generated with script 
    samples_to_markers.py
:param references: the reference genomes
:param secondary_samples: the folder containing the markers to secondary samples generated with 
    script samples_to_markers.py
:param secondary_references: the secondary reference genomes
:param clade: the clade to investigate
:param output_dir: the output directory
:param trim_sequences: the number of bases to remove when trimming markers
:param samples_with_n_markers: threshold defining the minimun number of markers 
    to keep a main sample
:param marker_in_n_samples: threshold defining the minimum percentage of samples
    to keep a marker
:param secondary_samples_with_n_markers: threshold defining the minimun number of markers 
    to keep a secondary sample
:param phylophlan_mode: the precision of the phylogenetic analysis
:param phylophlan_configuration: the PhyloPhlAn configuration file
:param mutation_rates: whether get  the mutation rates for the markers
:param print_clades_only: whether print only the potential clades and stop
:param debug: wether to save the tmp folders
:param nprocs: the threads used for execution
"""
def strainphlan(database, clade_markers, samples, references, secondary_samples, 
    secondary_references, clade, output_dir, trim_sequences, samples_with_n_markers, 
    marker_in_n_samples, secondary_samples_with_n_markers, phylophlan_mode, 
    phylophlan_configuration, mutation_rates, print_clades_only, debug, nprocs):
    if print_clades_only:
        print_clades(database, samples, samples_with_n_markers, marker_in_n_samples)
    else:
        info("Creating temporary directory...", init_new_line=True)
        tmp_dir = tempfile.mkdtemp(dir=output_dir) + "/"
        info("Done.", init_new_line=True)
        info("Getting markers from main sample files...", init_new_line=True)
        markers_matrix, clade_markers_file = get_markers_matrix(database, clade_markers, 
            samples, clade, tmp_dir, nprocs)
        info("Done.", init_new_line=True)    
        info("Getting markers from main reference files...", init_new_line=True)
        markers_matrix = get_markers_from_references(tmp_dir, clade_markers_file, 
            markers_matrix, references, nprocs)
        info("Done.", init_new_line=True)    
        info("Removing bad markers / samples...", init_new_line=True)
        cleaned_markers_matrix = clean_markers_matrix(markers_matrix, samples_with_n_markers, 
            marker_in_n_samples)
        info("Done.", init_new_line=True)
        cleaned_markers_matrix = add_secondary_samples_and_references(secondary_samples, 
            secondary_references, cleaned_markers_matrix, secondary_samples_with_n_markers, 
            clade_markers_file, tmp_dir, nprocs)
        info("Writing samples as markers' FASTA files...", init_new_line=True)
        samples_as_markers_dir = matrix_markers_to_fasta(cleaned_markers_matrix, clade,
            samples+secondary_samples, references+secondary_references, trim_sequences, 
            tmp_dir, nprocs)
        info("Done.", init_new_line=True)
        info("Writing filtered clade markers as FASTA file...", init_new_line=True)
        cleaned_clade_markers_to_fasta(cleaned_markers_matrix, tmp_dir, clade, clade_markers_file)
        info("Done.", init_new_line=True)
        info("Calculating polymorphic rates...", init_new_line=True)
        num_markers_for_clade = calculate_polimorfic_rates(samples+secondary_samples, clade_markers_file, 
            clade, output_dir)
        info("Done.", init_new_line=True)   
        info("Executing PhyloPhlAn 3.0...", init_new_line=True)
        compute_phylogeny(samples_as_markers_dir, len(cleaned_markers_matrix), tmp_dir, 
            output_dir, clade, marker_in_n_samples, phylophlan_mode, phylophlan_configuration,
            mutation_rates, nprocs)
        info("Done.", init_new_line=True)     
        info("Writing information file...", init_new_line=True)
        write_info(cleaned_markers_matrix, num_markers_for_clade, clade, output_dir,
            samples, secondary_samples, references, secondary_references, trim_sequences, 
            samples_with_n_markers, marker_in_n_samples, secondary_samples_with_n_markers, 
            phylophlan_mode, nprocs)
        info("Done.", init_new_line=True)
        if not debug:
            info("Removing temporary files...", init_new_line=True)
            rmtree(tmp_dir, ignore_errors=False, onerror=None)
            info("Done.", init_new_line=True)


"""
Main function

:param database: the MetaPhlan2 markers database
:param clade_markers: a FASTA containing the markers of a specific clade
:param samples: the folder containing the markers generated with script samples_to_markers.py
:param references: the reference genomes
:param secondary_samples: the folder containing the markers to secondary samples generated with 
    script samples_to_markers.py
:param secondary_references: the secondary reference genomes
:param clade: the clade to investigate
:param output_dir: the output directory
:param trim_sequences: the number of bases to remove when trimming markers
:param samples_with_n_markers: threshold defining the minimun number of markers to keep a main 
    sample
:param marker_in_n_samples: threshold defining the minimum percentage of samples to keep a marker
:param secondary_samples_with_n_markers: threshold defining the minimun number of markers to keep 
    a secondary sample
:param phylophlan_mode: the precision of the phylogenetic analysis
:param phylophlan_configuration: the PhyloPhlAn configuration file
:param mutation_rates: whether get  the mutation rates for the markers
:param print_clades_only: whether print only the potential clades and stop
:param debug: wether to save the tmp folders
:param nprocs: the threads used for execution
"""
def main():
    t0 = time.time()
    args = read_params()
    args = check_params(args)
    info("Start StrainPhlAn " + __version__ + " execution")
    strainphlan(args.database, args.clade_markers, args.samples, args.references, 
        args.secondary_samples, args.secondary_references,  args.clade, args.output_dir, 
        args.trim_sequences, args.sample_with_n_markers, args.marker_in_n_samples,
        args.secondary_sample_with_n_markers, args.phylophlan_mode, args.phylophlan_configuration, 
        args.mutation_rates, args.print_clades_only, args.debug, args.nprocs)
    exec_time = time.time() - t0
    if not args.print_clades_only:
        info("Finish StrainPhlAn " + __version__ + " execution ("+str(round(exec_time, 2))+
            " seconds): Results are stored at \""+args.output_dir+"\"\n",
             init_new_line=True)
    else:
        info("Finish StrainPhlAn " + __version__ + " execution ("+str(round(exec_time, 2))+
            " seconds)\n", init_new_line=True)

if __name__ == '__main__':
    main()
