#!/usr/bin/env python
# Author: Duy Tin Truong (duytin.truong@unitn.it)
#		at CIBIO, University of Trento, Italy

__author__ = 'Duy Tin Truong (duytin.truong@unitn.it)'
__version__ = '0.1'
__date__ = '1st Sep 2014'

import sys
import os
ABS_PATH = os.path.abspath(sys.argv[0])
MAIN_DIR = os.path.dirname(ABS_PATH)
os.environ['PATH'] += ':' + MAIN_DIR
os.environ['PATH'] += ':' + os.path.join(MAIN_DIR, 'mpa3src')
sys.path.append(os.path.join(MAIN_DIR, 'mpa3src'))

import which
import argparse as ap
import cPickle as pickle
import glob
import ooSubprocess
from ooSubprocess import statistics
from ooSubprocess import trace_unhandled_exceptions
import bz2
import gzip
from collections import defaultdict
from tempfile import SpooledTemporaryFile
from Bio import SeqIO, Seq, SeqRecord
from Bio.Alphabet import IUPAC
import pandas
import logging
import logging.config
import sample2markers
import copy
import threading
import numpy
import random
import gc
#import ipdb

shared_variables = type('shared_variables', (object,), {})

# logging config
ifn_logging_config = '%s/logging.ini'%MAIN_DIR
logging.config.fileConfig(ifn_logging_config, disable_existing_loggers=False)
logger = logging.getLogger(__name__)


# functions
def read_params():
    p = ap.ArgumentParser()
    p.add_argument(
        '--ifn_samples',
        nargs='+',
        required=True,
        default=None,
        type=str,
        help='The list of sample files (space separated).'\
                'The wildcard can also be used.')
    p.add_argument(
        '--mpa_pkl', 
        required=True, 
        default=None, 
        type=str, 
        help='The database of metaphlan3.py.')
    p.add_argument(
        '--output_dir', 
        required=True, 
        default='strainer_output', 
        type=str,
        help='The output directory.')
    p.add_argument(
        '--ifn_markers', 
        required=False, 
        default=None, 
        type=str,
        help='The marker file in fasta format.')
    p.add_argument(
        '--nprocs_main', 
        required=False, 
        default=1, 
        type=int,
        help='The number of processors are used for the main threads. '\
             'Default 1.')
    p.add_argument(
        '--nprocs_load_samples', 
        required=False, 
        default=None, 
        type=int,
        help='The number of processors are used for loading samples. '\
             'Default nprocs_main.')
    p.add_argument(
        '--nprocs_align_clean', 
        required=False, 
        default=None, 
        type=int,
        help='The number of processors are used for aligning and cleaning markers. '\
             'Default nprocs_main.')
    p.add_argument(
        '--nprocs_raxml', 
        required=False, 
        default=None, 
        type=int,
        help='The number of processors are used for running raxml. '\
             'Default nprocs_main.')
    p.add_argument(
        '--ifn_ref_genomes',
        nargs='+',
        required=False,
        default=None,
        type=str,
        help='The reference genome file names. They are separated by spaces.')
    p.add_argument(
        '--N_in_marker',
        required=False,
        default=0.2,
        type=float,
        help='The consensus markers with the percentage of N nucleotides greater than '\
                'this threshold are removed. Default 0.2.')
    p.add_argument(
        '--marker_strip_length',
        required=False,
        default=50,
        type=int,
        help='The number of nucleotides will be deleted from each of two ends '\
                'of a marker. Default 50.')
    p.add_argument(
        '--marker_in_clade',
        required=False,
        default=0.8,
        type=float,
        help='In each sample, the clades with the percentage of present markers less than '\
                'this threshold are removed. Default 0.8.')
    p.add_argument(
        '--sample_in_clade',
        required=False,
        default=2,
        type=int,
        help='Only clades present in at least sample_in_clade samples '\
             'are kept. Default 2.')
    p.add_argument(
        '--sample_in_marker',
        required=False,
        default=0.8,
        type=float,
        help='If the percentage of samples that a marker present in is '\
             'less than this threhold, that marker is removed. Default 0.8.')
    p.add_argument(
        '--gap_in_trailing_col',
        required=False,
        default=0.2,
        type=float,
        help='If the number of the trailing nucleotide columns in aligned '\
              'markers with the percentage of gaps greater than '\
              'gap_in_trailing_col is less than gap_trailing_col_limit, '\
              'these columns will be removed. '\
              'Default 0.2.')
    p.add_argument(
        '--gap_trailing_col_limit',
        required=False,
        default=101,
        type=float,
        help='If the number of the trailing nucleotide columns in aligned '\
              'markers with the percentage of gaps greater than '\
              'gap_in_trailing_col is less than gap_trailing_col_limit, '\
              'these columns will be removed. '\
              'Default 101.')
    p.add_argument(
        '--gap_in_internal_col',
        required=False,
        default=0.3,
        type=float,
        help='The internal nucleotide columns in aligned '\
              'markers with the percentage of gaps greater than '\
              'gap_in_internal_col will be removed. '\
              'Default 0.3.')
    p.add_argument(
        '--gap_in_sample',
        required=False,
        default=0.2,
        type=float,
        help='The samples with full sequences from all markers '\
            'and having the percentage of gaps greater than this threshold '\
            'will be removed. Default 0.2.')
    p.add_argument(
        '--N_col',
        required=False,
        default=0.8,
        type=float,
        help='In aligned markers, if the percentage of nucleotide columns '\
              'containing more than N_count Ns '\
              'less than this threshold, these columns will be removed. '
              'Default 0.8.')
    p.add_argument(
        '--N_count',
        required=False,
        default=0,
        type=int,
        help='In aligned markers, if the percentage of nucleotide columns '\
              'containing more than N_count Ns '\
              'less than N_col threshold, these columns will be removed. '\
              'Default 0.')
    p.add_argument(
        '--long_gap_length',
        required=False,
        default=2,
        type=int,
        help='In each concatenated sequence of a sample, sequential '\
                'gap positions is a gap group. '\
                'A gap group with length greater than this '\
                'threshold is considered as '\
                'a long gap group. If the ratio between the number of unique '\
                'positions in all long gap groups and the concatenated sequence '\
                'length is less than long_gap_percentage, these positions '\
                'will be removed from all concatenated sequences. '\
                'Default 2.')
    p.add_argument(
        '--long_gap_percentage',
        required=False,
        default=0.8,
        type=float,
        help='Combining this threshold with long_gap_length to removed long '\
             'gaps. Default 0.8.')
    p.add_argument(
        '--p_value',
        required=False,
        default=0.05,
        type=float,
        help='The p_value to reject a non-polymorphic site.'\
             'Default 0.05.')
    p.add_argument(
        '--clades', 
        nargs='+',
        required=False, 
        default=['all'], 
        type=str,
        help='The clades (space seperated) for which the script will compute '\
                'the marker alignments in fasta format and the phylogenetic '\
                'trees. Default "all".')
    p.add_argument(
        '--print_clades_only', 
        required=False, 
        dest='print_clades_only',
        action='store_true',
        help='Only print the potential clades and stop without building any '\
             'tree. This option is useful when you want to check quickly '\
             'all possible clades and rerun only for some specific ones. '\
             'Default "False".')
    p.set_defaults(print_clades_only=False)
    p.add_argument(
        '--reduce_memory', 
        required=False, 
        dest='reduce_memory', 
        action='store_true',
        help='Do not run many threads in parallel to reduce memory. '\
             'Default "not set".')
    p.set_defaults(reduce_memory=False)
    p.add_argument(
        '--use_threads', 
        required=False, 
        action='store_true',
        dest='use_threads',
        help='Use multithreading. Default "Use multiprocessing".')
    p.set_defaults(use_threads=False)

    return vars(p.parse_args())



def filter_sequence(marker2seq, marker_strip_length, N_in_marker):
    '''
    Filter markers with percentage of N-bases greater than a threshold.

    :param marker2seq: a dictionary containing sequences of a sample.
    marker2seq[marker]['seq'] should return the sequence of the marker.
    :returns: a dictionary containing filtered sequences of samples.
    '''
    remove_markers = [marker for marker in marker2seq if
                      float(marker2seq[marker]['seq'].count('N')) /
                      len(marker2seq[marker]['seq']) > N_in_marker]
    for marker in remove_markers:
        del marker2seq[marker]
    logger.debug('number of markers after N_in_marker: %d'\
                    %(len(marker2seq)))

    remove_markers = []
    for marker in marker2seq:
        if marker_strip_length > 0:
            marker2seq[marker]['seq'] = \
                marker2seq[marker]['seq'][marker_strip_length:-marker_strip_length]
            marker2seq[marker]['freq'] = \
                marker2seq[marker]['freq'][marker_strip_length:-marker_strip_length]
        #marker2seq[marker]['seq'] = marker2seq[marker]['seq'].strip('N')
        if len(marker2seq[marker]['seq']) == 0:
            remove_markers.append(marker)
    for marker in remove_markers:
        del marker2seq[marker]
    logger.debug('number of markers after marker_strip_length: %d'\
                    %(len(marker2seq)))

    return marker2seq




def get_db_clades(db):
    # find singleton clades
    sing_clades = []
    clade2subclades = defaultdict(set)
    for tax in db['taxonomy']:
        tax_clades = tax.split('|')
        for i, clade in enumerate(tax_clades):
            if 't__' not in clade and 's__' not in clade:
                if i < len(tax_clades)-1:
                    if 't__' in tax_clades[-1]:
                        clade2subclades[clade].add('|'.join(tax_clades[i+1:-1]))
                    else:
                        clade2subclades[clade].add('|'.join(tax_clades[i+1:]))
    sing_clades = [clade for clade in clade2subclades if
                             len(clade2subclades[clade]) == 1]

    # extract species
    clade2num_markers = defaultdict(int)
    level = 's__'
    for marker in db['markers']:
        clade = db['markers'][marker]['taxon'].split('|')[-1]
        if level in clade or clade in sing_clades:
            clade2num_markers[clade] = clade2num_markers[clade] + 1
    clade2num_markers = dict(clade2num_markers)

    return sing_clades, clade2num_markers




def align(marker_file):
    oosp = ooSubprocess.ooSubprocess()
    alignment_file = oosp.ex(
        'muscle',
        args=['-quiet', '-in', '-', '-out', '-'],
        in_pipe=marker_file,
        get_out_pipe=True,
        verbose=False)
    return alignment_file




def clean_alignment(
        samples,
        sample2seq,
        sample2freq,
        gap_in_trailing_col,
        gap_trailing_col_limit,
        gap_in_internal_col,
        N_count,
        N_col):

    length = len(sample2seq[sample2seq.keys()[0]])
    logger.debug('marker length: %d', length)
    aligned_samples = sample2seq.keys()
    for sample in samples:
        if sample not in aligned_samples:
            sample2seq[sample] = ['-' for i in range(length)]
            sample2freq[sample] = [(0.0, 0.0, 0.0) for i in range(length)]

    df_seq = pandas.DataFrame.from_dict(sample2seq, orient='index')
    df_freq = pandas.DataFrame.from_dict(sample2freq, orient='index')

    # remove trailing gap columns
    del_cols = []
    for i in range(len(df_seq.columns)):
        if float(list(df_seq[df_seq.columns[i]]).count('-')) / len(samples) <= gap_in_trailing_col:
            break
        else:
            del_cols.append(df_seq.columns[i])
    for i in reversed(range(len(df_seq.columns))):
        if float(list(df_seq[df_seq.columns[i]]).count('-')) / len(samples) <= gap_in_trailing_col:
            break
        else:
            del_cols.append(df_seq.columns[i])
    if len(del_cols) < gap_trailing_col_limit:
        df_seq.drop(del_cols, axis=1, inplace=True)
        df_freq.drop(del_cols, axis=1, inplace=True)
        logger.debug('length after gap_in_trailing_col: %d', len(df_seq.columns))
    else:
        logger.debug('do not use gap_in_trailing_col as the number of del_cols is %d'%len(del_cols))

    # remove internal gap columns
    del_cols = []
    for i in range(len(df_seq.columns)):
        if float(list(df_seq[df_seq.columns[i]]).count('-')) / len(samples) > gap_in_internal_col:
            del_cols.append(df_seq.columns[i])
    df_seq.drop(del_cols, axis=1, inplace=True)
    df_freq.drop(del_cols, axis=1, inplace=True)
    logger.debug('length after gap_in_internal_col: %d', len(df_seq.columns))

    # remove N columns
    if len(df_seq.columns) > 0:
        del_cols = []
        remove_N_col = False
        for i in range(len(df_seq.columns)):
            if list(df_seq[df_seq.columns[i]]).count('N') > N_count:
                del_cols.append(df_seq.columns[i])
        if float(len(del_cols)) / len(df_seq.columns) < N_col:
            remove_N_col = True
            df_seq.drop(del_cols, axis=1, inplace=True)   
            df_freq.drop(del_cols, axis=1, inplace=True)   
            logger.debug('length after N_col: %d', len(df_seq.columns))
            
        if N_count > 0 or not remove_N_col:
            logger.debug('replace Ns by gaps for all samples')
            for sample in samples:
                seq = ''.join(df_seq.loc[sample])
                logger.debug('sample %s, number of Ns: %d'\
                                %(sample, seq.count('N')))
                sample2seq[sample] = list(seq.replace('N', '-'))
        else:
            for sample in samples:
                sample2seq[sample] = df_seq.loc[sample].tolist()
        for sample in samples:
            sample2freq[sample] = df_freq.loc[sample].tolist()
    else:
        sample2seq = {}
        sample2freq = {}

    return sample2seq, sample2freq




def add_ref_genomes(sample2marker, marker_records, ifn_ref_genomes, tmp_dir):
    logger.debug('add %d reference genomes'%len(ifn_ref_genomes))
    logger.debug('Number of samples: %d'%len(sample2marker))

    # marker list
    if len(sample2marker) == 0:
        unique_markers = set(marker_records.keys())
    else:
        unique_markers = set([])
        for sample in sample2marker:
            for marker in sample2marker[sample]:
                if marker not in unique_markers:
                    unique_markers.add(marker)
    logger.debug('Number of unique markers: %d'%len(unique_markers))

    # add ifn_ref_genomes
    oosp = ooSubprocess.ooSubprocess(tmp_dir=tmp_dir)
    logger.debug('load genome contigs')
    p1 = SpooledTemporaryFile(dir=tmp_dir)
    contigs = defaultdict(dict)
    for ifn_genome in ifn_ref_genomes:
        genome = ooSubprocess.splitext2(ifn_genome)[0]
        if ifn_genome[-4:] == '.bz2':
            ifile_genome = bz2.BZ2File(ifn_genome, 'r')
        elif ifn_genome[-3:] == '.gz':
            ifile_genome = gzip.GzipFile(ifn_genome, 'r')
        elif ifn_genome[-4:] == '.fna':
            ifile_genome = open(ifn_genome, 'r')
        else:
            logger.error('Unknown file type of %s. '%ifn_genome +\
                            'It should be .fna.bz2, .fna.gz, or .fna!')
            exit(1)

        # extract genome contigs
        for rec in SeqIO.parse(ifile_genome, 'fasta'):
            #rec.name = genome + '___' + rec.name
            if rec.name in contigs:
                logger.error(
                            'Error: Contig %s in genome%s'\
                            %(rec.name.split('___')[-1], genome)\
                            + ' are not unique!')
                exit(1)
            contigs[rec.name]['seq'] = str(rec.seq)
            contigs[rec.name]['genome'] = genome
            SeqIO.write(rec, p1, 'fasta')

        ifile_genome.close()
    p1.seek(0)
                        
    # build blastdb
    logger.debug('build blastdb')
    blastdb_prefix = oosp.ftmp('genome_blastn_db_%s'%(random.random()))
    if len(glob.glob('%s*'%blastdb_prefix)):
        logger.error('blastdb exists! Please remove it or rerun!')
        exit(1)
    oosp.ex('makeblastdb', 
                args=[
                        '-dbtype', 'nucl', 
                        '-title', 'genome_db',
                        '-out', blastdb_prefix],
                in_pipe=p1,
                verbose=True)

    # blast markers against contigs
    logger.debug('blast markers against contigs')
    p1 = SpooledTemporaryFile(dir=tmp_dir)
    for marker in unique_markers:
        SeqIO.write(marker_records[marker], p1, 'fasta')
    p1.seek(0)
    blastn_args = [
                    '-db', blastdb_prefix,
                    '-outfmt', '6',
                    '-evalue', '1e-10',
                    '-max_target_seqs', '1000000000']
    if args['nprocs_main'] > 1:
        blastn_args += ['-num_threads', str(args['nprocs_main'])]
    output = oosp.ex(
                        'blastn', 
                        args=blastn_args,
                        in_pipe=p1,
                        get_out_pipe=True,
                        verbose=True)
    
    #output = output.split('\n')
    for line in output:
        if line.strip() == '':
            break
        line = line.strip().split()
        query = line[0]
        target = line[1]
        pstart = int(line[8])-1
        pend = int(line[9])-1
        genome = contigs[target]['genome']
        if query not in sample2marker[genome]:
            sample2marker[genome][query] = {}
            if pstart < pend:
                sample2marker[genome][query]['seq'] = contigs[target]['seq'][pstart:pend+1]
            else:
                sample2marker[genome][query]['seq'] = \
                     str(Seq.Seq(
                                contigs[target]['seq'][pend:pstart+1],
                                IUPAC.unambiguous_dna).reverse_complement())
            sample2marker[genome][query]['freq'] = [(0.0, 0.0, 0.0) for i in \
                                                    range(len(sample2marker[genome][query]['seq']))]
            sample2marker[genome][query]['seq'] = sample2marker[genome][query]['seq'].upper()

    # remove database
    for fn in glob.glob('%s*'%blastdb_prefix):
        os.remove(fn)

    logger.debug('Number of samples and genomes: %d'%len(sample2marker))
    return sample2marker




@trace_unhandled_exceptions
def align_clean(args):
    marker = args['marker']
    sample2marker = shared_variables.sample2marker #args['sample2marker']
    clade = args['clade']
    gap_in_trailing_col = args['gap_in_trailing_col']
    gap_trailing_col_limit = args['gap_trailing_col_limit']
    gap_in_internal_col = args['gap_in_internal_col']
    N_col = args['N_col']
    N_count = args['N_count']
    sample_in_marker = args['sample_in_marker']
    tmp_dir = args['tmp_dir']

    logger.debug('align and clean for marker: %s'%marker)
    marker_file = SpooledTemporaryFile(dir=tmp_dir)
    sample_count = 0
    for sample in iter(sample2marker.keys()):
        if marker in iter(sample2marker[sample].keys()):
            sample_count += 1
            SeqIO.write(
                SeqRecord.SeqRecord(
                    id=sample,
                    description='',
                    seq=Seq.Seq(sample2marker[sample][marker]['seq'])),
                marker_file,
                'fasta')
    ratio = float(sample_count) / len(sample2marker)
    if  ratio < sample_in_marker:
        logger.debug('skip this marker because percentage of samples '\
                     'it present is %f < sample_in_marker'%ratio)
        return {}, {}

    marker_file.seek(0)
    alignment_file = align(marker_file)

    sample2seq = {}
    sample2freq = {}
    for rec in SeqIO.parse(alignment_file, 'fasta'):
        sample = rec.name
        sample2seq[sample] = list(str(rec.seq))
        sample2freq[sample] = list(sample2marker[sample][marker]['freq'])
        for i, c in enumerate(sample2seq[sample]):
            if c == '-':
                sample2freq[sample].insert(i, (0.0, 0.0, 0.0))
        logger.debug('alignment length of sample %s is %d, %d'%(
                        sample,
                        len(sample2seq[sample]),
                        len(sample2freq[sample])))
    alignment_file.close()
    logger.debug('alignment for marker %s is done'%marker)


    if len(sample2seq) == 0:
        logger.error('Fatal error in alignment step!')
        exit(1)

    sample2seq, sample2freq = clean_alignment(
                                    sample2marker.keys(),
                                    sample2seq, 
                                    sample2freq,
                                    gap_in_trailing_col,
                                    gap_trailing_col_limit,
                                    gap_in_internal_col,
                                    N_count,
                                    N_col)
    logger.debug('cleaning for marker %s is done'%marker)

    return sample2seq, sample2freq




def build_tree(
        clade,
        sample2marker, 
        clade2num_markers,
        sample_in_clade,
        sample_in_marker,
        gap_in_trailing_col,
        gap_trailing_col_limit,
        gap_in_internal_col,
        N_count,
        N_col,
        gap_in_sample,
        long_gap_length,
        long_gap_percentage,
        p_value,
        output_dir,
        nprocs_align_clean,
        nprocs_raxml,
        use_threads):

    # build the tree for each clade
    if len(sample2marker) < sample_in_clade:
        logger.debug(
                        'skip clade %s because number of present samples '
                        'is %d'%(clade, len(sample2marker)))
        return

    ofn_cladeinfo = os.path.join(output_dir, '%s.info'%clade)
    ofile_cladeinfo = open(ofn_cladeinfo, 'w')

    logger.debug('clade: %s', clade)
    ofile_cladeinfo.write('clade: %s\n'%clade)
    logger.debug('number of samples: %d', len(sample2marker))
    ofile_cladeinfo.write('number of samples: %d\n'\
                          %len(sample2marker))
    logger.debug('number of markers of the clade in db: %d'\
                 %clade2num_markers[clade])
    ofile_cladeinfo.write('number of markers of the clade in db: %d\n'\
                          %clade2num_markers[clade])

    # align sequences in each marker
    markers = set([]) 
    for sample in sample2marker:
        for marker in sample2marker[sample]:
            if marker not in markers:
                markers.add(marker)
    markers = sorted(list(markers))

    logger.debug('number of used markers: %d'%len(markers))
    ofile_cladeinfo.write('number of used markers: %d\n'%len(markers))
    logger.debug('fraction of used markers: %f'\
                %(float(len(markers)) / clade2num_markers[clade]))
    ofile_cladeinfo.write('fraction of used markers: %f\n'\
                %(float(len(markers)) / clade2num_markers[clade]))

    logger.debug('align and clean')
    args_list = []

    # parallelize
    for i in range(len(markers)):
        args_list.append({})
        args_list[i]['marker'] = markers[i]
        args_list[i]['clade'] = clade
        args_list[i]['gap_in_trailing_col'] = gap_in_trailing_col
        args_list[i]['gap_trailing_col_limit'] = gap_trailing_col_limit
        args_list[i]['gap_in_internal_col'] = gap_in_internal_col
        args_list[i]['N_count'] = N_count
        args_list[i]['N_col'] = N_col
        args_list[i]['sample_in_marker'] = sample_in_marker
        args_list[i]['tmp_dir'] = output_dir

    logger.debug('start to align_clean for all markers')
    results = ooSubprocess.parallelize(
                                       align_clean, 
                                       args_list, 
                                       nprocs_align_clean,
                                       use_threads=use_threads)

    sample2seqs, sample2freqs = zip(*results)
    sample2fullseq = defaultdict(list)
    sample2fullfreq = defaultdict(list)
    empty_markers = []
    pos = 0
    marker_pos = []
    for i in range(len(sample2seqs)):
        #logger.debug('marker_name: %s, seq: %s'%(markers[i], sample2seqs[i]))
        if len(sample2seqs[i]):
            for sample in sample2seqs[i]:
                sample2fullseq[sample] += sample2seqs[i][sample]
                sample2fullfreq[sample] += sample2freqs[i][sample]
            marker_pos.append([markers[i], pos])
            pos += len(sample2seqs[i][sample])
        else:
            empty_markers.append(markers[i])

    logger.debug(
            'number of markers after deleting empty markers: %d',
            len(markers) - len(empty_markers))
    ofile_cladeinfo.write('number of markers after deleting '\
                          'empty markers: %d\n'%
                          (len(markers) - len(empty_markers)))
    logger.debug('fraction of used markers after deleting empty markers: '\
            '%f'%(float(len(markers) - len(empty_markers)) / clade2num_markers[clade]))
    ofile_cladeinfo.write('fraction of used markers after deleting empty '\
            'markers: %f\n'\
            %(float(len(markers) - len(empty_markers)) / clade2num_markers[clade]))


    if len(sample2fullseq) == 0:
        logger.debug('all markers were removed, skip this clade!')
        ofile_cladeinfo.write('all markers were removed, skip this clade!\n')
        return
    
    # remove long gaps
    logger.debug('full sequence length before long_gap_length: %d'\
                    %(len(sample2fullseq[sample2fullseq.keys()[0]])))
    ofile_cladeinfo.write(
                    'full sequence length before long_gap_length: %d\n'\
                    %(len(sample2fullseq[sample2fullseq.keys()[0]])))

    df_seq = pandas.DataFrame.from_dict(sample2fullseq, orient='index')
    df_freq = pandas.DataFrame.from_dict(sample2fullfreq, orient='index')
    del_cols = []
    del_pos = []
    for sample in sample2fullseq:
        row = df_seq.loc[sample]
        gap_in_cols = []
        gap_in_pos = []
        for i in range(len(row)):
            if row[df_seq.columns[i]] == '-':
                gap_in_cols.append(df_seq.columns[i])
                gap_in_pos.append(i)
            else:
                if len(gap_in_cols) > long_gap_length:
                    del_cols += gap_in_cols
                    del_pos += gap_in_pos
                gap_in_cols = []
                gap_in_pos = []

        if len(gap_in_cols) > long_gap_length:
            del_cols += gap_in_cols
            del_pos += gap_in_pos

    del_cols = list(set(del_cols))
    del_pos = sorted(list(set(del_pos)))
    del_ratio = float(len(del_cols)) / len(sample2fullseq[sample])
    if del_ratio < long_gap_percentage:
        df_seq.drop(del_cols, axis=1, inplace=True)
        df_freq.drop(del_cols, axis=1, inplace=True)
        for sample in sample2fullseq:
            sample2fullseq[sample] = df_seq.loc[sample].tolist()
            sample2fullfreq[sample] = df_freq.loc[sample].tolist()
        logger.debug('full sequence length after long_gap_length: %d'\
                        %(len(sample2fullseq[sample2fullseq.keys()[0]])))
        ofile_cladeinfo.write(
                        'full sequence length after long_gap_length: %d\n'\
                        %(len(sample2fullseq[sample2fullseq.keys()[0]])))

        for i in range(len(marker_pos)):
            num_del = 0 
            for p in del_pos:
                if marker_pos[i][1] > p:
                    num_del += 1
            marker_pos[i][1] -= num_del
    else:
        logger.debug('do not apply long_gap_length because '\
                        'long_gap_percentage is not satisfied. '\
                        'del_ratio: %f'%del_ratio)
        ofile_cladeinfo.write('do not apply long_gap_length because '\
                                'long_gap_percentage is not satisfied. '\
                                'del_ratio: %f\n'%del_ratio)

    ofn_clademarker = os.path.join(output_dir, '%s.marker_pos'%clade)
    with open(ofn_clademarker, 'w') as ofile_clademarker:
        for m, p in marker_pos:
            ofile_clademarker.write('%s\t%d\n'%(m, p))

    # remove samples with more than 10% of gaps
    logger.debug(
            'number of samples before gap_in_sample: %d'\
            %len(sample2fullseq))
    ofile_cladeinfo.write(
            'number of samples before gap_in_sample: %d\n'\
            %len(sample2fullseq))
    for sample in sample2marker:
        if float(sample2fullseq[sample].count('-')) / \
            len(sample2fullseq[sample]) > gap_in_sample:
            del sample2fullseq[sample]
    logger.debug(
            'number of samples after gap_in_sample: %d'\
            %len(sample2fullseq))
    ofile_cladeinfo.write(
            'number of samples after gap_in_sample: %d\n'\
            %len(sample2fullseq))

    # log gaps
    sequential_gaps = []
    all_gaps = []
    for sample in sample2fullseq:
        agap = 0
        sgap = 0
        row2 = sample2fullseq[sample]
        for i in range(len(row2)):
            if row2[i] == '-':
                sgap += 1
                agap += 1
            elif sgap > 0:
                sequential_gaps.append(sgap)
                sgap = 0
        all_gaps.append(agap)
        
    ofile_cladeinfo.write('all_gaps:\n' + statistics(all_gaps)[1])
    if sequential_gaps == []:
        sequential_gaps = [0]
    ofile_cladeinfo.write('sequential_gaps:\n' + \
                          statistics(sequential_gaps)[1])
    ofile_cladeinfo.close()

    # compute ppercentage of polymorphic sites
    ofn_pol = os.path.join(output_dir, '%s.polymorphic'%clade)
    logger.debug('polymorphic file: %s'%ofn_pol)
    with open(ofn_pol, 'w') as ofile:
        ofile.write('#sample\tpercentage_of_polymorphic_sites\tavg_freq\tmedian_freq\tstd_freq\tmin_freq\tmax_freq\tq90_freq\tq10_freq\tavg_coverage\tmedian_coverage\tstd_coverage\tmin_coverage\tmax_coverage\tq90_coverage\tq10_coverage\n')
        for sample in sample2fullfreq:
            freqs = [x[0] for x in sample2fullfreq[sample] if x[0] > 0 and x[0] < 1 and x[2] < p_value]
            coverages = [x[1] for x in sample2fullfreq[sample] if x[0] > 0 and x[0] < 1 and x[2] < p_value]
            ofile.write('%s\t%f'%(sample, float(len(freqs)) / len(sample2fullfreq[sample])))
            for vals in [freqs, coverages]:
                if len(vals):
                    ofile.write('\t%f\t%f\t%f\t%f\t%f\t%f\t%f'%(\
                                numpy.average(vals),
                                numpy.percentile(vals,50),
                                numpy.std(vals),
                                numpy.min(vals),
                                numpy.max(vals),
                                numpy.percentile(vals,90),
                                numpy.percentile(vals,10),
                               ))
                else:
                    ofile.write('\t0\t0\t0\t0\t0\t0\t0')
            ofile.write('\n')


    # save merged alignment
    ofn_align = os.path.join(output_dir, '%s.fasta'%clade)
    logger.debug('alignment file: %s'%ofn_align)
    with open(ofn_align, 'w') as ofile:
        for sample in sample2fullseq:
            SeqIO.write(
                SeqRecord.SeqRecord(
                    id=sample,
                    description='',
                    seq=Seq.Seq(''.join(sample2fullseq[sample]))),
                ofile, 
                'fasta')

    # produce tree
    oosp = ooSubprocess.ooSubprocess() 
    #ofn_tree = os.path.join(output_dir, '%s.tree'%clade)
    #oosp.ex('FastTree', args=['-quiet', '-nt', ofn_align], out_fn=ofn_tree)
    ofn_tree = clade + '.tree'
    logger.debug('tree file: %s'%ofn_tree)
    try:
        for fn in glob.glob('%s/RAxML_*%s'
                    %(os.path.abspath(output_dir), ofn_tree)):
            os.remove(fn)
        raxml_args = [
                            '-m', 'GTRCAT', 
                            '-s', os.path.abspath(ofn_align), 
                            '-w', os.path.abspath(output_dir),
                            '-n', ofn_tree,
                            '-p', '1234'
                        ]

        if nprocs_raxml > 1:
            raxml_args += ['-T', str(nprocs_raxml)]
            raxml_prog = 'raxmlHPC-PTHREADS-SSE3'
        else:
            raxml_prog = 'raxmlHPC'
        oosp.ex(
                raxml_prog,
                args=raxml_args
                )
    except:
        logger.info('Cannot build the tree! The number of samples is too few '\
                    'or there is some error with raxmlHMP')
        pass




@trace_unhandled_exceptions
def load_sample(args):
    ifn_sample = args['ifn_sample']
    logger.debug('load %s'%ifn_sample)
    output_dir = args['output_dir'] 
    ifn_markers = args['ifn_markers']
    clades = args['clades']
    kept_clade = args['kept_clade']
    db = shared_variables.db
    sing_clades = shared_variables.sing_clades
    clade2num_markers = shared_variables.clade2num_markers
    marker_in_clade = args['marker_in_clade']

    sample = ooSubprocess.splitext2(ifn_sample)[0]
    #gc.disable()
    marker2seq = pickle.load(open(ifn_sample, 'rb'))
    #gc.enable()

    if kept_clade:
        # remove redundant clades and markers
        nmarkers = 0
        for marker in marker2seq.keys():
            clade = db['markers'][marker]['taxon'].split('|')[-1]
            if clade == kept_clade:
                nmarkers += 1
            else:
                del marker2seq[marker]
        if float(nmarkers) / clade2num_markers[kept_clade] < marker_in_clade:
            marker2seq = {}

        # reformat 'pileup'
        for m in marker2seq:
            freq = marker2seq[m]['freq']
            marker2seq[m]['freq'] = [(0.0, 0.0, 0.0) for i in \
                                     range(len(marker2seq[m]['seq']))]
            for p in freq:
                marker2seq[m]['freq'][p] = freq[p]
            marker2seq[m]['seq'] = marker2seq[m]['seq'].replace('-', 'N') # make sure we have no gaps in the sequence

        return marker2seq
    else:
        # remove redundant clades and markers
        clade2n_markers = defaultdict(int)
        remove_clade = []
        for marker in marker2seq:
            clade = db['markers'][marker]['taxon'].split('|')[-1]
            if 's__' in clade or clade in sing_clades:
                clade2n_markers[clade] = clade2n_markers[clade] + 1
            else:
                remove_clade.append(clade)
        remove_clade += [clade for clade in clade2n_markers if
                         float(clade2n_markers[clade]) \
                         / float(clade2num_markers[clade]) < marker_in_clade]
        remove_marker = [marker for marker in marker2seq if
                         db['markers'][marker]['taxon'].split('|')[-1] in
                         remove_clade]
        for marker in remove_marker:
            del marker2seq[marker]
     
        sample_clades = set([])
        for marker in marker2seq:
            clade = db['markers'][marker]['taxon'].split('|')[-1]
            sample_clades.add(clade)
        return sample_clades




def load_all_samples(args, kept_clade):
    ifn_samples = args['ifn_samples']
    if ifn_samples == ['None']:
        ifn_samples = None

    if ifn_samples == None:
        return None
    else:
        args_list = []
        for ifn_sample in ifn_samples:
            func_args = {}
            func_args['ifn_sample'] = ifn_sample
            func_args['kept_clade'] = kept_clade
            for k in [ 
                      'output_dir',
                      'ifn_markers', 'nprocs_load_samples', 
                      'clades',
                      'marker_in_clade',
                      'mpa_pkl',
                      ]:
                func_args[k] = args[k]
            args_list.append(func_args)

        results = ooSubprocess.parallelize(
                                           load_sample,
                                           args_list,
                                           args['nprocs_load_samples'],
                                           use_threads=args['use_threads'])
        if kept_clade:
            sample2marker = {}
            for i in range(len(ifn_samples)):
                sample = ooSubprocess.splitext2(ifn_samples[i])[0]
                if len(results[i]): # skip samples with no markers
                    sample2marker[sample] = results[i] 
            return sample2marker
        else:
            all_clades = set([])
            for r in results:
                for c in r:
                    all_clades.add(c)
            all_clades = sorted(list(all_clades))
            return all_clades




def strainer(args):
    # check conditions
    ooSubprocess.mkdir(args['output_dir'])
    with open(os.path.join(args['output_dir'], 'arguments.txt'), 'w') as ofile:
        #for para in args:
        #    ofile.write('%s\n'%para)
        #    ofile.write('%s\n'%args[para])
        ofile.write('%s\n'%' '.join(sys.argv))
        ofile.write('%s'%args)

    if args['ifn_markers'] == None and args['ifn_ref_genomes'] != None:
        logger.error('ifn_ref_genomes is set but ifn_markers is not set!')
        exit(1)

    if args['ifn_markers'] != None and args['ifn_ref_genomes'] != None:
        if len(args['clades']) != 1 or args['clades'] == ['all']:
            logger.error('Only one clade can be specified when adding '\
                         'reference genomes')
            exit(1)

    if args['ifn_markers'] == None and args['clades'] == ['singleton']:
        logger.error('clades is set to singleton but ifn_markers is not set!')
        exit(1)

    if args['ifn_samples'] == None:
        args['clades'] = ['singleton']

    if args['nprocs_load_samples'] == None:
        args['nprocs_load_samples'] = args['nprocs_main']

    if args['nprocs_align_clean'] == None:
        args['nprocs_align_clean'] = args['nprocs_main']

    if args['nprocs_raxml'] == None:
        args['nprocs_raxml'] = args['nprocs_main']

    # logging config
    # create a file handler
    handler = logging.FileHandler(
                                    os.path.join(
                                                 args['output_dir'], 
                                                 'log.txt'),
                                'w')
    handler.setLevel(logger.getEffectiveLevel())

    # create a logging format
    formatter = None
    with open(ifn_logging_config, 'r') as ifile:
        for line in ifile:
            line=line.strip()
            if 'format=' in line:
                formatter = logging.Formatter(line.replace('format=', ''))
                break
    if formatter == None:
        logger.error('Cannot find the FORMAT line in %s'%ifn_logging_config)
        exit(1)
    handler.setFormatter(formatter)

    # add the handlers to the logger
    logger.addHandler(handler)

    # load mpa_pkl
    logger.info('Load mpa_pkl')
    db = pickle.load(bz2.BZ2File(args['mpa_pkl']))
    shared_variables.db = db

    # reduce and convert to shared memory
    #logger.debug('converting db')
    db['taxonomy'] = db['taxonomy'].keys()
    for m in db['markers']:
        del db['markers'][m]['clade']
        del db['markers'][m]['ext']
        del db['markers'][m]['len']
        del db['markers'][m]['score']
    gc.collect()
    #logger.debug('converted db')
    
    # get clades from db
    logger.info('Get clades from db')
    sing_clades, clade2num_markers = get_db_clades(db)
    shared_variables.sing_clades = sing_clades
    shared_variables.clade2num_markers = clade2num_markers


    # get clades from samples
    if args['clades'] == ['all']:
        logger.info('Get clades from samples')
        args['clades'] = load_all_samples(args, kept_clade=None)

    if args['print_clades_only']:
        for c in args['clades']:
            print c
        return

    # add reference genomes
    ref2marker = defaultdict(dict)
    if args['ifn_markers'] != None and args['ifn_ref_genomes'] != None:
        logger.info('Add reference genomes')
        marker_records = {}
        for rec in SeqIO.parse(open(args['ifn_markers'], 'r'), 'fasta'):
            marker_records[rec.name] = rec
        add_ref_genomes(
                        ref2marker, 
                        marker_records, 
                        args['ifn_ref_genomes'],
                        args['output_dir'])

        # remove bad reference genomes
        nmarkers = 0
        for rec in SeqIO.parse(open(args['ifn_markers'], 'r'), 'fasta'):
            nmarkers += 1

        remove_ref = []
        for ref in ref2marker:
            if float(len(ref2marker[ref])) / nmarkers < args['marker_in_clade']:
                remove_ref.append(ref)
        for ref in remove_ref:
            del ref2marker[ref]
    ref2marker = dict(ref2marker)

    # build tree for each clade
    for clade in args['clades']:
        logger.info('Build the tree for %s'%clade)

        # load samples and reference genomes
        if clade == 'singleton':
            sample2marker = {}
        else:
            sample2marker = load_all_samples(args, kept_clade=clade)

        for r in ref2marker:
            sample2marker[r] = ref2marker[r]
        logger.debug('number of samples and reference genomes: %d'%len(sample2marker))

        for sample in sample2marker:
            logger.debug('number of markers in sample %s: %d'%(
                          sample,
                          len(sample2marker[sample])))

        # Filter sequences
        logger.debug('Filter consensus marker sequences')
        for sample in sample2marker:
            sample2marker[sample] = filter_sequence(
                                                    sample2marker[sample],
                                                    args['marker_strip_length'],
                                                    args['N_in_marker'])

        # remove samples with percentage of markers less than marker_in_clade
        logger.debug('remove samples with percentage of markers '\
                     'less than marker_in_clade')
        for sample in sample2marker.keys():
            if len(sample2marker[sample]):
                marker = sample2marker[sample].keys()[0]
                clade = db['markers'][marker]['taxon'].split('|')[-1]
                if len(sample2marker[sample]) / \
                    float(clade2num_markers[clade]) < args['marker_in_clade']:
                        del sample2marker[sample]
            else:
                del sample2marker[sample]

        # build trees
        logger.info('Build trees')
        shared_variables.sample2marker = sample2marker
        build_tree(
            clade=clade,
            sample2marker=sample2marker, 
            clade2num_markers=clade2num_markers,
            sample_in_clade=args['sample_in_clade'],
            sample_in_marker=args['sample_in_marker'],
            gap_in_trailing_col=args['gap_in_trailing_col'],
            gap_trailing_col_limit=args['gap_trailing_col_limit'],
            gap_in_internal_col=args['gap_in_internal_col'],
            N_count=args['N_count'],
            N_col=args['N_col'],
            gap_in_sample=args['gap_in_sample'],
            long_gap_length=args['long_gap_length'],
            long_gap_percentage=args['long_gap_percentage'],
            p_value=args['p_value'],
            output_dir=args['output_dir'],
            nprocs_align_clean=args['nprocs_align_clean'],
            nprocs_raxml=args['nprocs_raxml'],
            use_threads=args['use_threads'])
        del shared_variables.sample2marker
        del sample2marker
        #gc.collect()

    logger.info('Finished!')




def check_dependencies(args):
    programs = ['muscle']

    if args['ifn_markers'] != None or args['ifn_ref_genomes'] != None:
            programs += ['blastn', 'makeblastdb']

    if args['nprocs_main'] > 1:
        programs += ['raxmlHPC-PTHREADS-SSE3']
    else:
        programs += ['raxmlHPC']

    for prog in programs:
        if not which.is_exe(prog):
            logger.error('Cannot find %s in the executable path!'%prog)
            exit(1)




if __name__ == "__main__":
    args = read_params()
    check_dependencies(args)
    strainer(args)
