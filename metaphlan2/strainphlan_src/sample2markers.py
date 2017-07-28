#!/usr/bin/env python
# Author: Duy Tin Truong (duytin.truong@unitn.it)
#		at CIBIO, University of Trento, Italy


import sys
import os
ABS_PATH = os.path.abspath(sys.argv[0])
MAIN_DIR = os.path.dirname(ABS_PATH)
os.environ['PATH'] += ':%s'%MAIN_DIR
os.environ['PATH'] += ':%s'%os.path.join(MAIN_DIR, 'strainphlan_src')
sys.path.append(MAIN_DIR)
sys.path.append(os.path.join(MAIN_DIR, 'strainphlan_src'))

import argparse as ap
import glob
import ooSubprocess
from ooSubprocess import print_stderr, trace_unhandled_exceptions
import ConfigParser
from Bio import SeqIO, Seq, SeqRecord
import cStringIO
import msgpack
import random
import subprocess
import bz2
import gzip
import logging
import logging.config
import tarfile
import threading
import multiprocessing
import pysam
from collections import defaultdict
from scipy import stats
import numpy

logging.basicConfig(level=logging.DEBUG, stream=sys.stderr,
                    disable_existing_loggers=False,
                    format='%(asctime)s | %(levelname)s | %(name)s | %(funcName)s | %(lineno)d | %(message)s')
logger = logging.getLogger(__name__)



def read_params():
    p = ap.ArgumentParser()
    p.add_argument('--ifn_samples', nargs='+', required=True, default=None, type=str)
    p.add_argument('--ifn_markers', required=False, default=None, type=str)
    p.add_argument('--output_dir', required=True, default=None, type=str)
    p.add_argument('--nprocs', required=False, default=1, type=int)
    p.add_argument('--min_read_len', required=False, default=90, type=int)
    p.add_argument('--min_align_score', required=False, default=None, type=int)
    p.add_argument('--min_base_quality', required=False, default=30, type=float)
    p.add_argument('--min_read_depth', required=False, default=0, type=int)
    p.add_argument('--error_rate', required=False, default=0.01, type=float)
    p.add_argument('--marker2file_ext', required=False, default='.markers', type=str)
    p.add_argument('--sam2file_ext', required=False, default='.sam.bz2', type=str)
    p.add_argument('--samtools_exe', required=False, default='samtools', type=str)
    p.add_argument('--bcftools_exe', required=False, default='bcftools', type=str)
    p.add_argument(
        '--verbose', 
        required=False, 
        dest='quiet',
        action='store_false',
        help='Show all information. Default "not set".')
    p.set_defaults(quiet=True)
    '''
    p.add_argument(
        '--use_processes', 
        required=False, 
        default=False, 
        action='store_false',
        dest='use_threads',
        help='Use multiprocessing. Default "Use multithreading".')
    p.set_defaults(use_threads=True)
    '''
    p.add_argument(
        '--input_type',
        required=True,
        default=None,
        type=str,
        choices=['fastq', 'sam'],
        help='The input type:'\
                ' fastq, sam. Sam '\
                ' files can be obtained from the previous run of'\
                ' this script or strainphlan.py).')

    return vars(p.parse_args())


def build_bowtie2db(ifn_markers, tmp_dir, error_pipe=None):
    # build bowtie2-index
    if not os.path.isfile(ifn_markers):
        error = 'ifn_markers %s does not exist!'%ifn_markers
        logger.error(error)
        raise Exception(error)

    if not os.path.isdir(tmp_dir):
        ooSubprocess.mkdir(tmp_dir)
    bt2_base = ooSubprocess.splitext(ifn_markers)[0]
    index_fns = glob.glob('%s/%s.*'%(tmp_dir, bt2_base))
    index_path = os.path.join(tmp_dir, bt2_base)
    oosp = ooSubprocess.ooSubprocess(tmp_dir)
    if index_fns == []:
        oosp.ex(
                'bowtie2-build', 
                ['--quiet', ifn_markers, index_path],
                stderr=error_pipe)
    else:
        logger.warning('bowtie2-indexes of %s are ready, skip rebuilding!'
                        %(bt2_base))
    return index_path



def sample2markers(
        ifn_sample, 
        min_read_len,
        min_align_score,
        min_base_quality,
        min_read_depth,
        error_rate,
        ifn_markers, 
        index_path,
        nprocs=1, 
        sam2file=None,
        marker2file=None, 
        tmp_dir='tmp',
        quiet=False,
        samtools_exe='samtools',
        bcftools_exe='bcftools'):

    '''
    Compute the consensus markers in a sample file ifn_sample.

    :param ifn_sample: the sample file in fastq format.
    :param marker2file: the file name to store the consensus markers.
    :param sam2file: the file name to store the sam content.
    :returns: if marker2file==None, return the dictionary of the consensus
    markers. Otherwise, save the result in fasta format to the file declared by 
    marker2file
    '''

    if quiet:
        error_pipe = open(os.devnull, 'w')
    else:
        error_pipe = None
    oosp = ooSubprocess.ooSubprocess(tmp_dir)

    # sample to sam
    sample_pipe = oosp.chain(
                             'dump_file.py', 
                             args=['--input_file', ifn_sample],
                             stderr=error_pipe
                             )
    filter_length_pipe = oosp.chain(
                                    'fastx_len_filter.py',
                                    args=['--min_len', str(min_read_len)],
                                    in_pipe=sample_pipe,
                                    stderr=error_pipe
                                    )
    bowtie2_pipe = oosp.chain(
                                'bowtie2', 
                                args=[
                                        '-U', '-',
                                        '-x', index_path,
                                        '--very-sensitive',
                                        '--no-unal',
                                        '-p', str(nprocs)],
                                in_pipe=filter_length_pipe,
                                stderr=error_pipe)
    if sam2file == None:
        sam_pipe = bowtie2_pipe
    else:
        oosp.chain(
                    'compress_file.py',
                    args=['--output_file', sam2file], 
                    in_pipe=bowtie2_pipe,
                    stderr=error_pipe,
                    stop=True)

        sam_pipe = oosp.chain(
                                'dump_file.py', 
                                args=['--input_file', sam2file],
                                stderr=error_pipe)
    ofn_bam_sorted_prefix = os.path.join(
                                tmp_dir,
                                os.path.basename(ifn_sample) + '.bam.sorted')

    return sam2markers(
               sam_file=sam_pipe, 
               ofn_bam_sorted_prefix=ofn_bam_sorted_prefix,
               min_align_score=min_align_score,
               min_base_quality=min_base_quality,
               min_read_depth=min_read_depth,
               error_rate=error_rate,
               marker2file=marker2file, 
               oosp=oosp, 
               tmp_dir=tmp_dir,
               quiet=quiet,
               samtools_exe=samtools_exe,
               bcftools_exe=bcftools_exe)


def save2file(tmp_file, ofn):
    logger.debug('save %s'%ofn)
    with open(ofn, 'w') as ofile:
        for line in tmp_file:
            ofile.write(line)
    tmp_file.seek(0)



def sam2markers(
                sam_file, 
                ofn_bam_sorted_prefix,
                min_align_score=None,
                min_base_quality=30,
                min_read_depth=0,
                error_rate=0.01,
                marker2file=None, 
                oosp=None, 
                tmp_dir='tmp',
                quiet=False,
                samtools_exe='samtools',
                bcftools_exe='bcftools'):
    '''
    Compute the consensus markers in a sample from a sam content.

    :param sam_file: a file name, a file object or subprocess.Popen object
    containing the content of a sam file.
    :param marker2file: the file name to store the consensus genomes.
    :param ofn_bam_sorted_prefix: the bam sorted file prefix
    :param oosp: an instance of ooSubprocess for running a pipe
    :returns: if marker2file==None, return the dictionary of the consensus
    genomes. Otherwise, save the result in fasta format to the file declared by
    marker2file
    '''

    if quiet:
        error_pipe = open(os.devnull, 'w')
    else:
        error_pipe = None

    # sam content to file object
    if oosp is None:
        oosp = ooSubprocess.ooSubprocess()

    if type(sam_file) == str:
        p1 = oosp.chain(
                        'dump_file.py', 
                        args=['--input_file', sam_file],
                        stderr=error_pipe)
    else:
        p1 = sam_file
    
    # filter sam
    if min_align_score == None:
        p1_filtered = p1
    else:
        p1_filtered = oosp.chain('sam_filter.py',
                                 args=['--min_align_score',
                                       str(min_align_score)],
                                 in_pipe=p1,
                                 stderr=error_pipe)
    # sam to bam
    p2 = oosp.chain(
                    samtools_exe, 
                    args=['view', '-bS', '-'], 
                    in_pipe=p1_filtered,
                    stderr=error_pipe)

    # sort bam
    tmp_fns = glob.glob('%s*'%ofn_bam_sorted_prefix)
    for tmp_fn in tmp_fns:
        os.remove(tmp_fn)
    p3 = oosp.chain(
                    samtools_exe, 
                    args=['sort', '-o', '-', ofn_bam_sorted_prefix], 
                    in_pipe=p2,
                    stderr=error_pipe)

    # extract polimorphic information 
    marker2seq = defaultdict(dict)
    pysam.index(p3.name)
    samfile = pysam.AlignmentFile(p3.name)
    for pileupcolumn in samfile.pileup():
        rname = samfile.getrname(pileupcolumn.reference_id)
        pileup = defaultdict(int)
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:  # query position is None if is_del or is_refskip is set.
                b = pileupread.alignment.query_sequence[pileupread.query_position]
                q = pileupread.alignment.query_qualities[pileupread.query_position]
                if q >= min_base_quality:
                    pileup[b] += 1
        if len(pileup) and sum(pileup.values()) > min_read_depth:
            f = float(max(pileup.values())) / sum(pileup.values())
            p = stats.binom.cdf(max(pileup.values()), sum(pileup.values()), 1.0 - error_rate)
            freq = (f, sum(pileup.values()), p)
        else:
            freq = (0.0, 0.0, 0.0)
        if 'freq' not in marker2seq[rname]:
            marker2seq[rname]['freq'] = {}
        marker2seq[rname]['freq'][pileupcolumn.pos] = freq
    samfile.close()
    os.remove(p3.name + '.bai')

    # bam to mpileup
    p3.seek(0)
    p4 = oosp.chain(
                    samtools_exe, 
                    args=['mpileup', '-u', '-'], 
                    in_pipe=p3,
                    stderr=error_pipe)

    # mpileup to vcf
    p5 = oosp.chain(
                    bcftools_exe, 
                    args=['view', '-c', '-g', '-p', '1.1', '-'], 
                    in_pipe=p4,
                    stderr=error_pipe)
                    #stderr=open(os.devnull, 'w'))

    # fix AF1=0
    p6 = oosp.chain(
                    'fix_AF1.py', 
                    args=['--input_file', '-'], 
                    in_pipe=p5,
                    stderr=error_pipe)

    # vcf to fastq
    p7 = oosp.chain(
                      'vcfutils.pl', 
                      args=['vcf2fq'], 
                      in_pipe=p6,
                      get_out_pipe=True,
                      stderr=error_pipe,
                      stop=True)

    try:
        for rec in SeqIO.parse(p7, 'fastq'):
            marker2seq[rec.name]['seq'] = str(rec.seq).upper()
        marker2seq = dict(marker2seq)
    except Exception as e:
        logger.error("sam2markers failed on file " + sam_file)
        raise 

    if type(p1) == file:
        p1.close()

    if marker2file:
        with open(marker2file, 'wb') as ofile:
            msgpack.dump(marker2seq, ofile)

    return marker2seq




@trace_unhandled_exceptions
def run_sample(args_list):
    ifn_sample = args_list[0]
    args = args_list[1]
    base_name = ooSubprocess.splitext(ifn_sample)[0]
    output_prefix = os.path.join(args['output_dir'], base_name)
    if args['sam2file_ext'] != None:
        sam2file = output_prefix + args['sam2file_ext']
    else:
        sam2file = None
    marker2file = output_prefix + args['marker2file_ext']
    if args['input_type'] == 'fastq':
        sample2markers(
                    ifn_sample=ifn_sample, 
                    min_read_len=args['min_read_len'],
                    min_align_score=args['min_align_score'],
                    min_base_quality=args['min_base_quality'],
                    min_read_depth=args['min_read_depth'],
                    error_rate=args['error_rate'],
                    ifn_markers=args['ifn_markers'], 
                    index_path=args['index_path'],
                    nprocs=args['nprocs'],
                    sam2file=sam2file,
                    marker2file=marker2file, 
                    tmp_dir=args['output_dir'],
                    quiet=args['quiet'],
                    samtools_exe=args['samtools_exe'],
                    bcftools_exe=args['bcftools_exe'])
    else:
        ofn_bam_sorted_prefix = os.path.join(
                            args['output_dir'],
                            os.path.basename(ifn_sample) + '.bam.sorted')
        sam2markers(
                    sam_file=ifn_sample, 
                    ofn_bam_sorted_prefix=ofn_bam_sorted_prefix,
                    min_align_score=args['min_align_score'],
                    min_base_quality=args['min_base_quality'],
                    min_read_depth=args['min_read_depth'],
                    error_rate=args['error_rate'],
                    marker2file=marker2file,
                    quiet=args['quiet'],
                    samtools_exe=args['samtools_exe'],
                    bcftools_exe=args['bcftools_exe'])
    return 0




def compute_polymorphic_sites(sample2pileup, ifn_alignment):
    return




def main(args):
    ooSubprocess.mkdir(args['output_dir'])
    manager = multiprocessing.Manager()

    if args['input_type'] == 'fastq':
        index_path = build_bowtie2db(args['ifn_markers'], args['output_dir'])
        args['index_path'] = index_path

    args_list = []
    for ifn_sample in args['ifn_samples']:
        args_list.append([ifn_sample, args])

    #ooSubprocess.parallelize(run_sample, args_list, args['nprocs'])
    pool = multiprocessing.Pool(args['nprocs'])
    results = []
    for a in args_list:
        r = pool.apply_async(run_sample, [a])
        results.append(r)
    for r in results:
        try:
            r.get()
        except Exception as e:
            print e


        
if __name__ == "__main__":
    args = read_params()
    main(args)
