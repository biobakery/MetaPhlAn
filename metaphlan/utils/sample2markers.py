#!/usr/bin/env python
__author__ = ('Aitor Blanco Miguez (aitor.blancomiguez@unitn.it), '
              'Duy Tin Truong (duytin.truong@unitn.it), '
              'Francesco Asnicar (f.asnicar@unitn.it), '
              'Moreno Zolfo (moreno.zolfo@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)')
__version__ = '3.0'
__date__ = '21 Feb 2020'

import sys
try:
    from .util_fun import error
except ImportError:
    from util_fun import error

if sys.version_info[0] < 3:
    error("StrainPhlAn 3.0 requires Python 3, your current Python version is {}.{}.{}"
                    .format(sys.version_info[0], sys.version_info[1], 
                        sys.version_info[2]), exit=True)

import os, time, shutil, pickle
import subprocess as sb
import argparse as ap
from cmseq import cmseq
try:
    from .external_exec import samtools_sam_to_bam, samtools_sort_bam_v1, decompress_bz2
    from .util_fun import info, optimized_dump, get_breath
    from .parallelisation import execute_pool
except ImportError:
    from external_exec import samtools_sam_to_bam, samtools_sort_bam_v1, decompress_bz2
    from util_fun import info, optimized_dump, get_breath
    from parallelisation import execute_pool


"""
Reads and parses the command line arguments of the script.

:returns: the parsed arguments
"""
def read_params():
    p = ap.ArgumentParser(description="")
    p.add_argument('-i', '--input', type=str, 
                   nargs='+', default=[],
                   help="The input samples as SAM or BAM files")
    p.add_argument('--sorted', action='store_true', default=False,
                   help="Whether the BAM input files are sorted. Default false")
    p.add_argument('-f', '--input_format', type=str, default="bz2",
                   help="The input samples format {bam, sam, bz2}. Default bz2")
    p.add_argument('-o', '--output_dir', type=str, default=None,
                   help="The output directory")
    p.add_argument('-b', '--breadth_threshold', type=int, default=80,
                   help="The breadth of coverage threshold for the consensus markers. Default 80 (%%)")
    p.add_argument('-n', '--nprocs', type=int, default=1,
                   help="The number of threads to execute the script")
    
    return p.parse_args()


"""
Checks the mandatory command line arguments of the script.

:param args: the arguments of the script
:returns: the checked args
"""
def check_params(args):
    if not args.input:
        error('-i (or --input) must be specified', exit=True, 
            init_new_line=True)
    elif not args.input_format:
        error('-f (or --input_format) must be specified', exit=True, 
            init_new_line=True)
    elif not args.output_dir:
        error('-o (or --output_dir) must be specified', exit=True, 
            init_new_line=True)
    elif args.input_format.lower() != "bam" and args.input_format.lower() != "sam" and args.input_format.lower() != "bz2":
        error('The input format must be SAM, BAM, or compressed in BZ2 format', 
            exit=True, init_new_line=True)
    else:
        check_input_files(args.input, args.input_format)
    if not args.output_dir.endswith('/'):
        args.output_dir += '/'    
    return args


"""
Checks the input sample files

:param input: the input files
:param input_format: the format of the input files
"""
def check_input_files(input, input_format):
    for s in input:
        _, extension = os.path.splitext(s)
        if not os.path.exists(s):
            error('The input file \"'+s+'\" does not exist', exit=True, 
                init_new_line=True)
        elif not input_format.lower() == extension[1:].lower():
            error('The the input file \"'+s+'\" must be in \"'+
                input_format.upper()+'\" format',
                exit=True, init_new_line=True)
    return True


#ToDo: Check CMSeq and samtools
"""
Checks the mandatory programs to execute of the script.

"""
def check_dependencies():
    try:
        # sb.check_call(["samtools", "tview"], stdout=sb.DEVNULL, stderr=sb.DEVNULL)
        sb.check_call(["bzip2", "--help"], stdout=sb.DEVNULL, stderr=sb.DEVNULL)
    except Exception as e:
        error('Program not installed or not present in the system path\n'+str(e), 
            init_new_line=True, exit=True)


"""
Decompressed SAM.BZ2 files

:param input: the list of samples as BZ2 files
:param tmp_dir: the temporal output directory
:param nprocs: the number of threads to use
:returns: the list of decompressed files
"""
def decompress_from_bz2(input, tmp_dir, nprocs):
    decompressed = []
    decompressed_format = []
    results = execute_pool(((decompress_bz2_file, i, tmp_dir) for i in input), 
        nprocs)
    for r in results:
        decompressed.append(r[0])
        decompressed_format.append(r[1])   

    if decompressed_format[1:] == decompressed_format[:-1]:
        if decompressed_format[0][1:].lower() == "sam":
            return decompressed, "sam"
        elif decompressed_format[0][1:].lower() == "bam":
            return decompressed, "bam"
        else:
            error("Decompressed files are not in SAM or BAM format",
                exit=True, init_new_line=True)
    else:
        error("Decompressed files have different formats",
            exit=True, init_new_line=True)


"""
Decompress a BZ2 file and returns the decompressed file and the 
decompressed file format

:param input: the input BZ2 file
:param output_dir: the output directory
:returns: the decompressed file and the decompressed file format
"""
def decompress_bz2_file(input, output_dir):
    decompressed_file = decompress_bz2(input, output_dir)
    _, e = os.path.splitext(decompressed_file)
    return decompressed_file, e


"""
Convert input sample files to sorted BAMs

:param input: the samples as SAM or BAM files
:param sorted: whether the BAM files are sorted
:param input_format: format of the sample files [bam or sam]
:param tmp_dir: the temporal output directory
:param nprocs: number of threads to use in the execution
:returns: the new list of input BAM files
"""
def convert_inputs(input, sorted, input_format, tmp_dir, nprocs):
    if input_format.lower() == "bz2":
        info("Decompressing samples...\n", init_new_line=True)
        input, input_format = decompress_from_bz2(input, tmp_dir, nprocs)
        info("Done.")
    if input_format.lower() == "sam":
        info("Converting samples to BAM format...\n", init_new_line=True)
        input = execute_pool(((samtools_sam_to_bam, i, tmp_dir) for i in input), 
            nprocs)
        info("Done.")
        info("Sorting BAM samples...\n", init_new_line=True)
        input = execute_pool(((samtools_sort_bam_v1, i, tmp_dir) for i in input), 
            nprocs)
        info("Done.")
    elif sorted == False:        
        info("Sorting BAM samples...\n", init_new_line=True)
        input = execute_pool(((samtools_sort_bam_v1, i, tmp_dir) for i in input), 
            nprocs)
        info("Done.")
    
    return input


#ToDo: minimumReadsAlignning as command line parameter
"""
Gets the markers for each sample and writes the Pickle files

:param input: the samples as sorted BAM files
:param output_dir: the output directory
:param breath_threshold: the breath threshold for the consensus markers
:param nprocs: number of threads to use in the execution
"""
def execute_cmseq(input, output_dir, breath_threshold, nprocs):
    info("Getting consensus markers from samples...", init_new_line=True)
    for i in input:
        info("\tProcessing sample: "+i, init_new_line=True)
        n, _ = os.path.splitext(os.path.basename(i))
        consensus = []
        collection = cmseq.BamFile(i, index=True, minimumReadsAligning=8)
        results = collection.parallel_reference_free_consensus(ncores=nprocs, consensus_rule=cmseq.BamContig.majority_rule_polymorphicLoci)
        for c, seq in results:
            breath = get_breath(seq)
            if breath > breath_threshold:
                consensus.append({"marker":c, "breath":breath, "sequence":seq})
        markers_pkl = open(output_dir+n+'.pkl', 'wb')
        optimized_dump(markers_pkl, consensus)
        info("\tDone.", init_new_line=True)
    info("Done.", init_new_line=True)


"""
Gets the clade-specific markers from a list of aligned samples in 
SAM or BAM format and writes the results in Pickles files in the
user-selected output directory

:param input: the samples as SAM or BAM files
:param sorted: whether the BAM files are sorted
:param input_format: format of the sample files [bam or sam]
:param output_dir: the output directory
:param breath_threshold: the breath threshold for the consensus markers
:param nprocs: number of threads to use in the execution
"""
def samples_to_markers(input, sorted, input_format, output_dir, breath_threshold, nprocs):
    tmp_dir = output_dir+'tmp/'
    try:
        os.mkdir(tmp_dir)
    except Exception as e:
        error('Folder \"'+tmp_dir+'\" already exists!\n'+str(e), exit=True,
            init_new_line=True)
    
    input = convert_inputs(input, sorted, input_format, tmp_dir, nprocs)
    execute_cmseq(input, output_dir, breath_threshold, nprocs)        
    
    shutil.rmtree(tmp_dir, ignore_errors=False, onerror=None)


"""
Main function

:param input: the samples as SAM or BAM files
:param sorted: whether the BAM files are sorted
:param input_format: format of the sample files [bam or sam]
:param output_dir: the output directory
:param breadth_threshold: the breadth threshold for the consensus markers
:param nprocs: number of threads to use in the execution
"""
def main():
    t0 = time.time()
    args = read_params()
    info("Start samples to markers execution")
    check_dependencies()
    args = check_params(args)
    samples_to_markers(args.input, args.sorted, args.input_format, args.output_dir, 
        args.breadth_threshold, args.nprocs)
    exec_time = time.time() - t0
    info("Finish samples to markers execution ("+str(round(exec_time, 2))+
        " seconds): Results are stored at \""+args.output_dir+"\"\n",
         init_new_line=True)

if __name__ == '__main__':
    main()
