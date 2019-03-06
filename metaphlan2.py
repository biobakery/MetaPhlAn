#!/usr/bin/env python3
from __future__ import with_statement

# ==============================================================================
# MetaPhlAn v2.x: METAgenomic PHyLogenetic ANalysis for taxonomic classification
#                 of metagenomic data
#
# Authors: Nicola Segata (nicola.segata@unitn.it),
#          Duy Tin Truong,
#          Francesco Asnicar (f.asnicar@unitn.it)
#
# Please type "./metaphlan2.py -h" for usage help
#
# ==============================================================================

__author__ = ('Nicola Segata (nicola.segata@unitn.it), '
              'Duy Tin Truong, '
              'Francesco Asnicar (f.asnicar@unitn.it)'
              'Francesco Beghini (francesco.beghini@unitn.it)')
__version__ = '2.8'
__date__ = '20 Feb 2019'


import sys
import os
import stat
import re
import time
import tarfile
# from binascii import b2a_uu
try:
    import numpy as np
except ImportError:
    sys.stderr.write("Error! numpy python library not detected!!\n")
    sys.exit(1)
import tempfile as tf
import argparse as ap
import subprocess as subp
try:
    from subprocess import DEVNULL  # py3k
except ImportError:
    DEVNULL = open(os.devnull, 'wb')
# import multiprocessing as mp
from collections import defaultdict as defdict
import bz2
import itertools
from distutils.version import LooseVersion
try:
    import cPickle as pickle
except ImportError:
    import pickle
# try to import urllib.request.urlretrieve for python3
try:
    from urllib.request import urlretrieve
except ImportError:
    from urllib import urlretrieve
from glob import glob
import hashlib


# set the location of the database download url
DATABASE_DOWNLOAD = "https://bitbucket.org/biobakery/metaphlan2/downloads/"
# get the directory that contains this script
metaphlan2_script_install_folder = os.path.dirname(os.path.abspath(__file__))
# get the default database folder
DEFAULT_DB_FOLDER = os.path.join(metaphlan2_script_install_folder, "metaphlan_databases")


#**********************************************************************************************
#  Modification of Code :                                                                     *
#  Modified the code so instead of using the current clade IDs, which are numbers, we will    *
#      use the clade_names                                                                    *
#      Users reported the biom output is invalid and also the IDs were changing from run to   *
#      run.                                                                                   *
#  George Weingart    05/22/2017   george.weingart@mail.com                                   *
#**********************************************************************************************



#*************************************************************
#*  Imports related to biom file generation                  *
#*************************************************************
try:
    import biom
    import biom.table
    # import numpy as np  # numpy already imported above
except ImportError:
    sys.stderr.write("Warning! Biom python library not detected!"
                     "\n Exporting to biom format will not work!\n")
try:
    import json
except ImportError:
    sys.stderr.write("Warning! json python library not detected!"
                     "\n Exporting to biom format will not work!\n")


tax_units = "kpcofgst"

if float(sys.version_info[0]) < 3.0:
    def read_and_split(ofn):
        return (l.strip().split('\t') for l in ofn)

    def read_and_split_line(line):
        return line.strip().split('\t')
else:
    def read_and_split(ofn):
        return (l.decode('utf-8').strip().split('\t') for l in ofn)

    def read_and_split_line(line):
        return line.decode('utf-8').strip().split('\t')


def plain_read_and_split(ofn):
    return (l.strip().split('\t') for l in ofn)


def plain_read_and_split_line(l):
    return l.strip().split('\t')


if float(sys.version_info[0]) < 3.0:
    def mybytes(val):
        return val
else:
    def mybytes(val):
        return bytes(val, encoding='utf-8')


def read_params(args):
    p = ap.ArgumentParser( description=
            "DESCRIPTION\n"
            " MetaPhlAn version "+__version__+" ("+__date__+"): \n"
            " METAgenomic PHyLogenetic ANalysis for metagenomic taxonomic profiling.\n\n"
            "AUTHORS: "+__author__+"\n\n"
            "COMMON COMMANDS\n\n"
            " We assume here that metaphlan2.py is in the system path and that mpa_dir bash variable contains the\n"
            " main MetaPhlAn folder. Also BowTie2 should be in the system path with execution and read\n"
            " permissions, and Perl should be installed)\n\n"

            "\n========== MetaPhlAn 2 clade-abundance estimation ================= \n\n"
            "The basic usage of MetaPhlAn 2 consists in the identification of the clades (from phyla to species and \n"
            "strains in particular cases) present in the metagenome obtained from a microbiome sample and their \n"
            "relative abundance. This correspond to the default analysis type (--analysis_type rel_ab).\n\n"

            "*  Profiling a metagenome from raw reads:\n"
            "$ metaphlan2.py metagenome.fastq --input_type fastq\n\n"

            "*  You can take advantage of multiple CPUs and save the intermediate BowTie2 output for re-running\n"
            "   MetaPhlAn extremely quickly:\n"
            "$ metaphlan2.py metagenome.fastq --bowtie2out metagenome.bowtie2.bz2 --nproc 5 --input_type fastq\n\n"

            "*  If you already mapped your metagenome against the marker DB (using a previous MetaPhlAn run), you\n"
            "   can obtain the results in few seconds by using the previously saved --bowtie2out file and \n"
            "   specifying the input (--input_type bowtie2out):\n"
            "$ metaphlan2.py metagenome.bowtie2.bz2 --nproc 5 --input_type bowtie2out\n\n"

            "*  You can also provide an externally BowTie2-mapped SAM if you specify this format with \n"
            "   --input_type. Two steps: first apply BowTie2 and then feed MetaPhlAn2 with the obtained sam:\n"
            "$ bowtie2 --sam-no-hd --sam-no-sq --no-unal --very-sensitive -S metagenome.sam -x ${mpa_dir}/db_v20/mpa_v20_m200 -U metagenome.fastq\n"
            "$ metaphlan2.py metagenome.sam --input_type sam > profiled_metagenome.txt\n\n"

            # "*  Multiple alternative ways to pass the input are also available:\n"
            # "$ cat metagenome.fastq | metaphlan2.py --input_type fastq \n"
            # "$ tar xjf metagenome.tar.bz2 --to-stdout | metaphlan2.py --input_type fastq \n"
            # "$ metaphlan2.py --input_type fastq < metagenome.fastq\n"
            # "$ metaphlan2.py --input_type fastq <(bzcat metagenome.fastq.bz2)\n"
            # "$ metaphlan2.py --input_type fastq <(zcat metagenome_1.fastq.gz metagenome_2.fastq.gz)\n\n"

            "*  We can also natively handle paired-end metagenomes, and, more generally, metagenomes stored in \n"
            "  multiple files (but you need to specify the --bowtie2out parameter):\n"
            "$ metaphlan2.py metagenome_1.fastq,metagenome_2.fastq --bowtie2out metagenome.bowtie2.bz2 --nproc 5 --input_type fastq\n\n"
            "\n------------------------------------------------------------------- \n \n\n"


            "\n========== Marker level analysis ============================ \n\n"
            "MetaPhlAn 2 introduces the capability of charachterizing organisms at the strain level using non\n"
            "aggregated marker information. Such capability comes with several slightly different flavours and \n"
            "are a way to perform strain tracking and comparison across multiple samples.\n"
            "Usually, MetaPhlAn 2 is first ran with the default --analysis_type to profile the species present in\n"
            "the community, and then a strain-level profiling can be performed to zoom-in into specific species\n"
            "of interest. This operation can be performed quickly as it exploits the --bowtie2out intermediate \n"
            "file saved during the execution of the default analysis type.\n\n"

            "*  The following command will output the abundance of each marker with a RPK (reads per kil-base) \n"
            "   higher 0.0. (we are assuming that metagenome_outfmt.bz2 has been generated before as \n"
            "   shown above).\n"
            "$ metaphlan2.py -t marker_ab_table metagenome_outfmt.bz2 --input_type bowtie2out > marker_abundance_table.txt\n"
            "   The obtained RPK can be optionally normalized by the total number of reads in the metagenome \n"
            "   to guarantee fair comparisons of abundances across samples. The number of reads in the metagenome\n"
            "   needs to be passed with the '--nreads' argument\n\n"

            "*  The list of markers present in the sample can be obtained with '-t marker_pres_table'\n"
            "$ metaphlan2.py -t marker_pres_table metagenome_outfmt.bz2 --input_type bowtie2out > marker_abundance_table.txt\n"
            "   The --pres_th argument (default 1.0) set the minimum RPK value to consider a marker present\n\n"

            "*  The list '-t clade_profiles' analysis type reports the same information of '-t marker_ab_table'\n"
            "   but the markers are reported on a clade-by-clade basis.\n"
            "$ metaphlan2.py -t clade_profiles metagenome_outfmt.bz2 --input_type bowtie2out > marker_abundance_table.txt\n\n"

            "*  Finally, to obtain all markers present for a specific clade and all its subclades, the \n"
            "   '-t clade_specific_strain_tracker' should be used. For example, the following command\n"
            "   is reporting the presence/absence of the markers for the B. fragulis species and its strains\n"
            "   the optional argument --min_ab specifies the minimum clade abundance for reporting the markers\n\n"
            "$ metaphlan2.py -t clade_specific_strain_tracker --clade s__Bacteroides_fragilis metagenome_outfmt.bz2 --input_type bowtie2out > marker_abundance_table.txt\n"

            "\n------------------------------------------------------------------- \n\n"
            "",
            formatter_class=ap.RawTextHelpFormatter,
            add_help=False )
    arg = p.add_argument

    arg( 'inp', metavar='INPUT_FILE', type=str, nargs='?', default=None, help=
         "the input file can be:\n"
         "* a fastq file containing metagenomic reads\n"
         "OR\n"
         "* a BowTie2 produced SAM file. \n"
         "OR\n"
         "* an intermediary mapping file of the metagenome generated by a previous MetaPhlAn run \n"
         "If the input file is missing, the script assumes that the input is provided using the standard \n"
         "input, or named pipes.\n"
         "IMPORTANT: the type of input needs to be specified with --input_type" )

    arg( 'output', metavar='OUTPUT_FILE', type=str, nargs='?', default=None,
         help= "the tab-separated output file of the predicted taxon relative abundances \n"
               "[stdout if not present]")


    g = p.add_argument_group('Required arguments')
    arg = g.add_argument
    input_type_choices = ['fastq','fasta','multifasta','multifastq','bowtie2out','sam']
    arg( '--input_type', choices=input_type_choices, required = '--install' not in args, help =
         "set whether the input is the multifasta file of metagenomic reads or \n"
         "the SAM file of the mapping of the reads against the MetaPhlAn db.\n"
         "[default 'automatic', i.e. the script will try to guess the input format]\n" )

    g = p.add_argument_group('Mapping arguments')
    arg = g.add_argument
    arg('--mpa_pkl', type=str, default=None,
        help="The metadata pickled MetaPhlAn file [deprecated]")

    arg('--bowtie2db', metavar="METAPHLAN_BOWTIE2_DB", type=str, default=DEFAULT_DB_FOLDER,
        help=("The BowTie2 database file of the MetaPhlAn database. Used if "
              "--input_type is fastq, fasta, multifasta, or multifastq [default "+DEFAULT_DB_FOLDER+"]\n"))

    INDEX = 'v25_CHOCOPhlAn_0.2'
    arg('-x', '--index', type=str, default='v25_CHOCOPhlAn_0.2',
        help=("Specify the id of the database version to use. If the database\n"
              "files are not found on the local MetaPhlAn2 installation they\n"
              "will be automatically downloaded [default "+INDEX+"]\n"))

    bt2ps = ['sensitive', 'very-sensitive', 'sensitive-local', 'very-sensitive-local']
    arg('--bt2_ps', metavar="BowTie2 presets", default='very-sensitive',
        choices=bt2ps, help="Presets options for BowTie2 (applied only when a "
                            "multifasta file is provided)\n"
                            "The choices enabled in MetaPhlAn are:\n"
                            " * sensitive\n"
                            " * very-sensitive\n"
                            " * sensitive-local\n"
                            " * very-sensitive-local\n"
                            "[default very-sensitive]\n")
    arg('--bowtie2_exe', type=str, default=None,
        help='Full path and name of the BowTie2 executable. This option allows'
             'MetaPhlAn to reach the executable even when it is not in the '
             'system PATH or the system PATH is unreachable')
    arg('--bowtie2_build', type=str, default='bowtie2-build',
        help="Full path to the bowtie2-build command to use, deafult assumes "
             "that 'bowtie2-build is present in the system path")
    arg('--bowtie2out', metavar="FILE_NAME", type=str, default=None,
        help="The file for saving the output of BowTie2")
    arg('--no_map', action='store_true',
        help="Avoid storing the --bowtie2out map file")
    arg('--tmp_dir', metavar="", default=None, type=str,
        help="The folder used to store temporary files [default is the OS "
             "dependent tmp dir]")

    g = p.add_argument_group('Post-mapping arguments')
    arg = g.add_argument
    stat_choices = ['avg_g','avg_l','tavg_g','tavg_l','wavg_g','wavg_l','med']
    arg( '--tax_lev', metavar='TAXONOMIC_LEVEL', type=str,
         choices='a'+tax_units, default='a', help =
         "The taxonomic level for the relative abundance output:\n"
         "'a' : all taxonomic levels\n"
         "'k' : kingdoms\n"
         "'p' : phyla only\n"
         "'c' : classes only\n"
         "'o' : orders only\n"
         "'f' : families only\n"
         "'g' : genera only\n"
         "'s' : species only\n"
         "[default 'a']" )
    arg( '--min_cu_len', metavar="", default="2000", type=int, help =
         "minimum total nucleotide length for the markers in a clade for\n"
         "estimating the abundance without considering sub-clade abundances\n"
         "[default 2000]\n"   )
    arg( '--min_alignment_len', metavar="", default=None, type=int, help =
         "The sam records for aligned reads with the longest subalignment\n"
         "length smaller than this threshold will be discarded.\n"
         "[default None]\n"   )
    arg( '--ignore_viruses', action='store_true', help=
         "Do not profile viral organisms" )
    arg( '--ignore_eukaryotes', action='store_true', help=
         "Do not profile eukaryotic organisms" )
    arg( '--ignore_bacteria', action='store_true', help=
         "Do not profile bacterial organisms" )
    arg( '--ignore_archaea', action='store_true', help=
         "Do not profile archeal organisms" )
    arg( '--stat_q', metavar="", type = float, default=0.1, help =
         "Quantile value for the robust average\n"
         "[default 0.1]"   )
    arg( '--ignore_markers', type=str, default = None, help =
         "File containing a list of markers to ignore. \n")
    arg( '--avoid_disqm', action="store_true", help =
         "Deactivate the procedure of disambiguating the quasi-markers based on the \n"
         "marker abundance pattern found in the sample. It is generally recommended \n"
         "to keep the disambiguation procedure in order to minimize false positives\n")
    arg( '--stat', metavar="", choices=stat_choices, default="tavg_g", type=str, help =
         "EXPERIMENTAL! Statistical approach for converting marker abundances into clade abundances\n"
         "'avg_g'  : clade global (i.e. normalizing all markers together) average\n"
         "'avg_l'  : average of length-normalized marker counts\n"
         "'tavg_g' : truncated clade global average at --stat_q quantile\n"
         "'tavg_l' : trunated average of length-normalized marker counts (at --stat_q)\n"
         "'wavg_g' : winsorized clade global average (at --stat_q)\n"
         "'wavg_l' : winsorized average of length-normalized marker counts (at --stat_q)\n"
         "'med'    : median of length-normalized marker counts\n"
         "[default tavg_g]"   )

    arg = p.add_argument



    g = p.add_argument_group('Additional analysis types and arguments')
    arg = g.add_argument
    analysis_types = ['rel_ab', 'rel_ab_w_read_stats', 'reads_map', 'clade_profiles', 'marker_ab_table', 'marker_counts', 'marker_pres_table', 'clade_specific_strain_tracker']
    arg( '-t', metavar='ANALYSIS TYPE', type=str, choices = analysis_types,
         default='rel_ab', help =
         "Type of analysis to perform: \n"
         " * rel_ab: profiling a metagenomes in terms of relative abundances\n"
         " * rel_ab_w_read_stats: profiling a metagenomes in terms of relative abundances and estimate the number of reads comming from each clade.\n"
         " * reads_map: mapping from reads to clades (only reads hitting a marker)\n"
         " * clade_profiles: normalized marker counts for clades with at least a non-null marker\n"
         " * marker_ab_table: normalized marker counts (only when > 0.0 and normalized by metagenome size if --nreads is specified)\n"
         " * marker_counts: non-normalized marker counts [use with extreme caution]\n"
         " * marker_pres_table: list of markers present in the sample (threshold at 1.0 if not differently specified with --pres_th\n"
         "[default 'rel_ab']" )
    arg( '--nreads', metavar="NUMBER_OF_READS", type=int, default = None, help =
         "The total number of reads in the original metagenome. It is used only when \n"
         "-t marker_table is specified for normalizing the length-normalized counts \n"
         "with the metagenome size as well. No normalization applied if --nreads is not \n"
         "specified" )
    arg( '--pres_th', metavar="PRESENCE_THRESHOLD", type=int, default = 1.0, help =
         'Threshold for calling a marker present by the -t marker_pres_table option' )
    arg( '--clade', metavar="", default=None, type=str, help =
         "The clade for clade_specific_strain_tracker analysis\n"  )
    arg( '--min_ab', metavar="", default=0.1, type=float, help =
         "The minimum percentage abundace for the clade in the clade_specific_strain_tracker analysis\n"  )

    g = p.add_argument_group('Output arguments')
    arg = g.add_argument
    arg( '-o', '--output_file',  metavar="output file", type=str, default=None, help =
         "The output file (if not specified as positional argument)\n")
    arg('--sample_id_key',  metavar="name", type=str, default="#SampleID",
        help =("Specify the sample ID key for this analysis."
               " Defaults to '#SampleID'."))
    arg('--sample_id',  metavar="value", type=str,
        default="Metaphlan2_Analysis",
        help =("Specify the sample ID for this analysis."
               " Defaults to 'Metaphlan2_Analysis'."))
    arg( '-s', '--samout', metavar="sam_output_file",
        type=str, default=None, help="The sam output file\n")

    arg( '--legacy-output', action='store_true', help="Old two columns output\n")
    #*************************************************************
    #* Parameters related to biom file generation                *
    #*************************************************************
    arg( '--biom', '--biom_output_file',  metavar="biom_output", type=str, default=None, help =
         "If requesting biom file output: The name of the output file in biom format \n")

    arg( '--mdelim', '--metadata_delimiter_char',  metavar="mdelim", type=str, default="|", help =
         "Delimiter for bug metadata: - defaults to pipe. e.g. the pipe in k__Bacteria|p__Proteobacteria \n")
    #*************************************************************
    #* End parameters related to biom file generation            *
    #*************************************************************

    g = p.add_argument_group('Other arguments')
    arg = g.add_argument
    arg('--nproc', metavar="N", type=int, default=4,
        help="The number of CPUs to use for parallelizing the mapping [default 4]")
    arg('--install', action='store_true',
        help="Only checks if the MetaPhlAn2 DB is installed and installs it if not. All other parameters are ignored.")
    arg('--read_min_len', type=int, default=70,
        help="Specify the minimum length of the reads to be considered when parsing the input file with "
             "'read_fastx.py' script, default value is 70")
    arg('-v', '--version', action='version',
        version="MetaPhlAn version {} ({})".format(__version__, __date__),
        help="Prints the current MetaPhlAn version and exit")
    arg("-h", "--help", action="help", help="show this help message and exit")

    return vars(p.parse_args())


def byte_to_megabyte(byte):
    """
    Convert byte value to megabyte
    """

    return byte / (1024.0**2)


class ReportHook():
    def __init__(self):
        self.start_time = time.time()

    def report(self, blocknum, block_size, total_size):
        """
        Print download progress message
        """

        if blocknum == 0:
            self.start_time = time.time()
            if total_size > 0:
                sys.stderr.write("Downloading file of size: {:.2f} MB\n"
                                 .format(byte_to_megabyte(total_size)))
        else:
            total_downloaded = blocknum * block_size
            status = "{:3.2f} MB ".format(byte_to_megabyte(total_downloaded))

            if total_size > 0:
                percent_downloaded = total_downloaded * 100.0 / total_size
                # use carriage return plus sys.stderr to overwrite stderr
                download_rate = total_downloaded / (time.time() - self.start_time)
                estimated_time = (total_size - total_downloaded) / download_rate
                estimated_minutes = int(estimated_time / 60.0)
                estimated_seconds = estimated_time - estimated_minutes * 60.0
                status += ("{:3.2f} %  {:5.2f} MB/sec {:2.0f} min {:2.0f} sec "
                           .format(percent_downloaded,
                                   byte_to_megabyte(download_rate),
                                   estimated_minutes, estimated_seconds))

            status += "        \r"
            sys.stderr.write(status)


def download(url, download_file):
    """
    Download a file from a url
    """

    if not os.path.isfile(download_file):
        try:
            sys.stderr.write("\nDownloading " + url + "\n")
            file, headers = urlretrieve(url, download_file,
                                        reporthook=ReportHook().report)
        except EnvironmentError:
            sys.stderr.write("\nWarning: Unable to download " + url + "\n")
    else:
        sys.stderr.write("\nFile {} already present!\n".format(download_file))


def download_unpack_tar(url, download_file_name, folder, bowtie2_build, nproc):
    """
    Download the url to the file and decompress into the folder
    """

    # Create the folder if it does not already exist
    if not os.path.isdir(folder):
        try:
            os.makedirs(folder)
        except EnvironmentError:
            sys.exit("ERROR: Unable to create folder for database install: " + folder)

    # Check the directory permissions
    if not os.access(folder, os.W_OK):
        sys.exit("ERROR: The directory is not writeable: " + folder + ". "
                 "Please modify the permissions.")

    tar_file = os.path.join(folder, "mpa_" + download_file_name + ".tar")
    url_tar_file = os.path.join(url, "mpa_" + download_file_name + ".tar")
    download(url_tar_file, tar_file)

    # download MD5 checksum
    md5_file = os.path.join(folder, "mpa_" + download_file_name + ".md5")
    url_md5_file = os.path.join(url, "mpa_" + download_file_name + ".md5")
    download(url_md5_file, md5_file)

    md5_md5 = None
    md5_tar = None

    if os.path.isfile(md5_file):
        with open(md5_file) as f:
            for row in f:
                md5_md5 = row.strip().split(' ')[0]
    else:
        sys.stderr.write('File "{}" not found!\n'.format(md5_file))

    # compute MD5 of .tar.bz2
    if os.path.isfile(tar_file):
        hash_md5 = hashlib.md5()

        with open(tar_file, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)

        md5_tar = hash_md5.hexdigest()[:32]
    else:
        sys.stderr.write('File "{}" not found!\n'.format(tar_file))

    if (md5_tar is None) or (md5_md5 is None):
        sys.exit("MD5 checksums not found, something went wrong!")

    # compare checksums
    if md5_tar != md5_md5:
        sys.exit("MD5 checksums do not correspond! If this happens again, you should remove the database files and "
                 "rerun MetaPhlAn2 so they are re-downloaded")

    # untar
    try:
        tarfile_handle = tarfile.open(tar_file)
        tarfile_handle.extractall(path=folder)
        tarfile_handle.close()
    except EnvironmentError:
        sys.stderr.write("Warning: Unable to extract {}.\n".format(tar_file))

    # uncompress sequences
    bz2_file = os.path.join(folder, "mpa_" + download_file_name + ".fna.bz2")
    fna_file = os.path.join(folder, "mpa_" + download_file_name + ".fna")

    if not os.path.isfile(fna_file):
        sys.stderr.write('\n\nDecompressing {} into {}\n'.format(bz2_file, fna_file))

        with open(fna_file, 'wb') as fna_h, bz2.BZ2File(bz2_file, 'rb') as bz2_h:
            for data in iter(lambda: bz2_h.read(100 * 1024), b''):
                fna_h.write(data)

    # build bowtie2 indexes
    if not glob(os.path.join(folder, "mpa_" + download_file_name + "*.bt2")):
        bt2_base = os.path.join(folder, "mpa_" + download_file_name)
        bt2_cmd = [bowtie2_build, '--quiet']

        if nproc > 1:
            bt2_build_output = subp.check_output([bowtie2_build, '--usage'], stderr=subp.STDOUT)

            if 'threads' in str(bt2_build_output):
                bt2_cmd += ['--threads', str(nproc)]

        bt2_cmd += ['-f', fna_file, bt2_base]

        sys.stderr.write('\nBuilding Bowtie2 indexes\n')

        try:
            subp.check_call(bt2_cmd)
        except Exception as e:
            sys.stderr.write("Fatal error running '{}'\nError message: '{}'\n\n".format(' '.join(bt2_cmd), e))
            sys.exit(1)

    sys.stderr.write('Removing uncompress database {}\n'.format(fna_file))
    os.remove(fna_file)


def check_and_install_database(index, bowtie2_db, bowtie2_build, nproc):
    """ Check if the database is installed, if not download and install """

    if len(glob(os.path.join(bowtie2_db, "mpa_{}*".format(index)))) >= 7:
        return

    # download the tar archive and decompress
    sys.stderr.write("\nDownloading MetaPhlAn2 database\nPlease note due to "
                     "the size this might take a few minutes\n")
    download_unpack_tar(DATABASE_DOWNLOAD, index, bowtie2_db, bowtie2_build, nproc)
    sys.stderr.write("\nDownload complete\n")


def set_mapping_arguments(index, bowtie2_db):
    mpa_pkl = 'mpa_pkl'
    bowtie2db = 'bowtie2db'

    if os.path.isfile(os.path.join(bowtie2_db, "mpa_{}.pkl".format(index))):
        mpa_pkl = os.path.join(bowtie2_db, "mpa_{}.pkl".format(index))

    if glob(os.path.join(bowtie2_db, "mpa_{}*.bt2".format(index))):
        bowtie2db = os.path.join(bowtie2_db, "mpa_{}".format(index))

    return (mpa_pkl, bowtie2db)


def run_bowtie2(fna_in, outfmt6_out, bowtie2_db, preset, nproc, file_format="multifasta",
                exe=None, samout=None, min_alignment_len=None, read_min_len=0):
    # checking read_fastx.py
    read_fastx = "read_fastx.py"

    try:
        subp.check_call([read_fastx, "-h"], stdout=DEVNULL)
    except Exception as e:
        try:
            read_fastx = os.path.join(os.path.join(os.path.dirname(__file__), "utils"), read_fastx)
            subp.check_call([read_fastx, "-h"], stdout=DEVNULL)
        except Exception as e:
            sys.stderr.write("OSError: fatal error running '{}'. Is it in the system path?\n".format(read_fastx))
            sys.exit(1)

    # checking bowtie2
    try:
        subp.check_call([exe if exe else 'bowtie2', "-h"], stdout=DEVNULL)
    except Exception as e:
        sys.stderr.write('OSError: "{}"\nFatal error running BowTie2. Is BowTie2 in the system path?\n'.format(e))
        sys.exit(1)

    try:
        if fna_in:
            readin = subp.Popen([read_fastx, '-l', str(read_min_len), fna_in], stdout=subp.PIPE, stderr=subp.PIPE)

        else:
            readin = subp.Popen([read_fastx, '-l', str(read_min_len)], stdin=sys.stdin, stdout=subp.PIPE, stderr=subp.PIPE)

        bowtie2_cmd = [exe if exe else 'bowtie2', "--quiet", "--no-unal", "--{}".format(preset),
                       "-S", "-", "-x", bowtie2_db]

        if int(nproc) > 1:
            bowtie2_cmd += ["-p", str(nproc)]

        bowtie2_cmd += ["-U", "-"]  # if not stat.S_ISFIFO(os.stat(fna_in).st_mode) else []

        if file_format == "multifasta":
            bowtie2_cmd += ["-f"]

        p = subp.Popen(bowtie2_cmd, stdout=subp.PIPE, stdin=readin.stdout)
        readin.stdout.close()
        lmybytes, outf = (mybytes, bz2.BZ2File(outfmt6_out, "w")) if outfmt6_out.endswith(".bz2") else (str, open(outfmt6_out, "w"))
        try:
            if samout:
                if samout[-4:] == '.bz2':
                    sam_file = bz2.BZ2File(samout, 'w')
                else:
                    sam_file = open(samout, 'wb')
        except IOError as e:
            sys.stderr.write('IOError: "{}"\nUnable to open sam output file.\n'.format(e))
            sys.exit(1)
        for line in p.stdout:
            if samout:
                sam_file.write(line)

            o = read_and_split_line(line)

            if not o[0].startswith('@'):
                if not o[2].endswith('*'):
                    if ((min_alignment_len is None) or
                            (max([int(x.strip('M')) for x in re.findall(r'(\d*M)', o[5]) if x]) >= min_alignment_len)):
                        outf.write(lmybytes("\t".join([o[0], o[2]]) + "\n"))


        if samout:
            sam_file.close()

        p.communicate()
        
        n_metagenome_reads = int(readin.stderr.readline())
        outf.write(lmybytes('#nreads\t{}'.format(n_metagenome_reads)))
        outf.close()

    except OSError as e:
        sys.stderr.write('OSError: "{}"\nFatal error running BowTie2.\n'.format(e))
        sys.exit(1)
    except ValueError as e:
        sys.stderr.write('ValueError: "{}"\nFatal error running BowTie2.\n'.format(e))
        sys.exit(1)
    except IOError as e:
        sys.stderr.write('IOError: "{}"\nFatal error running BowTie2.\n'.format(e))
        sys.exit(1)

    if p.returncode == 13:
        sys.stderr.write("Permission Denied Error: fatal error running BowTie2."
                         "Is the BowTie2 file in the path with execution and read permissions?\n")
        sys.exit(1)
    elif p.returncode != 0:
        sys.stderr.write("Error while running bowtie2.\n")
        sys.exit(1)
#def guess_input_format( inp_file ):
#    if "," in inp_file:
#        sys.stderr.write( "Sorry, I cannot guess the format of the input, when "
#                          "more than one file is specified. Please set the --input_type parameter \n" )
#        sys.exit(1)
#
#    with open( inp_file ) as inpf:
#        for i,l in enumerate(inpf):
#            line = l.strip()
#            if line[0] == '#': continue
#            if line[0] == '>': return 'multifasta'
#            if line[0] == '@': return 'multifastq'
#            if len(l.split('\t')) == 2: return 'bowtie2out'
#            if i > 20: break
#    return None

class TaxClade:
    min_cu_len = -1
    markers2lens = None
    stat = None
    quantile = None
    avoid_disqm = False

    def __init__( self, name, tax_id, uncl = False):
        self.children, self.markers2nreads = {}, {}
        self.name, self.father = name, None
        self.uncl, self.subcl_uncl = uncl, False
        self.abundance, self.uncl_abundance = None, 0
        self.nreads, self.uncl_nreads = 0, 0
        self.tax_id = tax_id

    def add_child( self, name, tax_id ):
        new_clade = TaxClade( name, tax_id )
        self.children[name] = new_clade
        new_clade.father = self
        return new_clade


    def get_terminals( self ):
        terms = []
        if not self.children:
            return [self]
        for c in self.children.values():
            terms += c.get_terminals()
        return terms

    def get_full_taxids( self ):
        fullname = [self.tax_id]
        cl = self.father
        while cl:
            fullname = [cl.tax_id] + fullname
            cl = cl.father
        return "|".join(fullname[1:])

    def get_full_name( self ):
        fullname = [self.name]
        cl = self.father
        while cl:
            fullname = [cl.name] + fullname
            cl = cl.father
        return "|".join(fullname[1:])

    def get_normalized_counts( self ):
        return [(m,float(n)*1000.0/self.markers2lens[m])
                    for m,n in self.markers2nreads.items()]

    def compute_mapped_reads( self ):
        if self.name.startswith('s__'):
            return self.nreads
        for c in self.children.values():
            self.nreads += c.compute_mapped_reads()
        return self.nreads
        
    def compute_abundance( self ):
        if self.abundance is not None: return self.abundance

        sum_ab = sum([c.compute_abundance() for c in self.children.values()])

        # rat_nreads = sorted([(self.markers2lens[marker], n_reads)
        #                             for marker,n_reads in self.markers2nreads.items()],
        #                                     key = lambda x: x[1])

        rat_nreads, removed = [], []
        for marker, n_reads in sorted(self.markers2nreads.items(),key=lambda x:x[0]):
            misidentified = False

            if not self.avoid_disqm:
                for ext in self.markers2exts[marker]:
                    ext_clade = self.taxa2clades[ext]
                    m2nr = ext_clade.markers2nreads
                    
                    tocladetmp = ext_clade
                    while len(tocladetmp.children) == 1:
                        tocladetmp = list(tocladetmp.children.values())[0]
                        m2nr = tocladetmp.markers2nreads

                    nonzeros = sum([v>0 for v in m2nr.values()])
                    if len(m2nr):
                        if float(nonzeros) / len(m2nr) > 0.33:
                            misidentified = True
                            removed.append( (self.markers2lens[marker],n_reads) )
                            break
            if not misidentified:
                rat_nreads.append( (self.markers2lens[marker],n_reads) )

        if not self.avoid_disqm and len(removed):
            n_rat_nreads = float(len(rat_nreads))
            n_removed = float(len(removed))
            n_tot = n_rat_nreads + n_removed
            n_ripr = 10

            if len(self.get_terminals()) < 2:
                n_ripr = 0

            if "k__Viruses" in self.get_full_name():
                n_ripr = 0

            if n_rat_nreads < n_ripr and n_tot > n_rat_nreads:
                rat_nreads += removed[:n_ripr-int(n_rat_nreads)]


        rat_nreads = sorted(rat_nreads, key = lambda x: x[1])

        rat_v,nreads_v = zip(*rat_nreads) if rat_nreads else ([],[])
        rat, nrawreads, loc_ab = float(sum(rat_v)) or -1.0, sum(nreads_v), 0.0
        quant = int(self.quantile*len(rat_nreads))
        ql,qr,qn = (quant,-quant,quant) if quant else (None,None,0)

        if self.name[0] == 't' and (len(self.father.children) > 1 or "_sp" in self.father.name or "k__Viruses" in self.get_full_name()):
            non_zeros = float(len([n for r,n in rat_nreads if n > 0]))
            nreads = float(len(rat_nreads))
            if nreads == 0.0 or non_zeros / nreads < 0.7:
                self.abundance = 0.0
                return 0.0

        if rat < 0.0:
            pass
        elif self.stat == 'avg_g' or (not qn and self.stat in ['wavg_g','tavg_g']):
            loc_ab = nrawreads / rat if rat >= 0 else 0.0
        elif self.stat == 'avg_l' or (not qn and self.stat in ['wavg_l','tavg_l']):
            loc_ab = np.mean([float(n)/r for r,n in rat_nreads])
        elif self.stat == 'tavg_g':
            wnreads = sorted([(float(n)/r,r,n) for r,n in rat_nreads], key=lambda x:x[0])
            den,num = zip(*[v[1:] for v in wnreads[ql:qr]])
            loc_ab = float(sum(num))/float(sum(den)) if any(den) else 0.0
        elif self.stat == 'tavg_l':
            loc_ab = np.mean(sorted([float(n)/r for r,n in rat_nreads])[ql:qr])
        elif self.stat == 'wavg_g':
            vmin, vmax = nreads_v[ql], nreads_v[qr]
            wnreads = [vmin]*qn+list(nreads_v[ql:qr])+[vmax]*qn
            loc_ab = float(sum(wnreads)) / rat
        elif self.stat == 'wavg_l':
            wnreads = sorted([float(n)/r for r,n in rat_nreads])
            vmin, vmax = wnreads[ql], wnreads[qr]
            wnreads = [vmin]*qn+list(wnreads[ql:qr])+[vmax]*qn
            loc_ab = np.mean(wnreads)
        elif self.stat == 'med':
            loc_ab = np.median(sorted([float(n)/r for r,n in rat_nreads])[ql:qr])

        self.abundance = loc_ab
        if rat < self.min_cu_len and self.children:
            self.abundance = sum_ab
        elif loc_ab < sum_ab:
            self.abundance = sum_ab

        if self.abundance > sum_ab and self.children: # *1.1??
            self.uncl_abundance = self.abundance - sum_ab
        self.subcl_uncl = not self.children and self.name[0] not in tax_units[-2:]

        return self.abundance

    def get_all_abundances( self ):
        ret = [(self.name, self.tax_id, self.abundance)]
        if self.uncl_abundance > 0.0:
            lchild = list(self.children.values())[0].name[:3]
            ret += [(lchild+self.name[3:]+"_unclassified", "", self.uncl_abundance)]
        if self.subcl_uncl and self.name[0] != tax_units[-2]:
            cind = tax_units.index( self.name[0] )
            ret += [(   tax_units[cind+1]+self.name[1:]+"_unclassified","",
                        self.abundance)]
        for c in self.children.values():
            ret += c.get_all_abundances()
        return ret


class TaxTree:
    def __init__( self, mpa, markers_to_ignore = None ): #, min_cu_len ):
        self.root = TaxClade( "root", 0)
        self.all_clades, self.markers2lens, self.markers2clades, self.taxa2clades, self.markers2exts = {}, {}, {}, {}, {}
        TaxClade.markers2lens = self.markers2lens
        TaxClade.markers2exts = self.markers2exts
        TaxClade.taxa2clades = self.taxa2clades

        # clades_txt = ((l.strip().split("|"),n) for l,n in mpa_pkl['taxonomy'].items())
        clades_txt = ((l.strip().split("|"), t.strip().split("|"), n) for l, (t, n) in mpa['taxonomy'].items())

        for clade, taxids, lenc in clades_txt:
            father = self.root
            for i in range(len(clade)):
                clade_lev = clade[i]
                clade_taxid = taxids[i] if i < 7 else clade_lev[3:]
                if not clade_lev in father.children:
                    father.add_child(clade_lev, tax_id=clade_taxid)              
                    self.all_clades[clade_lev] = father.children[clade_lev]
                if clade_lev[0] == "t":
                    self.taxa2clades[clade_lev[3:]] = father
                father = father.children[clade_lev]
                if clade_lev[0] == "t":
                    father.glen = lenc

        def add_lens( node ):
            if not node.children:
                return node.glen
            lens = []
            for c in node.children.values():
                lens.append( add_lens( c ) )
            node.glen = sum(lens) / len(lens)
            return node.glen
        add_lens(self.root)

        # for k,p in mpa_pkl['markers'].items():
        for k, p in mpa['markers'].items():
            if k in markers_to_ignore:
                continue

            self.markers2lens[k] = p['len']
            self.markers2clades[k] = p['clade']
            self.add_reads(k, 0)
            self.markers2exts[k] = p['ext']

    def set_min_cu_len( self, min_cu_len ):
        TaxClade.min_cu_len = min_cu_len

    def set_stat( self, stat, quantile, avoid_disqm = False):
        TaxClade.stat = stat
        TaxClade.quantile = quantile
        TaxClade.avoid_disqm = avoid_disqm

    def add_reads(  self, marker, n,
                    ignore_viruses = False, ignore_eukaryotes = False,
                    ignore_bacteria = False, ignore_archaea = False  ):
        clade = self.markers2clades[marker]
        cl = self.all_clades[clade]
        if ignore_viruses or ignore_eukaryotes or ignore_bacteria or ignore_archaea:
            cn = cl.get_full_name()
            if ignore_viruses and cn.startswith("k__Viruses"):
                return ""
            if ignore_eukaryotes and cn.startswith("k__Eukaryota"):
                return ""
            if ignore_archaea and cn.startswith("k__Archaea"):
                return ""
            if ignore_bacteria and cn.startswith("k__Bacteria"):
                return ""
        # while len(cl.children) == 1:
            # cl = list(cl.children.values())[0]
        cl.markers2nreads[marker] = n
        return (cl.get_full_name(), cl.get_full_taxids(), )


    def markers2counts( self ):
        m2c = {}
        for k,v in self.all_clades.items():
            for m,c in v.markers2nreads.items():
                m2c[m] = c
        return m2c

    def clade_profiles( self, tax_lev, get_all = False  ):
        cl2pr = {}
        for k,v in self.all_clades.items():
            if tax_lev and not k.startswith(tax_lev):
                continue
            prof = v.get_normalized_counts()
            if not get_all and ( len(prof) < 1 or not sum([p[1] for p in prof]) > 0.0 ):
                continue
            cl2pr[v.get_full_name()] = prof
        return cl2pr

    def relative_abundances( self, tax_lev  ):
        clade2abundance_n = dict([(tax_label, clade) for tax_label, clade in self.all_clades.items()
                    if tax_label.startswith("k__") and not clade.uncl])

        clade2abundance, clade2genomelen, clade2est_nreads, tot_ab, tot_reads = {}, {}, {}, 0.0, 0

        for tax_label, clade in clade2abundance_n.items():
            tot_ab += clade.compute_abundance()

        for tax_label, clade in clade2abundance_n.items():
            for clade_label, tax_id, abundance in sorted(clade.get_all_abundances(), key=lambda pars:pars[0]):
                if clade_label[:3] != 't__':
                    if not tax_lev:
                        if clade_label not in self.all_clades:
                            to = tax_units.index(clade_label[0])
                            t = tax_units[to-1]
                            clade_label = t + clade_label.split("_unclassified")[0][1:]
                            tax_id = self.all_clades[clade_label].get_full_taxids()
                            clade_label = self.all_clades[clade_label].get_full_name()
                            spl = clade_label.split("|")
                            clade_label = "|".join(spl+[tax_units[to]+spl[-1][1:]+"_unclassified"])
                            glen = self.all_clades[spl[-1]].glen
                        else:
                            glen = self.all_clades[clade_label].glen
                            tax_id = self.all_clades[clade_label].get_full_taxids()
                            
                            if 's__' in clade_label and abundance > 0:
                                self.all_clades[clade_label].nreads = int(round(abundance*glen,0))

                            clade_label = self.all_clades[clade_label].get_full_name()
                    elif not clade_label.startswith(tax_lev):
                        if clade_label in self.all_clades:
                            glen = self.all_clades[clade_label].glen
                        else:
                            glen = 1.0
                        continue
                    clade2abundance[(clade_label, tax_id)] = abundance
                    clade2genomelen[(clade_label, tax_id)] = glen
        
        for tax_label, clade in clade2abundance_n.items():
            tot_reads += clade.compute_mapped_reads()

        for clade_label, clade in self.all_clades.items():
            nreads = clade.nreads
            clade_label = clade.get_full_name()
            tax_id = clade.get_full_taxids()
            clade2est_nreads[(clade_label, tax_id)] = nreads

        ret_d = dict([( tax, float(abundance) / tot_ab if tot_ab else 0.0) for tax, abundance in clade2abundance.items()])

        ret_r = dict([( tax, (abundance, clade2genomelen[tax], clade2est_nreads[tax] )) for tax, abundance in clade2abundance.items()])

        if tax_lev:
            ret_d[tax_lev+"unclassified"] = 1.0 - sum(ret_d.values())
        return ret_d, ret_r, tot_reads


def map2bbh(mapping_f, input_type='bowtie2out', min_alignment_len=None):
    if not mapping_f:
        ras, ras_line, inpf = plain_read_and_split, plain_read_and_split_line, sys.stdin
    else:
        if mapping_f.endswith(".bz2"):
            ras, ras_line, inpf = read_and_split, read_and_split_line, bz2.BZ2File(mapping_f, "r")
        else:
            ras, ras_line, inpf = plain_read_and_split, plain_read_and_split_line, open(mapping_f)

    reads2markers = {}
    n_metagenoges_reads = None

    if input_type == 'bowtie2out':
        for r, c in ras(inpf):
            if r.startswith('#') and 'nreads' in r:
                n_metagenoges_reads = int(c)
            else:
                reads2markers[r] = c
    elif input_type == 'sam':
        for line in inpf:
            o = ras_line(line)

            if ((o[0][0] != '@') and
                (o[2][-1] != '*') and
                ((min_alignment_len is None) or
                 (max([int(x.strip('M')) for x in re.findall(r'(\d*M)', o[5]) if x]) >= min_alignment_len))):
                    reads2markers[o[0]] = o[2]

    inpf.close()
    markers2reads = defdict(set)

    for r, m in reads2markers.items():
        markers2reads[m].add(r)

    return (markers2reads, n_metagenoges_reads)


def maybe_generate_biom_file(tree, pars, abundance_predictions):
    json_key = "MetaPhlAn2"

    if not pars['biom']:
        return None
    if not abundance_predictions:
        biom_table = biom.Table([], [], [])  # create empty BIOM table

        with open(pars['biom'], 'w') as outfile:
            biom_table.to_json(json_key, direct_io=outfile)

        return True

    delimiter = "|" if len(pars['mdelim']) > 1 else pars['mdelim']

    def istip(clade_name):
        end_name = clade_name.split(delimiter)[-1]
        return end_name.startswith("t__") or end_name.endswith("_unclassified")

    def findclade(clade_name):
        if clade_name.endswith('_unclassified'):
            name = clade_name.split(delimiter)[-2]
        else:
            name = clade_name.split(delimiter)[-1]
        return tree.all_clades[name]

    def to_biomformat(clade_name):
        return {'taxonomy': clade_name.split(delimiter)}

    clades = iter((abundance, findclade(name))
                  for (name, abundance) in abundance_predictions if istip(name))
    packed = iter(([abundance], clade.get_full_name(), clade.tax_id)
                  for (abundance, clade) in clades)

    # unpack that tuple here to stay under 80 chars on a line
    data, clade_names, clade_ids = zip(*packed)
    # biom likes column vectors, so we give it an array like this:
    # np.array([a],[b],[c])
    data = np.array(data)
    sample_ids = [pars['sample_id']]
    table_id = 'MetaPhlAn2_Analysis'




    #**********************************************************************************************
    #  Modification of Code :                                                                     *
    #  Modified the code so instead of using the current clade IDs, which are numbers, we will    *
    #      use the clade_names                                                                    *
    #      Users reported the biom output is invalid and also the IDs were changing from run to   *
    #      run.                                                                                   *
    #  George Weingart    05/22/2017   george.weingart@mail.com                                   *
    #**********************************************************************************************
    if LooseVersion(biom.__version__) < LooseVersion("2.0.0"):
        biom_table = biom.table.table_factory(
            data,
        sample_ids,
            ######## clade_ids,     #Modified by George Weingart 5/22/2017 - We will use instead the clade_names
            clade_names,            #Modified by George Weingart 5/22/2017 - We will use instead the clade_names
            sample_metadata      = None,
            observation_metadata = map(to_biomformat, clade_names),
            table_id             = table_id,
            constructor          = biom.table.DenseOTUTable
        )
        with open(pars['biom'], 'w') as outfile:
            json.dump( biom_table.getBiomFormatObject(json_key),
                           outfile )
    else:  # Below is the biom2 compatible code
        biom_table = biom.table.Table(
            data,
            #clade_ids,           #Modified by George Weingart 5/22/2017 - We will use instead the clade_names
            clade_names,          #Modified by George Weingart 5/22/2017 - We will use instead the clade_names
            sample_ids,
            sample_metadata      = None,
            observation_metadata = map(to_biomformat, clade_names),
            table_id             = table_id,
            input_is_dense       = True
        )

        with open(pars['biom'], 'w') as outfile:
            biom_table.to_json( json_key,
                                direct_io = outfile )

    return True


def metaphlan2():
    pars = read_params(sys.argv)
    
    # check if the database is installed, if not then install
    check_and_install_database(pars['index'], pars['bowtie2db'], pars['bowtie2_build'], pars['nproc'])

    if pars['install']:
        sys.stderr.write('The database is installed\n')
        return

    # set correct map_pkl and bowtie2db variables
    pars['mpa_pkl'], pars['bowtie2db'] = set_mapping_arguments(pars['index'], pars['bowtie2db'])

    if (pars['bt2_ps'] in ["sensitive-local", "very-sensitive-local"]) and (pars['min_alignment_len'] is None):
            pars['min_alignment_len'] = 100
            sys.stderr.write('Warning! bt2_ps is set to local mode, and min_alignment_len is None, I automatically '
                             'set min_alignment_len to 100! If you do not like, rerun the command and set '
                             'min_alignment_len to a specific value.\n')

    if pars['input_type'] == 'fastq':
        pars['input_type'] = 'multifastq'
    if pars['input_type'] == 'fasta':
        pars['input_type'] = 'multifasta'

    # check for the mpa_pkl file
    if not os.path.isfile(pars['mpa_pkl']):
        sys.stderr.write("Error: Unable to find the mpa_pkl file at: " + pars['mpa_pkl'] +
                         "\nExpecting location ${mpa_dir}/db_v20/map_v20_m200.pkl "
                         "\nSelect the file location with the option --mpa_pkl.\n"
                         "Exiting...\n\n")
        sys.exit(1)

    if pars['ignore_markers']:
        with open(pars['ignore_markers']) as ignv:
            ignore_markers = set([l.strip() for l in ignv])
    else:
        ignore_markers = set()

    no_map = False
    if pars['input_type'] == 'multifasta' or pars['input_type'] == 'multifastq':
        bow = pars['bowtie2db'] is not None

        if not bow:
            sys.stderr.write( "No MetaPhlAn BowTie2 database provided\n "
                              "[--bowtie2db options]!\n"
                              "Exiting...\n\n" )
            sys.exit(1)

        if pars['no_map']:
            pars['bowtie2out'] = tf.NamedTemporaryFile(dir=pars['tmp_dir']).name
            no_map = True
        else:
            if bow and not pars['bowtie2out']:
                if pars['inp'] and "," in  pars['inp']:
                    sys.stderr.write("Error! --bowtie2out needs to be specified when multiple "
                                     "fastq or fasta files (comma separated) are provided\n")
                    sys.exit(1)
                fname = pars['inp']
                if fname is None:
                    fname = "stdin_map"
                elif stat.S_ISFIFO(os.stat(fname).st_mode):
                    fname = "fifo_map"
                pars['bowtie2out'] = fname + ".bowtie2out.txt"

            if os.path.exists( pars['bowtie2out'] ):
                sys.stderr.write(
                    "BowTie2 output file detected: " + pars['bowtie2out'] + "\n"
                    "Please use it as input or remove it if you want to "
                    "re-perform the BowTie2 run.\n"
                    "Exiting...\n\n" )
                sys.exit(1)

        if bow and not all([os.path.exists(".".join([str(pars['bowtie2db']), p]))
                            for p in ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]]):
            sys.stderr.write("No MetaPhlAn BowTie2 database found (--index "
                             "option)!\nExpecting location {}\nExiting..."
                             .format(pars['bowtie2db']))
            sys.exit(1)

        if bow:
            run_bowtie2(pars['inp'], pars['bowtie2out'], pars['bowtie2db'],
                                pars['bt2_ps'], pars['nproc'], file_format=pars['input_type'],
                                exe=pars['bowtie2_exe'], samout=pars['samout'],
                                min_alignment_len=pars['min_alignment_len'], read_min_len=pars['read_min_len'])
            pars['input_type'] = 'bowtie2out'
        pars['inp'] = pars['bowtie2out'] # !!!

    with open( pars['mpa_pkl'], 'rb' ) as a:
        mpa_pkl = pickle.loads( bz2.decompress( a.read() ) )

    tree = TaxTree( mpa_pkl, ignore_markers )
    tree.set_min_cu_len( pars['min_cu_len'] )
    tree.set_stat( pars['stat'], pars['stat_q'], pars['avoid_disqm'])

    markers2reads, n_metagenome_reads = map2bbh(pars['inp'], pars['input_type'],
                            pars['min_alignment_len'])
    if no_map:
        os.remove( pars['inp'] )

    if not n_metagenome_reads and not pars['nreads']:
        sys.stderr.write(
                "Please provide the size of the metagenome using the"
                "--nreads parameter when running MetaPhlAn2"
                "Exiting...\n\n" )
        sys.exit(1)


    map_out = []
    for marker,reads in sorted(markers2reads.items(), key=lambda pars: pars[0]):
        if marker not in tree.markers2lens:
            continue
        tax_seq, ids_seq = tree.add_reads( marker, len(reads),
                                  ignore_viruses = pars['ignore_viruses'],
                                  ignore_eukaryotes = pars['ignore_eukaryotes'],
                                  ignore_bacteria = pars['ignore_bacteria'],
                                  ignore_archaea = pars['ignore_archaea'],
                                  )
        if tax_seq:
            map_out +=["\t".join([r,tax_seq, ids_seq]) for r in sorted(reads)]

    if pars['output'] is None and pars['output_file'] is not None:
        pars['output'] = pars['output_file']

    with (open(pars['output'],"w") if pars['output'] else sys.stdout) as outf:
        if not pars['legacy_output']:
            outf.write('#{}\n'.format(pars['index']))

        outf.write('\t'.join((pars["sample_id_key"], pars["sample_id"])) + '\n')

        if pars['t'] == 'reads_map':
            if not pars['legacy_output']:
               outf.write('#read_id\tNCBI_taxlineage_str\tNCBI_taxlineage_ids\n')
            outf.write( "\n".join( map_out ) + "\n" )

        elif pars['t'] == 'rel_ab':
            if not pars['legacy_output']:
                outf.write('#clade_name\tNCBI_tax_id\trelative_abundance\n')

            cl2ab, _, tot_nreads = tree.relative_abundances(
                        pars['tax_lev']+"__" if pars['tax_lev'] != 'a' else None )

            fraction_mapped_reads = tot_nreads/n_metagenome_reads
            unmapped_reads = n_metagenome_reads - tot_nreads

            outpred = [(taxstr, taxid,round(relab*100.0,5)) for (taxstr, taxid),relab in cl2ab.items() if relab > 0.0]

            if outpred:
                outf.write( "\t".join( [    "UNKNOWN",
                                            "-1",
                                            str((1-fraction_mapped_reads)*100) ]) + "\n" )
                                                
                for clade, taxid, relab in sorted(  outpred, reverse=True,
                                    key=lambda x:x[2 if not pars['legacy_output'] else 1]+(100.0*(8-(x[0].count("|"))))):
                    if not pars['legacy_output']:
                        outf.write( "\t".join( [clade, 
                                                taxid, 
                                                str(relab*fraction_mapped_reads)] ) + "\n" )
                    else:
                        outf.write( "\t".join( [clade, 
                                                str(relab*fraction_mapped_reads)] ) + "\n" )
            else:
                if not pars['legacy_output']:
                    outf.write( "unclassified\t\-1\t100.0\n" )
                else:
                    outf.write( "unclassified\t100.0\n" )
            maybe_generate_biom_file(tree, pars, outpred)

        elif pars['t'] == 'rel_ab_w_read_stats':
            cl2ab, rr, tot_nreads = tree.relative_abundances(
                        pars['tax_lev']+"__" if pars['tax_lev'] != 'a' else None )
            fraction_mapped_reads = tot_nreads/n_metagenome_reads
            unmapped_reads = n_metagenome_reads - tot_nreads

            outpred = [(taxstr, taxid,round(relab*100.0,5)) for (taxstr, taxid),relab in cl2ab.items() if relab > 0.0]

            if outpred:
                outf.write( "#estimated_reads_mapped_to_known_clades:{}\n".format(tot_nreads) )
                outf.write( "\t".join( [    "#clade_name",
                                            "clade_taxid",
                                            "relative_abundance",
                                            "coverage",
                                            "average_genome_length_in_the_clade",
                                            "estimated_number_of_reads_from_the_clade" ]) +"\n" )
                outf.write( "\t".join( [    "UNKNOWN",
                                            "-1",
                                            str((1-fraction_mapped_reads)*100),
                                            "-",
                                            "-",
                                            str(unmapped_reads) ]) + "\n" )
                                                
                for taxstr, taxid, relab in sorted(  outpred, reverse=True,
                                    key=lambda x:x[2 if not pars['legacy_output'] else 1]+(100.0*(8-(x[0].count("|"))))):
                    outf.write( "\t".join( [    taxstr,
                                                taxid,
                                                str(relab*fraction_mapped_reads),
                                                str(rr[(taxstr, taxid)][0]) if (taxstr, taxid) in rr else "-",          #coverage
                                                str(rr[(taxstr, taxid)][1]) if (taxstr, taxid) in rr else "-",          #avg genome length in clade
                                                str(int(round(rr[(taxstr, taxid)][2],0)) if (taxstr, taxid) in rr else "-")      #estimated_number_of_reads_from_the_clade
                                                ] ) + "\n" )
            else:
                if not pars['legacy_output']:
                    outf.write( "unclassified\t-1\t100.0\n" )
                else:
                    outf.write( "unclassified\t100.0\n" )
            maybe_generate_biom_file(tree, pars, outpred)

        elif pars['t'] == 'clade_profiles':
            cl2pr = tree.clade_profiles( pars['tax_lev']+"__" if pars['tax_lev'] != 'a' else None  )
            for c,p in cl2pr.items():
                mn,n = zip(*p)
                outf.write( "\t".join( [""]+[str(s) for s in mn] ) + "\n" )
                outf.write( "\t".join( [c]+[str(s) for s in n] ) + "\n" )

        elif pars['t'] == 'marker_ab_table':
            cl2pr = tree.clade_profiles( pars['tax_lev']+"__" if pars['tax_lev'] != 'a' else None  )
            for v in cl2pr.values():
                outf.write( "\n".join(["\t".join([str(a),str(b/float(pars['nreads'])) if pars['nreads'] else str(b)])
                                for a,b in v if b > 0.0]) + "\n" )

        elif pars['t'] == 'marker_pres_table':
            cl2pr = tree.clade_profiles( pars['tax_lev']+"__" if pars['tax_lev'] != 'a' else None  )
            for v in cl2pr.values():
                strout = ["\t".join([str(a),"1"]) for a,b in v if b > pars['pres_th']]
                if strout:
                    outf.write( "\n".join(strout) + "\n" )

        elif pars['t'] == 'marker_counts':
            outf.write( "\n".join( ["\t".join([m,str(c)]) for m,c in tree.markers2counts().items() ]) +"\n" )

        elif pars['t'] == 'clade_specific_strain_tracker':
            cl2pr = tree.clade_profiles( None, get_all = True  )
            cl2ab, _ = tree.relative_abundances( None )
            strout = []
            for clade,v in cl2pr.items():
                if clade.endswith(pars['clade']) and cl2ab[clade]*100.0 < pars['min_ab']:
                    strout = []
                    break
                if pars['clade'] in clade:
                    strout += ["\t".join([str(a),str(int(b > pars['pres_th']))]) for a,b in v]
            if strout:
                strout = sorted(strout,key=lambda x:x[0])
                outf.write( "\n".join(strout) + "\n" )
            else:
                sys.stderr.write("Clade "+pars['clade']+" not present at an abundance >"+str(round(pars['min_ab'],2))+"%, "
                                 "so no clade specific markers are reported\n")


if __name__ == '__main__':
    metaphlan2()
