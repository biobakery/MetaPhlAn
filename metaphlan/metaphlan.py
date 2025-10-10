#!/usr/bin/env python
__author__ = ('Aitor Blanco-Miguez (aitor.blancomiguez@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it), '
              'Moreno Zolfo (moreno.zolfo@unitn.it), '
              'Nicola Segata (nicola.segata@unitn.it), '
              'Duy Tin Truong, '
              'Francesco Asnicar (f.asnicar@unitn.it), ' 
              'Claudia Mengoni (claudia.mengoni@unitn.it)')
__version__ = '4.1.1'
__date__ = '11 Mar 2024'

import sys
try:
    from metaphlan import mybytes, plain_read_and_split, plain_read_and_split_line, read_and_split, read_and_split_line, check_and_install_database, remove_prefix
except ImportError:
    sys.exit("CRITICAL ERROR: Unable to find the MetaPhlAn python package. Please check your install.")

if float(sys.version_info[0]) < 3.0:
    sys.stderr.write("MetaPhlAn requires Python 3, your current Python version is {}.{}.{}\n"
                    .format(sys.version_info[0], sys.version_info[1], sys.version_info[2]))
    sys.exit(1)
import os
import stat
import re
import time
import random
from collections import defaultdict as defdict
from distutils.version import LooseVersion
from glob import glob
from subprocess import DEVNULL
import argparse as ap
import bz2
import gzip
import pickle
import subprocess as subp
import tempfile as tf

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import Counter
try:
    from .utils.parallelisation import execute_pool
except ImportError:
    from utils.parallelisation import execute_pool
try:
    import pandas as pd
    import numpy as np
    import pysam
except Exception as e:
    sys.stderr.write("Error! python library not detected\n{}\n".format(e))
    sys.exit(1)

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
except ImportError:
    sys.stderr.write("Warning! Biom python library not detected!"
                     "\n Exporting to biom format will not work!\n")
import json

# get the directory that contains this script
metaphlan_script_install_folder = os.path.dirname(os.path.abspath(__file__))
# get the default database folder
DEFAULT_DB_FOLDER = os.path.join(metaphlan_script_install_folder, "metaphlan_databases")
DEFAULT_DB_FOLDER= os.environ.get('METAPHLAN_DB_DIR', DEFAULT_DB_FOLDER)
#Wether to execute a SGB-based analysis
SGB_ANALYSIS = True
INDEX = 'latest'
tax_units = "kpcofgst"

def read_params(args):
    p = ap.ArgumentParser( description =
            "DESCRIPTION\n"
            " MetaPhlAn version "+__version__+" ("+__date__+"): \n"
            " METAgenomic PHyLogenetic ANalysis for metagenomic taxonomic profiling.\n\n"
            "AUTHORS: "+__author__+"\n\n"
            "COMMON COMMANDS\n\n"
            " We assume here that MetaPhlAn is installed using the several options available (pip, conda, PyPi)\n"
            " Also BowTie2 should be in the system path with execution and read permissions, and Perl should be installed)\n\n"

            "\n========== MetaPhlAn clade-abundance estimation ================= \n\n"
            "The basic usage of MetaPhlAn consists in the identification of the clades (from phyla to species ) \n"
            "present in the metagenome obtained from a microbiome sample and their \n"
            "relative abundance. This correspond to the default analysis type (-t rel_ab).\n\n"

            "*  Profiling a metagenome from raw reads:\n"
            "$ metaphlan metagenome.fastq --input_type fastq -o profiled_metagenome.txt\n\n"

            "*  You can take advantage of multiple CPUs and save the intermediate BowTie2 output for re-running\n"
            "   MetaPhlAn extremely quickly:\n"
            "$ metaphlan metagenome.fastq --bowtie2out metagenome.bowtie2.bz2 --nproc 5 --input_type fastq -o profiled_metagenome.txt\n\n"

            "*  If you already mapped your metagenome against the marker DB (using a previous MetaPhlAn run), you\n"
            "   can obtain the results in few seconds by using the previously saved --bowtie2out file and \n"
            "   specifying the input (--input_type bowtie2out):\n"
            "$ metaphlan metagenome.bowtie2.bz2 --nproc 5 --input_type bowtie2out -o profiled_metagenome.txt\n\n"
            
            "*  bowtie2out files generated with MetaPhlAn versions below 3 are not compatibile.\n"
            "   Starting from MetaPhlAn 3.0, the BowTie2 ouput now includes the size of the profiled metagenome and the average read length.\n"
            "   If you want to re-run MetaPhlAn using these file you should provide the metagenome size via --nreads:\n"
            "$ metaphlan metagenome.bowtie2.bz2 --nproc 5 --input_type bowtie2out --nreads 520000 -o profiled_metagenome.txt\n\n"

            "*  You can also provide an externally BowTie2-mapped SAM if you specify this format with \n"
            "   --input_type. Two steps: first apply BowTie2 and then feed MetaPhlAn with the obtained sam:\n"
            "$ bowtie2 --sam-no-hd --sam-no-sq --no-unal --very-sensitive -S metagenome.sam -x ${mpa_dir}/metaphlan_databases/mpa_v30_CHOCOPhlAn_201901 -U metagenome.fastq\n"
            "$ metaphlan metagenome.sam --input_type sam -o profiled_metagenome.txt\n\n"

            "*  We can also natively handle paired-end metagenomes, and, more generally, metagenomes stored in \n"
            "  multiple files (but you need to specify the --bowtie2out parameter):\n"
            "$ metaphlan metagenome_1.fastq,metagenome_2.fastq --bowtie2out metagenome.bowtie2.bz2 --nproc 5 --input_type fastq\n\n"
            "\n------------------------------------------------------------------- \n \n\n"

            "\n========== Marker level analysis ============================ \n\n"
            "MetaPhlAn introduces the capability of characterizing organisms at the strain level using non\n"
            "aggregated marker information. Such capability comes with several slightly different flavours and \n"
            "are a way to perform strain tracking and comparison across multiple samples.\n"
            "Usually, MetaPhlAn is first ran with the default -t to profile the species present in\n"
            "the community, and then a strain-level profiling can be performed to zoom-in into specific species\n"
            "of interest. This operation can be performed quickly as it exploits the --bowtie2out intermediate \n"
            "file saved during the execution of the default analysis type.\n\n"

            "*  The following command will output the abundance of each marker with a RPK (reads per kilo-base) \n"
            "   higher 0.0. (we are assuming that metagenome_outfmt.bz2 has been generated before as \n"
            "   shown above).\n"
            "$ metaphlan -t marker_ab_table metagenome_outfmt.bz2 --input_type bowtie2out -o marker_abundance_table.txt\n"
            "   The obtained RPK can be optionally normalized by the total number of reads in the metagenome \n"
            "   to guarantee fair comparisons of abundances across samples. The number of reads in the metagenome\n"
            "   needs to be passed with the '--nreads' argument\n\n"

            "*  The list of markers present in the sample can be obtained with '-t marker_pres_table'\n"
            "$ metaphlan -t marker_pres_table metagenome_outfmt.bz2 --input_type bowtie2out -o marker_abundance_table.txt\n"
            "   The --pres_th argument (default 1.0) set the minimum RPK value to consider a marker present\n\n"

            "*  The list '-t clade_profiles' analysis type reports the same information of '-t marker_ab_table'\n"
            "   but the markers are reported on a clade-by-clade basis.\n"
            "$ metaphlan -t clade_profiles metagenome_outfmt.bz2 --input_type bowtie2out -o marker_abundance_table.txt\n\n"

            "*  Finally, to obtain all markers present for a specific clade and all its subclades, the \n"
            "   '-t clade_specific_strain_tracker' should be used. For example, the following command\n"
            "   is reporting the presence/absence of the markers for the B. fragilis species and its strains\n"
            "   the optional argument --min_ab specifies the minimum clade abundance for reporting the markers\n\n"
            "$ metaphlan -t clade_specific_strain_tracker --clade s__Bacteroides_fragilis metagenome_outfmt.bz2 --input_type bowtie2out -o marker_abundance_table.txt\n"

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
    input_type_choices = ['fastq','fasta','bowtie2out','sam']
    arg( '--input_type', choices=input_type_choices, required = '--install' not in args, help =
         "set whether the input is the FASTA file of metagenomic reads or \n"
         "the SAM file of the mapping of the reads against the MetaPhlAn db.\n"
        )

    g = p.add_argument_group('Mapping arguments')
    arg = g.add_argument
    arg('--force', action='store_true', help="Force profiling of the input file by removing the bowtie2out file")
    arg('--bowtie2db', metavar="METAPHLAN_BOWTIE2_DB", type=str, default=DEFAULT_DB_FOLDER,
        help=("Folder containing the MetaPhlAn database. You can specify the location by exporting the DEFAULT_DB_FOLDER variable in the shell."
              "[default "+DEFAULT_DB_FOLDER+"]\n"))

    arg('-x', '--index', type=str, default=INDEX,
        help=("Specify the id of the database version to use. "
              "If \"latest\", MetaPhlAn will get the latest version.\n"
              "If an index name is provided, MetaPhlAn will try to use it, if available, and skip the online check.\n"
              "If the database files are not found on the local MetaPhlAn installation they\n"
              "will be automatically downloaded [default "+INDEX+"]\n"))

    bt2ps = ['sensitive', 'very-sensitive', 'sensitive-local', 'very-sensitive-local']
    arg('--bt2_ps', metavar="BowTie2 presets", default='very-sensitive',
        choices=bt2ps, help="Presets options for BowTie2 (applied only when a "
                            "FASTA file is provided)\n"
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
    arg('--min_mapq_val', type=int, default=5,
        help="Minimum mapping quality value (MAPQ) [default 5]")
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
         "'t' : SGBs only\n"
         "[default 'a']" )
    arg( '--min_cu_len', metavar="", default="2000", type=int, help =
         "minimum total nucleotide length for the markers in a clade for\n"
         "estimating the abundance without considering sub-clade abundances\n"
         "[default 2000]\n"   )
    arg( '--min_alignment_len', metavar="", default=None, type=int, help =
         "The sam records for aligned reads with the longest subalignment\n"
         "length smaller than this threshold will be discarded.\n"
         "[default None]\n"   )
    arg( '--add_viruses', action='store_true', help=
         "Together with --mpa3, allow the profiling of viral organisms" )
    arg( '--ignore_eukaryotes', action='store_true', help=
         "Do not profile eukaryotic organisms" )
    arg( '--ignore_bacteria', action='store_true', help=
         "Do not profile bacterial organisms" )
    arg( '--ignore_archaea', action='store_true', help=
         "Do not profile archeal organisms" )
    arg( '--ignore_ksgbs', action='store_true', help=
         "Do not profile known SGBs (together with --sgb option)" )
    arg( '--ignore_usgbs', action='store_true', help=
         "Do not profile unknown SGBs (together with --sgb option)" )
    arg( '--stat_q', metavar="", type = float, default=0.2, help =
         "Quantile value for the robust average\n"
         "[default 0.2]" )
    arg( '--perc_nonzero', metavar="", type = float, default=0.33, help =
         "Percentage of markers with a non zero relative abundance for misidentify a species\n"
         "[default 0.33]"   )
    arg( '--ignore_markers', type=str, default = None, help =
         "File containing a list of markers to ignore. \n")
    arg( '--avoid_disqm', action="store_true", help =
         "Deactivate the procedure of disambiguating the quasi-markers based on the \n"
         "marker abundance pattern found in the sample. It is generally recommended \n"
         "to keep the disambiguation procedure in order to minimize false positives\n")
    arg( '--stat', metavar="", choices=stat_choices, default="tavg_g", type=str, help =
         "Statistical approach for converting marker abundances into clade abundances\n"
         "'avg_g'  : clade global (i.e. normalizing all markers together) average\n"
         "'avg_l'  : average of length-normalized marker counts\n"
         "'tavg_g' : truncated clade global average at --stat_q quantile\n"
         "'tavg_l' : truncated average of length-normalized marker counts (at --stat_q)\n"
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
         " * rel_ab_w_read_stats: profiling a metagenomes in terms of relative abundances and estimate the number of reads coming from each clade.\n"
         " * reads_map: mapping from reads to clades (only reads hitting a marker)\n"
         " * clade_profiles: normalized marker counts for clades with at least a non-null marker\n"
         " * marker_ab_table: normalized marker counts (only when > 0.0 and normalized by metagenome size if --nreads is specified)\n"
         " * marker_counts: non-normalized marker counts [use with extreme caution]\n"
         " * marker_pres_table: list of markers present in the sample (threshold at 1.0 if not differently specified with --pres_th\n"
         " * clade_specific_strain_tracker: list of markers present for a specific clade, specified with --clade, and all its subclades\n"
         "[default 'rel_ab']" )
    arg( '--nreads', metavar="NUMBER_OF_READS", type=int, default = None, help =
         "The total number of reads in the original metagenome. It is mandatory when the --input_type is a SAM file." )
    arg( '--pres_th', metavar="PRESENCE_THRESHOLD", type=int, default = 1.0, help =
         'Threshold for calling a marker present by the -t marker_pres_table option' )
    arg( '--clade', metavar="", default=None, type=str, help =
         "The clade for clade_specific_strain_tracker analysis\n"  )
    arg( '--min_ab', metavar="", default=0.1, type=float, help =
         "The minimum percentage abundance for the clade in the clade_specific_strain_tracker analysis\n"  )

    g = p.add_argument_group('Viral Sequence Clusters Analisys')
    arg = g.add_argument
    arg("--profile_vsc", action="store_true",help="Add this parameter to profile Viruses with VSCs approach.")
    arg("--vsc_out", help="Path to the VSCs breadth-of-coverage output file", default="mp3_viruses.csv")
    arg("--vsc_breadth", help="Minimum Breadth of Coverage for a Viral Group to be reported.\n"
        "Default is 0.75 (at least 75 percent breadth to report)", default=0.75,type=float)

    g = p.add_argument_group('Output arguments')
    arg = g.add_argument
    arg( '-o', '--output_file',  metavar="output file", type=str, default=None, help =
         "The output file (if not specified as positional argument)\n")
    arg('--sample_id_key',  metavar="name", type=str, default="SampleID",
        help =("Specify the sample ID key for this analysis."
               " Defaults to 'SampleID'."))
    arg('--use_group_representative', action='store_true',  help =("Use a species as representative for species groups."))
    arg('--sample_id',  metavar="value", type=str,
        default="Metaphlan_Analysis",
        help =("Specify the sample ID for this analysis."
               " Defaults to 'Metaphlan_Analysis'."))
    arg( '-s', '--samout', metavar="sam_output_file",
        type=str, default=None, help="The sam output file\n")

    arg( '--legacy-output', action='store_true', help="Old MetaPhlAn2 two columns output\n")
    arg( '--CAMI_format_output', action='store_true', help="Report the profiling using the CAMI output format\n")
    arg( '--unclassified_estimation', action='store_true', help="Scale relative abundances to the number of reads mapping to identified clades in order to estimate unclassified taxa\n")
    arg( '--mpa3', action='store_true', help="Perform the analysis using the MetaPhlAn 3 algorithm\n")

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
    arg('--subsampling', type=int, default=None,
        help="Specify the number of reads to be considered from the input metagenomes [default None]")
    arg('--subsampling_output', type=str, default=None,
        help="The output file for the subsampled reads. If --subsampling_paired is specified two files are created with suffixes R1 and R2. If not specified the subsampled reads will not be saved.")
    arg('--subsampling_paired',  type=int, default=None,
        help="Specify the number of paired reads to be considered from the input metagenomes [default None]")
    arg('-1', type=str, default=None, metavar='FORWARD_READS',
        help="Specify the fastq file with forward reads of the input metagenomes. Reads are assumed to be in the same order in the forward and reverse files! [default None]")
    arg('-2', type=str, default=None, metavar='REVERSE_READS',
        help="Specify the fastq file with reverse reads of the input metagenomes. Reads are assumed to be in the same order in the forward and reverse files! [default None]")
    arg('--mapping_subsampling', action='store_true',
        help="If used, the subsamping will be done on the mapping results instead of on the reads.")
    arg('--subsampling_seed', type=str, default='1992',
        help="Random seed to use in the selection of the subsampled reads. Choose \"random\r for a random behaviour") 
    arg('--install', action='store_true',
        help="Only checks if the MetaPhlAn DB is installed and installs it if not. All other parameters are ignored.")
    arg('--offline', action='store_true',
        help="If used, MetaPhlAn will not check for new database updates.")
    arg('--force_download', action='store_true',
        help="Force the re-download of the latest MetaPhlAn database.")
    arg('--read_min_len', type=int, default=70,
        help="Specify the minimum length of the reads to be considered when parsing the input file with "
             "'read_fastx.py' script, default value is 70")
    arg('-v', '--version', action='store_true',
        help="Prints the current MetaPhlAn version and exit")
    arg("-h", "--help", action="help", help="show this help message and exit")

    return vars(p.parse_args())

def get_installed_db_version(db_dir: str = None):
    """
    Return a sorted list of DB base names (e.g., mpa_vJan24_CHOCOPhlAnSGB_202401)
    for which both .pkl and .fna exist in db_dir.
    """
    if db_dir:
        db_folder = db_dir
    else:
        db_folder = DEFAULT_DB_FOLDER
    if not os.path.exists(db_folder):
        return None
    files = os.listdir(db_folder)
    pkl_bases = {os.path.splitext(f)[0] for f in files if f.endswith(".pkl")}
    bt2l_bases = {os.path.basename(f).split('.')[0] for f in files if f.endswith(".bt2l")}
    common = sorted(pkl_bases & bt2l_bases)
    return common or None

def get_cli_db_dir(argv=None):
    """
    Robustly detect database folder from CLI (--db_dir or --bowtie2db),
    supporting both '--flag value' and '--flag=value'.
    """
    if argv is None:
        argv = sys.argv[1:]

    p = ap.ArgumentParser(add_help=False)
    p.add_argument("--db_dir", dest="db_dir")
    p.add_argument("--bowtie2db", dest="bowtie2db")
    # Ignore everything else for now
    args, _ = p.parse_known_args(argv)

    val = args.db_dir or args.bowtie2db
    if not val:
        return None
    # normalize: expand ~ and $VARS, strip trailing slash, absolutize
    val = os.path.abspath(os.path.expanduser(os.path.expandvars(val).rstrip("/")))
    return val

def set_mapping_arguments(index, bowtie2_db):
    mpa_pkl = 'mpa_pkl'
    bowtie2db = 'bowtie2db'
    bt2_ext = 'bt2l' if SGB_ANALYSIS else 'bt2'

    if os.path.isfile(os.path.join(bowtie2_db, "{}.pkl".format(index))):
        mpa_pkl = os.path.join(bowtie2_db, "{}.pkl".format(index))

    if glob(os.path.join(bowtie2_db, "{}*.{}".format(index, bt2_ext))):
        bowtie2db = os.path.join(bowtie2_db, "{}".format(index))

    return (mpa_pkl, bowtie2db)


def set_vsc_parameters(index, bowtie2_db):
    vsc_fna = os.path.join(bowtie2_db, "{}_VSG.fna".format(index))
    vsc_vinfo = os.path.join(bowtie2_db, "{}_VINFO.csv".format(index))

    if not os.path.isfile(vsc_fna):
        sys.stderr.write("Error:\n {} file not found in bowtie2db folder ({}). Re-download MetaPhlAn database.".format("{}_VSG.fna".format(index, bowtie2_db)))
        sys.exit(1)

    if not os.path.isfile(vsc_vinfo):
        sys.stderr.write("Error:\n {} file not found in bowtie2db folder ({}). Re-download MetaPhlAn database.".format("{}_VINFO.csv".format(index, bowtie2_db)))
        sys.exit(1)

    return(vsc_fna, vsc_vinfo)

def vsc_bowtie2(profile_vsc_folder, nproc, file_format="fasta",
                exe=None,bt2build_exe=None, min_alignment_len=None, read_min_len=0):

    try:
        subp.check_call([exe if exe else 'bowtie2', "-h"], stdout=DEVNULL)
        subp.check_call([bt2build_exe if bt2build_exe else 'bowtie2-build', "-h"], stdout=DEVNULL)
    except Exception as e:
        sys.stderr.write('OSError: "{}"\nFatal error running BowTie2 and bowtie2-build. Is BowTie2 in the system path?\n'.format(e))
        sys.exit(1)

    # checking files
    markerfile= profile_vsc_folder+'/v_mks.fa' 
    readsfile= profile_vsc_folder+'/v_reads.fq' 
    vscBamFile= profile_vsc_folder+'/v_align.bam'

    if not os.path.exists(markerfile) or not os.path.exists(readsfile):
        
        sys.stderr.write('Error:\nThere was an error in the VSCs file lookup.\n \
            It may be that there are not enough reads to profile (not enough depth).\n \
            Passing without reporting any virus.\n')
        return []

    if os.stat(markerfile).st_size == 0:
        sys.stderr.write('Warning:\nNo reads aligning to VSC markers in this file.\n')        
        #return an empty list. MetaPhlAn will continue as normal, without VSCs
        return []

    #build the db

    bt2build_call = bt2build_exe if bt2build_exe else 'bowtie2-build'
    dbpath = profile_vsc_folder+'/v_mks'
    subp.check_call( [bt2build_call, markerfile, dbpath, '-q','--threads', str(nproc)] )


    if not all([os.path.exists(".".join([str(dbpath), p]))
                            for p in ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]]):
        sys.stderr.write('Error:\nThere was an error in running bowtie2-map and not all files were generated. Stopping.\n')
        sys.exit(1)

    try:

        bt2_call = exe if exe else 'bowtie2'
        bt2_command = [bt2_call,'--very-sensitive','--quiet','-p',str(nproc),'-x',dbpath,'--no-unal','-U',readsfile]
        stv_command = ['samtools','view','-bS','-','-@',str(nproc)]
        sts_command = ['samtools','sort','-','-@',str(nproc),'-o',vscBamFile]
        sti_command = ['samtools','index',vscBamFile]

        p3 = subp.Popen(sts_command, stdin=subp.PIPE)
        p2 = subp.Popen(stv_command, stdin=subp.PIPE, stdout=p3.stdin)
        p1 = subp.Popen(bt2_command, stdout=p2.stdin)

        p1.communicate() 
        p2.communicate()
        p3.communicate()

        #index
        subp.check_call(sti_command)

    except Exception as e:
        sys.stderr.write('Error: "{}"\nFatal error running BowTie2 for Viruses\n'.format(e))
        sys.exit(1)

    if not os.path.exists(vscBamFile):
        sys.stderr.write('Error:\nUnable to create BAM FILE for viruses\n')
        sys.exit(1)

    try:
        bamHandle = pysam.AlignmentFile(vscBamFile, "rb")
    except Exception as e:
        sys.stderr.write('Error: "{}"\nCheck PySam is correctly working\n'.format(e))
        sys.exit(1)

    VSC_report=[]

    for c, length in zip(bamHandle.references,bamHandle.lengths):
        coverage_positions = {}

        for pileupcolumn in bamHandle.pileup(c):

            tCoverage = 0
            for pileupread in pileupcolumn.pileups:

                if not pileupread.is_del and not pileupread.is_refskip \
                        and pileupread.alignment.query_qualities[pileupread.query_position] >= 20 \
                        and pileupread.alignment.query_sequence[pileupread.query_position].upper() in ('A','T','C','G'):
                        tCoverage +=1

            if tCoverage >= 1:
                coverage_positions[pileupcolumn.pos] = tCoverage
        
        breadth = float(len(coverage_positions.keys()))/float(length)
        
        if breadth > 0:
            cvals=list(coverage_positions.values())
            VSC_report.append({'M-Group/Cluster':c.split('|')[2].split('-')[0], 'genomeName':c, 'len':length, 'breadth_of_coverage':breadth, 'depth_of_coverage_mean': np.mean(cvals), 'depth_of_coverage_median': np.median(cvals)})

    if bamHandle:
        bamHandle.close()

    return VSC_report


def run_bowtie2(fna_in, outfmt6_out, bowtie2_db, preset, nproc, min_mapq_val, file_format="fasta",
                exe=None, samout=None, min_alignment_len=None, read_min_len=0, profile_vsc_folder=False):
    # checking read_fastx.py
    read_fastx = "read_fastx.py"

    try:
        subp.check_call([read_fastx, "-h"], stdout=DEVNULL, stderr=DEVNULL)
    except Exception as e:
        try:
            read_fastx = os.path.join(os.path.join(os.path.dirname(__file__), "utils"), read_fastx)
            subp.check_call([read_fastx, "-h"], stdout=DEVNULL, stderr=DEVNULL)
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

        bowtie2_cmd = [exe if exe else 'bowtie2', "--seed", "1992", "--quiet", "--no-unal", "--{}".format(preset),
                       "-S", "-", "-x", bowtie2_db]

        if int(nproc) > 1:
            bowtie2_cmd += ["-p", str(nproc)]

        bowtie2_cmd += ["-U", "-"]  # if not stat.S_ISFIFO(os.stat(fna_in).st_mode) else []

        if file_format == "fasta":
            bowtie2_cmd += ["-f"]

        p = subp.Popen(bowtie2_cmd, stdout=subp.PIPE, stdin=readin.stdout)
        readin.stdout.close()
        lmybytes, outf = (mybytes, bz2.BZ2File(outfmt6_out, "w")) if outfmt6_out.endswith(".bz2") else (str, open(outfmt6_out, "w"))

        if profile_vsc_folder:
            CREAD=[]
            list_of_viral_markers = open(profile_vsc_folder+'/viralmk.txt','w')

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
                    if (hex(int(o[1]) & 0x100) == '0x0'): #no secondary
                        if mapq_filter(o[2], int(o[4]), min_mapq_val) :  # filter low mapq reads
                            if ((min_alignment_len is None) or
                                    (max([int(x.strip('M')) for x in re.findall(r'(\d*M)', o[5]) if x]) >= min_alignment_len)):
                                                                # Profile viral markers in a different way
                                if profile_vsc_folder and o[2].startswith('VDB|'):

                                    mCluster = o[2]
                                    mGroup = o[2].split('|')[2].split('-')[0]

                                    list_of_viral_markers.write(mGroup+'\t'+mCluster+'\n')

                                    if (hex(int(o[1]) & 0x10) == '0x0'): #front read
                                        rr=SeqRecord(Seq(o[9]),letter_annotations={'phred_quality':[ord(_)-33 for _ in o[10][::-1]]}, id=o[0])
                                    else:
                                        rr=SeqRecord(Seq(o[9]).reverse_complement(),letter_annotations={'phred_quality':[ord(_)-33 for _ in o[10][::-1]]}, id=o[0])

                                    CREAD.append(rr)

                                # normal route for non-viral markers
                                outf.write(lmybytes("\t".join([ o[0], o[2].split('/')[0] ]) + "\n"))

        if profile_vsc_folder and os.path.isdir(profile_vsc_folder):
            SeqIO.write(CREAD,profile_vsc_folder+'/v_reads.fq','fastq')
            list_of_viral_markers.close()

        if samout:  
            sam_file.close()

        p.communicate()
        read_fastx_stderr = readin.stderr.readlines()
        nreads = None
        avg_read_length = None
        try:
            nreads, avg_read_length = list(map(float, read_fastx_stderr[0].decode().split()))
            if not nreads:
                sys.stderr.write('Fatal error running MetaPhlAn. Total metagenome size was not estimated.\nPlease check your input files.\n')
                sys.exit(1)
            if not avg_read_length:
                sys.stderr.write('Fatal error running MetaPhlAn. The average read length was not estimated.\nPlease check your input files.\n')
                sys.exit(1)
            outf.write(lmybytes('#nreads\t{}\n'.format(int(nreads))))
            outf.write(lmybytes('#avg_read_length\t{}'.format(avg_read_length)))
            outf.close()
        except ValueError:
            sys.stderr.write(b''.join(read_fastx_stderr).decode())
            outf.close()
            os.unlink(outfmt6_out)
            sys.exit(1)

    except OSError as e:
        sys.stderr.write('OSError: "{}"\nFatal error running BowTie2.\n'.format(e))
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

class TaxClade:
    min_cu_len = -1
    markers2lens = None
    stat = None
    perc_nonzero = None
    quantile = None
    avoid_disqm = False
    avg_read_length = 1

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
        fullname = ['']
        if self.tax_id:
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
        return [(m,float(n)*1000.0/(np.absolute(self.markers2lens[m] - self.avg_read_length) +1) )
                    for m,n in self.markers2nreads.items()]

    def compute_mapped_reads( self ):    
        tax_level = 't__' if SGB_ANALYSIS else 's__'
        if self.nreads != 0 or self.name.startswith(tax_level):
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
                        if float(nonzeros) / len(m2nr) > self.perc_nonzero:
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

        if not SGB_ANALYSIS and self.name[0] == 't' and (len(self.father.children) > 1 or "_sp" in self.father.name or "k__Viruses" in self.get_full_name()):
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
            loc_ab = np.mean([float(n)/(np.absolute(r - self.avg_read_length) + 1) for r,n in rat_nreads])
        elif self.stat == 'tavg_g':
            wnreads = sorted([(float(n)/(np.absolute(r-self.avg_read_length)+1),(np.absolute(r - self.avg_read_length)+1) ,n) for r,n in rat_nreads], key=lambda x:x[0])
            den,num = zip(*[v[1:] for v in wnreads[ql:qr]])
            loc_ab = float(sum(num))/float(sum(den)) if any(den) else 0.0
        elif self.stat == 'tavg_l':
            loc_ab = np.mean(sorted([float(n)/(np.absolute(r - self.avg_read_length) + 1) for r,n in rat_nreads])[ql:qr])
        elif self.stat == 'wavg_g':
            vmin, vmax = nreads_v[ql], nreads_v[qr]
            wnreads = [vmin]*qn+list(nreads_v[ql:qr])+[vmax]*qn
            loc_ab = float(sum(wnreads)) / rat
        elif self.stat == 'wavg_l':
            wnreads = sorted([float(n)/(np.absolute(r - self.avg_read_length) + 1) for r,n in rat_nreads])
            vmin, vmax = wnreads[ql], wnreads[qr]
            wnreads = [vmin]*qn+list(wnreads[ql:qr])+[vmax]*qn
            loc_ab = np.mean(wnreads)
        elif self.stat == 'med':
            loc_ab = np.median(sorted([float(n)/(np.absolute(r - self.avg_read_length) +1) for r,n in rat_nreads])[ql:qr])

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
        self.avg_read_length = 1

        for clade, value in mpa['taxonomy'].items():
            clade = clade.strip().split("|")
            if isinstance(value,tuple):
                taxids, lenc = value
                taxids = taxids.strip().split("|")
            if isinstance(value,int):
                lenc = value
                taxids = None

            father = self.root
            for i in range(len(clade)):
                clade_lev = clade[i]
                if SGB_ANALYSIS:
                    clade_taxid = taxids[i] if i < 8 and taxids is not None else None
                else:
                    clade_taxid = taxids[i] if i < 7 and taxids is not None else None
                if not clade_lev in father.children:
                    father.add_child(clade_lev, tax_id=clade_taxid)
                    self.all_clades[clade_lev] = father.children[clade_lev]
                if SGB_ANALYSIS: father = father.children[clade_lev]
                if clade_lev[0] == "t":
                    self.taxa2clades[clade_lev[3:]] = father
                if not SGB_ANALYSIS: father = father.children[clade_lev]
                if clade_lev[0] == "t":
                    father.glen = lenc

        def add_lens( node ):
            if not node.children:
                return node.glen
            lens = []
            for c in node.children.values():
                lens.append( add_lens( c ) )
            node.glen = min(np.mean(lens), np.median(lens))
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

    def set_stat( self, stat, quantile, perc_nonzero, avg_read_length, avoid_disqm = False):
        TaxClade.stat = stat
        TaxClade.perc_nonzero = perc_nonzero
        TaxClade.quantile = quantile
        TaxClade.avoid_disqm = avoid_disqm
        TaxClade.avg_read_length = avg_read_length

    def add_reads(  self, marker, n,
                    add_viruses = False,
                    ignore_eukaryotes = False,
                    ignore_bacteria = False, ignore_archaea = False, 
                    ignore_ksgbs = False, ignore_usgbs = False  ):
        clade = self.markers2clades[marker]
        cl = self.all_clades[clade]
        if ignore_bacteria or ignore_archaea or ignore_eukaryotes:
            cn = cl.get_full_name()
            if ignore_archaea and cn.startswith("k__Archaea"):
                return (None, None)
            if ignore_bacteria and cn.startswith("k__Bacteria"):
                return (None, None)
            if ignore_eukaryotes and cn.startswith("k__Eukaryota"):
                return (None, None)
        if not SGB_ANALYSIS and not add_viruses:
            cn = cl.get_full_name()
            if not add_viruses and cn.startswith("k__Vir"):
                return (None, None)
        if SGB_ANALYSIS and (ignore_ksgbs or ignore_usgbs):
            cn = cl.get_full_name()
            if ignore_ksgbs and not '_SGB' in cn.split('|')[-2]:
                return (None, None)
            if ignore_usgbs and '_SGB' in cn.split('|')[-2]:
                return (None, None)
        # while len(cl.children) == 1:
            # cl = list(cl.children.values())[0]
        cl.markers2nreads[marker] = n
        return (cl.get_full_name(), cl.get_full_taxids(), )


    def markers2counts( self ):
        m2c = {}
        for _ ,v in self.all_clades.items():
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

        clade2abundance, clade2est_nreads, tot_ab, tot_reads = {}, {}, 0.0, 0

        for tax_label, clade in clade2abundance_n.items():
            tot_ab += clade.compute_abundance()

        for tax_label, clade in clade2abundance_n.items():
            for clade_label, tax_id, abundance in sorted(clade.get_all_abundances(), key=lambda pars:pars[0]):
                if SGB_ANALYSIS or clade_label[:3] != 't__':
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
                            tax_level = 't__' if SGB_ANALYSIS else 's__' 
                            if tax_level in clade_label and abundance > 0:
                                self.all_clades[clade_label].nreads = int(np.floor(abundance*glen))

                            clade_label = self.all_clades[clade_label].get_full_name()
                    elif not clade_label.startswith(tax_lev):
                        if clade_label in self.all_clades:
                            glen = self.all_clades[clade_label].glen
                        else:
                            glen = 1.0
                        continue
                    clade2abundance[(clade_label, tax_id)] = abundance
        
        for tax_label, clade in clade2abundance_n.items():
            tot_reads += clade.compute_mapped_reads()

        for clade_label, clade in self.all_clades.items():
            if SGB_ANALYSIS or clade.name[:3] != 't__':
                nreads = clade.nreads
                clade_label = clade.get_full_name()
                tax_id = clade.get_full_taxids()
                clade2est_nreads[(clade_label, tax_id)] = nreads

        ret_d = dict([( tax, float(abundance) / tot_ab if tot_ab else 0.0) for tax, abundance in clade2abundance.items()])

        ret_r = dict([( tax, (abundance, clade2est_nreads[tax] )) for tax, abundance in clade2abundance.items() if tax in clade2est_nreads])

        if tax_lev:
            ret_d[("UNCLASSIFIED", '-1')] = 1.0 - sum(ret_d.values())
        return ret_d, ret_r, tot_reads

def mapq_filter(marker_name, mapq_value, min_mapq_val):
    ##if 'GeneID:' in marker_name:
    if 'GeneID:' in marker_name or 'VDB' in marker_name:
        return True
    else:
        if mapq_value > min_mapq_val:
            return True
    return False


def separate_reads2markers(reads2markers):
    if not SGB_ANALYSIS:
        return reads2markers, {}
    else:
        return {r: m for r, m in reads2markers.items() if ('SGB' in m or 'EUK' in m) and not 'VDB' in m}, {r: m for r, m in reads2markers.items() if 'VDB' in m and not ('SGB' in m or 'EUK' in m)}

def map2bbh(mapping_f, min_mapq_val, input_type='bowtie2out', min_alignment_len=None, nreads=None, mapping_subsampling=False, subsampling=None, subsampling_seed='1992', remove_input=False):
    if not mapping_f:
        ras, ras_line, inpf = plain_read_and_split, plain_read_and_split_line, sys.stdin
    else:
        if mapping_f.endswith(".bz2"):
            ras, ras_line, inpf = read_and_split, read_and_split_line, bz2.BZ2File(mapping_f, "r")
        else:
            ras, ras_line, inpf = plain_read_and_split, plain_read_and_split_line, open(mapping_f)

    reads2markers = {}
    n_metagenome_reads = None
    avg_read_length = 1 #Set to 1 if it is not calculated from read_fastx

    if input_type == 'bowtie2out':
        for r, c in ras(inpf):
            if r.startswith('#') and 'nreads' in r:
                n_metagenome_reads = int(c)
            if r.startswith('#') and 'avg_read_length' in r:
                avg_read_length = float(c)
            else:
                reads2markers[r] = c
    elif input_type == 'sam':
        n_metagenome_reads = nreads 
        for line in inpf:
            o = ras_line(line)
            if ((o[0][0] != '@') and #no header
                (o[2][-1] != '*') and # no unmapped reads
                (hex(int(o[1]) & 0x100) == '0x0') and #no secondary
                mapq_filter(o[2], int(o[4]), min_mapq_val) and # filter low mapq reads
                ( (min_alignment_len is None) or ( max(int(x.strip('M')) for x in re.findall(r'(\d*M)', o[5]) if x) >= min_alignment_len ) )
            ):
                    reads2markers[o[0]] = o[2].split('/')[0]
    inpf.close()
    
    if subsampling is not None and mapping_subsampling:
        if subsampling >= n_metagenome_reads:
            sys.stderr.write("WARNING: The specified subsampling ({}) is equal or higher than the original number of reads ({}). Subsampling will be skipped.\n".format(subsampling, n_metagenome_reads))
        else:
            reads2markers =  dict(sorted(reads2markers.items()))
            if subsampling_seed.lower() != 'random':
                random.seed(int(subsampling_seed))
            reads2filtmarkers = {}
            sgb_reads2markers, viral_reads2markers = separate_reads2markers(reads2markers)            
            n_sgb_mapped_reads = int((len(sgb_reads2markers) * subsampling) / n_metagenome_reads)
            reads2filtmarkers = { r:sgb_reads2markers[r] for r in random.sample(list(sgb_reads2markers.keys()), n_sgb_mapped_reads) }     
            if SGB_ANALYSIS:       
                n_viral_mapped_reads = int((len(viral_reads2markers) * subsampling) / n_metagenome_reads)
                reads2filtmarkers.update({ r:viral_reads2markers[r] for r in random.sample(list(viral_reads2markers.keys()), n_viral_mapped_reads) })            
            reads2markers = reads2filtmarkers
            sgb_reads2markers.clear()
            viral_reads2markers.clear()
            n_metagenome_reads = subsampling
    elif subsampling is None and n_metagenome_reads < 10000:
        sys.stderr.write("WARNING: The number of reads in the sample ({}) is below the recommended minimum of 10,000 reads.\n".format(n_metagenome_reads))
        
    markers2reads = defdict(set)   
    for r, m in reads2markers.items():
        markers2reads[m].add(r)

    return (markers2reads, n_metagenome_reads, avg_read_length)

def maybe_generate_biom_file(tree, pars, abundance_predictions):
    json_key = "MetaPhlAn"

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
        return end_name.startswith("s__") or end_name.endswith("_unclassified")

    def findclade(clade_name):
        if clade_name.endswith('_unclassified'):
            name = clade_name.split(delimiter)[-2]
        else:
            name = clade_name.split(delimiter)[-1]
        return tree.all_clades[name]

    def to_biomformat(clade_name):
        return {'taxonomy': clade_name.split(delimiter)}

    clades = iter((abundance, findclade(name))
                  for (name, taxid, abundance) in abundance_predictions if istip(name))
    packed = iter(([abundance], clade.get_full_name(), clade.tax_id)
                  for (abundance, clade) in clades)

    # unpack that tuple here to stay under 80 chars on a line
    data, clade_names, _ = zip(*packed)
    # biom likes column vectors, so we give it an array like this:
    # np.array([a],[b],[c])
    data = np.array(data)
    sample_ids = [pars['sample_id']]
    table_id = 'MetaPhlAn_Analysis'




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
            observation_metadata = list(map(to_biomformat, clade_names)),
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
            observation_metadata = list(map(to_biomformat, clade_names)),
            table_id             = table_id,
            input_is_dense       = True
        )

        with open(pars['biom'], 'w') as outfile:
            biom_table.to_json( json_key,
                                direct_io = outfile )

    return True

def _make_gen_fastq(reader):
    b = reader(1024 * 1024) 
    while (b):
        yield b
        b = reader(1024 * 1024)

def rawpycount(filename):
    f = bz2.BZ2File(filename, 'rb') if filename.endswith(".bz2") else gzip.open(filename, 'rb') if filename.endswith('.gz') else open(filename, 'rb')
    f_gen = _make_gen_fastq(f.read)
    n_metagenome_reads = sum( buf.count(b'\n') for buf in f_gen)
    f.close()
    return n_metagenome_reads

def subsample_file(inp, paired, subsampling, out, sample):         
    length, read_num = 0, -1
    counter=0
    out_l = out
    for inp_f in inp.split(','):
        if paired:
            length, read_num = 0, -1  
            out=out_l[counter]
            counter+=1
        line = 1
        reader = bz2.open(inp_f, 'rt') if inp_f.endswith(".bz2") else gzip.open(inp_f, 'rt') if inp_f.endswith('.gz') else open(inp_f, 'r')
        while (line and length < subsampling):
            line = reader.readline()
            read_num += 1
            if read_num in sample:
                length += 1
                out.write(line)       
                for _ in range(3):
                    out.write(reader.readline())
            else:
                for _ in range(3):
                    reader.readline()
        reader.close()
    out.close()
    
def subsample_reads(inp, subsampling, subsampling_seed, subsampling_output, tmp_dir, paired):

    n_metagenome_reads = execute_pool(((rawpycount, inp_f) for inp_f in inp.split(',')), 2)
    n_metagenome_reads = [int(n/4) for n in n_metagenome_reads]

    if paired:
        subsampling //= 2
        if n_metagenome_reads[0] != n_metagenome_reads[1]:
            sys.stderr.write("Error: The specified reads file are not the same length! Make sure the forward and reverse reads are files are not damaged and reads are in the same order. Exiting ...\n\n") 
            sys.exit(1)
        n_metagenome_reads = n_metagenome_reads[0]
    else:
        n_metagenome_reads = sum(n_metagenome_reads)

    if subsampling >= n_metagenome_reads:
        sys.stderr.write("WARNING: The specified subsampling ({}) is equal or higher than the original number of reads ({}). Subsampling will be skipped.\n".format(subsampling, sum(n_metagenome_reads))) 
        return inp, None

    if subsampling_seed.lower() != 'random':
        random.seed(int(subsampling_seed))
        
    sample = set(random.sample(range(n_metagenome_reads), subsampling))
    
    if subsampling_output is None:
        if paired:
            subsampling_output = list()
            out=list()
            for _ in inp.split(','):
                out_f = tf.NamedTemporaryFile(dir=tmp_dir, mode='w', delete=False)
                out.append(out_f)
                subsampling_output.append(out_f.name)
        else:
            out = tf.NamedTemporaryFile(dir=tmp_dir, mode='w', delete=False)
            subsampling_output=out.name
    else:
        if paired:
            r, ext = os.path.splitext(subsampling_output)

            subsampling_output = list()
            out = list()

            subsampling_output.append('.'.join([r, 'R1'+ext]))
            subsampling_output.append('.'.join([r, 'R2'+ext]))

            for s in subsampling_output:
                out.append(bz2.open(s, 'wt') if s.endswith(".bz2") else gzip.open(s, 'wt') if s.endswith('.gz') else open(s, 'w'))
        else:
            out = bz2.open(subsampling_output, 'wt') if subsampling_output.endswith(".bz2") else gzip.open(subsampling_output, 'wt') if subsampling_output.endswith('.gz') else open(subsampling_output, 'w')
            
    subsample_file(inp, paired, subsampling, out, sample)

    if isinstance(subsampling_output, list):
        subsampling_output = ','.join(subsampling_output)
        
    return subsampling_output, subsampling

def main():
    ranks2code = { 'k' : 'superkingdom', 'p' : 'phylum', 'c':'class',
                   'o' : 'order', 'f' : 'family', 'g' : 'genus', 's' : 'species'}
    pars = read_params(sys.argv)

    #Set SGB- / species- analysis
    global SGB_ANALYSIS
    SGB_ANALYSIS = not pars['mpa3']

    ESTIMATE_UNK = pars['unclassified_estimation']

    if pars['subsampling'] and pars['subsampling_paired']:
        sys.stderr.write("Error: You specified both --subsampling and --subsampling_paired. Choose only one of the two options. Exiting...")
        sys.exit(1)

    if pars['subsampling']:
        sys.stderr.write("WARNING: If you use --subsampling the reads will be subsampled NOT taking into account paired information. Use --subsampling_paired if you have paired ends reads.\n".format(pars['subsampling']))
    
    if pars['mapping_subsampling'] and pars['subsampling'] is None:
        sys.stderr.write("Error: The --mapping_subsampling parameter should be used together with the --subsampling parameter. Exiting...\n\n")
        sys.exit(1)

    if pars['subsampling_paired']:
        subsampling_paired=True
        pars['subsampling']=pars['subsampling_paired']
    else:
        subsampling_paired=False

    if subsampling_paired and (pars['1'] is None or pars['2'] is None):
        sys.stderr.write("Error: If you specify --subsampling_paired you have to provide forward and reverse reads as -1 and -2 respectively. Reads are assumed to be in the same order in the two files.\n Exiting...\n\n".format(pars['subsampling']))
        sys.exit(1)
    
    if pars['subsampling'] is not None and pars['subsampling'] < 10000:
        sys.stderr.write("WARNING: The specified subsampling ({}) is below the recommended minimum of 10,000 reads.\n".format(pars['subsampling']))
    
    if not (pars['subsampling_seed'].lower() == 'random' or pars['subsampling_seed'].isdigit()):
        sys.stderr.write("Error: The --subsampling_seed parameter is not accepted. It should contain an integer number or \"random\". Exiting...\n\n")
        sys.exit(1)
        
    if not pars['mapping_subsampling'] and pars['subsampling'] is not None:
        if pars['input_type'] != 'fastq':
            sys.stderr.write("Error: The reads' subsampling procedure requires fastq input! Exiting...\n\n")
            sys.exit(1)
        elif (not pars['inp']) and (not subsampling_paired):
            sys.stderr.write("Error: Input reads for the subsampling must be provided as parameter. Stdin input is not allowed! Exiting...\n\n")
            sys.exit(1)
        
        if subsampling_paired:
            if not os.path.exists(pars['2']) or not os.path.exists(pars['1']):
                sys.stderr.write("Error: Files passed with -1 ({}) or -2 ({}) not found. Exiting...\n\n".format(pars['1'],pars['2']))
                sys.exit(1)

            if pars['inp']:
                sys.stderr.write("WARNING: since --subsampling_paired has been specified, reads are taken from -1 ({}) and -2 ({}), not from -inp.\n".format(pars['1'],pars['2']))
            pars['inp'] = pars['1']+','+pars['2']

        pars['inp'], pars['subsampling'] = subsample_reads(pars['inp'], pars['subsampling'], pars['subsampling_seed'], pars['subsampling_output'], pars['tmp_dir'], paired=subsampling_paired)
        
    # check if the database is installed, if not then install
    pars['index'] = check_and_install_database(pars['index'], pars['bowtie2db'], pars['bowtie2_build'], pars['nproc'], pars['force_download'], pars['offline'])

    if pars['install']:
        sys.stderr.write('The database is installed\n')
        return

    #if we are to profile viral clusters:

    if pars['profile_vsc']:
        #TODO: comply to other input types
        if pars['input_type'] != 'fastq':
            sys.stderr.write("The Viral Sequence Clusters mode requires fastq input!\n")
            sys.exit(1)
        
        if not pars['vsc_out']:
            sys.stderr.write("The Viral Sequence Clusters mode requires to specify an output file with the --vsc_out parameter\n")
            sys.exit(1)

        tmpvirdir=tf.TemporaryDirectory(dir=pars['tmp_dir'])
        viralTempFolder=tmpvirdir.name
        
        # Uncomment to preserve tmp dir
        #tmpvirdir=tf.mkdtemp(dir=pars['tmp_dir'])
        #viralTempFolder=tmpvirdir
        #print("VF: ",viralTempFolder)

        vsc_fna, vsc_vinfo = set_vsc_parameters(pars['index'], pars['bowtie2db'])

        # implies SGBs (uses the same database with SGBs, so this is needed to avoid
        # the need to specify --sgb when using --profile_vscs)
        SGB_ANALYSIS = True


    else:
        viralTempFolder = None


    # set correct map_pkl and bowtie2db variables
    pars['mpa_pkl'], pars['bowtie2db'] = set_mapping_arguments(pars['index'], pars['bowtie2db'])

    if (pars['bt2_ps'] in ["sensitive-local", "very-sensitive-local"]) and (pars['min_alignment_len'] is None):
            pars['min_alignment_len'] = 100
            sys.stderr.write('Warning! bt2_ps is set to local mode, and min_alignment_len is None, I automatically '
                             'set min_alignment_len to 100! If you do not like, rerun the command and set '
                             'min_alignment_len to a specific value.\n')

    # check for the mpa_pkl file
    if not os.path.isfile(pars['mpa_pkl']):
        sys.stderr.write("Error: Unable to find the mpa_pkl file at: " + pars['mpa_pkl'] +
                         "Exiting...\n\n")
        sys.exit(1)

    if pars['ignore_markers']:
        with open(pars['ignore_markers']) as ignv:
            ignore_markers = set([l.strip() for l in ignv])
    else:
        ignore_markers = set()

    no_map = False
    if pars['input_type'] == 'fasta' or pars['input_type'] == 'fastq':
        bow = pars['bowtie2db'] is not None

        if not bow:
            sys.stderr.write( "No MetaPhlAn BowTie2 database provided\n "
                              "[--bowtie2db and --index options]!\n"
                              "Exiting...\n\n" )
            sys.exit(1)

        if pars['no_map']:
            pars['bowtie2out'] = tf.NamedTemporaryFile(dir=pars['tmp_dir']).name
            no_map = True
        else:
            if bow and not pars['bowtie2out']:
                if pars['inp'] and "," in  pars['inp']:
                    sys.stderr.write("Error! --bowtie2out needs to be specified when multiple "
                                     "FASTQ or FASTA files (comma separated) are provided\n")
                    sys.exit(1)
                fname = pars['inp']
                if fname is None:
                    fname = "stdin_map"
                elif stat.S_ISFIFO(os.stat(fname).st_mode):
                    fname = "fifo_map"
                pars['bowtie2out'] = fname + ".bowtie2out.txt"

            if os.path.exists( pars['bowtie2out'] ) and not pars['force']:
                sys.stderr.write(
                    "BowTie2 output file detected: " + pars['bowtie2out'] + "\n"
                    "Please use it as input or remove it if you want to "
                    "re-perform the BowTie2 run.\n"
                    "Exiting...\n\n" )
                sys.exit(1)
            if pars['force']:
                if os.path.exists(pars['bowtie2out']):
                    os.remove( pars['bowtie2out'] )

        bt2_ext = 'bt2l' if SGB_ANALYSIS else 'bt2'            
        if bow and not all([os.path.exists(".".join([str(pars['bowtie2db']), p]))
                            for p in ["1." + bt2_ext, "2." + bt2_ext, "3." + bt2_ext, "4." + bt2_ext, "rev.1." + bt2_ext, "rev.2." + bt2_ext]]):
            sys.stderr.write("No MetaPhlAn BowTie2 database found (--index "
                             "option)!\nExpecting location {}\nExiting..."
                             .format(pars['bowtie2db']))
            sys.exit(1)
        if bow and not (abs(os.path.getsize(".".join([str(pars['bowtie2db']), "1." + bt2_ext])) - os.path.getsize(".".join([str(pars['bowtie2db']), "rev.1." + bt2_ext]))) <= 1000):
            sys.stderr.write("Partial MetaPhlAn BowTie2 database found at {}. "
                             "Please remove and rebuild the database.\nExiting..."
                             .format(pars['bowtie2db']))
            sys.exit(1)

        if bow:
            run_bowtie2(pars['inp'], pars['bowtie2out'], pars['bowtie2db'],
                                pars['bt2_ps'], pars['nproc'], file_format=pars['input_type'],
                                exe=pars['bowtie2_exe'], samout=pars['samout'],

                                min_alignment_len=pars['min_alignment_len'], read_min_len=pars['read_min_len'], min_mapq_val=pars['min_mapq_val'],profile_vsc_folder=viralTempFolder)
            if pars['subsampling_output'] is None and not pars['mapping_subsampling'] and pars['subsampling'] is not None:
                for inp_f in pars['inp'].split(','):
                    os.remove(inp_f)
            pars['input_type'] = 'bowtie2out'
        pars['inp'] = pars['bowtie2out'] # !!!
    with bz2.BZ2File( pars['mpa_pkl'], 'r' ) as a:
        mpa_pkl = pickle.load( a )

    REPORT_MERGED = mpa_pkl.get('merged_taxon',False)
    tree = TaxTree( mpa_pkl, ignore_markers )
    tree.set_min_cu_len( pars['min_cu_len'] )

    if pars['input_type'] == 'sam' and not pars['nreads']:
        sys.stderr.write(
                "Please provide the size of the metagenome using the "
                "--nreads parameter when running MetaPhlAn using SAM files as input"
                "\nExiting...\n\n" )
        sys.exit(1)

    markers2reads, n_metagenome_reads, avg_read_length = map2bbh(pars['inp'], pars['min_mapq_val'], pars['input_type'], pars['min_alignment_len'], pars['nreads'], pars['mapping_subsampling'], pars['subsampling'], pars['subsampling_seed'])

    if pars['profile_vsc']:
        
        try:
            VSCs_markers = SeqIO.index(vsc_fna, "fasta")
        except Exception as e:
            sys.stderr.write('Error: "{}"\nCould not access to {}.\n'.format(e,vsc_fna))
            sys.exit(1)

        with open(viralTempFolder+'/viralmk.txt') as viral_markers_to_remap:

            allViralMarkers = {}
            for vmarker in viral_markers_to_remap:
                group,marker = vmarker.strip().split()
                if group not in allViralMarkers:
                    allViralMarkers[group] = [marker]
                else:
                    allViralMarkers[group].append(marker)

            selectedMakrers=[]
            for grp,v in allViralMarkers.items():

                cv=Counter(v)
                if (len(cv) > 1):
                    normalized_occurrencies=sorted([(mk,mk_occurrencies,len(VSCs_markers[mk].seq)) for (mk,mk_occurrencies) in cv.items()],key=lambda x: x[1]/x[2],reverse=True)
                    topMarker= VSCs_markers[normalized_occurrencies[0][0]] #name of the viralGenome
                else:
                    topMarker=VSCs_markers[v[0]] # the only hitted viralGenome is the winning one

                selectedMakrers.append(topMarker)


            SeqIO.write( selectedMakrers,viralTempFolder+'/v_mks.fa','fasta' )

            VSC_report = vsc_bowtie2(viralTempFolder, pars['nproc'], file_format=pars['input_type'],
                        exe=pars['bowtie2_exe'], bt2build_exe=pars['bowtie2_build']) 


            if not VSC_report == []:
                with open(pars['vsc_out'],'w') as outf:

                    outf.write('#{}\n'.format(pars['index']))
                    outf.write('#{}\n'.format(' '.join(sys.argv)))
                    outf.write('#SampleID\t{}\n'.format(pars['sample_id']))
                    vsc_out_df = pd.DataFrame.from_dict(VSC_report).query('breadth_of_coverage >= {}'.format(pars['vsc_breadth']))
                    vsc_info_df = pd.read_table(vsc_vinfo, sep='\t')
                    vsc_out_df = vsc_out_df.merge(vsc_info_df, on='M-Group/Cluster').sort_values(by='breadth_of_coverage', ascending=False).set_index('M-Group/Cluster')

                    vsc_out_df.to_csv(outf,sep='\t',na_rep='-')


    tree.set_stat( pars['stat'], pars['stat_q'], pars['perc_nonzero'], avg_read_length, pars['avoid_disqm'])

    if no_map:
        os.remove( pars['inp'] )

    map_out = []
    for marker,reads in sorted(markers2reads.items(), key=lambda pars: pars[0]):
        if marker not in tree.markers2lens:
            continue
        tax_seq, ids_seq = tree.add_reads( marker, len(reads),
                                  add_viruses = pars['add_viruses'],
                                  ignore_eukaryotes = pars['ignore_eukaryotes'],
                                  ignore_bacteria = pars['ignore_bacteria'],
                                  ignore_archaea = pars['ignore_archaea'],
                                  ignore_ksgbs = pars['ignore_ksgbs'],
                                  ignore_usgbs = pars['ignore_usgbs']
                                  )
        if tax_seq:
            map_out +=["\t".join([r,tax_seq, ids_seq]) for r in sorted(reads)]

    if pars['output'] is None and pars['output_file'] is not None:
        pars['output'] = pars['output_file']

    out_stream = open(pars['output'],"w") if pars['output'] else sys.stdout
    MPA2_OUTPUT = pars['legacy_output']
    CAMI_OUTPUT = pars['CAMI_format_output']

    with out_stream as outf:
        if not MPA2_OUTPUT:
            outf.write('#{}\n'.format(pars['index']))
            outf.write('#{}\n'.format(' '.join(sys.argv)))
            outf.write('#{} reads processed\n'.format(n_metagenome_reads))
        
        if pars['t'] == 'rel_ab_w_read_stats':
            outf.write('#Average read length {}\n'.format(avg_read_length))           

        if not CAMI_OUTPUT:
            outf.write('#' + '\t'.join((pars["sample_id_key"], pars["sample_id"])) + '\n')

        if ESTIMATE_UNK:
            mapped_reads = 0
            cl2pr = tree.clade_profiles( pars['tax_lev']+"__" if pars['tax_lev'] != 'a' else None  )
            cl2ab, _, _ = tree.relative_abundances( pars['tax_lev']+"__" if pars['tax_lev'] != 'a' else None )
            confident_taxa = [taxstr for (taxstr, _),relab in cl2ab.items() if relab > 0.0]
            for c, m in cl2pr.items():
                if c in confident_taxa:
                    markers_cov = [a  / 1000 for _, a in m if a > 0]
                    mapped_reads += np.mean(markers_cov) * tree.all_clades[c.split('|')[-1]].glen
            # If the mapped reads are over-estimated, set the ratio at 1
            fraction_mapped_reads = min(mapped_reads/float(n_metagenome_reads), 1.0)
        else:
            fraction_mapped_reads = 1.0
      
        if pars['t'] == 'reads_map':
            if not MPA2_OUTPUT:
               outf.write('#read_id\tNCBI_taxlineage_str\tNCBI_taxlineage_ids\n')
            outf.write( "\n".join( map_out ) + "\n" )

        elif pars['t'] == 'rel_ab':
            if CAMI_OUTPUT:
                outf.write('''@SampleID:{}\n@Version:0.10.0\n@Ranks:superkingdom|phylum|class|order|family|genus|species|strain\n@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n'''.format(pars["sample_id"]))
            if not MPA2_OUTPUT and not CAMI_OUTPUT:
                if not pars['use_group_representative']:
                    outf.write('#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species\n')
                else:
                    outf.write('#clade_name\tNCBI_tax_id\trelative_abundance\n')

            cl2ab, _, tot_nreads = tree.relative_abundances(
                        pars['tax_lev']+"__" if pars['tax_lev'] != 'a' else None )
            
            outpred = [(taxstr, taxid,round(relab*100.0,5)) for (taxstr, taxid), relab in cl2ab.items() if relab > 0.0]
            has_repr = False
            
            if outpred:
                if CAMI_OUTPUT:
                    for clade, taxid, relab in sorted(  outpred, reverse=True,
                                        key=lambda x:x[2]+(100.0*(8-(x[0].count("|"))))):
                        if taxid and clade.split('|')[-1][0] != 't':
                            rank = ranks2code[clade.split('|')[-1][0]]
                            leaf_taxid = taxid.split('|')[-1]
                            taxpathsh = '|'.join([remove_prefix(name) if '_unclassified' not in name else '' for name in clade.split('|')])
                            outf.write( '\t'.join( [ leaf_taxid, rank, taxid, taxpathsh, str(relab*fraction_mapped_reads) ] ) + '\n' )
                else:
                    if ESTIMATE_UNK:
                        outf.write( "\t".join( [    "UNCLASSIFIED",
                                                    "-1",
                                                    str(round((1-fraction_mapped_reads)*100,5)),""]) + "\n" )
                                                    
                    for clade, taxid, relab in sorted(  outpred, reverse=True,
                                        key=lambda x:x[2]+(100.0*(8-(x[0].count("|"))))):
                        add_repr = ''
                        if REPORT_MERGED and (clade, taxid) in mpa_pkl['merged_taxon']:
                            if pars['use_group_representative'] and not SGB_ANALYSIS:
                                if '_group' in clade:
                                    clade, taxid, _ = sorted(mpa_pkl['merged_taxon'][(clade, taxid)], key=lambda x:x[2], reverse=True)[0]
                            elif not pars['use_group_representative']:
                                add_repr = '{}'.format(','.join( [ n[0] for n in mpa_pkl['merged_taxon'][(clade, taxid)]] ))
                                has_repr = True
                        if not MPA2_OUTPUT:
                            outf.write( "\t".join( [clade, 
                                                    taxid, 
                                                    str(relab*fraction_mapped_reads), 
                                                    add_repr
                                                ] ) + "\n" )
                        else:
                            outf.write( "\t".join( [clade, 
                                                    str(relab*fraction_mapped_reads)] ) + "\n" )
                if REPORT_MERGED and has_repr:
                    sys.stderr.write("WARNING: The metagenome profile contains clades that represent multiple species merged into a single representant.\n"
                                     "An additional column listing the merged species is added to the MetaPhlAn output.\n"
                                    )
            else:
                if not MPA2_OUTPUT:
                    outf.write( "UNCLASSIFIED\t-1\t100.0\t\n" )
                else:
                    outf.write( "UNCLASSIFIED\t100.0\n" )
                sys.stderr.write("WARNING: MetaPhlAn did not detect any microbial taxa in the sample.\n")
            maybe_generate_biom_file(tree, pars, outpred)

        elif pars['t'] == 'rel_ab_w_read_stats':
            cl2ab, rr, tot_nreads = tree.relative_abundances(
                        pars['tax_lev']+"__" if pars['tax_lev'] != 'a' else None )

            unmapped_reads = max(n_metagenome_reads - tot_nreads, 0)

            outpred = [(taxstr, taxid,round(relab*100.0*fraction_mapped_reads,5)) for (taxstr, taxid),relab in cl2ab.items() if relab > 0.0]

            if outpred:
                outf.write( "#estimated_reads_mapped_to_known_clades:{}\n".format(round(tot_nreads)) )
                outf.write( "\t".join( [    "#clade_name",
                                            "clade_taxid",
                                            "relative_abundance",
                                            "coverage",
                                            "estimated_number_of_reads_from_the_clade" ]) +"\n" )
                if ESTIMATE_UNK:
                    outf.write( "\t".join( [    "UNCLASSIFIED",
                                                "-1",
                                                str(round((1-fraction_mapped_reads)*100,5)),
                                                "-",
                                                str(round(unmapped_reads)) ]) + "\n" )
                                                
                for taxstr, taxid, relab in sorted(  outpred, reverse=True,
                                    key=lambda x:x[2]+(100.0*(8-(x[0].count("|"))))):
                    outf.write( "\t".join( [    taxstr,
                                                taxid,
                                                str(relab),
                                                str(round(rr[(taxstr, taxid)][0],5)) if (taxstr, taxid) in rr else '-',          #coverage
                                                str( int( round( rr[(taxstr, taxid)][1], 0) )  if (taxstr, taxid) in rr else '-')       #estimated_number_of_reads_from_the_clade
                                                ] ) + "\n" )
            else:
                if not MPA2_OUTPUT:
                    outf.write( "#estimated_reads_mapped_to_known_clades:0\n")
                    outf.write( "\t".join( [    "#clade_name",
                                                "clade_taxid",
                                                "relative_abundance",
                                                "coverage",
                                                "estimated_number_of_reads_from_the_clade" ]) +"\n" )
                    outf.write( "unclassified\t-1\t100.0\t0\t0\n" )
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
            cl2ab, _, _ = tree.relative_abundances( None )
            strout = []
            for (taxstr, taxid), relab in cl2ab.items():
                clade = taxstr
                if clade.endswith(pars['clade']) and relab*100.0 < pars['min_ab']:
                    strout = []
                    break
                if pars['clade'] in clade:
                    strout += ["\t".join([str(a),str(int(b > pars['pres_th']))]) for a,b in cl2pr[clade]]
            if strout:
                strout = sorted(strout,key=lambda x:x[0])
                outf.write( "\n".join(strout) + "\n" )
            else:
                sys.stderr.write("Clade "+pars['clade']+" not present at an abundance >"+str(round(pars['min_ab'],2))+"%, "
                                 "so no clade specific markers are reported\n")


if __name__ == '__main__':
    t0 = time.time()
    if '-v' in sys.argv or '--version' in sys.argv:
        print(f"MetaPhlAn version {__version__} ({__date__})")
        db_dir=get_cli_db_dir()
        db_version=get_installed_db_version(db_dir)
        if db_version:
            print("Installed databases:", ", ".join(db_version))
        else:
            print("No complete MetaPhlAn Bowtie2 database found")
    else: 
        main()
        sys.stderr.write('Elapsed time to run MetaPhlAn: {} s\n'.format( (time.time()-t0) ) )

