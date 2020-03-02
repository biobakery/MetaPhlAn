#!/usr/bin/env python3


__author__ = ('Francesco Asnicar (f.asnicar@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it), '
              'Claudia Mengoni (claudia.mengoni@studenti.unitn.it), '
              'Mattia Bolzan (mattia.bolzan@unitn.it), '
              'Nicola Segata (nicola.segata@unitn.it)')
__version__ = '0.41'
__date__ = '31 December 2019'


import os
import sys
import glob
import shutil
import argparse as ap
import configparser as cp
import subprocess as sb
import multiprocessing as mp
from Bio import SeqIO  # Biopython requires NumPy
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from collections import Counter
import bz2
import math
import re
import time
import pickle
from itertools import combinations
import dendropy  # DendroPy
from urllib.request import urlretrieve
import tarfile
import hashlib
import gzip
import random as lib_random


if sys.version_info[0] < 3:
    raise Exception("PhyloPhlAn2 requires Python 3, your current Python version is {}.{}.{}"
                    .format(sys.version_info[0], sys.version_info[1], sys.version_info[2]))

CONFIG_SECTIONS_MANDATORY = [['map_dna', 'map_aa'], ['msa'], ['tree1']]
CONFIG_SECTIONS_ALL = ['map_dna', 'map_aa', 'msa', 'trim', 'gene_tree1', 'gene_tree2', 'tree1', 'tree2']
CONFIG_OPTIONS_MANDATORY = [['program_name'], ['command_line']]
CONFIG_OPTIONS_ALL = ['program_name', 'params', 'threads', 'input', 'database',
                      'output_path', 'output', 'version', 'environment', 'command_line']
CONFIG_OPTIONS_TO_EXCLUDE = ['version', 'environment']
INPUT_FOLDER = 'input/'
DATABASES_FOLDER = 'phylophlan2_databases/'
SUBMAT_FOLDER = 'phylophlan2_substitution_matrices/'
SUBMOD_FOLDER = 'phylophlan2_substitution_models/'
CONFIGS_FOLDER = 'phylophlan2_configs/'
OUTPUT_FOLDER = ''
MIN_NUM_PROTEINS = 1
MIN_LEN_PROTEIN = 50
MIN_NUM_MARKERS = 1
DB_TYPE_CHOICES = ['n', 'a']
TRIM_CHOICES = ['gap_trim', 'gap_perc', 'not_variant', 'greedy']
SUBSAMPLE_CHOICES = ['phylophlan', 'onethousand', 'sevenhundred', 'fivehundred', 'threehundred', 'onehundred',
                     'fifty', 'twentyfive', 'tenpercent', 'twentyfivepercent', 'fiftypercent']
SCORING_FUNCTION_CHOICES = ['trident', 'muscle', 'random']
DIVERSITY_CHOICES = ['low', 'medium', 'high']
MIN_NUM_ENTRIES = 4
GENOME_EXTENSION = '.fna'
PROTEOME_EXTENSION = '.faa'
NOT_VARIANT_THRESHOLD = 0.99
GAP_PERC_THRESHOLD = 0.67
FRAGMENTARY_THRESHOLD = 0.85
UNKNOWN_FRACTION = 0.3
DATABASE_DOWNLOAD_URL = "https://bitbucket.org/nsegata/phylophlan/downloads/"


def info(s, init_new_line=False, exit=False, exit_value=0):
    if init_new_line:
        sys.stdout.write('\n')

    sys.stdout.write('{}'.format(s))
    sys.stdout.flush()

    if exit:
        sys.exit(exit_value)


def error(s, init_new_line=False, exit=False, exit_value=1):
    if init_new_line:
        sys.stderr.write('\n')

    sys.stderr.write('[e] {}\n'.format(s))
    sys.stderr.flush()

    if exit:
        sys.exit(exit_value)


def read_params():
    p = ap.ArgumentParser(description=("PhyloPhlAn2 is an accurate, rapid, and easy-to-use method for large-scale microbial genome "
                                       "characterization and phylogenetic analysis at multiple levels of resolution. PhyloPhlAn2 can assign "
                                       "finished, draft, or metagenome-assembled genomes (MAGs) to species-level genome bins (SGBs). For "
                                       "individual clades of interest (e.g. newly sequenced genome sets), PhyloPhlAn2 reconstructs strain-level "
                                       "phylogenies from among the closest species using clade-specific maximally informative markers. At the "
                                       "other extreme of resolution, PhyloPhlAn2 scales to very-large phylogenies comprising >17,000 microbial "
                                       "species"),
                          formatter_class=ap.ArgumentDefaultsHelpFormatter)

    group = p.add_mutually_exclusive_group()
    group.add_argument('-i', '--input', metavar='PROJECT_NAME', type=str, default=None,
                       help="")
    group.add_argument('-c', '--clean', type=str, default=None,
                       help="Clean the output and partial data produced for the specified project")

    p.add_argument('-o', '--output', type=str, default=None,
                   help=("Output folder name, otherwise it will be the name of the input folder concatenated with the name of "
                         "the database used"))
    p.add_argument('-d', '--database', type=str, default=None, help="The name of the database of markers to use.")
    p.add_argument('-t', '--db_type', default=None, choices=DB_TYPE_CHOICES,
                   help=('Specify the type of the database of markers, where "n" stands for nucleotides and "a" for amino '
                         'acids. If not specified, PhyloPhlAn2 will automatically detect the type of database'))
    p.add_argument('-f', '--config_file', type=str, default=None,
                   help=('The configuration file to load, four ready-to-use configuration files can be generated using the '
                         '"write_default_configs.sh" script present in the "configs" folder'))
    p.add_argument('--diversity', default=None, choices=DIVERSITY_CHOICES, required=True,
                   help=('Specify the expected diversity of the phylogeny, automatically adjust some parameters: '
                         '"low": for genus-/species-/strain-level phylogenies; '
                         '"medium": for class-/order-level phylogenies; '
                         '"high": for phylum-/tree-of-life size phylogenies'))

    group = p.add_mutually_exclusive_group()
    group.add_argument('--accurate', action='store_true', default=False,
                       help=('Use more phylogenetic signal which can result in more accurate phylogeny; '
                             'affected parameters depend on the "--diversity" level'))
    group.add_argument('--fast', action='store_true', default=False,
                       help=('Perform more a faster phylogeny reconstruction by reducing the phylogenetic positions to use; '
                             'affected parameters depend on the "--diversity" level'))

    p.add_argument('--clean_all', action='store_true', default=False,
                   help="Remove all installation and database files automatically generated")
    p.add_argument('--database_list', action='store_true', default=False,
                   help="List of all the available databases that can be specified with the -d/--database option")
    p.add_argument('-s', '--submat', type=str, default=None,
                   help='Specify the substitution matrix to use, available substitution matrices can be listed with "--submat_list"')
    p.add_argument('--submat_list', action='store_true', default=False,
                   help="List of all the available substitution matrices that can be specified with the -s/--submat option")
    p.add_argument('--submod_list', action='store_true', default=False,
                   help="List of all the available substitution models that can be specified with the --maas option")
    p.add_argument('--nproc', type=int, default=1, help="The number of cores to use")
    p.add_argument('--min_num_proteins', type=int, default=MIN_NUM_PROTEINS,
                   help=("Proteomes with less than this number of proteins will be discarded"))
    p.add_argument('--min_len_protein', type=int, default=MIN_LEN_PROTEIN,
                   help="Proteins in proteomes shorter than this value will be discarded")
    p.add_argument('--min_num_markers', type=int, default=MIN_NUM_MARKERS,
                   help="Input genomes or proteomes that map to less than the specified number of markers will be discarded")
    p.add_argument('--trim', default=None, choices=TRIM_CHOICES,
                   help=('Specify which type of trimming to perform: '
                         '"gap_trim": execute what specified in the "trim" section of the configuration file; '
                         '"gap_perc": remove columns with a percentage of gaps above a certain threshold '
                                      '(see "--gap_perc_threshold" parameter)'
                         '"not_variant": remove columns with at least one nucleotide/amino acid appearing above a certain threshold '
                                        '(see "--not_variant_threshold" parameter); '
                         '"greedy": performs all the above trimming steps; '
                         'If not specified, no trimming will be performed'))
    p.add_argument('--gap_perc_threshold', type=float, default=GAP_PERC_THRESHOLD,
                   help='Specify the value used to consider a column not variant when "--trim not_variant" is specified')
    p.add_argument('--not_variant_threshold', type=float, default=NOT_VARIANT_THRESHOLD,
                   help='Specify the value used to consider a column not variant when "--trim not_variant" is specified')
    p.add_argument('--subsample', default=None, choices=SUBSAMPLE_CHOICES,
                   help=('The number of positions to retain from each single marker, available option are: '
                         '"phylophlan": specific number of positions for each PhyloPhlAn marker (only when "--database phylophlan"); '
                         '"onethousand": return the top 1000 positions; '
                         '"sevenhundred": return the top 700; '
                         '"fivehundred": return the top 500; '
                         '"threehundred" return the top 300; '
                         '"onehundred": return the top 100 positions; '
                         '"fifty": return the top 50 positions; '
                         '"twentyfive": return the top 25 positions; '
                         '"fiftypercent": return the top 50 percent positions; '
                         '"twentyfivepercent": return the top 25% positions; '
                         '"tenpercent": return the top 10% positions; '
                         'If not specified, the complete alignment will be used').replace(r"%", r"%%"))
    p.add_argument('--unknown_fraction', type=float, default=UNKNOWN_FRACTION,
                   help='Define the amount of unknowns ("X" and "-") allowed in each column of the MSA of the markers')
    p.add_argument('--scoring_function', default=None, choices=SCORING_FUNCTION_CHOICES,
                   help="Specify which scoring function to use to evaluate columns in the MSA results")
    p.add_argument('--sort', action='store_true', default=False,
                   help=('If specified, the markers will be ordered, when using the '
                         'PhyloPhlAn database, it will be automatically set to "True"'))
    p.add_argument('--remove_fragmentary_entries', action='store_true', default=False,
                   help=("If specified the MSAs will be checked and cleaned from fragmentary entries. See --fragmentary_threshold "
                         "for the threshold values above which an entry will be considered fragmentary"))
    p.add_argument('--fragmentary_threshold', type=float, default=FRAGMENTARY_THRESHOLD,
                   help="The fraction of gaps in the MSA to be considered fragmentary and hence discarded")
    p.add_argument('--min_num_entries', type=int, default=MIN_NUM_ENTRIES,
                   help="The minimum number of entries to be present for each of the markers in the database")
    p.add_argument('--maas', type=str, default=None,
                   help=("Select a mapping file that specifies the substitution model of amino acid to use for each of the markers "
                         "for the gene tree reconstruction. File must be tab-separated"))
    p.add_argument('--remove_only_gaps_entries', action='store_true', default=False,
                   help=('If specified, entries in the MSAs composed only of gaps ("-") will be removed. This is equivalent to '
                         'specify "--remove_fragmentary_entries --fragmentary_threshold 1"'))
    p.add_argument('--mutation_rates', action='store_true', default=False,
                   help=("If specified will produced a mutation rates table for each of the aligned markers and a summary table "
                         "for the concatenated MSA. This operation can take a long time to finish"))
    p.add_argument('--force_nucleotides', action='store_true', default=False,
                   help=("If specified force PhyloPhlAn2 to use nucleotide sequences for the phylogenetic analysis,  "
                         "even in the case of a database of amino acids"))

    group = p.add_argument_group(title="Folder paths", description="Parameters for setting the folder locations")
    group.add_argument('--input_folder', type=str, default=INPUT_FOLDER, help="Path to the folder containing the input data")
    group.add_argument('--data_folder', type=str, default=None,
                       help=('Path to the folder where to store the intermediate files, '
                             'default is "tmp" inside the project\'s output folder'))
    group.add_argument('--databases_folder', type=str, default=DATABASES_FOLDER, help="Path to the folder containing the database files")
    group.add_argument('--submat_folder', type=str, default=SUBMAT_FOLDER,
                       help=("Path to the folder containing the substitution matrices to "
                             "use to compute the column score for the subsampling step"))
    group.add_argument('--submod_folder', type=str, default=SUBMOD_FOLDER,
                       help=("Path to the folder containing the mapping file with substitution models for each marker for the "
                             "gene tree building"))
    group.add_argument('--configs_folder', type=str, default=CONFIGS_FOLDER, help="Path to the folder containing the configuration files")
    group.add_argument('--output_folder', type=str, default=OUTPUT_FOLDER, help="Path to the output folder where to save the results")

    group = p.add_argument_group(title="Filename extensions",
                                 description="Parameters for setting the extensions of the input files")
    group.add_argument('--genome_extension', type=str, default=GENOME_EXTENSION, help="Extension for input genomes")
    group.add_argument('--proteome_extension', type=str, default=PROTEOME_EXTENSION, help="Extension for input proteomes")

    p.add_argument('--verbose', action='store_true', default=False, help="Makes PhyloPhlAn2 verbose")
    p.add_argument('-v', '--version', action='version', version='PhyloPhlAn2 version {} ({})'.format(__version__, __date__),
                   help="Prints the current PhyloPhlAn2 version and exit")

    return p.parse_args()


def read_configs(config_file, verbose=False):
    configs = {}
    config = cp.ConfigParser()
    config.read(config_file)

    if verbose:
        info('Loading configuration file "{}"\n'.format(config_file))

    for section in config.sections():  # "DEFAULT" section not included!
        configs[section.lower()] = {}

        for option in config[section]:
            configs[section.lower()][option.lower()] = config[section][option]

    return configs


def check_args(args, command_line_arguments, verbose=True):
    if args.database:
        if os.path.isdir(args.database):
            if '--databases_folder' not in command_line_arguments:
                args.databases_folder = os.path.dirname(os.path.abspath(args.database))
                args.database = os.path.basename(os.path.abspath(args.database))

                if verbose:
                    info('Automatically setting "database={}" and "databases_folder={}"\n'
                         .format(args.database, args.databases_folder))
            elif not glob.glob(os.path.join(args.databases_folder, args.database) + '*'):
                error('no database found at "{}"'.format(os.path.join(args.databases_folder, args.database)), exit=True)

    args.databases_folder = check_and_create_folder(os.path.join(args.databases_folder),
                                                    try_local=True, create=True, exit=True, verbose=verbose)
    args.submat_folder = check_and_create_folder(os.path.join(args.submat_folder),
                                                 try_local=True, create=True, exit=True, verbose=verbose)
    args.submod_folder = check_and_create_folder(os.path.join(args.submod_folder),
                                                 try_local=True, create=True, exit=True, verbose=verbose)

    if args.clean_all:
        return None
    elif args.database_list:
        database_list(args.databases_folder, exit=True)
    elif args.submat_list:
        submat_list(args.submat_folder, exit=True)
    elif args.submod_list:
        submod_list(args.submod_folder, exit=True)
    elif (not args.input) and (not args.clean):
        error('either -i/--input or -c/--clean must be specified', exit=True)
    elif not args.database:
        error('-d/--database must be specified')
        database_list(args.databases_folder, exit=True)

    input_folder_setted = False

    if args.input:
        if os.path.isdir(args.input):
            if '--input_folder' not in command_line_arguments:
                input_folder_setted = True
                args.input_folder = os.path.dirname(os.path.abspath(args.input))
                args.input = os.path.basename(os.path.abspath(args.input))

                if verbose:
                    info('Automatically setting "input={}" and "input_folder={}"\n'
                         .format(args.input, args.input_folder))
            else:
                error('--input_folder is specified and the -i/--input is a path to a folder', exit=True)

        project_name = args.input
    elif args.clean:
        if os.path.isdir(args.clean):
            if '--input_folder' not in command_line_arguments:
                input_folder_setted = True
                args.input_folder = os.path.dirname(os.path.abspath(args.clean))
                args.clean = os.path.basename(os.path.abspath(args.clean))

                if verbose:
                    info('Automatically setting "clean={}" and "input_folder={}"\n'
                         .format(args.input, args.input_folder))
            else:
                error('--input_folder is specified and the -c/--clean is a path to a folder', exit=True)

        project_name = args.clean

    if input_folder_setted:
        args.input_folder = os.path.join(args.input_folder, project_name)

    if not args.output:
        args.output = os.path.join(args.output_folder, project_name + '_' + args.database)
    elif args.output_folder:
        args.output = os.path.join(args.output_folder, args.output)

    if not args.data_folder:
        args.data_folder = os.path.join(args.output, 'tmp')

    if args.clean:
        check_and_create_folder(args.data_folder, exit=True, verbose=verbose)
        check_and_create_folder(args.output, exit=True, verbose=verbose)
        return None

    check_and_create_folder(args.input_folder, exit=True, verbose=verbose)
    args.configs_folder = check_and_create_folder(args.configs_folder, try_local=True, exit=True, verbose=verbose)
    check_and_create_folder(args.output, create=True, exit=True, verbose=verbose)
    check_and_create_folder(args.data_folder, create=True, exit=True, verbose=verbose)

    if not args.genome_extension.startswith('.'):
        args.genome_extension = '.' + args.genome_extension

    if args.genome_extension.endswith('.'):
        args.genome_extension = args.genome_extension[:-1]

    if not args.proteome_extension.startswith('.'):
        args.proteome_extension = '.' + args.proteome_extension

    if args.proteome_extension.endswith('.'):
        args.proteome_extension = args.proteome_extension[:-1]

    # if both --accurate and --fast are not specified, then go with accurate
    if (not args.accurate) and (not args.fast):
        args.accurate = True

    # set pre-config params
    trim = None
    submat = 'pfasum60'
    subsample = None
    scoring_function = None
    fragmentary_threshold = None
    gap_perc_threshold = None
    not_variant_threshold = None
    remove_fragmentary_entries = True

    if args.diversity == 'low':
        # accurate
        trim = 'not_variant'
        not_variant_threshold = 0.99

        if args.fast:
            trim = 'greedy'
            subsample = 'fivehundred'
            scoring_function = 'trident'
            gap_perc_threshold = 0.67
            fragmentary_threshold = 0.85

    if args.diversity == 'medium':
        scoring_function = 'trident'

        if args.accurate:
            trim = 'gap_trim'
            subsample = 'onehundred'
            fragmentary_threshold = 0.85

        if args.fast:
            trim = 'greedy'
            subsample = 'fifty'
            fragmentary_threshold = 0.75
            gap_perc_threshold = 0.75
            not_variant_threshold = 0.97

    if args.diversity == 'high':
        trim = 'greedy'
        scoring_function = 'trident'
        gap_perc_threshold = 0.85

        if args.accurate:
            subsample = 'twentyfive'
            fragmentary_threshold = 0.75
            not_variant_threshold = 0.95

        if args.fast:
            subsample = 'phylophlan' if args.database == 'phylophlan' else 'tenpercent'
            fragmentary_threshold = 0.67
            not_variant_threshold = 0.9

    if verbose:
        info('"{}" preset\n'.format('{}-{}'.format(args.diversity, 'accurate' if args.accurate else 'fast' if args.fast else 'none')))

    # check substitution model
    if args.maas:
        if not os.path.isfile(args.maas):
            if os.path.isfile(os.path.join(args.submod_folder, args.maas)):
                args.maas = os.path.join(args.submod_folder, args.maas)
            else:
                error('file "{}" not found in "{}"'.format(args.maas, args.submod_folder))
                submod_list(args.submod_folder, exit=True)

    # check subsample settings with pre-config ones
    if args.subsample:
        if '--subsample' in command_line_arguments:
            if (args.database != 'phylophlan') and (args.subsample == 'phylophlan'):
                error('scoring function "phylophlan" is compatible only with "phylophlan" database', exit=True)
    elif subsample:
        args.subsample = subsample

    # check database settings
    if (args.database == 'phylophlan') and (not args.sort):
        args.sort = True

        if verbose:
            info('Setting "sort={}" because "database={}"\n'.format(args.sort, args.database))

    # check min_num_markers settings
    if not args.min_num_markers:
        if args.database == 'phylophlan':
            args.min_num_markers = 100

            if verbose:
                info('Setting "min_num_markers={}" since no value has been specified and the "database={}"\n'
                     .format(args.min_num_markers, args.database))
        elif args.database == 'amphora2':
            args.min_num_markers = 34

            if verbose:
                info('Setting "min_num_markers={}" since no value has been specified and the "database={}"\n'
                     .format(args.min_num_markers, args.database))
        else:
            args.min_num_markers = MIN_NUM_MARKERS

    # check not_variant_threshold settings
    if not_variant_threshold and ('--not_variant_threshold' not in command_line_arguments):
        args.not_variant_threshold = not_variant_threshold

    # check fragmentary_threshold settings
    if fragmentary_threshold and ('--fragmentary_threshold' not in command_line_arguments):
        args.fragmentary_threshold = fragmentary_threshold

    # check remove_only_gaps_entries settings
    if args.remove_only_gaps_entries:
        args.fragmentary_threshold = 1.0

        if args.remove_fragmentary_entries:
            args.fragmentary_threshold = min(FRAGMENTARY_THRESHOLD, args.fragmentary_threshold)
            info('[w] both "--remove_only_gaps_entries" and '
                 '"--remove_fragmentary_entries" have been specified, setting '
                 '"fragmentary_threshold={}"\n'.format(args.fragmentary_threshold))

    # check trim settings
    if trim and ('--trim' not in command_line_arguments):
        args.trim = trim

    # check submat settings
    if submat and (('-s' not in command_line_arguments) or ('--submat' not in command_line_arguments)):
        args.submat = submat

    # checking configuration file
    if not args.config_file:
        error('-f/--config_file must be specified')
        config_list(args.configs_folder, exit=True, exit_value=1)

    if os.path.isfile(args.config_file):
        pass
    elif os.path.isfile(os.path.join(args.configs_folder, args.config_file)):
        args.config_file = os.path.join(args.configs_folder, args.config_file)
    else:
        error('configuration file "{}" not found'.format(args.config_file))
        config_list(args.configs_folder, exit=True, exit_value=1)

    # checking substitution matrix
    if args.subsample and (not args.submat):
        error('-s/--submat must be specified')
        submat_list(args.submat_folder, exit=True, exit_value=1)
    elif (args.submat and (not os.path.isfile(os.path.join(args.submat_folder, args.submat + '.pkl')))):
        error('substitution matrix "{}" not found in "{}"'.format(args.submat, args.submat_folder))
        submat_list(args.submat_folder, exit=True, exit_value=1)

    # get scoring function
    score_function = (args.scoring_function
                      if args.scoring_function and ('--scoring_function' in command_line_arguments)
                      else scoring_function)

    if score_function:
        if score_function in globals():
            args.scoring_function = globals().get(score_function)
        elif (not score_function) and (score_function in locals()):
            args.scoring_function = locals().get(score_function)
        else:
            error('cannot find scoring function "{}"'.format(score_function), exit=True)

    # get position function
    if args.subsample:
        npos_function = None

        if args.subsample in globals():
            npos_function = globals().get(args.subsample)
        elif (not npos_function) and (args.subsample in locals()):
            npos_function = locals().get(args.subsample)
        else:
            error('cannot find subsampling function "{}"'.format(args.subsample), exit=True)

        args.subsample = npos_function

    if verbose:
        info('Arguments: {}\n'.format(vars(args)))

    return project_name


def check_configs(configs, verbose=False):
    if verbose:
        info('Checking configuration file\n')

    # checking whether mandatory sections and options are present
    for sections in CONFIG_SECTIONS_MANDATORY:
        mandatory_sections = False

        for section in sections:
            if section in configs:
                mandatory_sections = True

                for options in CONFIG_OPTIONS_MANDATORY:
                    mandatory_options = False

                    for option in options:
                        if (option in configs[section]) and configs[section][option]:
                            mandatory_options = True
                            break

                    if not mandatory_options:
                        error('could not find "{}" mandatory option in section "{}"'.format(option, section), exit=True)

        if not mandatory_sections:
            error('could not find "{}" section'.format(section), exit=True)

    # all the options (but 'version') must be defined in 'command_line'
    for section, options in configs.items():
        actual_options = []
        mandatory_options = []

        for option in options:
            if option in ['command_line']:
                actual_options = [a.strip() for a in configs[section][option].split('#') if a.strip()]
            elif option not in CONFIG_OPTIONS_TO_EXCLUDE:
                mandatory_options.append(option)

        if mandatory_options and actual_options:
            for option in mandatory_options:
                if option not in actual_options:
                    error('option "{}" not defined in section "{}" in your configuration file'.format(option, section),
                          exit=True)
        else:
            error('wrongly formatted configuration file?', exit=True)


def check_and_create_folder(folder, try_local=False, create=False, exit=False, verbose=False):
    folders_to_test = [folder]
    msg_err = ''

    if try_local:
        folders_to_test.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), folder))

    for f in folders_to_test:
        if not os.path.isdir(f):
            if not create:
                msg_err = '"{}" folder does not exists'.format(f)
        else:
            return f

    if msg_err:
        error(msg_err, exit=exit)
        return None

    if create:
        if not os.path.isdir(f):
            os.mkdir(folder, mode=0o775)

            if verbose:
                info('Creating folder "{}"\n'.format(folder))
        elif verbose:
            info('Folder "{}" exists\n'.format(folder))

        return folder


def check_dependencies(configs, nproc, verbose=False):
    cmds_tested = []

    for params in configs:
        cmd = compose_command(configs[params], check=True, nproc=nproc)

        if cmd not in cmds_tested:
            cmds_tested.append(cmd)

            if verbose:
                info('Checking "{}"\n'.format(' '.join(cmd['command_line'])))

            inp_f = None
            out_f = sb.DEVNULL

            if cmd['stdin']:
                inp_f = open(cmd['stdin'], 'r')

            if cmd['stdout']:
                out_f = open(cmd['stdout'], 'w')

            try:
                if 'diamond' in cmd['command_line']:
                    stdout = sb.check_output(cmd['command_line'], env=cmd['env'])
                    wrong = False

                    # diamond before version 0.8.30 is not able to map with the blastx algorithm as they are missing
                    # the --max-hsps param. According to the changelog:
                    # 0.8.31 (1 Jan 2017) - removed --single-domain option and replaced by --max-hsps
                    # phylophlan-users@googlegroups.com, subject: "phylophlan.fna is empty"
                    for a, b in zip(stdout.decode().strip().split(' ')[-1].split('.'), [0, 8, 31]):
                        if a != str(b):
                            wrong = wrong or (int(a) < b)
                            break

                    if wrong:
                        error('the "{}" of diamond is not supported, please update it to at least version 0.8.31'.format(stdout),
                              init_new_line=True, exit=True)
                else:
                    sb.check_call(cmd['command_line'], stdin=inp_f, stdout=out_f, stderr=sb.DEVNULL, env=cmd['env'])
            except Exception as e:
                error(str(e), init_new_line=True)
                error('program not installed or not present in the system path\n'
                      '    {}'.format('\n    '.join(['{:>12}: {}'.format(a, ' '.join(cmd[a]) if type(cmd[a]) is list else cmd[a])
                                                     for a in ['command_line', 'stdin', 'stdout', 'env']])),
                      init_new_line=True, exit=True)

            if cmd['stdin']:
                inp_f.close()

            if cmd['stdout']:
                out_f.close()


def check_database(db_name, databases_folder, verbose=False):
    is_dir = os.path.isdir(os.path.join(databases_folder, db_name))
    is_faa = os.path.isfile(os.path.join(databases_folder, db_name + '.faa'))
    is_faa_bz2 = os.path.isfile(os.path.join(databases_folder, db_name + '.faa.bz2'))
    is_tar = (os.path.isfile(os.path.join(databases_folder, db_name + '.tar')) and
              os.path.isfile(os.path.join(databases_folder, db_name + '.md5')))

    if not (is_dir or is_faa or is_faa_bz2 or is_tar):
        error('database "{}" not found in "{}"'.format(db_name, databases_folder))
        database_list(databases_folder, exit=True)


def database_list(databases_folder, exit=False, exit_value=0):
    if not os.path.isdir(databases_folder):
        error('folder "{}" does not exists'.format(databases_folder), exit=True, exit_value=exit_value)

    info('Available databases in "{}":\n    '.format(databases_folder))
    info('\n    '.join([a for a in sorted(os.listdir(databases_folder)) if os.path.isdir(os.path.join(databases_folder, a))]) + '\n',
         exit=exit, exit_value=exit_value)


def submat_list(submat_folder, exit=False, exit_value=0):
    if not os.path.isdir(submat_folder):
        error('folder "{}" does not exists'.format(submat_folder), exit=True, exit_value=exit_value)

    info('Available substitution matrices in "{}":\n    '.format(submat_folder))
    info('\n    '.join(sorted([os.path.splitext(os.path.basename(a))[0]
                               for a in glob.iglob(os.path.join(submat_folder, '*.pkl'))])) + '\n', exit=exit, exit_value=exit_value)


def submod_list(submod_folder, exit=False, exit_value=0):
    if not os.path.isdir(submod_folder):
        error('folder "{}" does not exists'.format(submod_folder), exit=True, exit_value=exit_value)

    info('Available substitution models in "{}":\n    '.format(submod_folder))
    info('\n    '.join(sorted([os.path.basename(a) for a in glob.iglob(os.path.join(submod_folder, '*.tsv'))])) + '\n',
         exit=exit, exit_value=exit_value)


def config_list(config_folder, exit=False, exit_value=0):
    if not os.path.isdir(config_folder):
        error('folder "{}" does not exists'.format(config_folder), exit=True, exit_value=exit_value)

    info('Available configuration files in "{}":\n    '.format(config_folder))
    info('\n    '.join(sorted(glob.iglob(os.path.join(config_folder, '*.cfg')))) + '\n', exit=exit, exit_value=exit_value)


def compose_command(params, check=False, sub_mod=None, input_file=None, database=None, output_path=None, output_file=None, nproc=1):
    program_name = None
    stdin = None
    stdout = None
    environment = os.environ.copy()
    r_output_path = None
    r_output_file = None
    command_line = params['command_line']

    if 'program_name' in params.keys():
        command_line = command_line.replace('#program_name#', params['program_name'])
        program_name = params['program_name']
    else:
        error('something wrong... "program_name" not found!', exit=True)

    if check:
        command_line = program_name

        if 'version' in params:
            command_line = '{} {}'.format(program_name, params['version'])
    else:
        if 'params' in params:
            command_line = command_line.replace('#params#', params['params'])

        if 'threads' in params:
            command_line = command_line.replace('#threads#', '{} {}'.format(params['threads'], nproc))

        if output_path:
            r_output_path = output_path

            if 'output_path' in params:
                command_line = command_line.replace('#output_path#', '{} {}'.format(params['output_path'], output_path))
            else:
                output_file = os.path.join(output_path, output_file)

        if sub_mod:
            mod = sub_mod

            if 'model' in params:
                mod = '{} {}'.format(params['model'], sub_mod)

            command_line = command_line.replace('#model#', mod)

        if input_file:
            inp = input_file

            if 'input' in params:
                inp = '{} {}'.format(params['input'], input_file)

            if '<' in command_line:
                command_line = command_line.replace('<', '')
                command_line = command_line.replace('#input#', '')
                stdin = inp
            else:
                command_line = command_line.replace('#input#', inp)

        if database and ('database' in params):
            command_line = command_line.replace('#database#', '{} {}'.format(params['database'], database))

        if output_file:
            out = output_file
            r_output_file = output_file

            if 'output' in params:
                out = '{} {}'.format(params['output'], output_file)

            if '>' in command_line:
                command_line = command_line.replace('>', '')
                command_line = command_line.replace('#output#', '')
                stdout = out
            else:
                command_line = command_line.replace('#output#', out)

        if 'environment' in params:
            new_environment = dict([(var.strip(), val.strip())
                                    for var, val in [a.strip().split('=') for a in params['environment'].split(',')]])
            environment.update(new_environment)

    # find string sourrunded with " and make them as one string
    quotes = [j for j, e in enumerate(command_line) if e == '"']

    for s, e in zip(quotes[0::2], quotes[1::2]):
        command_line = command_line.replace(command_line[s + 1:e], command_line[s + 1:e])

    return {'command_line': [str(a) for a in re.sub(' +', ' ', command_line.replace('"', '')).split(' ') if a],
            'stdin': stdin, 'stdout': stdout, 'env': environment, 'output_path': r_output_path, 'output_file': r_output_file}


def remove_file(filename, path=None, verbose=False):
    to_remove = None

    if os.path.isfile(filename):
        to_remove = filename
    elif os.path.isfile(os.path.join(path, filename)):
        to_remove = os.path.join(path, filename)

    if to_remove:
        if os.path.isfile(to_remove):
            if verbose:
                info('Removing "{}"\n'.format(to_remove))

            os.remove(to_remove)
        elif verbose:
            error('cannot remove "{}", file not found'.format(to_remove))


def remove_files(file_list, path=None, verbose=False):
    for f in file_list:
        remove_file(f, path=path, verbose=verbose)


def init_database(database, databases_folder, db_type, params, key_dna, key_aa, verbose=False):
    db_dna, db_aa = None, None

    if not db_type:
        fna = os.path.join(databases_folder, database, database + '.fna')
        fna_bz2 = os.path.join(databases_folder, database, database + '.fna.bz2')
        faa = os.path.join(databases_folder, database, database + '.faa')
        faa_bz2 = os.path.join(databases_folder, database, database + '.faa.bz2')
        folder = os.path.join(databases_folder, database)
        d = None

        if os.path.isfile(fna) or os.path.isfile(fna_bz2):
            d = Counter([len(set(seq))
                         for _, seq in SimpleFastaParser(open(fna) if os.path.isfile(fna) else bz2.open(fna_bz2, 'rt'))])
        elif os.path.isfile(faa) or os.path.isfile(faa_bz2):
            d = Counter([len(set(seq))
                         for _, seq in SimpleFastaParser(open(faa) if os.path.isfile(faa) else bz2.open(faa_bz2, 'rt'))])
        elif os.path.isdir(folder):
            d = Counter([len(set(seq))
                         for f in glob.iglob(os.path.join(folder, '*'))
                         for _, seq in SimpleFastaParser(bz2.open(f, 'rt') if f.endswith('.bz2') else open(f))])
        else:
            error("-t/--db_type not specified and could not automatically detect the input database file(s)", exit=True)

        # nucleotides --> A, T, G, C, S (G or C), N, -
        if len([f for f, _ in d.items() if f > 7]) > (len(d) / 2):  # if >50% are >7
            db_type = 'a'  # amino acids
        elif len([f for f, _ in d.items() if f < 8]) > (len(d) / 2):  # if >50% are <=7
            db_type = 'n'  # nucleotides
        else:
            error("-t/--db_type not specified and I could not automatically infer the db type", exit=True)

    if db_type == 'a':
        db_dna, db_aa = init_database_aa(database, databases_folder, params, key_dna, key_aa, verbose=verbose)
    elif db_type == 'n':
        db_dna, db_aa = init_database_nt(database, databases_folder, params, key_dna, key_aa, verbose=verbose)
    else:
        error('\n    '.join(['{}: {}'.format(lbl, val) for lbl, val in
                             zip(['init_database', 'database', 'databases_folder', 'db_type', 'params', 'key_dna', 'key_aa', 'verbose'],
                                 [init_database, database, databases_folder, db_type, params, key_dna, key_aa, verbose])]), exit=True)

    if (not db_dna) and (not db_aa):
        error('both db_dna and db_aa are None!', exit=True)

    return (db_type, db_dna, db_aa)


def init_database_aa(database, databases_folder, params, key_dna, key_aa, verbose=False):
    db_fasta = os.path.join(databases_folder, database, database + '.faa')
    db_folder = os.path.join(databases_folder, database)
    db_dna = None
    db_aa = None

    # assumed to be a fasta file containing the markers
    if os.path.isfile(db_fasta):
        markers = None
    elif os.path.isfile(os.path.join(db_folder, database + '.faa.bz2')):
        markers = [os.path.join(db_folder, database + '.faa.bz2')]
    elif os.path.isdir(db_folder):  # assumed to be a folder with a fasta file for each marker
        markers = glob.iglob(os.path.join(db_folder, '*.faa*'))
    else:  # what's that??
        error('database format ("{}", "{}", or "{}") not recognize'
              .format(db_fasta, os.path.join(db_folder, database + '.faa.bz2'), db_folder), exit=True)

    if key_dna in params:
        if 'diamond' in params[key_dna]['program_name']:
            db_dna = database + '.dmnd'

            if not os.path.isfile(os.path.join(db_folder, db_dna)):
                make_database(params[key_dna], db_fasta, markers, db_folder, database, key_dna, verbose=verbose)

                if not os.path.isfile(os.path.join(db_folder, db_dna)):
                    error('database "{}" has not been created... something went wrong!'
                          .format(os.path.join(db_folder, db_dna)), exit=True)
            elif verbose:
                info('"{}" database "{}" present\n'.format(key_dna, db_dna))
            else:
                db_dna = None

            db_dna = os.path.join(db_folder, db_dna)
    else:
        db_dna = None

    if key_aa in params:
        if 'usearch' in params[key_aa]['program_name']:
            db_aa = os.path.join(db_folder, database + '.udb')

            if not os.path.isfile(db_aa):
                make_database(params[key_aa], db_fasta, markers, db_folder, database + '.udb', key_aa, verbose=verbose)

                if not os.path.isfile(db_aa):
                    error('database "{}" has not been created... something went wrong!'.format(db_aa), exit=True)
            elif verbose:
                info('"{}" database "{}" present\n'.format(key_aa, db_aa))
            else:
                db_aa = None
        elif 'diamond' in params[key_aa]['program_name']:
            db_aa = os.path.join(db_folder, database + '.dmnd')

            if not os.path.isfile(db_aa):
                make_database(params[key_aa], db_fasta, markers, db_folder, database, key_aa, verbose=verbose)

                if not os.path.isfile(db_aa):
                    error('database "{}" has not been created... something went wrong!'.format(db_aa), exit=True)
            elif verbose:
                info('"{}" database "{}" present\n'.format(key_aa, db_aa))

            db_dna = db_aa  # if there are genomes the database is the same
        else:
            error('program "{}" not recognize'
                  .format(params[key_aa]['program_name'] if 'program_name' in params[key_aa] else 'None'), exit=True)
    else:
        db_aa = None

    return (db_dna, db_aa)


def init_database_nt(database, databases_folder, params, key_dna, key_aa, verbose=False):
    db_fasta = os.path.join(databases_folder, database, database + '.fna')
    db_folder = os.path.join(databases_folder, database)
    db_aa = None
    db_dna = None
    makeblastdb_exts = ['.nhr', '.nin', '.nog', '.nsd', '.nsi', '.nsq']

    # assumed to be a fasta file containing the markers
    if os.path.isfile(db_fasta):
        markers = None
    elif os.path.isfile(os.path.join(db_folder, database + '.fna.bz2')):  # is it a compressed fasta file?
        markers = [os.path.join(db_folder, database + '.fna.bz2')]
    elif os.path.isdir(db_folder):  # assumed to be a folder with a fasta file for each marker
        markers = glob.iglob(os.path.join(db_folder, '*.fna*'))
    else:  # what's that??
        error('database format ("{}", "{}", or "{}") not recognize'
              .format(db_fasta, os.path.join(db_folder, database + '.fna.bz2'), db_folder), exit=True)

    if key_dna in params:
        if 'makeblastdb' in params[key_dna]['program_name']:
            db_dna = os.path.join(databases_folder, database, database)

            if [None for ext in makeblastdb_exts if not os.path.isfile(os.path.join(db_folder, database + ext))]:
                make_database(params[key_dna], db_fasta, markers, db_folder, database, key_dna,
                              output_exts=makeblastdb_exts, verbose=verbose)

                if [None for ext in makeblastdb_exts if not os.path.isfile(os.path.join(db_folder, database + ext))]:
                    error('database "{}" ({}) has not been created... something went wrong!'
                          .format(key_dna, os.path.join(databases_folder, database), ', '.join(makeblastdb_exts)), exit=True)
            elif verbose:
                info('"{}" database "{}" ({}) present\n'
                     .format(key_dna, os.path.join(databases_folder, database), ', '.join(makeblastdb_exts)))
        else:
            error('program "{}" not recognize'
                  .format(params[key_dna]['program_name'] if 'program_name' in params[key_dna] else 'None'), exit=True)
    else:
        db_dna = None
        info('[w] cannot create database "{}", section "{}" not present in configurations\n'.format(database, key_dna))

    return (db_dna, db_aa)


def make_database(command, fasta, markers, db_folder, db, label, output_exts=[''], verbose=False):
    if fasta and (not os.path.isfile(fasta)):
        with open(fasta, 'w') as f:
            for i in markers:
                g = bz2.open(i, 'rt') if i.endswith('.bz2') else open(i)
                f.write(g.read())
                g.close()
    elif verbose:
        info('File "{}" present\n'.format(fasta))

    info('Generating "{}" indexed database "{}"\n'.format(label, db))
    cmd = compose_command(command, input_file=fasta, output_path=db_folder, output_file=db)
    inp_f = None
    out_f = sb.DEVNULL

    if cmd['stdin']:
        inp_f = open(cmd['stdin'], 'r')

    if cmd['stdout']:
        out_f = open(cmd['stdout'], 'w')

    try:
        sb.check_call(cmd['command_line'], stdin=inp_f, stdout=out_f, stderr=sb.DEVNULL, env=cmd['env'])
    except Exception as e:
        remove_files([os.path.join(db_folder, db + ext) for ext in output_exts], path=cmd['output_path'], verbose=verbose)
        error(str(e), init_new_line=True)
        error('cannot execute command\n    {}'
              .format('\n    '.join(['{:>12}: {}'.format(a, ' '.join(cmd[a]) if type(cmd[a]) is list else cmd[a])
                                     for a in ['command_line', 'stdin', 'stdout', 'env']])), init_new_line=True, exit=True)

    if cmd['stdin']:
        inp_f.close()

    if cmd['stdout']:
        out_f.close()

    info('"{}" ("{}") generated\n'.format(label, db))


def clean_all(databases_folder, verbose=False):
    for f in glob.glob(os.path.join(databases_folder, '*.udb')) + glob.glob(os.path.join(databases_folder, '*.dmnd')):
        if verbose:
            info('Removing "{}"\n'.format(f))

        os.remove(f)
        f_clean, _ = os.path.splitext(f)

        if (os.path.isfile(f_clean + '.faa') and os.path.isfile(f_clean + '.faa.bz2')):
            if verbose:
                info('Removing "{}"\n'.format(f_clean + '.faa'))

            os.remove(f_clean + '.faa')

    for database in os.listdir(databases_folder):
        for f in (glob.glob(os.path.join(databases_folder, database, database + '.faa')) +
                  glob.glob(os.path.join(databases_folder, database, database + '.udb')) +
                  glob.glob(os.path.join(databases_folder, database, database + '.dmnd'))):
            if verbose:
                info('Removing "{}"\n'.format(f))

            os.remove(f)

    sys.exit(0)


def clean_project(data_folder, output_folder, verbose=False):
    if os.path.exists(data_folder):
        if verbose:
            info('Removing folder "{}"\n'.format(data_folder))

        shutil.rmtree(data_folder)

    if os.path.exists(output_folder):
        if verbose:
            info('Removing folder "{}"\n'.format(output_folder))

        shutil.rmtree(output_folder)

    info('Folders "{}" and "{}" removed\n'.format(data_folder, output_folder))
    sys.exit(0)


def load_input_files(input_folder, tmp_folder, extension, verbose=False):
    inputs = {}
    done = False

    if os.path.isdir(input_folder):
        info('Loading files from "{}"\n'.format(input_folder))
        files = glob.iglob(os.path.join(input_folder, '*' + extension + '*'))

        for f in files:
            if f.endswith('.bz2'):
                if verbose and (not done):
                    info('Decompressing input files\n')
                    done = True

                check_and_create_folder(tmp_folder, create=True, exit=True, verbose=verbose)
                file_clean = os.path.splitext(os.path.basename(f))[0]

                if not os.path.isfile(os.path.join(tmp_folder, file_clean)):
                    with open(os.path.join(tmp_folder, file_clean), 'w') as g, bz2.open(f, 'rt') as h:
                        SeqIO.write((SeqRecord(Seq(seq), id=header.split(' ')[0]) for header, seq in SimpleFastaParser(h)), g, "fasta")
                elif verbose:
                    info('File "{}" already decompressed\n'.format(os.path.join(tmp_folder, file_clean)))

                inputs[file_clean] = tmp_folder
            elif f.endswith('.gz'):
                if verbose and (not done):
                    info('Decompressing input files\n')
                    done = True

                check_and_create_folder(tmp_folder, create=True, exit=True, verbose=verbose)
                file_clean = os.path.splitext(os.path.basename(f))[0]

                if not os.path.isfile(os.path.join(tmp_folder, file_clean)):
                    with open(os.path.join(tmp_folder, file_clean), 'w') as g, gzip.open(f, 'rt') as h:
                        SeqIO.write((SeqRecord(Seq(seq), id=header.split(' ')[0]) for header, seq in SimpleFastaParser(h)), g, "fasta")
                elif verbose:
                    info('File "{}" already decompressed\n'.format(os.path.join(tmp_folder, file_clean)))

                inputs[file_clean] = tmp_folder
            elif f.endswith(extension):
                inputs[os.path.basename(f)] = input_folder
            elif verbose:
                info('Input file "{}" not recognized\n'.format(f))
    elif verbose:
        info('Folder "{}" does not exists\n'.format(input_folder))

    return inputs


def initt(terminating_):
    # This places terminating in the global namespace of the worker subprocesses.
    # This allows the worker function to access `terminating` even though it is
    # not passed as an argument to the function.
    global terminating
    terminating = terminating_


def check_input_proteomes(inputs, min_num_proteins, min_len_protein, data_folder, nproc=1, verbose=False):
    good_inputs = []

    if os.path.isfile(os.path.join(data_folder, 'checked_inputs.pkl')):
        info('Inputs already checked\n')

        if verbose:
            info('Loading checked inputs from "{}"\n'.format(os.path.join(data_folder, 'checked_inputs.pkl')))

        with open(os.path.join(data_folder, 'checked_inputs.pkl'), 'rb') as f:
            good_inputs = pickle.load(f)
    else:
        info('Checking {} inputs\n'.format(len(inputs)))
        terminating = mp.Event()

        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                good_inputs = [a for a in
                               pool.imap_unordered(check_input_proteomes_rec,
                                                   ((os.path.join(inp_fol, inp), min_len_protein, min_num_proteins, verbose)
                                                    for inp, inp_fol in inputs.items()), chunksize=1) if a]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('check_input_proteomes crashed', init_new_line=True, exit=True)

        with open(os.path.join(data_folder, 'checked_inputs.pkl'), 'wb') as f:
            pickle.dump(good_inputs, f, protocol=pickle.HIGHEST_PROTOCOL)

    return good_inputs


def check_input_proteomes_rec(x):
    if not terminating.is_set():
        try:
            inp, min_len_protein, min_num_proteins, verbose = x
            info('Checking "{}"\n'.format(inp))
            num_proteins = len([1 for _, seq in SimpleFastaParser(open(inp)) if len(seq) >= min_len_protein])

            if num_proteins >= min_num_proteins:
                return inp
            elif verbose:
                info('"{}" discarded, not enough proteins ({}/{}) of at least {} AAs\n'
                     .format(inp, num_proteins, min_num_proteins, min_len_protein))

            return None
        except Exception as e:
            terminating.set()
            error(str(e), init_new_line=True)
            error('error while checking\n    {}'.format('\n    '.join([str(a) for a in x])), init_new_line=True)
            raise
    else:
        terminating.set()


def clean_input_proteomes(inputs, output_folder, nproc=1, verbose=False):
    commands = []
    check_and_create_folder(output_folder, create=True, verbose=verbose)
    commands = [(inp, os.path.join(output_folder, os.path.basename(inp)))
                for inp in inputs if not os.path.isfile(os.path.join(output_folder, os.path.basename(inp)))]

    if commands:
        info('Cleaning {} inputs\n'.format(len(commands)))
        terminating = mp.Event()

        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                [_ for _ in pool.imap_unordered(clean_input_proteomes_rec, commands, chunksize=1)]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('clean_input_proteomes crashed', init_new_line=True, exit=True)
    else:
        info('Inputs already cleaned\n')


def clean_input_proteomes_rec(x):
    if not terminating.is_set():
        try:
            t0 = time.time()
            inp, out = x
            inp_clean, _ = os.path.splitext(os.path.basename(inp))
            info('Cleaning "{}"\n'.format(inp))

            # http://biopython.org/DIST/docs/api/Bio.Alphabet.IUPAC.ExtendedIUPACProtein-class.html
            # B = "Asx"; Aspartic acid (R) or Asparagine (N)
            # X = "Xxx"; Unknown or 'other' amino acid
            # Z = "Glx"; Glutamic acid (E) or Glutamine (Q)
            # J = "Xle"; Leucine (L) or Isoleucine (I), used in mass-spec (NMR)
            # U = "Sec"; Selenocysteine
            # O = "Pyl"; Pyrrolysine
            output = (SeqRecord(Seq(e[1].replace('B', 'X').replace('Z', 'X').replace('J', 'X').replace('U', 'X').replace('O', 'X')),
                                id='{}_{}'.format(inp_clean, counter), description='')
                      for counter, e in enumerate(SimpleFastaParser(open(inp))) if e[1])

            with open(out, 'w') as f:
                SeqIO.write(output, f, "fasta")

            t1 = time.time()
            info('"{}" generated in {}s\n'.format(out, int(t1 - t0)))
        except Exception as e:
            terminating.set()
            remove_file(out)
            error(str(e), init_new_line=True)
            error('error while cleaning\n    {}'.format('\n    '.join([str(a) for a in x])), init_new_line=True)
            raise
    else:
        terminating.set()


def gene_markers_identification(configs, key, inputs, output_folder, database_name, database, min_num_proteins,
                                nproc=1, verbose=False):
    commands = []
    check_and_create_folder(output_folder, create=True, verbose=verbose)

    for inp, inp_fol in inputs.items():
        out = os.path.splitext(inp)[0] + '.b6o.bkp'

        if not os.path.isfile(os.path.join(output_folder, out)):
            commands.append((configs[key], os.path.join(inp_fol, inp), database, output_folder, out, min_num_proteins, verbose))

    if commands:
        info('Mapping "{}" on {} inputs (key: "{}")\n'.format(database_name, len(commands), key))
        terminating = mp.Event()

        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                [_ for _ in pool.imap_unordered(gene_markers_identification_rec, commands, chunksize=1)]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('gene_markers_identification crashed', init_new_line=True, exit=True)
    else:
        info('"{}" markers already mapped (key: "{}")\n'.format(database_name, key))


def gene_markers_identification_rec(x):
    if not terminating.is_set():
        try:
            t0 = time.time()
            params, inp, db, out_fld, out, min_num_proteins, verbose = x
            info('Mapping "{}"\n'.format(inp))
            cmd = compose_command(params, input_file=inp, database=db, output_path=out_fld, output_file=out)
            inp_f = None
            out_f = sb.DEVNULL

            if cmd['stdin']:
                inp_f = open(cmd['stdin'], 'r')

            if cmd['stdout']:
                out_f = open(cmd['stdout'], 'w')

            try:
                sb.check_call(cmd['command_line'], stdin=inp_f, stdout=out_f, stderr=sb.DEVNULL, env=cmd['env'])
            except Exception as e:
                terminating.set()
                remove_file(cmd['output_file'], path=cmd['output_path'], verbose=verbose)
                error(str(e), init_new_line=True)
                error('cannot execute command\n    {}'
                      .format('\n    '.join(['{:>12}: {}'.format(a, ' '.join(cmd[a]) if type(cmd[a]) is list else cmd[a])
                                             for a in ['command_line', 'stdin', 'stdout', 'env']])), init_new_line=True)
                raise

            if cmd['stdin']:
                inp_f.close()

            if cmd['stdout']:
                out_f.close()

            t1 = time.time()
            info('"{}" generated in {}s\n'.format(out, int(t1 - t0)))
        except Exception as e:
            terminating.set()
            remove_file(out, out_fld, verbose=verbose)
            error(str(e), init_new_line=True)
            error('error while mapping\n    {}'.format('\n    '.join([str(a) for a in x])), init_new_line=True)
            raise
    else:
        terminating.set()


def gene_markers_selection(input_folder, function, min_num_proteins, nucleotides, nproc=1, verbose=False):
    commands = [(f, f[:-4], function, min_num_proteins, nucleotides)
                for f in glob.iglob(os.path.join(input_folder, '*.b6o.bkp')) if not os.path.isfile(f[:-4])]

    if commands:
        info('Selecting {} markers from "{}"\n'.format(len(commands), input_folder))
        terminating = mp.Event()

        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                [_ for _ in pool.imap_unordered(gene_markers_selection_rec, commands, chunksize=1)]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('gene_markers_selection crashed', init_new_line=True, exit=True)
    else:
        info('Markers already selected\n')


def gene_markers_selection_rec(x):
    if not terminating.is_set():
        try:
            t0 = time.time()
            inp, out, function, min_num_proteins, nt = x
            info('Selecting "{}"\n'.format(inp))
            matches = function(inp, nt)

            # there should be at least min_num_proteins mapped
            if len(matches) >= min_num_proteins:
                with open(out, 'w') as f:
                    f.write('{}\n'.format('\n'.join(['\t'.join(m) for m in matches])))

                t1 = time.time()
                info('"{}" generated in {}s\n'.format(out, int(t1 - t0)))
            else:
                info('Not enough markers mapped ({}/{}) in "{}"\n'.format(len(matches), min_num_proteins, inp))
        except Exception as e:
            terminating.set()
            remove_file(out)
            error(str(e), init_new_line=True)
            error('error while selecting\n    {}'.format('\n    '.join([str(a) for a in x])), init_new_line=True)
            raise
    else:
        terminating.set()


def best_hit(f, nucleotides):
    tab = (ll.strip().split('\t') for ll in open(f))
    best_matches = {}

    for entry in tab:
        c = entry[0].split(' ')[0]
        m = entry[1].split('_')[1]
        cs = entry[6]
        ce = entry[7]
        ms = entry[8]
        me = entry[9]
        b = entry[-1]
        cl = int(ce) - (int(cs) - 1) if int(cs) < int(ce) else int(cs) - (int(ce) - 1)
        ml = int(me) - (int(ms) - 1) if int(ms) < int(me) else int(ms) - (int(me) - 1)

        if nucleotides:
            cl = cl / 3

        ratio = cl / ml if ml > cl else ml / cl

        if ratio < 0.34:  # skip too short hits (less than 34% of the length of the matching marker)
            continue

        rev = '0' if (int(cs) < int(ce)) and (int(ms) < int(me)) else '1'

        if m in best_matches:
            if float(b) > float(best_matches[m][-1]):
                best_matches[m] = [c, m, cs, ce, rev, b]
        else:
            best_matches[m] = [c, m, cs, ce, rev, b]

    return [v for _, v in best_matches.items()]


def largest_cluster(f, nucleotides):
    tab = (ll.strip().split('\t') for ll in open(f))
    clusters = {}
    largest_clusters = []

    for entry in tab:
        c = entry[0].split(' ')[0]
        m = entry[1].split('_')[1]
        pid = float(entry[2])

        if pid > 50.0:  # consider only matches with at least 50% of identity
            if (c, m) in clusters:
                clusters[(c, m)].append(entry)
            else:
                clusters[(c, m)] = [entry]

    for (c, m), entries in clusters.items():
        cs = [int(s) for _, _, _, _, _, _, s, _, _, _, _, _ in entries]
        ce = [int(e) for _, _, _, _, _, _, _, e, _, _, _, _ in entries]
        ms = (int(s) for _, _, _, _, _, _, _, _, s, _, _, _ in entries)
        me = (int(e) for _, _, _, _, _, _, _, _, _, e, _, _ in entries)
        b = max((float(b) for _, _, _, _, _, _, _, _, _, _, _, b in entries))
        rev = 0

        for s, e in zip(cs, ce):
            if s > e:  # check if the contig is reverse
                rev = 1
                break

        if not rev:  # if contig position are forward
            for s, e in zip(ms, me):
                if s > e:  # check if the marker is reverse
                    rev = 1
                    break

        largest_clusters.append([c, m, str(min(cs + ce)), str(max(cs + ce)), str(rev), str(b)])

    return largest_clusters


def gene_markers_extraction(inputs, input_folder, output_folder, extension, min_num_markers,
                            frameshifts=False, nproc=1, verbose=False):
    commands = []
    check_and_create_folder(output_folder, create=True, verbose=verbose)

    for f in glob.iglob(os.path.join(input_folder, '*.b6o')):
        f_clean = os.path.basename(f).replace('.b6o', extension)
        src_file = os.path.join(inputs[f_clean], f_clean)
        out_file = os.path.join(output_folder, f_clean)

        if os.path.isfile(src_file) and (not os.path.isfile(out_file)):
            commands.append((out_file, src_file, f, min_num_markers, frameshifts))

    if commands:
        info('Extracting markers from {} inputs\n'.format(len(commands)))
        terminating = mp.Event()

        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                [_ for _ in pool.imap_unordered(gene_markers_extraction_rec, commands, chunksize=1)]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('gene_markers_extraction crashed', init_new_line=True, exit=True)
    else:
        info('Markers already extracted\n')


def gene_markers_extraction_rec(x):
    if not terminating.is_set():
        try:
            t0 = time.time()
            out_file, src_file, b6o_file, min_num_markers, frameshifts = x
            out_file_seq = []
            contig2marker2b6o = {}
            info('Extracting "{}"\n'.format(b6o_file))

            for l in open(b6o_file):
                row = l.strip().split('\t')
                contig = row[0]
                marker = row[1]
                start = int(row[2])
                end = int(row[3])
                rev = bool(int(row[4]))

                if (contig in contig2marker2b6o) and (marker in contig2marker2b6o[contig]):
                    error('contig: {} and marker: {} already present into contig2marker2b6o'.format(contig, marker))

                if contig not in contig2marker2b6o:
                    contig2marker2b6o[contig] = {}

                contig2marker2b6o[contig][marker] = (end, start, rev) if end < start else (start, end, rev)

            for record in SimpleFastaParser(open(src_file)):
                fid = record[0].split(' ')[0]

                if fid not in contig2marker2b6o:
                    continue

                for marker in contig2marker2b6o[fid]:
                    s, e, rev = contig2marker2b6o[fid][marker]
                    idd = '{}_{}:'.format(fid, marker)
                    seq = Seq(record[1][s - 1:e]) if (s - 1 >= 0) else None

                    if not seq:  # skip empty sequences
                        continue

                    if rev:
                        idd += 'c'
                        seq = seq.reverse_complement()

                    out_file_seq.append(SeqRecord(seq, id='{}{}-{}'.format(idd, s, e), description=''))

                    if frameshifts:
                        if not rev:
                            if record[1][s:e]:  # skip empty sequences
                                out_file_seq.append(SeqRecord(Seq(record[1][s:e]),
                                                              id='{}{}-{}'.format(idd, s + 1, e), description=''))

                            if record[1][s + 1:e]:  # skip empty sequences
                                out_file_seq.append(SeqRecord(Seq(record[1][s + 1:e]),
                                                              id='{}{}-{}'.format(idd, s + 2, e), description=''))
                        else:
                            if (s - 1 >= 0) and (e - 1 >= 0) and (record[1][s - 1:e - 1]):  # skip empty sequences
                                out_file_seq.append(SeqRecord(Seq(record[1][s - 1:e - 1]).reverse_complement(),
                                                              id='{}{}-{}'.format(idd, s, e - 1), description=''))

                            if (s - 1 >= 0) and (e - 2 >= 0) and (record[1][s - 1:e - 2]):  # skip empty sequences
                                out_file_seq.append(SeqRecord(Seq(record[1][s - 1:e - 2]).reverse_complement(),
                                                              id='{}{}-{}'.format(idd, s, e - 2), description=''))

            len_out_file_seq = int(len(out_file_seq) / 3) if frameshifts else len(out_file_seq)

            if out_file_seq and (len_out_file_seq >= min_num_markers):
                with open(out_file, 'w') as f:
                    SeqIO.write(out_file_seq, f, 'fasta')

                t1 = time.time()
                info('"{}" generated in {}s\n'.format(out_file, int(t1 - t0)))
            else:
                info('Not enough markers ({}/{}) found in "{}"\n'.format(len_out_file_seq, min_num_markers, b6o_file))
        except Exception as e:
            terminating.set()
            remove_file(out_file)
            error(str(e), init_new_line=True)
            error('error while extracting\n    {}'.format('\n    '.join([str(a) for a in x])), init_new_line=True)
            raise
    else:
        terminating.set()


def fake_proteome(input_folder, output_folder, in_extension, out_extension, min_len_protein, nproc=1, verbose=False):
    commands = []
    check_and_create_folder(output_folder, create=True, verbose=verbose)

    for f in glob.iglob(os.path.join(input_folder, '*' + in_extension)):
        out = os.path.join(output_folder, os.path.splitext(os.path.basename(f))[0] + out_extension)

        if not os.path.isfile(out):
            commands.append((f, out, min_len_protein))

    if commands:
        info('Generate proteomes from {} genomes\n'.format(len(commands)))
        terminating = mp.Event()

        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                [_ for _ in pool.imap_unordered(fake_proteome_rec, commands, chunksize=1)]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('fake_proteomes crashed', init_new_line=True, exit=True)
    else:
        info('Fake proteomes already generated\n')


def fake_proteome_rec(x):
    if not terminating.is_set():
        try:
            t0 = time.time()
            inp, out, min_len_protein = x
            proteome = []
            info('Generating "{}"\n'.format(inp))

            for idd, seq in SimpleFastaParser(open(inp)):
                idd = idd.split(' ')[0]
                s, e = idd.split(':')[-1].split('-')

                if s.startswith('c'):
                    s = s[1:]

                while (len(seq) % 3) != 0:
                    seq = seq[:-1]

                seq_t = Seq.translate(Seq(seq), to_stop=True)

                if len(seq_t) >= min_len_protein:
                    proteome.append(SeqRecord(seq_t, id=idd, description=''))

            with open(out, 'w') as f:
                SeqIO.write(proteome, f, 'fasta')

            t1 = time.time()
            info('"{}" generated in {}s\n'.format(out, int(t1 - t0)))
        except Exception as e:
            terminating.set()
            remove_file(out)
            error(str(e), init_new_line=True)
            error('error while generating\n    {}'.format('\n    '.join([str(a) for a in x])), init_new_line=True)
            raise
    else:
        terminating.set()


def inputs2markers(input_folder, output_folder, min_num_entries, extension, verbose=False):
    markers2inputs = {}
    check_and_create_folder(output_folder, create=True, verbose=verbose)

    if os.path.isdir(output_folder):
        for f in glob.iglob(os.path.join(output_folder, '*' + extension)):
            info('Inputs already translated into markers\n')
            return

    for f in glob.iglob(os.path.join(input_folder, '*' + extension)):
        inp, _ = os.path.splitext(os.path.basename(f))

        for idd, seq in SimpleFastaParser(open(f)):
            marker = idd.split(' ')[0].split(':')[0].split('_')[-1]

            if marker in markers2inputs:
                markers2inputs[marker].append(SeqRecord(Seq(seq), id=inp, description=''))
            else:
                markers2inputs[marker] = [SeqRecord(Seq(seq), id=inp, description='')]

    for marker, sequences in markers2inputs.items():
        if len(sequences) >= min_num_entries:
            with open(os.path.join(output_folder, marker + extension), 'w') as f:
                SeqIO.write(sequences, f, 'fasta')
        elif verbose:
            info('"{}" discarded, not enough inputs ({}/{})\n'.format(marker, len(sequences), min_num_entries))


def msas(configs, key, input_folder, extension, output_folder, nproc=1, verbose=False):
    commands = []
    check_and_create_folder(output_folder, create=True, verbose=verbose)

    for inp in glob.iglob(os.path.join(input_folder, '*' + extension)):
        out = os.path.splitext(os.path.basename(inp))[0] + '.aln'

        if (not os.path.isfile(os.path.join(output_folder, out))) and \
           (not os.path.isfile(os.path.join(output_folder, out + '_alignment_masked.fasta'))) and \
           (not os.path.isfile(os.path.join(output_folder, out + '_insertion_columns.txt'))) and \
           (not os.path.isfile(os.path.join(output_folder, out + '_alignment.fasta'))) and \
           (not os.path.isfile(os.path.join(output_folder, out + '_pasta.fasta'))) and \
           (not os.path.isfile(os.path.join(output_folder, out + '_pasta.fasttree'))):
            commands.append((configs[key], inp, os.path.abspath(output_folder), out))

    if commands:
        info('Aligning {} markers (key: "{}")\n'.format(len(commands), key))
        terminating = mp.Event()

        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                [_ for _ in pool.imap_unordered(msas_rec, commands, chunksize=1)]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('msas crashed', init_new_line=True, exit=True)

    else:
        info('Markers already aligned (key: "{}")\n'.format(key))


def msas_rec(x):
    if not terminating.is_set():
        try:
            t0 = time.time()
            params, inp, out_p, out = x
            info('Aligning "{}"\n'.format(inp))
            cmd = compose_command(params, input_file=inp, output_path=out_p, output_file=out)
            inp_f = None
            out_f = sb.DEVNULL

            if cmd['stdin']:
                inp_f = open(cmd['stdin'], 'r')

            if cmd['stdout']:
                out_f = open(cmd['stdout'], 'w')

            try:
                sb.check_call(cmd['command_line'], stdin=inp_f, stdout=out_f,
                              stderr=sb.DEVNULL, env=cmd['env'])
            except Exception as e:
                terminating.set()
                remove_files([cmd['output_file'], cmd['output_file'] + '_alignment_masked.fasta',
                              cmd['output_file'] + '_insertion_columns.txt',
                              cmd['output_file'] + '_alignment.fasta',
                              cmd['output_file'] + '_pasta.fasta',
                              cmd['output_file'] + '_pasta.fasttree'], path=cmd['output_path'])
                error(str(e), init_new_line=True)
                error('error while aligning\n    {}'
                      .format('\n    '.join(['{:>12}: {}'.format(a, ' '.join(cmd[a]) if type(cmd[a]) is list else cmd[a])
                                             for a in ['command_line', 'stdin', 'stdout', 'env']])), init_new_line=True)
                raise

            if cmd['stdin']:
                inp_f.close()

            if cmd['stdout']:
                out_f.close()

            t1 = time.time()
            info('"{}" generated in {}s\n'.format(out, int(t1 - t0)))
        except Exception as e:
            terminating.set()
            error(str(e), init_new_line=True)
            error('error while aligning\n    {}'.format('\n    '.join([str(a) for a in x])), init_new_line=True)
            raise
    else:
        terminating.set()


def is_msa_empty(msa, path=None):
    msa_path = os.path.join(path, msa) if path else msa

    if os.path.isfile(msa_path):
        if [True for aln in AlignIO.read(msa_path, "fasta") if len(aln.seq) <= 0]:
            return True  # there is at least an empty sequence that shouldn't be there
    else:
        error('file "{}" not found'.format(msa_path))

    return False  # all entries checked and are not empty or file not found


def trim_gap_trim(configs, key, inputt, output_folder, nproc=1, verbose=False):
    commands = []
    check_and_create_folder(output_folder, create=True, verbose=verbose)

    if os.path.isdir(inputt):
        for inp in glob.iglob(os.path.join(inputt, '*.aln')):
            out = os.path.basename(inp)

            if not os.path.isfile(os.path.join(output_folder, out)):
                commands.append((configs[key], inp, output_folder, out))

        if not commands:  # aligned with UPP?
            for inp in glob.iglob(os.path.join(inputt, '*.aln_alignment_masked.fasta')):
                out = os.path.basename(inp).replace('_alignment_masked.fasta', '')

                if not os.path.isfile(os.path.join(output_folder, out)):
                    commands.append((configs[key], inp, output_folder, out))

    elif os.path.isfile(inputt):
        base, ext = os.path.splitext(os.path.basename(inputt))
        out = base + '.trim_gap_trim' + ext

        if not os.path.isfile(os.path.join(output_folder, out)):
            commands.append((configs[key], inputt, output_folder, out))
    else:
        error('unrecognized input "{}" is not a folder nor a file'.format(inputt), exit=True)

    if commands:
        info('Trimming gappy regions for {} markers (key: "{}")\n'.format(len(commands), key))
        terminating = mp.Event()

        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                [_ for _ in pool.imap_unordered(trim_gap_trim_rec, commands, chunksize=1)]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('trim_gap_trim crashed', init_new_line=True, exit=True)
    else:
        info('Markers already trimmed (key: "{}")\n'.format(key))


def trim_gap_trim_rec(x):
    if not terminating.is_set():
        try:
            t0 = time.time()
            params, inp, out_fld, out = x
            info('Trimming gappy regions "{}"\n'.format(inp))
            cmd = compose_command(params, input_file=inp, output_path=out_fld, output_file=out)
            inp_f = None
            out_f = sb.DEVNULL

            if cmd['stdin']:
                inp_f = open(cmd['stdin'], 'r')

            if cmd['stdout']:
                out_f = open(cmd['stdout'], 'w')

            try:
                sb.check_call(cmd['command_line'], stdin=inp_f, stdout=out_f, stderr=sb.DEVNULL, env=cmd['env'])
            except Exception as e:
                terminating.set()
                remove_file(cmd['output_file'], path=cmd['output_path'])
                error(str(e), init_new_line=True)
                error('error while trimming\n    {}'
                      .format('\n    '.join(['{:>12}: {}'.format(a, ' '.join(cmd[a]) if type(cmd[a]) is list else cmd[a])
                                             for a in ['command_line', 'stdin', 'stdout', 'env']])), init_new_line=True)
                raise

            if cmd['stdin']:
                inp_f.close()

            if cmd['stdout']:
                out_f.close()

            outt = cmd['output_file']

            if cmd['output_path'] and (not outt.startswith(cmd['output_path'])):
                outt = os.path.join(cmd['output_path'], cmd['output_file'])

            if is_msa_empty(outt):  # sanity check
                info('Removing generated empty MSAs "{}"\n'.format(outt))
                os.remove(outt)
            else:
                t1 = time.time()
                info('"{}" generated in {}s\n'.format(out, int(t1 - t0)))
        except Exception as e:
            terminating.set()
            error(str(e), init_new_line=True)
            error('error while trimming gappy regions\n    {}'.format('\n    '.join([str(a) for a in x])), init_new_line=True)
            raise
    else:
        terminating.set()


def trim_gap_perc(inputt, output_folder, gap_perc_threshold, nproc=1, verbose=False):
    commands = []
    check_and_create_folder(output_folder, create=True, verbose=verbose)

    if os.path.isdir(inputt):
        for inp in glob.iglob(os.path.join(inputt, '*.aln')):
            out = os.path.join(output_folder, os.path.basename(inp))

            if not os.path.isfile(out):
                commands.append((inp, out, gap_perc_threshold, verbose))

        if not commands:  # aligned with UPP?
            for inp in glob.iglob(os.path.join(inputt, '*.aln_alignment_masked.fasta')):
                out = os.path.join(output_folder, os.path.basename(inp).replace('_alignment_masked.fasta', ''))

                if not os.path.isfile(out):
                    commands.append((inp, out, gap_perc_threshold, verbose))
    elif os.path.isfile(inputt):
        base, ext = os.path.splitext(inputt)
        out = base + '.trim_gap_perc' + ext

        if not os.path.isfile(out):
            commands.append((inputt, out, gap_perc_threshold, verbose))
    else:
        error('unrecognized input "{}" is not a folder nor a file'.format(inputt), exit=True)

    if commands:
        info('Trimming gappy columns from {} markers\n'.format(len(commands)))
        terminating = mp.Event()

        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                [_ for _ in pool.imap_unordered(trim_gap_perc_rec, commands, chunksize=1)]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('trim_gap_perc crashed', init_new_line=True, exit=True)
    else:
        info('Markers already trimmed\n')


def trim_gap_perc_rec(x):
    if not terminating.is_set():
        try:
            t0 = time.time()
            inp, out, thr, verbose = x
            info('Trimming gappy columns "{}"\n'.format(inp))
            inp_aln = AlignIO.read(inp, "fasta")
            cols_to_remove = []
            sub_aln = []

            for i in range(len(inp_aln[0])):
                if gap_cost(inp_aln[:, i], norm=True) >= thr:
                        cols_to_remove.append(i)

            for aln in inp_aln:
                seq = ''.join([c for i, c in enumerate(aln.seq) if i not in cols_to_remove])

                if seq:
                    sub_aln.append(SeqRecord(Seq(seq), id=aln.id, description=''))

            if sub_aln:
                with open(out, 'w') as f:
                    AlignIO.write(MultipleSeqAlignment(sub_aln), f, "fasta")

                if is_msa_empty(out):  # sanity check
                    info('Removing generated empty MSAs "{}"\n'.format(out))
                    os.remove(out)
                else:
                    t1 = time.time()
                    info('"{}" generated in {}s\n'.format(out, int(t1 - t0)))
            elif verbose:
                    info('"{}" discarded because no columns retained while removing not variant sites (thr: {})\n'.format(inp, thr))

        except Exception as e:
            terminating.set()
            remove_file(out)
            error(str(e), init_new_line=True)
            error('error while trimming gappy columns\n    {}'.format('\n    '.join([str(a) for a in x])), init_new_line=True)
            raise
    else:
        terminating.set()


def trim_not_variant(inputt, output_folder, not_variant_threshold, nproc=1, verbose=False):
    commands = []
    check_and_create_folder(output_folder, create=True, verbose=verbose)

    if os.path.isdir(inputt):
        for inp in glob.iglob(os.path.join(inputt, '*.aln')):
            out = os.path.join(output_folder, os.path.basename(inp))

            if not os.path.isfile(out):
                commands.append((inp, out, not_variant_threshold, verbose))

        if not commands:  # aligned with UPP?
            for inp in glob.iglob(os.path.join(inputt, '*.aln_alignment_masked.fasta')):
                out = os.path.join(output_folder, os.path.basename(inp).replace('_alignment_masked.fasta', ''))

                if not os.path.isfile(out):
                    commands.append((inp, out, not_variant_threshold, verbose))
    elif os.path.isfile(inputt):
        base, ext = os.path.splitext(inputt)
        out = base + '.trim_not_variant' + ext

        if not os.path.isfile(out):
            commands.append((inputt, out, not_variant_threshold, verbose))
    else:
        error('unrecognized input "{}" is not a folder nor a file'.format(inputt), exit=True)

    if commands:
        info('Trimming not variant from {} markers\n'.format(len(commands)))
        terminating = mp.Event()

        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                [_ for _ in pool.imap_unordered(trim_not_variant_rec, commands, chunksize=1)]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('trim_not_variant crashed', init_new_line=True, exit=True)
    else:
        info('Markers already trimmed\n')


def trim_not_variant_rec(x):
    if not terminating.is_set():
        try:
            t0 = time.time()
            inp, out, thr, verbose = x
            info('Trimming not variant "{}"\n'.format(inp))
            inp_aln = AlignIO.read(inp, "fasta")
            cols_to_remove = []
            sub_aln = []

            for i in range(len(inp_aln[0])):
                nrows = len(inp_aln) - inp_aln[:, i].count('-')

                for aa, fq in Counter(inp_aln[:, i].replace('-', '')).items():
                    if (fq / nrows) >= thr:
                        cols_to_remove.append(i)
                        break

            for aln in inp_aln:
                seq = ''.join([c for i, c in enumerate(aln.seq) if i not in cols_to_remove])

                if seq:
                    sub_aln.append(SeqRecord(Seq(seq), id=aln.id, description=''))

            if sub_aln:
                with open(out, 'w') as f:
                    AlignIO.write(MultipleSeqAlignment(sub_aln), f, "fasta")

                if is_msa_empty(out):  # sanity check
                    info('Removing generated empty MSAs "{}"\n'.format(out))
                    os.remove(out)
                else:
                    t1 = time.time()
                    info('"{}" generated in {}s\n'.format(out, int(t1 - t0)))
            elif verbose:
                    info('"{}" discarded because no columns retained while removing not variant sites (thr: {})\n'.format(inp, thr))
        except Exception as e:
            terminating.set()
            remove_file(out)
            error(str(e), init_new_line=True)
            error('error while trimming not variant\n    {}'.format('\n    '.join([str(a) for a in x])), init_new_line=True)
            raise
    else:
        terminating.set()


def remove_fragmentary_entries(input_folder, data_folder, output_folder, fragmentary_threshold, min_num_entries, nproc=1, verbose=False):
    commands = []
    frag_entries = []

    if os.path.isfile(os.path.join(data_folder, 'fragmentary_entries.pkl')):
        info('Fragmentary entries already removed\n')

        if verbose:
            info('Loading fragmentary entries from "{}"\n'.format(os.path.join(data_folder, 'fragmentary_entries.pkl')))

        with open(os.path.join(data_folder, 'fragmentary_entries.pkl'), 'rb') as f:
            frag_entries = pickle.load(f)

        return frag_entries

    check_and_create_folder(output_folder, create=True, verbose=verbose)

    for inp in glob.iglob(os.path.join(input_folder, '*.aln')):
        out = os.path.join(output_folder, os.path.basename(inp))

        if not os.path.isfile(out):
            commands.append((inp, out, fragmentary_threshold, min_num_entries, verbose))

        if not commands:  # aligned with UPP?
            for inp in glob.iglob(os.path.join(input_folder, '*.aln_alignment_masked.fasta')):
                out = os.path.join(output_folder, os.path.basename(inp).replace('_alignment_masked.fasta', ''))

                if not os.path.isfile(out):
                    commands.append((inp, out, fragmentary_threshold, min_num_entries, verbose))

    if commands:
        info('Checking {} alignments for fragmentary entries (thr: {})\n'.format(len(commands), fragmentary_threshold))
        terminating = mp.Event()

        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                frag_entries = [a for a in pool.imap_unordered(remove_fragmentary_entries_rec, commands, chunksize=1) if a]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('remove_fragmentary_entries crashed', init_new_line=True, exit=True)

            if frag_entries:
                with open(os.path.join(data_folder, 'fragmentary_entries.pkl'), 'wb') as f:
                    pickle.dump(frag_entries, f, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        info('Fragmentary entries already removed\n')

    return frag_entries


def remove_fragmentary_entries_rec(x):
    if not terminating.is_set():
        try:
            t0 = time.time()
            inp, out, frag_thr, min_num_entries, verbose = x
            info('Fragmentary "{}"\n'.format(inp))
            inp_aln = AlignIO.read(inp, "fasta")
            out_aln = []

            for aln in inp_aln:
                if gap_cost(aln.seq) < frag_thr:
                    out_aln.append(aln)

            if len(out_aln) >= min_num_entries:
                with open(out, 'w') as f:
                    AlignIO.write(MultipleSeqAlignment(out_aln), f, 'fasta')

                t1 = time.time()
                info('"{}" generated in {}s\n'.format(out, int(t1 - t0)))
            elif verbose:
                info('"{}" discarded, not enough inputs ({}/{})\n'
                     .format(inp, len(out_aln), min_num_entries))
                return inp

            return None
        except Exception as e:
            terminating.set()
            remove_file(out)
            error(str(e), init_new_line=True)
            error('error while removing fragmentary\n    {}'.format('\n    '.join([str(a) for a in x])),
                  init_new_line=True)
            raise
    else:
        terminating.set()


def inputs_list(input_folder, extension, output_file, nproc=1, verbose=False):
    file_input_list = []

    if os.path.isfile(output_file):
        with open(output_file, 'rb') as f:
            file_input_list = pickle.load(f)
    else:
        commands = glob.glob(os.path.join(input_folder, '*' + extension))

        if commands:
            tmp_input_list = None
            terminating = mp.Event()

            with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
                try:
                    tmp_input_list = [a for a in pool.imap_unordered(inputs_list_rec, commands, chunksize=1) if a]
                except Exception as e:
                    error(str(e), init_new_line=True)
                    error('inputs_list crashed', init_new_line=True, exit=True)

            if tmp_input_list:
                for inputs in tmp_input_list:
                    file_input_list += inputs

                file_input_list = list(set(file_input_list))

                with open(output_file, 'wb') as f:
                    pickle.dump(file_input_list, f, protocol=pickle.HIGHEST_PROTOCOL)

    return file_input_list


def inputs_list_rec(x):
    if not terminating.is_set():
        try:
            return [aln.id for aln in AlignIO.read(x, "fasta")]
        except Exception as e:
            terminating.set()
            error(str(e), init_new_line=True)
            error('error while computing inputs list\n    input_file: {}'.format(x),
                  init_new_line=True)
            raise
    else:
        terminating.set()


def subsample(input_folder, output_folder, positions_function, scoring_function, submat, unknown_fraction=0.3, nproc=1, verbose=False):
    commands = []
    mat = None

    if not os.path.isfile(submat):
        error('could not find substitution matrix "{}"'.format(submat), exit=True)
    else:
        with open(submat, 'rb') as f:
            mat = pickle.load(f)

        if verbose:
            info('Substitution matrix "{}" loaded\n'.format(submat))

    check_and_create_folder(output_folder, create=True, verbose=verbose)

    for inp in glob.iglob(os.path.join(input_folder, '*.aln')):
        out = os.path.join(output_folder, os.path.basename(inp))

        if not os.path.isfile(out):
            commands.append((inp, out, positions_function, scoring_function, unknown_fraction, mat))

    if not commands:  # aligned with UPP?
        for inp in glob.iglob(os.path.join(input_folder, '*_alignment_masked.fasta')):
            out = os.path.join(output_folder, os.path.basename(inp).replace('_alignment_masked.fasta', '.aln'))

            if not os.path.isfile(out):
                commands.append((inp, out, positions_function, scoring_function, unknown_fraction, mat))

    if commands:
        info('Subsampling {} markers\n'.format(len(commands)))
        terminating = mp.Event()

        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                [_ for _ in pool.imap_unordered(subsample_rec, commands, chunksize=1)]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('subsample crashed', init_new_line=True, exit=True)
    else:
        info('Markers already subsampled\n')


def subsample_rec(x):
    if not terminating.is_set():
        try:
            t0 = time.time()
            inp, out, npos_function, score_function, unknown_fraction, mat = x
            info('Subsampling "{}"\n'.format(inp))
            inp_aln = AlignIO.read(inp, "fasta")
            len_seq = inp_aln.get_alignment_length()
            scores = []
            out_aln = []

            for i in range(len(inp_aln[0])):
                col = Counter(inp_aln[:, i].upper())

                if (len(col) == 1) or \
                   (len(col) == 2 and ("-" in col or "X" in col)) or \
                   (len(col) == 3 and "X" in col and "-" in col):
                    continue

                if (col["-"] + col["X"]) > (len(inp_aln[:, i]) * unknown_fraction):
                    continue

                scores.append((score_function(inp_aln[:, i], mat), i))

            try:
                marker = os.path.splitext(os.path.basename(inp))[0][1:]
                marker = int(marker)
            except Exception as _:
                marker = None

            npos = npos_function(marker, len_seq)
            best_npos = [p for _, p in sorted(scores)[-npos:]]

            for aln in inp_aln:
                seq = ''.join([c for i, c in enumerate(aln.seq) if i in best_npos])
                out_aln.append(SeqRecord(Seq(seq), id=aln.id, description=''))

            with open(out, 'w') as f:
                AlignIO.write(MultipleSeqAlignment(out_aln), f, 'fasta')

            t1 = time.time()
            info('"{}" generated in {}s\n'.format(out, int(t1 - t0)))
        except Exception as e:
            terminating.set()
            remove_file(out)
            error(str(e), init_new_line=True)
            error('error while subsampling\n    {}'.format('\n    '.join([str(a) for a in x])),
                  init_new_line=True)
            raise
    else:
        terminating.set()


def phylophlan(marker, len_seq):
    # return max(int(max(int((400-marker)*30/400.0), 1)**2/30.0), 3)  # ~4k AAs (original PhyloPhlAn formulae)
    return max(int(math.ceil((max(int(math.ceil(((400 - marker) * 30) / 400.0)), 1)**2) / 30.0)), 3)  # ~4.6k AAs


def onethousand(marker, len_seq):
    return 1000 if 1000 < len_seq else len_seq


def sevenhundred(marker, len_seq):
    return 700 if 700 < len_seq else len_seq


def fivehundred(marker, len_seq):
    return 500 if 500 < len_seq else len_seq


def threehundred(marker, len_seq):
    return 300 if 300 < len_seq else len_seq


def onehundred(marker, len_seq):
    return 100 if 100 < len_seq else len_seq


def fifty(marker, len_seq):
    return 50 if 50 < len_seq else len_seq


def twentyfive(marker, len_seq):
    return 25 if 25 < len_seq else len_seq


def tenpercent(marker, len_seq):
    return int(math.ceil(len_seq * 0.1))


def twentyfivepercent(marker, len_seq):
    return int(math.ceil(len_seq * 0.25))


def fiftypercent(marker, len_seq):
    return int(math.ceil(len_seq * 0.5))


def trident(seq, submat, alpha=1, beta=0.5, gamma=3):
    return ((1 - symbol_diversity(seq))**alpha *
            (1 - stereochemical_diversity(seq, submat))**beta *
            (1 - gap_cost(seq))**gamma)


def muscle(seq, submat):
    combos = [submat[a, b] for a, b in combinations(seq.upper().replace('-', ''), 2) if (a, b) in submat]
    retval = 0.0

    if len(combos):
        # average score over pairs of letters in the column
        retval = sum(combos) / len(combos)

    return retval


def random(seq, submat):
    return lib_random.random()


def symbol_diversity(seq, log_base=21):
    """
    Sander C, Schneider R. Database of homology-derived protein
    structures and the structural meaning of sequence alignment.
    Proteins 1991;9:56 - 68.
    """
    sh = 0.0

    for aa, abs_freq in Counter(seq.upper()).items():
        rel_freq = abs_freq / len(seq)
        sh -= rel_freq * math.log(rel_freq)

    sh /= math.log(min(len(seq), log_base))
    return sh if (sh > 0.15) and (sh < 0.85) else 0.99


def stereochemical_diversity(seq, submat):
    """
    Valdar W. Scoring Residue Conservation. PROTEINS: Structure,
    Function, and Genetics 48:227-241 (2002)
    """
    set_seq = set(seq.upper())
    aa_avg = sum([normalized_submat_scores(aa, submat) for aa in set_seq])
    aa_avg /= len(set_seq)
    r = sum([abs(aa_avg - normalized_submat_scores(aa, submat)) for aa in set_seq])
    r /= len(set_seq)
    r /= math.sqrt(20 * (max(submat.values()) - min(submat.values()))**2)
    return r


def normalized_submat_scores(aa, submat):
    """
    Karlin S, Brocchieri L. Evolutionary conservation of RecA genes in
    relation to protein structure and function. J Bacteriol 1996;178:
    1881-1894.
    """
    aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
           'T', 'V', 'W', 'Y']
    aa = aa.upper()
    m = 0.0

    if (aa != '-') and (aa != 'X') and (aa != 'Z'):
        for bb in aas:
            try:
                m += submat[(aa, bb)] / math.sqrt(submat[(aa, aa)] * submat[(bb, bb)])
            except Exception as e:
                error('{}\n    {}'
                      .format(str(e), '\n'.join(['    {}: {}'.format(l, v) for l, v in
                                                 zip(['aa', 'bb', 'submat[(aa, bb)]', 'submat[(aa, aa)]', 'submat[(bb, bb)]'],
                                                     [aa, bb, submat[(aa, bb)], submat[(aa, aa)], submat[(bb, bb)]])])))
                error('error while normalizing submat scores', exit=True)

    return m


def gap_cost(seq, norm=True):
    gaps = seq.count('-')

    if norm:
        if len(seq) != 0:
            gaps /= len(seq)
        else:
            gaps = 1  # in the range [0, 1], 1 plays the role of infinity

    return gaps


def load_substitution_model(input_file):
    if not input_file:
        error('mapping file (--maas) not specified', exit=True)
    elif not os.path.isfile(input_file):
        error('file "{}" not found'.format(input_file), exit=True)

    return dict(((line.strip().split('\t')[0], line.strip().split('\t')[1])
                 for line in open(input_file) if not line.startswith('#')))


def concatenate(all_inputs, input_folder, output_file, sort=False, verbose=False):
    t0 = time.time()

    if os.path.isfile(output_file):
        info('Alignments already merged "{}"\n'.format(output_file))
        return

    info('Concatenating alignments\n')
    all_inputs = set(all_inputs)
    inputs2alignments = dict(((inp, SeqRecord(Seq(''), id=inp, description='')) for inp in all_inputs))
    markers = glob.glob(os.path.join(input_folder, '*.aln'))

    if not len(markers):  # aligned with UPP?
        markers = glob.glob(os.path.join(input_folder, '*.aln_alignment_masked.fasta'))

    if sort:
        markers = sorted(markers)

    if len(markers):
        for a in markers:
            alignment_length = None
            current_inputs = []

            for fid, seq in SimpleFastaParser(open(a)):
                idd = fid.split(' ')[0]
                current_inputs.append(idd)
                inputs2alignments[idd].seq += seq

                if not alignment_length:
                    alignment_length = len(seq)
                elif alignment_length != len(seq):
                    error('wrong alignment length ({} != {})... Something is wrong'.format(alignment_length, len(seq)))

            current_inputs = set(current_inputs)

            for inp in all_inputs - current_inputs:
                inputs2alignments[inp].seq += Seq('-' * alignment_length)

        with open(output_file, 'w') as f:
            # discard inputs with only gaps in alignment (RAxML gives error)
            SeqIO.write([v for _, v in inputs2alignments.items() if gap_cost(v.seq, norm=True) < 1], f, "fasta")

        t1 = time.time()
        info('Alignments concatenated "{}" in {}s\n'.format(output_file, int(t1 - t0)))
    else:
        error('No alignments found to concatenate', init_new_line=True, exit=True)


def build_gene_tree(configs, key, sub_mod, input_folder, output_folder, nproc=1, verbose=False):
    commands = []
    check_and_create_folder(output_folder, create=True, verbose=verbose)

    for inp in glob.iglob(os.path.join(input_folder, '*.aln')):
        marker, _ = os.path.splitext(os.path.basename(inp))
        out = marker + '.tre'
        model = sub_mod[marker]

        if (not os.path.isfile(os.path.join(output_folder, out))) and \
           (not os.path.isfile(os.path.join(output_folder, 'RAxML_bestTree.' + out))) and \
           (not os.path.isfile(os.path.join(output_folder, 'RAxML_info.' + out))) and \
           (not os.path.isfile(os.path.join(output_folder, 'RAxML_log.' + out))) and \
           (not os.path.isfile(os.path.join(output_folder, 'RAxML_result.' + out))):
            commands.append((configs[key], model, inp, os.path.abspath(output_folder), out))

    if not commands:  # aligned with UPP?
        for inp in glob.iglob(os.path.join(input_folder, '*.aln_alignment_masked.fasta')):
            marker, _ = os.path.basename(inp).replace('.aln_alignment_masked.fasta', '')
            out = marker + '.tre'
            model = sub_mod[marker]

        if (not os.path.isfile(os.path.join(output_folder, out))) and \
           (not os.path.isfile(os.path.join(output_folder, 'RAxML_bestTree.' + out))) and \
           (not os.path.isfile(os.path.join(output_folder, 'RAxML_info.' + out))) and \
           (not os.path.isfile(os.path.join(output_folder, 'RAxML_log.' + out))) and \
           (not os.path.isfile(os.path.join(output_folder, 'RAxML_result.' + out))):
            commands.append((configs[key], model, inp, os.path.abspath(output_folder), out))

    if commands:
        info('Building {} gene trees\n'.format(len(commands)))
        terminating = mp.Event()

        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                [_ for _ in pool.imap_unordered(build_gene_tree_rec, commands, chunksize=1)]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('build_gene_tree crashed', init_new_line=True, exit=True)
    else:
        info('Gene trees already built\n')


def build_gene_tree_rec(x):
    if not terminating.is_set():
        try:
            t0 = time.time()
            params, model, inp, abs_wf, out = x
            info('Building gene tree "{}"\n'.format(inp))
            cmd = compose_command(params, sub_mod=model, input_file=inp,
                                  output_path=abs_wf, output_file=out)
            inp_f = None
            out_f = sb.DEVNULL

            if cmd['stdin']:
                inp_f = open(cmd['stdin'], 'r')

            if cmd['stdout']:
                out_f = open(cmd['stdout'], 'w')

            try:
                sb.check_call(cmd['command_line'], stdin=inp_f, stdout=out_f,
                              stderr=sb.DEVNULL, env=cmd['env'])
            except Exception as e:
                terminating.set()
                remove_files([cmd['output_file'], 'RAxML_bestTree.' + cmd['output_file'],
                              'RAxML_info.' + cmd['output_file'], 'RAxML_log.' + cmd['output_file'],
                              'RAxML_result.' + cmd['output_file']], path=cmd['output_path'])
                error(str(e), init_new_line=True)
                error('error while building gene tree\n    {}'
                      .format('\n    '.join(['{:>12}: {}'.format(a, ' '.join(cmd[a]) if type(cmd[a]) is list else cmd[a])
                                             for a in ['command_line', 'stdin', 'stdout', 'env']])), init_new_line=True)
                raise

            if cmd['stdin']:
                inp_f.close()

            if cmd['stdout']:
                out_f.close()

            t1 = time.time()
            info('"{}" generated in {}s\n'.format(out, int(t1 - t0)))
        except Exception as e:
            terminating.set()
            error(str(e), init_new_line=True)
            error('error while building gene tree\n    {}'.format('\n    '.join([str(a) for a in x])), init_new_line=True)
            raise
    else:
        terminating.set()


def resolve_polytomies(inputt, output, nproc=1, verbose=False):
    commands = []

    if os.path.isfile(inputt):
        output_path = os.path.dirname(output)
        output_file = os.path.basename(output)

        if (not os.path.isfile(output)) and \
           (not os.path.isfile(os.path.join(output_path, "RAxML_bestTree." + output_file))):
            commands.append((inputt, output))

    elif os.path.isdir(inputt):
        check_and_create_folder(output, create=True, verbose=verbose)

        for inp in glob.iglob(os.path.join(inputt, '*.tre')):
            out = os.path.basename(inp)

            if (not os.path.isfile(os.path.join(output, out))) and \
               (not os.path.isfile(os.path.join(output, "RAxML_bestTree." + out))):
                commands.append((inp, os.path.join(output, out)))

    if commands:
        info('Resolving {} polytomies\n'.format(len(commands)))
        terminating = mp.Event()

        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                [_ for _ in pool.imap_unordered(resolve_polytomies_rec, commands, chunksize=1)]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('resolve_polytomies crashed', init_new_line=True, exit=True)
    else:
        info('Polytomies already resolved\n')


def resolve_polytomies_rec(x):
    if not terminating.is_set():
        try:
            t0 = time.time()
            inp_f, out_f = x
            info('Resolving polytomies for "{}"\n'.format(inp_f))
            tree = dendropy.Tree.get(path=inp_f, schema="newick", preserve_underscores=True)
            tree.resolve_polytomies()
            tree.write(path=out_f, schema="newick")
            t1 = time.time()
            info('"{}" generated in {}s\n'.format(out_f, int(t1 - t0)))
        except Exception as e:
            terminating.set()
            error(str(e), init_new_line=True)
            error('error while resolving polytomies\n    {}'.format('\n    '.join([str(a) for a in x])),
                  init_new_line=True)
            raise
    else:
        terminating.set()


def refine_gene_tree(configs, key, sub_mod, input_alns, input_trees, output_folder, nproc=1, verbose=False):
    commands = []
    check_and_create_folder(output_folder, create=True, verbose=verbose)

    for inp in glob.iglob(os.path.join(input_alns, '*.aln')):
        marker, _ = os.path.splitext(os.path.basename(inp))
        out = marker + '.tre'
        starting_tree = os.path.join(input_trees, out)
        model = sub_mod[marker]

        if os.path.isfile(starting_tree):
            if (not os.path.isfile(os.path.join(output_folder, out))) and \
               (not os.path.isfile(os.path.join(output_folder, 'RAxML_bestTree.' + out))) and \
               (not os.path.isfile(os.path.join(output_folder, 'RAxML_info.' + out))) and \
               (not os.path.isfile(os.path.join(output_folder, 'RAxML_log.' + out))) and \
               (not os.path.isfile(os.path.join(output_folder, 'RAxML_result.' + out))):
                commands.append((configs[key], model, inp, starting_tree, os.path.abspath(output_folder), out))
        else:
            error('starting tree "{}" not found in "{}", built from "{}"'.format(starting_tree, input_trees, inp))

    if commands:
        info('Refining {} gene trees\n'.format(len(commands)))
        terminating = mp.Event()

        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                [_ for _ in pool.imap_unordered(refine_gene_tree_rec, commands, chunksize=1)]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('refine_gene_tree crashed', init_new_line=True, exit=True)
    else:
        info('Gene trees already refined\n')


def refine_gene_tree_rec(x):
    if not terminating.is_set():
        try:
            t0 = time.time()
            params, model, inp, st, abs_wf, out = x
            info('Refining gene tree "{}"\n'.format(inp))
            cmd = compose_command(params, sub_mod=model, input_file=inp, database=st, output_path=abs_wf, output_file=out)
            inp_f = None
            out_f = sb.DEVNULL

            if cmd['stdin']:
                inp_f = open(cmd['stdin'], 'r')

            if cmd['stdout']:
                out_f = open(cmd['stdout'], 'w')

            try:
                sb.check_call(cmd['command_line'], stdin=inp_f, stdout=out_f,
                              stderr=sb.DEVNULL, env=cmd['env'])
            except Exception as e:
                terminating.set()
                remove_files([cmd['output_file'], 'RAxML_bestTree.' + cmd['output_file'],
                              'RAxML_info.' + cmd['output_file'], 'RAxML_log.' + cmd['output_file'],
                              'RAxML_result.' + cmd['output_file']], path=cmd['output_path'])
                error(str(e), init_new_line=True)
                error('error while executing\n    {}'
                      .format('\n    '.join(['{:>12}: {}'.format(a, ' '.join(cmd[a]) if type(cmd[a]) is list else cmd[a])
                                             for a in ['command_line', 'stdin', 'stdout', 'env']])), init_new_line=True)
                raise

            if cmd['stdin']:
                inp_f.close()

            if cmd['stdout']:
                out_f.close()

            t1 = time.time()
            info('"{}" generated in {}s\n'.format(out, int(t1 - t0)))
        except Exception as e:
            terminating.set()
            error(str(e), init_new_line=True)
            error('error while refining gene tree\n    {}'
                  .format('\n    '.join([str(a) for a in x])), init_new_line=True)
            raise
    else:
        terminating.set()


def merging_gene_trees(trees_folder, output_file, verbose=False):
    if os.path.exists(output_file):
        info('Gene trees already merged "{}"\n'.format(output_file))
        return

    t0 = time.time()
    raxml = False
    info('Merging gene trees\n')

    with open(output_file, 'w') as f:
        # RAxML files
        for gtree in glob.iglob(os.path.join(trees_folder, "*bestTree*.tre")):
            raxml = True

            with open(gtree) as g:
                f.write(g.read())

        if not raxml:
            for gtree in glob.iglob(os.path.join(trees_folder, "*.tre")):
                with open(gtree) as g:
                    f.write(g.read())

    t1 = time.time()
    info('Gene trees merged into "{}" in {}s\n'
         .format(output_file, int(t1 - t0)))


def build_phylogeny(configs, key, inputt, output_path, output_tree, nproc=1,
                    verbose=False):
    if os.path.isfile(os.path.join(output_path, output_tree)) or \
       os.path.isfile(os.path.join(output_path, 'RAxML_bestTree.' + output_tree)) or \
       os.path.isfile(os.path.join(output_path, 'RAxML_info.' + output_tree)) or \
       os.path.isfile(os.path.join(output_path, 'RAxML_log.' + output_tree)) or \
       os.path.isfile(os.path.join(output_path, 'RAxML_result.' + output_tree)):
        info('Phylogeny "{}" (or "RAxML_*.{}") already built\n'
             .format(os.path.join(output_path, output_tree), output_tree))
        return

    nproc_loc = nproc

    # RAxML is slower when using more than 20 cores
    if (nproc > 20) and ('raxml' in configs[key]['program_name'].lower()):
        nproc_loc = 20

        if verbose:
            info('Reducing number of RAxML threads to {}, as it appears to underperform with more threads\n'.format(nproc_loc))

    # threaded RAxML complains if number of threads is 1!
    if (nproc_loc == 1) and ('raxml' in configs[key]['program_name'].lower()) and ('threads' in configs[key]):
        nproc_loc = 2

        if verbose:
            info('Setting RAxML threads to {}, cannot use 1 thread with RAxML parallel version\n'.format(nproc_loc))

    t0 = time.time()
    info('Building phylogeny "{}"\n'.format(inputt))
    cmd = compose_command(configs[key], input_file=inputt, output_path=output_path,
                          output_file=output_tree, nproc=nproc_loc)
    inp_f = None
    out_f = sb.DEVNULL

    if cmd['stdin']:
        inp_f = open(cmd['stdin'], 'r')

    if cmd['stdout']:
        out_f = open(cmd['stdout'], 'w')

    try:
        sb.check_call(cmd['command_line'], stdin=inp_f, stdout=out_f, stderr=sb.DEVNULL,
                      env=cmd['env'])
    except Exception as e:
        remove_files([cmd['output_file'], 'RAxML_bestTree.' + cmd['output_file'],
                      'RAxML_info.' + cmd['output_file'], 'RAxML_log.' + cmd['output_file'],
                      'RAxML_result.' + cmd['output_file']], path=cmd['output_path'])
        error(str(e), init_new_line=True)
        error('error while executing\n    {}'
              .format('\n    '.join(['{:>12}: {}'.format(a, ' '.join(cmd[a]) if type(cmd[a]) is list else cmd[a])
                                     for a in ['command_line', 'stdin', 'stdout', 'env']])), init_new_line=True, exit=True)

    if cmd['stdin']:
        inp_f.close()

    if cmd['stdout']:
        out_f.close()

    t1 = time.time()
    info('Phylogeny "{}" built in {}s\n'.format(output_tree, int(t1 - t0)))


def refine_phylogeny(configs, key, inputt, starting_tree, output_path, output_tree, nproc=1, verbose=False):
    if os.path.isfile(os.path.join(output_path, output_tree)) or \
       os.path.isfile(os.path.join(output_path, 'RAxML_bestTree.' + output_tree)) or \
       os.path.isfile(os.path.join(output_path, 'RAxML_info.' + output_tree)) or \
       os.path.isfile(os.path.join(output_path, 'RAxML_log.' + output_tree)) or \
       os.path.isfile(os.path.join(output_path, 'RAxML_result.' + output_tree)):
        info('Phylogeny "{}" already refined\n'.format(output_tree))
        return

    t0 = time.time()
    info('Refining phylogeny "{}"\n'.format(starting_tree))

    nproc_loc = nproc

    # RAxML is slower when using more than 20 cores
    if (nproc > 20) and ('raxml' in configs[key]['program_name'].lower()):
        nproc_loc = 20

        if verbose:
            info('Reducing number of RAxML threads to {}, as it appears to underperform with more threads\n'.format(nproc_loc))

    # threaded RAxML complains if number of threads is 1!
    if (nproc_loc == 1) and ('raxml' in configs[key]['program_name'].lower()) and ('threads' in configs[key]):
        nproc_loc = 2

        if verbose:
            info('Setting RAxML threads to {}, cannot use 1 thread with RAxML parallel version\n'.format(nproc_loc))

    cmd = compose_command(configs[key], input_file=inputt, database=starting_tree,
                          output_path=output_path, output_file=output_tree, nproc=nproc_loc)
    inp_f = None
    out_f = sb.DEVNULL

    if cmd['stdin']:
        inp_f = open(cmd['stdin'], 'r')

    if cmd['stdout']:
        out_f = open(cmd['stdout'], 'w')

    try:
        sb.check_call(cmd['command_line'], stdin=inp_f, stdout=out_f, stderr=sb.DEVNULL,
                      env=cmd['env'])
    except Exception as e:
        remove_files([cmd['output_file'], 'RAxML_bestTree.' + cmd['output_file'],
                      'RAxML_info.' + cmd['output_file'], 'RAxML_log.' + cmd['output_file'],
                      'RAxML_result.' + cmd['output_file']], path=cmd['output_path'])
        error(str(e), init_new_line=True)
        error('error while executing\n    {}'
              .format('\n    '.join(['{:>12}: {}'.format(a, ' '.join(cmd[a]) if type(cmd[a]) is list else cmd[a])
                                     for a in ['command_line', 'stdin', 'stdout', 'env']])), init_new_line=True, exit=True)

    if cmd['stdin']:
        inp_f.close()

    if cmd['stdout']:
        out_f.close()

    t1 = time.time()
    info('Phylogeny "{}" refined in {}s\n'.format(output_tree, int(t1 - t0)))


def compute_dists(s1, s2):
    idd = [s for s in range(len(s1)) if (s1[s] not in ['-', 'N']) and (s2[s] not in ['-', 'N'])]
    d = sum((a != b for a, b in zip((s1[s] for s in idd), (s2[s] for s in idd))))

    return (d, len(idd))


def mutation_rates(input_folder, output_folder, nproc=1, verbose=False):
    check_and_create_folder(output_folder, create=True, verbose=verbose)
    commands = [(inp, output_folder, os.path.splitext(os.path.basename(inp))[0], verbose)
                for inp in glob.iglob(os.path.join(input_folder, "*.aln"))
                if not os.path.exists(os.path.join(output_folder, os.path.splitext(os.path.basename(inp))[0] + '.tsv'))]

    if commands:
        info('Computing mutation rates for {} markers\n'.format(len(commands)))
        terminating = mp.Event()

        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                [_ for _ in pool.imap_unordered(mutation_rates_rec, commands, chunksize=1)]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('mutation_rates crashed', init_new_line=True, exit=True)
    else:
        info('Mutation rates already computed\n')


def mutation_rates_rec(x):
    if not terminating.is_set():
        try:
            inp_aln, out_fld, fn, verbose = x
            out_tsv = os.path.join(out_fld, fn + '.tsv')
            out_pkl = os.path.join(out_fld, fn + '.pkl')

            if verbose:
                info('Computing mutation rates of "{}"\n'.format(inp_aln))

            inp = AlignIO.read(inp_aln, "fasta")
            dists = dict([((inp[i].id, inp[j].id), compute_dists(inp[i], inp[j]))
                          for i in range(len(inp)) for j in range(i + 1, len(inp))])

            if verbose:
                info('Serializing to "{}"\n'.format(out_pkl))

            with open(out_pkl, 'wb') as f:
                pickle.dump(dists, f, protocol=pickle.HIGHEST_PROTOCOL)

            out_tbl = [['ids'] + [inp[i].id for i in range(len(inp))]]  # header

            for i in range(len(inp)):
                row = [inp[i].id]

                for j in range(len(inp)):
                    a, b = (i, j) if (inp[i].id, inp[j].id) in dists else (j, i)

                    if i < j:  # upper triangular
                        row.append(str(float(dists[(inp[a].id, inp[b].id)][0] / dists[(inp[a].id, inp[b].id)][1]))
                                   if dists[(inp[a].id, inp[b].id)][1] > 0 else 'NaN')
                    elif i > j:  # lower triangular
                        row.append('/'.join([str(e) for e in dists[(inp[a].id, inp[b].id)]]))
                    else:
                        row.append('0')

                out_tbl.append(row)

            if verbose:
                info('Writing output to "{}"\n'.format(out_tsv))

            with open(out_tsv, 'w') as f:
                f.write('\n'.join(['\t'.join(r) for r in out_tbl]))
        except Exception as e:
            terminating.set()
            remove_file(out_tsv)
            error(str(e), init_new_line=True)
            error('error while computing mutation rates\n    {}'.format('\n    '.join([str(a) for a in x])),
                  init_new_line=True)
            raise
    else:
        terminating.set()


def aggregate_mutation_rates(input_folder, output_file, verbose=False):
    mutation_rates = glob.glob(os.path.join(input_folder, "*.pkl"))
    aggregated = {}
    ids = set()

    if not mutation_rates:
        info('No mutation rates to aggregate\n')
        return

    if os.path.exists(output_file):
        info('Mutation rates already aggregated\n')
        return

    if verbose:
        info('Aggregating {} mutation rates\n'.format(len(mutation_rates)))

    for mr in mutation_rates:
        with open(mr, 'rb') as f:
            mr_dists = pickle.load(f)

        for k, d in mr_dists.items():
            a, b = k if k in aggregated else (k[1], k[0])
            ids.update(k)

            if (a, b) in aggregated:
                aggregated[(a, b)] = (d[0] + aggregated[(a, b)][0], d[1] + aggregated[(a, b)][1])
            else:
                aggregated[k] = d  # the assignment of (a, b) will swap the order if not in aggregated!

    s_ids = sorted(ids)
    out_tbl = [['ids'] + s_ids]  # header

    for i in s_ids:
        row = [i]

        for j in s_ids:
            a, b = (i, j) if (i, j) in aggregated else (j, i)

            if i < j:  # upper triangular
                row.append(str(float(aggregated[(a, b)][0] / aggregated[(a, b)][1]))
                           if aggregated[(a, b)][1] > 0 else 'NaN')
            elif i > j:  # lower triangular
                row.append('/'.join([str(e) for e in aggregated[(a, b)]]))
            else:
                row.append('0')

        out_tbl.append(row)

    if verbose:
        info('Writing output to "{}"\n'.format(output_file))

    with open(output_file, 'w') as f:
        f.write('\n'.join(['\t'.join(r) for r in out_tbl]))


def standard_phylogeny_reconstruction(project_name, configs, args, db_dna, db_aa):
    all_inputs = None
    input_faa = {}
    inp_bz2 = os.path.join(args.data_folder, 'uncompressed')
    input_fna = load_input_files(args.input_folder, inp_bz2, args.genome_extension, verbose=args.verbose)

    if input_fna:
        inp_f = os.path.join(args.data_folder, 'map_dna')
        gene_markers_identification(configs, 'map_dna', input_fna, inp_f, args.database, db_dna, args.min_num_proteins,
                                    nproc=args.nproc, verbose=args.verbose)
        gene_markers_selection(inp_f, largest_cluster if (args.db_type == 'a') and (not args.force_nucleotides) else best_hit,
                               args.min_num_proteins, nucleotides=True if args.db_type == 'a' else False,
                               nproc=args.nproc, verbose=args.verbose)
        out_f = os.path.join(args.data_folder, 'markers_dna')
        gene_markers_extraction(input_fna, inp_f, out_f, args.genome_extension, args.min_num_markers,
                                frameshifts=True if (args.db_type == 'a') and (not args.force_nucleotides) else False,
                                nproc=args.nproc, verbose=args.verbose)
        inp_f = out_f

    if (args.db_type == 'a') and (not args.force_nucleotides):
        if input_fna:
            out_f = os.path.join(args.data_folder, 'fake_proteomes')
            fake_proteome(inp_f, out_f, args.genome_extension, args.proteome_extension, args.min_len_protein,
                          nproc=args.nproc, verbose=args.verbose)
            inp_f = out_f
            input_faa = load_input_files(inp_f, inp_bz2, args.proteome_extension, verbose=args.verbose)

        faa = load_input_files(args.input_folder, inp_bz2, args.proteome_extension, verbose=args.verbose)

        if input_faa:  # if duplicates input keep the ones from 'faa'
            input_faa.update(faa)
        else:  # otherwise use only the ones from 'faa'
            input_faa = faa

        input_faa_checked = check_input_proteomes(input_faa, args.min_num_proteins, args.min_len_protein,
                                                  args.data_folder, nproc=args.nproc, verbose=args.verbose)

        if input_faa_checked:
            inp_f = os.path.join(args.data_folder, 'clean_aa')
            clean_input_proteomes(input_faa_checked, inp_f, nproc=args.nproc, verbose=args.verbose)
            input_faa_clean = load_input_files(inp_f, inp_bz2, args.proteome_extension, verbose=args.verbose)

            if input_faa_clean:
                inp_f = os.path.join(args.data_folder, 'map_aa')
                gene_markers_identification(configs, 'map_aa', input_faa_clean, inp_f, args.database, db_aa,
                                            args.min_num_proteins, nproc=args.nproc, verbose=args.verbose)
                gene_markers_selection(inp_f, best_hit, args.min_num_proteins, nucleotides=False,
                                       nproc=args.nproc, verbose=args.verbose)
                out_f = os.path.join(args.data_folder, 'markers_aa')
                gene_markers_extraction(input_faa_clean, inp_f, out_f, args.proteome_extension, args.min_num_markers,
                                        nproc=args.nproc, verbose=args.verbose)
                inp_f = out_f

    # check if inputs is empty
    if not (len(input_faa) + len(input_fna)):
        error('no inputs found, please check your params and inputs file extensions', exit=True)

    out_f = os.path.join(args.data_folder, 'markers')
    inputs2markers(inp_f, out_f, args.min_num_entries,
                   args.proteome_extension if (args.db_type == 'a') and (not args.force_nucleotides) else args.genome_extension,
                   verbose=args.verbose)
    inp_f = out_f
    out_f = os.path.join(args.data_folder, 'msas')
    msas(configs, 'msa', inp_f,
         args.proteome_extension if (args.db_type == 'a') and (not args.force_nucleotides) else args.genome_extension, out_f,
         nproc=args.nproc, verbose=args.verbose)
    inp_f = out_f

    if args.mutation_rates:
        mr_out_d = os.path.join(args.output, 'mutation_rates')
        mr_out_f = os.path.join(args.output, 'mutation_rates.tsv')
        mutation_rates(inp_f, mr_out_d, nproc=args.nproc, verbose=args.verbose)
        aggregate_mutation_rates(mr_out_d, mr_out_f, verbose=args.verbose)

    if args.trim:
        if ('trim' in configs) and ((args.trim == 'gap_trim') or (args.trim == 'greedy')):
            out_f = os.path.join(args.data_folder, 'trim_gap_trim')
            trim_gap_trim(configs, 'trim', inp_f, out_f, nproc=args.nproc, verbose=args.verbose)
            inp_f = out_f

        if (args.trim == 'gap_perc') or (args.trim == 'greedy'):
            out_f = os.path.join(args.data_folder, 'trim_gap_perc')
            trim_gap_perc(inp_f, out_f, args.gap_perc_threshold, nproc=args.nproc, verbose=args.verbose)
            inp_f = out_f

        if (args.trim == 'not_variant') or (args.trim == 'greedy'):
            out_f = os.path.join(args.data_folder, 'trim_not_variant')
            trim_not_variant(inp_f, out_f, args.not_variant_threshold, nproc=args.nproc, verbose=args.verbose)
            inp_f = out_f

    if args.remove_fragmentary_entries or args.remove_only_gaps_entries:
        out_f = os.path.join(args.data_folder, 'only_gaps')

        if args.remove_fragmentary_entries:
            out_f = os.path.join(args.data_folder, 'fragmentary')

        remove_fragmentary_entries(inp_f, args.data_folder, out_f, args.fragmentary_threshold,
                                   args.min_num_entries, nproc=args.nproc, verbose=args.verbose)
        inp_f = out_f

    # compute inputs list
    all_inputs = inputs_list(inp_f, '.aln', os.path.join(args.data_folder, project_name + '_input_list.pkl'),
                             nproc=args.nproc, verbose=args.verbose)

    if args.subsample:
        out_f = os.path.join(args.data_folder, 'sub')
        subsample(inp_f, out_f, args.subsample, args.scoring_function, os.path.join(args.submat_folder, args.submat + '.pkl'),
                  unknown_fraction=args.unknown_fraction, nproc=args.nproc, verbose=args.verbose)
        inp_f = out_f

    if 'gene_tree1' in configs:
        sub_mod = load_substitution_model(args.maas)
        out_f = os.path.join(args.data_folder, 'gene_tree1')
        build_gene_tree(configs, 'gene_tree1', sub_mod, inp_f, out_f, nproc=args.nproc, verbose=args.verbose)

        if 'gene_tree2' in configs:
            outt = os.path.join(args.data_folder, 'gene_tree1_polytomies')
            resolve_polytomies(out_f, outt, nproc=args.nproc, verbose=args.verbose)
            out_f = outt

            outt = os.path.join(args.data_folder, 'gene_tree2')
            refine_gene_tree(configs, 'gene_tree2', sub_mod, inp_f, out_f, outt, nproc=args.nproc, verbose=args.verbose)
            out_f = outt

        inp_f = out_f
        out_f = os.path.join(args.data_folder, 'gene_trees.tre')
        merging_gene_trees(inp_f, out_f, verbose=args.verbose)
        inp_f = out_f
    else:
        if not all_inputs:
            all_inputs = (os.path.splitext(os.path.basename(i))[0] for i in input_faa_clean)

        out_f = os.path.join(args.output, project_name + '_concatenated.aln')
        concatenate(all_inputs, inp_f, out_f, sort=args.sort, verbose=args.verbose)
        inp_f = out_f

    out_f = project_name + '.tre'
    build_phylogeny(configs, 'tree1', inp_f, os.path.abspath(args.output), out_f, nproc=args.nproc, verbose=args.verbose)

    if 'tree2' in configs:
        outt = project_name + '_resolved.tre'
        resolve_polytomies(os.path.join(args.output, out_f), os.path.join(args.output, outt), nproc=args.nproc, verbose=args.verbose)
        out_f = os.path.join(args.output, outt)

        refine_phylogeny(configs, 'tree2', inp_f, out_f, os.path.abspath(args.output), project_name + '_refined.tre',
                         nproc=args.nproc, verbose=args.verbose)


def byte_to_megabyte(byte):
    """
    Convert byte value to megabyte
    """

    return (byte / 1048576)


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
                info("Downloading file of size: {:.2f} MB\n".format(byte_to_megabyte(total_size)))
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
                           .format(percent_downloaded, byte_to_megabyte(download_rate), estimated_minutes, estimated_seconds))

            status += "        \r"
            info(status)


def download(url, download_file, verbose=False):
    """
    Download a file from a url
    """

    if not os.path.isfile(download_file):
        try:
            if verbose:
                info('Downloading "{}"\n'.format(url))

            urlretrieve(url, download_file, reporthook=ReportHook().report)
            info('\n')
        except EnvironmentError:
            error('unable to download "{}"'.format(url))
    elif verbose:
        info('File "{}" present\n'.format(download_file))


def download_and_unpack_db(url, db_name, folder, verbose=False):
    """
    Download the url to the file and decompress into the folder
    """

    if not os.path.isdir(folder):  # create the folder if it does not exist
        try:
            if verbose:
                info('Creating "{}" folder\n'.format(folder))

            os.makedirs(folder)
        except EnvironmentError:
            error('unable to create database folder "{}"'.format(folder), exit=True)

    if os.path.isdir(os.path.join(folder, db_name)):  # check if there already exists the folder
        if verbose:
            info('Database folder "{}" present\n'.format(os.path.join(folder, db_name)))

        return

    # check the directory permissions
    if not os.access(folder, os.W_OK):
        error('database directory "{}" is not writeable, please modify the permissions'.format(folder), exit=True)

    # download database
    tar_file = os.path.join(folder, db_name + ".tar")
    url_tar_file = os.path.join(url, db_name + ".tar")
    download(url_tar_file, tar_file, verbose=verbose)

    # download MD5 checksum
    md5_file = os.path.join(folder, db_name + ".md5")
    url_md5_file = os.path.join(url, db_name + ".md5")
    download(url_md5_file, md5_file, verbose=verbose)

    md5_md5 = None
    md5_tar = None

    if os.path.isfile(md5_file):
        with open(md5_file) as f:
            for row in f:
                md5_md5 = row.strip().split(' ')[0]
    else:
        error('file "{}" not found!'.format(md5_file))

    # compute MD5 of .tar.bz2
    if os.path.isfile(tar_file):
        hash_md5 = hashlib.md5()

        with open(tar_file, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)

        md5_tar = hash_md5.hexdigest()[:32]
    else:
        error('file "{}" not found!'.format(tar_file), exit=True)

    if (md5_tar is None) or (md5_md5 is None):
        error("MD5 checksums not found, something went wrong!", exit=True)

    # compare checksums
    if md5_tar != md5_md5:
        error('MD5 checksums do not correspond, if this happens again, please remove the database files '
              'in "{}" and re-run PhyloPhlAn2 so they will be re-downloaded'.format(folder), exit=True)

    # untar
    if verbose:
        info('Extracting "{}" into "{}"\n'.format(tar_file, folder))

    try:
        tarfile_handle = tarfile.open(tar_file)
        tarfile_handle.extractall(path=folder)
        tarfile_handle.close()
    except EnvironmentError:
        error('unable to extract "{}"'.format(tar_file))


def phylophlan2():
    args = read_params()

    if args.verbose:
        info('PhyloPhlAn2 version {} ({})\n'.format(__version__, __date__))
        info('Command line: {}\n\n'.format(' '.join(sys.argv)), init_new_line=True)

    project_name = check_args(args, sys.argv, verbose=args.verbose)

    if args.clean_all:
        clean_all(args.databases_folder, verbose=args.verbose)

    if args.clean:
        clean_project(args.data_folder, args.output, verbose=args.verbose)

    configs = read_configs(args.config_file, verbose=args.verbose)
    check_configs(configs, verbose=args.verbose)
    check_dependencies(configs, args.nproc, verbose=args.verbose)
    download_and_unpack_db(DATABASE_DOWNLOAD_URL, args.database, args.databases_folder, verbose=args.verbose)
    check_database(args.database, args.databases_folder, verbose=args.verbose)
    db_type, db_dna, db_aa = init_database(args.database, args.databases_folder, args.db_type,
                                           configs, 'db_dna', 'db_aa', verbose=args.verbose)
    if not args.db_type:
        args.db_type = db_type

    standard_phylogeny_reconstruction(project_name, configs, args, db_dna, db_aa)


if __name__ == '__main__':
    t0 = time.time()
    phylophlan2()
    t1 = time.time()
    info('\nTotal elapsed time {}s\n'.format(int(t1 - t0)))
    sys.exit(0)
