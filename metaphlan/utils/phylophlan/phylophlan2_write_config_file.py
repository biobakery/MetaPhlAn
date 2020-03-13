#!/usr/bin/env python3


__author__ = ('Francesco Asnicar (f.asnicar@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it), '
              'Claudia Mengoni (claudia.mengoni@studenti.unitn.it),'
              'Mattia Bolzan (mattia.bolzan@unitn.it), '
              'Nicola Segata (nicola.segata@unitn.it)')
__version__ = '0.11'
__date__ = '9 September 2019'


import os
import sys
import stat
import time
import argparse as ap
import configparser as cp


if sys.version_info[0] < 3:
    raise Exception("PhyloPhlAn2 requires Python 3, your current Python version is {}.{}.{}"
                    .format(sys.version_info[0], sys.version_info[1], sys.version_info[2]))

DB_TYPE_CHOICES = ['n', 'a']
DB_DNA_CHOICES = ['makeblastdb']
DB_AA_CHOICES = ['usearch', 'diamond']
MAP_DNA_CHOICES = ['blastn', 'tblastn', 'diamond']
MAP_AA_CHOICES = ['usearch', 'diamond']
MSA_CHOICES = ['muscle', 'mafft', 'opal', 'upp']
TRIM_CHOICES = ['trimal']
GENE_TREE1_CHOICES = ['fasttree', 'raxml', 'iqtree']
GENE_TREE2_CHOICES = ['raxml']
TREE1_CHOICES = ['fasttree', 'raxml', 'iqtree', 'astral', 'astrid']
TREE2_CHOICES = ['raxml']


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
    p = ap.ArgumentParser(description=("The phylophlan2_write_config_file.py script generates a configuration file to be used with the "
                                       "phylophlan2.py script. It implements some standard parameters for the software integrated, but if "
                                       "needed, the parameters of the selected software can be added/modified/removed by editing the "
                                       "generated configuration file using a text editor"),
                          formatter_class=ap.ArgumentDefaultsHelpFormatter)

    p.add_argument('-o', '--output', type=str, required=True,
                   help="Specify the output file where to write the configurations")
    p.add_argument('-d', '--db_type', required=True, choices=DB_TYPE_CHOICES,
                   help='Specify the type of the database, where "n" stands for nucleotides and "a" for amino acids')

    p_db = p.add_mutually_exclusive_group(required=True)
    p_db.add_argument('--db_dna', default=None, choices=DB_DNA_CHOICES,
                      help='Add the "db_dna" section of the selected software that will be used for building the indexed database')
    p_db.add_argument('--db_aa', default=None, choices=DB_AA_CHOICES,
                      help='Add the "db_aa" section of the selected software that will be used for building the indexed database')

    p.add_argument('--map_dna', default=None, choices=MAP_DNA_CHOICES,
                   help='Add the "map_dna" section of the selected software that will be used for mapping the database against '
                        'the input genomes')
    p.add_argument('--map_aa', default=None, choices=MAP_AA_CHOICES,
                   help='Add the "map_aa" section of the selected software that will be used for mapping the database against '
                        'the input proteomes')
    p.add_argument('--msa', required=True, default=None, choices=MSA_CHOICES,
                   help='Add the "msa" section of the selected software that will be used for producing the MSAs')
    p.add_argument('--trim', default=None, choices=TRIM_CHOICES,
                   help='Add the "trim" section of the selected software that will be used for the gappy regions removal of the MSAs')
    p.add_argument('--gene_tree1', default=None, choices=GENE_TREE1_CHOICES,
                   help='Add the "gene_tree1" section of the selected software that will be used for building the phylogenies for '
                        'the markers in the database')
    p.add_argument('--gene_tree2', default=None, choices=GENE_TREE2_CHOICES,
                   help='Add the "gene_tree2" section of the selected software that will be used for refining '
                        'the phylogenies previously built with what specified in the "gene_tree1" section')
    p.add_argument('--tree1', required=True, default=None, choices=TREE1_CHOICES,
                   help='Add the "tree1" section of the selected software that will be used for building the first phylogeny')
    p.add_argument('--tree2', default=None, choices=TREE2_CHOICES,
                   help=('Add the "tree2" section of the selected software that will be used for refining the '
                         'phylogeny previously built with what specified in the "tree1" section'))
    p.add_argument('-a', '--absolute_path', action='store_true', default=False,
                   help="Write the absolute path to the executable instead of the executable name as found in the system path environment")
    p.add_argument('--force_nucleotides', default=None, action='store_true',
                   help=('If specified sets parameters for phylogenetic analysis software so that they use '
                         'nucleotide sequences, even in the case of a database of amino acids'))
    p.add_argument('--overwrite', action='store_true', default=False, help="Overwrite output file if it exists")
    p.add_argument('--verbose', action='store_true', default=False, help="Prints more stuff")
    p.add_argument('-v', '--version', action='version',
                   version='phylophlan2_write_config_file.py version {} ({})'.format(__version__, __date__),
                   help="Prints the current phylophlan2_write_config_file.py version and exit")
    return p.parse_args()


def check_params(args, verbose=False):
    if (not args.map_dna) and (not args.map_aa):
        error('at least one of --map_dna and --map_aa must be specified', exit=True)

    if (os.path.isfile(args.output)) and (not args.overwrite):
        error('output file {} already exists, delete it or specify --overwrite'.format(args.output), exit=True)

    if (args.db_type == 'n') and (args.map_dna == 'diamond'):
        error('Incompatible choices.\n    Please, select "blastn" instead in --map_dna, as there is no "diamond" '
              'option for aligning against a nucleotide reference database', exit=True)

    if verbose:
        info('Arguments: {}\n'.format(vars(args)))


def find_executable(exe, absolute=False, rollback=False):
    current_env = os.environ.copy()

    if 'PATH' in current_env:
        for d in current_env['PATH'].split(os.pathsep):
            if not os.path.exists(d):
                continue

            if not os.path.isdir(d):
                continue

            for f in os.listdir(d):
                e = os.path.realpath(os.path.join(d, f))

                if (os.path.exists(e) and (not os.path.isdir(e)) and bool(os.stat(e)[stat.ST_MODE] & stat.S_IXUSR)):
                    if exe.lower() == f.lower():
                        return (f if not absolute else e,
                                True if rollback is None else False)

    if rollback:
        find_executable(rollback, absolute=absolute, rollback=None)

    error('could not find "{}" ("{}") executable in your PATH environment variable'.format(exe, rollback), exit=True)


# AVAILABLE OPTIONS:
# program_name: name of the executable to use
# params: params to use
# threads: specify the option to use to pass the number of threads
# input: specify the option to use for the input file
# database: specify the option to use for setting the database
# output_path: specify the to option to use to set the path of the folder that will
#              contains the output file
# output: specify the option to use for the output file
# version: specify the option to use to get the version of the software, used to verify
#          the software installation
# command_line: specify the command line to generate with the position of each argument,
#               '<' and '>' can be used to specify input/output redirection, respectively
# environment: specify variables and their values to be defined in the environment,
#              syntax VARIABLE1=VALUE1,VARIABLE2=VALUE2,...


def phylophlan2_write_config_file():
    args = read_params()

    if args.verbose:
        info('phylophlan2_write_config_file.py version {} ({})\n'.format(__version__, __date__))
        info('Command line: {}\n\n'.format(' '.join(sys.argv)), init_new_line=True)

    check_params(args, verbose=args.verbose)
    progs = {}

    # setting the program for building the nucleotides DB
    if args.db_dna:
        if 'makeblastdb' in args.db_dna:
            exe, _ = find_executable('makeblastdb', absolute=args.absolute_path)
            db_dna = {'program_name': exe,
                      'params': '-parse_seqids -dbtype nucl',
                      'input': '-in',
                      'output': '-out',
                      'version': '-version',
                      'command_line': '#program_name# #params# #input# #output#'}

        progs['db_dna'] = db_dna

    # setting the program for building the amino acids DB
    if args.db_aa:
        if 'usearch' in args.db_aa:
            exe, _ = find_executable('usearch', absolute=args.absolute_path)
            db_aa = {'program_name': exe,
                     'params': '-quiet',
                     'input': '-makeudb_ublast',
                     'output': '-output',
                     'version': '-version',
                     'command_line': '#program_name# #params# #input# #output#'}
        elif 'diamond' in args.db_aa:
            exe, _ = find_executable('diamond', absolute=args.absolute_path)
            db_aa = {'program_name': exe,
                     'params': 'makedb',
                     'threads': '--threads',
                     'input': '--in',
                     'output': '--db',
                     'version': 'version',
                     'command_line': '#program_name# #params# #threads# #input# #output#'}

        progs['db_aa'] = db_aa

    # setting the software for mapping the genomes (dna)
    if args.map_dna:
        if 'blastn' in args.map_dna:
            exe, _ = find_executable('blastn', absolute=args.absolute_path)
            map_dna = {'program_name': exe,
                       'params': '-outfmt 6 -max_target_seqs 1000000',
                       'input': '-query',
                       'database': '-db',
                       'output': '-out',
                       'version': '-version',
                       'command_line': '#program_name# #params# #input# #database# #output#'}
        elif 'tblastn' in args.map_dna:
            exe, _ = find_executable('tblastn', absolute=args.absolute_path)
            map_dna = {'program_name': exe,
                       'params': ('-outfmt "6 saccver qaccver pident length mismatch gapopen sstart send qstart '
                                  'qend evalue bitscore" -evalue 1e-50 -max_target_seqs 1000000'),
                       'input': '-subject',
                       'database': '-query',
                       'output': '-out',
                       'version': '-version',
                       'command_line': '#program_name# #params# #input# #database# #output#'}
        elif 'diamond' in args.map_dna:
            exe, _ = find_executable('diamond', absolute=args.absolute_path)
            map_dna = {'program_name': exe,
                       'params': 'blastx --quiet --threads 1 --outfmt 6 --more-sensitive --id 50 --max-hsps 35 -k 0',
                       'input': '--query',
                       'database': '--db',
                       'output': '--out',
                       'version': 'version',
                       'command_line': '#program_name# #params# #input# #database# #output#'}

        progs['map_dna'] = map_dna

    # setting the software for mapping the proteomes (aa)
    if args.map_aa:
        if 'usearch' in args.map_aa:
            exe, _ = find_executable('usearch', absolute=args.absolute_path)
            map_aa = {'program_name': exe,
                      'params': '-quiet -evalue 1e-10 -maxaccepts 8 -maxrejects 32',
                      'threads': '-threads',
                      'input': '-ublast',
                      'database': '-db',
                      'output': '-blast6out',
                      'version': '-version',
                      'command_line': '#program_name# #params# #threads# #input# #database# #output#'}
        elif 'diamond' in args.map_aa:
            exe, _ = find_executable('diamond', absolute=args.absolute_path)
            map_aa = {'program_name': exe,
                      'params': 'blastp --quiet --threads 1 --outfmt 6 --more-sensitive --id 50 --max-hsps 35 -k 0',
                      'input': '--query',
                      'database': '--db',
                      'output': '--out',
                      'version': 'version',
                      'command_line': '#program_name# #params# #input# #database# #output#'}

        progs['map_aa'] = map_aa

    # setting the MSA software
    if 'muscle' in args.msa:
        exe, _ = find_executable('muscle', absolute=args.absolute_path)
        msa = {'program_name': exe,
               'params': '-quiet -maxiters 2',
               'input': '-in',
               'output': '-out',
               'version': '-version',
               'command_line': '#program_name# #params# #input# #output#'}
    elif 'mafft' in args.msa:
        exe, _ = find_executable('mafft', absolute=args.absolute_path)
        msa = {'program_name': exe,
               'params': '--quiet --anysymbol --thread 1 --auto',
               'version': '--version',
               'command_line': '#program_name# #params# #input# > #output#'}

        for fld in ['/local-storage', '/tmp']:
            if os.path.isdir(fld) and bool(os.stat(fld)[stat.ST_MODE] & stat.S_IWUSR):
                if 'environment' in msa:
                    break

                msa['environment'] = 'TMPDIR={}'.format(fld)
                break
    elif 'opal' in args.msa:
        exe, _ = find_executable('opal', absolute=args.absolute_path)
        msa = {'program_name': exe,
               'input': '--in',
               'output': '--out',
               'params': '--quiet',
               'command_line': '#program_name# #params# #input# #output#'}

        if (args.db_type == 'a' and (not args.force_nucleotides)):
            msa['params'] += ' --protein'
    elif 'upp' in args.msa:
        exe, _ = find_executable('upp', rollback='run-upp.sh', absolute=args.absolute_path)
        msa = {'program_name': exe,
               'params': '-x 1 -M -1 -T 0.66 -B 999999999',
               'input': '-s',
               'output': '-o',
               'output_path': '-d',
               'version': '--version',
               'command_line': '#program_name# #params# #input# #output_path# #output#'}

        if (args.db_type == 'n' or args.force_nucleotides):
            msa['params'] += ' -m dna',
        elif args.db_type == 'a':
            msa['model'] += ' -m amino'

    progs['msa'] = msa

    # setting the trimming software
    if args.trim:
        if 'trimal' in args.trim:
            exe, _ = find_executable('trimal', absolute=args.absolute_path)
            trim = {'program_name': exe,
                    'params': '-gappyout',
                    'input': '-in',
                    'output': '-out',
                    'version': '--version',
                    'command_line': '#program_name# #params# #input# #output#'}

        progs['trim'] = trim

    # setting gene_tree1
    if args.gene_tree1:
        if 'fasttree' in args.gene_tree1:
            exe, _ = find_executable('FastTree', rollback='fasttree', absolute=args.absolute_path)
            gene_tree1 = {'program_name': exe,
                          'params': '-quiet -mlacc 2 -slownni -spr 4 -fastest -mlnni 4 -no2nd',
                          'output': '-out',
                          'command_line': '#program_name# #params# #output# #input#'}

            if (args.db_type == 'n' or args.force_nucleotides):
                gene_tree1['params'] += ' -gtr -nt'
            elif args.db_type == 'a':
                gene_tree1['params'] += ' -lg'
        elif 'raxml' in args.gene_tree1:
            exe, _ = find_executable('raxmlHPC', rollback='raxml', absolute=args.absolute_path)
            gene_tree1 = {'program_name': exe,
                          'params': '-p 1989',
                          'input': '-s',
                          'output_path': '-w',
                          'output': '-n',
                          'version': '-v'}

            if (args.db_type == 'n' or args.force_nucleotides):
                gene_tree1['params'] += ' -m GTRCAT',
                gene_tree1['command_line'] = '#program_name# #params# #output_path# #input# #output#'
            elif args.db_type == 'a':
                gene_tree1['model'] = '-m'
                gene_tree1['command_line'] = '#program_name# #model# #params# #output_path# #input# #output#'
        elif 'iqtree' in args.gene_tree1:
            exe, _ = find_executable('iqtree', absolute=args.absolute_path)
            gene_tree1 = {'program_name': exe,
                          'params': '-quiet -nt 1',
                          'input': '-s',
                          'output': '-pre',
                          'version': '-version',
                          'command_line': '#program_name# #params# #input# #output#'}

            if (args.db_type == 'n' or args.force_nucleotides):
                gene_tree1['params'] += ' -m GTR',
            elif args.db_type == 'a':
                gene_tree1['params'] = ' -m LG'

        progs['gene_tree1'] = gene_tree1

    # setting gene_tree2
    if args.gene_tree2:
        if 'raxml' in args.gene_tree2:
            exe, _ = find_executable('raxmlHPC', rollback='raxml', absolute=args.absolute_path)
            gene_tree2 = {'program_name': exe,
                          'params': '-p 1989',
                          'database': '-t',  # starting tree
                          'input': '-s',
                          'output_path': '-w',
                          'output': '-n',
                          'version': '-v'}

            if (args.db_type == 'n' or args.force_nucleotides):
                gene_tree2['params'] += '-p 1989 -m GTRCAT'
                gene_tree2['command_line'] = '#program_name# #params# #database# #output_path# #input# #output#'
            elif args.db_type == 'a':
                gene_tree2['model'] = '-m'
                gene_tree2['command_line'] = '#program_name# #model# #params# #database# #output_path# #input# #output#'

        progs['gene_tree2'] = gene_tree2

    # setting tree1
    if 'astral' in args.tree1:
        tree1 = {'program_name': 'java -jar /CM/tools/astral-4.11.1/astral.4.11.1.jar',
                 'input': '-i',
                 'output': '-o',
                 # 'version': '--help',
                 'version': '-i /CM/tools/astral-4.11.1/test_data/song_mammals.424.gene.tre',
                 'command_line': '#program_name# #input# #output#'}
    if 'astrid' in args.tree1:
        exe, _ = find_executable('ASTRID', absolute=args.absolute_path)
        tree1 = {'program_name': exe,
                 'input': '-i',
                 'params': '-m auto',
                 'output': '-o',
                 'version': '--help',
                 'command_line': '#program_name# #input# #params# #output#'}
    elif 'fasttree' in args.tree1:
        exe, rb = find_executable('FastTreeMP', rollback='fasttree', absolute=args.absolute_path)
        tree1 = {'program_name': exe,
                 'params': '-quiet -mlacc 2 -slownni -spr 4 -fastest -mlnni 4 -no2nd',
                 'output': '-out',
                 'command_line': '#program_name# #params# #output# #input#'}

        if not rb:
            tree1['environment'] = 'OMP_NUM_THREADS=3'

        if (args.db_type == 'n' or args.force_nucleotides):
            tree1['params'] += ' -gtr -nt'
        elif args.db_type == 'a':
            tree1['params'] += ' -lg'
    elif 'raxml' in args.tree1:
        exe, rb = find_executable('raxmlHPC-PTHREADS-SSE3', rollback='raxml', absolute=args.absolute_path)
        tree1 = {'program_name': exe,
                 'params': '-p 1989',
                 'input': '-s',
                 'output_path': '-w',
                 'output': '-n',
                 'version': '-v',
                 'command_line': '#program_name# #params# #threads# #database# #output_path# #input# #output#'}

        if not rb:
            tree1['threads'] = '-T'

        if (args.db_type == 'n' or args.force_nucleotides):
            tree1['params'] += ' -m GTRCAT'
        elif args.db_type == 'a':
            tree1['params'] += ' -m PROTCATLG'

    elif 'iqtree' in args.tree1:
        exe, _ = find_executable('iqtree', absolute=args.absolute_path)
        tree1 = {'program_name': exe,
                 'params': '-quiet -nt AUTO',
                 'input': '-s',
                 'output': '-pre',
                 'version': '-version',
                 'command_line': '#program_name# #params# #input# #output#'}

        if (args.db_type == 'n' or args.force_nucleotides):
            tree1['params'] += ' -m GTR'
        elif args.db_type == 'a':
            tree1['params'] += ' -m LG'

    progs['tree1'] = tree1

    # setting tree2
    if args.tree2:
        if 'raxml' in args.tree2:
            exe, rb = find_executable('raxmlHPC-PTHREADS-SSE3', rollback='raxml', absolute=args.absolute_path)
            tree2 = {'program_name': exe,
                     'params': '-p 1989',
                     'database': '-t',  # starting tree
                     'input': '-s',
                     'output_path': '-w',
                     'output': '-n',
                     'version': '-v',
                     'command_line': '#program_name# #params# #threads# #database# #output_path# #input# #output#'}

            if not rb:
                tree2['threads'] = '-T'

            if (args.db_type == 'n' or args.force_nucleotides):
                tree2['params'] += ' -m GTRCAT'
            elif args.db_type == 'a':
                tree2['params'] += ' -m PROTCATLG'

        progs['tree2'] = tree2

    config = cp.ConfigParser()

    for prog, options in progs.items():
        config[prog] = {}

        for option, value in options.items():
            config[prog][option] = value

    if (os.path.isfile(args.output)) and args.overwrite and args.verbose:
        info('Output file "{}" will be overwritten\n'.format(args.output))

    with open(args.output, 'w') as f:
        config.write(f)


if __name__ == '__main__':
    t0 = time.time()
    phylophlan2_write_config_file()
    t1 = time.time()
    info('Total elapsed time {}s\n'.format(int(t1 - t0)), init_new_line=True)
    sys.exit(0)
