__author__ = ('Aitor Blanco Miguez (aitor.blancomiguez@unitn.it), '
              'Duy Tin Truong (duytin.truong@unitn.it), '
              'Francesco Asnicar (f.asnicar@unitn.it), '
              'Moreno Zolfo (moreno.zolfo@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)')
__version__ = '4.0.3'
__date__ = '24 Oct 2022'


import os
import re
import shutil
import tempfile
import subprocess as sb
try:
    from .util_fun import info, error
except ImportError:
    from util_fun import info, error


def execute(cmd):
    """Runs a command line executable

    Args:
        cmd (dict): the dict with the command line information
    """
    inp_f = None
    out_f = sb.DEVNULL
    if cmd['stdin']:
        inp_f = open(cmd['stdin'], 'r')
    if cmd['stdout']:
        out_f = open(cmd['stdout'], 'w')
    exec_res = sb.run(cmd['command_line'], stdin=inp_f, stdout=out_f)
    if exec_res.returncode == 1:
        error("An error was ocurred executing a external tool, exiting...", exit=True)
    if cmd['stdin']:
        inp_f.close()
    if cmd['stdout']:
        out_f.close()


def create_blastn_db(output_dir, reference):
    """Creates the BLASTn database to align the reference genomes"""
    reference_name = os.path.splitext(os.path.basename(reference))[0]
    params = {
        "program_name": "makeblastdb",
        "params": "-parse_seqids -dbtype nucl",
        "input": "-in",
        "output": "-out",
        "command_line": "#program_name# #params# #input# #output#"
    }
    execute(compose_command(params, input_file=reference,
            output_file=os.path.join(output_dir, reference_name)))
    return os.path.join(output_dir, reference_name)


def execute_blastn(output_dir, clade_markers, blastn_db, nprocs=1):
    """Executes BLASTn"""
    db_name = os.path.splitext(os.path.basename(blastn_db))[0]
    params = {
        "program_name": "blastn",
        "params": "-outfmt \"6 qseqid sseqid qlen qstart qend sstart send\" -evalue 1e-10 -max_target_seqs 1",
        "input": "-query",
        "database": "-db",
        "output": "-out",
        "threads": "-num_threads",
        "command_line": "#program_name# #params# #threads# #database# #input# #output#"
    }
    execute(compose_command(params, input_file=clade_markers, database=blastn_db,
                            output_file=os.path.join(output_dir, "{}.blastn".format(db_name)), nproc=nprocs))
    return os.path.join(output_dir, "{}.blastn".format(db_name))


def create_phylophlan_db(output_dir, clade):
    """Creates the PhyloPhlAn database"""
    markers = os.path.join(output_dir, clade)
    params = {
        "program_name": "phylophlan_setup_database",
        "params": "-d {} -e fna -t n --overwrite".format(clade),
        "input": "-i",
        "command_line": "#program_name# #input# #params#"
    }
    execute(compose_command(params, input_file=markers))


def generate_phylophlan_config_file(output_dir, configuration):
    """Generates the PhyloPhlan configuration file"""
    conf_file = os.path.join(output_dir, "phylophlan.cfg")
    params = {
        "program_name": "phylophlan_write_config_file",
        "params": "-d n --db_dna makeblastdb --map_dna {} --msa {} --trim {} --tree1 {}".format(configuration['map'], configuration['aligner'], configuration['trim'], configuration['tree1']),
        "output": "-o",
        "command_line": "#program_name# #output# #params#"
    }
    execute(compose_command(params, output_file=conf_file))
    return conf_file


def execute_phylophlan(samples_markers_dir, conf_file, min_entries, min_markers, tmp_dir, output_dir, clade, phylogeny_conf, mutation_rates, nprocs):
    """Executes PhyloPhlAn"""
    advanced_params = "--"+phylogeny_conf
    if mutation_rates:
        advanced_params = advanced_params + " --mutation_rates"
    params = {
        "program_name": "phylophlan",
        "params": "-d {} --data_folder {} --databases_folder {} -t n -f {} --diversity low {} --genome_extension fna --force_nucleotides --min_num_entries {} --convert_N2gap --min_num_markers {}".format(clade[:30], tmp_dir, tmp_dir, conf_file, advanced_params, min_entries, min_markers),
        "input": "-i",
        "output_path": "--output_folder",
        "output": "-o",
        "threads": "--nproc",
        "command_line": "#program_name# #input# #output# #output_path# #params# #threads#"
    }
    execute(compose_command(params=params, input_file=samples_markers_dir, output_path=output_dir,
                            output_file=".", nproc=nprocs))


def execute_treeshrink(input_tree, output_dir, tmp=None, centroid=False, q_value=0.05, m_value='all-genes'):
    """Executes TreeShrink"""
    advanced_string = ''
    if tmp is not None:
        advanced_string += '-p {} '.format(tmp)
    if centroid:
        advanced_string += '-c '
    if not os.path.exists('{}/treeshrink'.format(output_dir)):
        os.mkdir('{}/treeshrink'.format(output_dir))
    params = {
        "program_name": "treeshrink.py",
        "params": "{}-q {} -m {} -f".format(advanced_string, q_value, m_value),
        "input": "-t",
        "output_path": '-o',
        "output": "-O",
        "command_line": "#program_name# #input# #output_path# #output# #params#"
    }
    execute(compose_command(params=params, input_file=input_tree, output_path='{}/treeshrink'.format(
        output_dir), output_file=input_tree.split('/')[-1].replace('.StrainPhlAn4.tre', '.TreeShrink')))


def decompress_bz2(input_file, output_dir):
    """Decompressed BZ2 files"""
    n, _ = os.path.splitext(os.path.basename(input_file))
    params = {
        "program_name": "bzip2",
        "input": "-cdk",
        "command_line": "#program_name# #input# > #output#"
    }
    decompressed_file = os.path.join(output_dir, n)
    execute(compose_command(params, input_file=input_file,
            output_file=decompressed_file))
    if decompressed_file.endswith('_sam'):
        os.rename(decompressed_file, decompressed_file[:-4] + '.sam')
        decompressed_file = decompressed_file[:-4] + '.sam'
    return decompressed_file


def samtools_sam_to_bam(input_file, output_dir):
    """Converts SAM files to sorted BAM files using samtools"""
    n, _ = os.path.splitext(os.path.basename(input_file))
    params = {
        "program_name": "samtools",
        "params": "view",
        "input": "-Sb",
        "command_line": "#program_name# #params# #input# > #output#"
    }
    execute(compose_command(params, input_file=input_file,
            output_file=os.path.join(output_dir, "{}.bam".format(n))))
    return os.path.join(output_dir, "{}.bam".format(n))


def samtools_sort_bam_v1(input_file, output_dir):
    """Sort BAM files using samtools"""
    n, _ = os.path.splitext(os.path.basename(input_file))
    params = {
        "program_name": "samtools",
        "params": "sort",
        "output": "-o",
        "command_line": "#program_name# #params# #input# #output#"
    }
    execute(compose_command(params, input_file=input_file,
            output_file=os.path.join(output_dir, "{}.sorted.bam".format(n))))
    shutil.move(os.path.join(output_dir, "{}.sorted.bam".format(n)),
                os.path.join(output_dir, "{}.bam".format(n)))
    return os.path.join(output_dir, "{}.bam".format(n))


def generate_markers_fasta(database, output_dir):
    """Generates a FASTA file with the markers form a MetaPhlAn database"""
    file_out, db_markers_faa = tempfile.mkstemp(
        dir=output_dir, prefix="db_markers_temp_", suffix=".fna")
    os.close(file_out)
    bowtie_database, _ = os.path.splitext(database)
    params = {
        "program_name": "bowtie2-inspect",
        "command_line": "#program_name# #input# > #output#"
    }
    execute(compose_command(params, input_file=bowtie_database,
            output_file=db_markers_faa))
    return db_markers_faa


def compose_command(params, check=False, input_file=None, database=None, output_path=None, output_file=None, nproc=1):
    """Compose a command for further executions. Copied from the PhyloPhlAn 3 code"""
    program_name = None
    stdin = None
    stdout = None
    environment = os.environ.copy()
    r_output_path = None
    r_output_file = None
    command_line = params['command_line']

    if 'program_name' in list(params):
        command_line = command_line.replace(
            '#program_name#', params['program_name'])
        program_name = params['program_name']
    else:
        error('Error: something wrong... ' +
              program_name+' not found!', exit=True)

    if check:
        command_line = program_name

        if 'version' in list(params):
            command_line = '{} {}'.format(program_name, params['version'])
    else:
        if 'params' in list(params):
            command_line = command_line.replace('#params#', params['params'])

        if 'threads' in list(params):
            command_line = command_line.replace(
                '#threads#', '{} {}'.format(params['threads'], nproc))

        if output_path:
            r_output_path = output_path

            if 'output_path' in list(params):
                command_line = command_line.replace(
                    '#output_path#', '{} {}'.format(params['output_path'], output_path))
            else:
                output_file = os.path.join(output_path, output_file)

        if input_file:
            inp = input_file

            if 'input' in list(params):
                inp = '{} {}'.format(params['input'], input_file)

            if '<' in command_line:
                command_line = command_line.replace('<', '')
                command_line = command_line.replace('#input#', '')
                stdin = inp
            else:
                command_line = command_line.replace('#input#', inp)

        if database and ('database' in list(params)):
            command_line = command_line.replace(
                '#database#', '{} {}'.format(params['database'], database))

        if output_file:
            out = output_file
            r_output_file = output_file

            if 'output' in list(params):
                out = '{} {}'.format(params['output'], output_file)

            if '>' in command_line:
                command_line = command_line.replace('>', '')
                command_line = command_line.replace('#output#', '')
                stdout = out
            else:
                command_line = command_line.replace('#output#', out)

        if 'environment' in list(params):
            new_environment = dict([(var.strip(), val.strip())
                                    for var, val in [a.strip().split('=') for a in params['environment'].split(',')]])
            environment.update(new_environment)

    # find string sourrunded with " and make them as one string
    quotes = [j for j, e in enumerate(command_line) if e == '"']

    for s, e in zip(quotes[0::2], quotes[1::2]):
        command_line = command_line.replace(
            command_line[s + 1:e], command_line[s + 1:e].replace(' ', '#'))

    return {'command_line': [str(a).replace('#', ' ') for a in re.sub(' +', ' ', command_line.replace('"', '')).split(' ') if a],
            'stdin': stdin, 'stdout': stdout, 'env': environment, 'output_path': r_output_path, 'output_file': r_output_file}
