__author__ = ('Aitor Blanco Miguez (aitor.blancomiguez@unitn.it), '
              'Duy Tin Truong (duytin.truong@unitn.it), '
              'Francesco Asnicar (f.asnicar@unitn.it), '
              'Moreno Zolfo (moreno.zolfo@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)')
__version__ = '4.1.1'
__date__ = '11 Mar 2024'


import os
import re
import shlex
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
    out_f = sb.PIPE
    if cmd['stdin']:
        inp_f = open(cmd['stdin'], 'r')
    if cmd['stdout']:
        out_f = open(cmd['stdout'], 'w')
    exec_res = sb.run(cmd['command_line'], stdin=inp_f, stdout=out_f)
    if exec_res.returncode == 1:
        error("An error was ocurred executing a external tool, exiting...", exit=True)
        print(exec_res.stdout)
        print('===')
        print(exec_res.stderr)
    if cmd['stdin']:
        inp_f.close()
    if cmd['stdout']:
        out_f.close()


def build_bowtie2_db(input_fasta, output_database):
    params = {
        "program_name": "bowtie2-build",
        "params": "--large-index {} {}".format(input_fasta, output_database),
        "command_line": "#program_name# #params#"
    }
    execute(compose_command(params, input_file=input_fasta, output_file=output_database))


def generate_phylophlan_config_file(output_dir, configuration):
    """Generates the PhyloPhlan configuration file"""
    conf_file = os.path.join(output_dir, "phylophlan.cfg")
    params = {
        "program_name": "phylophlan_write_config_file",
        "params": "-d n --db_dna makeblastdb --map_dna {} --msa {} --tree1 {}".format(configuration['map'], configuration['aligner'], configuration['tree1']),
        "output": "-o",
        "command_line": "#program_name# #output# #params#"
    }
    execute(compose_command(params, output_file=conf_file))
    return conf_file


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


def compress_bz2(input_file):
    """Compress BZ2 files"""
    params = {
        "program_name": "bzip2",
        "params": "{}".format(input_file),
        "command_line": "#program_name# #params#"
    }
    execute(compose_command(params, input_file=input_file))
    

def decompress_bz2(input_file, output_dir):
    """Decompresses BZ2 files"""
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
    stdin = None
    stdout = None
    environment = os.environ.copy()
    r_output_path = None
    r_output_file = None
    command_line = params['command_line']
    program_name = params['program_name']
    command_line = command_line.replace('#program_name#', program_name)

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
    # TODO: we should use shlex.split for this or drop this ugly function altogether (Michal)
    quotes = [j for j, e in enumerate(command_line) if e == '"']

    for s, e in zip(quotes[0::2], quotes[1::2]):
        command_line = command_line.replace(
            command_line[s + 1:e], command_line[s + 1:e].replace(' ', '#'))

    return {'command_line': [str(a).replace('#', ' ') for a in re.sub(' +', ' ', command_line.replace('"', '')).split(' ') if a],
            'stdin': stdin, 'stdout': stdout, 'env': environment, 'output_path': r_output_path, 'output_file': r_output_file}


def run_command(cmd, shell=False, **kwargs):
    """
    Runs a command and checks for exit code

    Args:
        cmd (str): A command to run
        shell (bool): Whether to invoke shell
        **kwargs: additional arguments passed to subprocess.run function

    Returns:

    """
    if not shell:
        cmd_s = shlex.split(cmd)
    else:
        cmd_s = cmd

    r = sb.run(cmd_s, shell=shell, capture_output=True, **kwargs)

    if r.returncode != 0:
        stdout = r.stdout
        stderr = r.stderr
        if isinstance(stdout, bytes):
            stdout = stdout.decode()
        if isinstance(stderr, bytes):
            stderr = stderr.decode()
        error(f'Execution failed for command {cmd}', exit=False)
        error(f'stdout:\n{stdout}', exit=False)
        error(f'stderr:\n{stderr}', exit=False)
        error('Exiting', exit=True)

    return r
