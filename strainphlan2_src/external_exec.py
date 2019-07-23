#!/usr/bin/env python

__author__ = ('Duy Tin Truong (duytin.truong@unitn.it), '
              'Francesco Asnicar (f.asnicar@unitn.it), '
              'Moreno Zolfo (moreno.zolfo@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it), '
              'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it)')
__version__ = '2.0.0'
__date__ = '17 Jul 2019'

import os, re
import subprocess as sb
from utils import info, error

#ToDo: We have to recover the info/errors from the std_out and std_err
"""
Executes a command

:param cmd: the command to execute
"""
def execute(cmd):
    inp_f = None
    out_f = sb.DEVNULL

    if cmd['stdin']:
        inp_f = open(cmd['stdin'], 'r')
    if cmd['stdout']:
        out_f = open(cmd['stdout'], 'w')
    info("\t"+str(cmd['command_line']), init_new_line=True)

    sb.run(cmd['command_line'], stdin=inp_f, stdout=out_f, stderr=sb.DEVNULL)

    if cmd['stdin']:
        inp_f.close()
    if cmd['stdout']:
        out_f.close()


#ToDo: Parametrize this function: default output_dir, remove the compressed file...
"""
Decompressed BZ2 files

:param reference: the reference FASTA as BZ2 files
:param tmp_dir: the temporal output directory
:returns: the decompressed file
"""
def decompress_bz2(reference, output_dir):
    params = {
        "program_name" : "bzip2",
        "input" : "-cdk",
        "command_line" : "#program_name# #input# > #output#"
    }      
    n, _ = os.path.splitext(os.path.basename(reference))
    execute(compose_command(params, input_file=reference, output_file=output_dir+n))
    return output_dir+n


"""
Compose a command for further executions

:param params: the params of the command
:param check: [default=False] check if program is available
:param sub_mod: [optional] the model
:param input_file: [optional] the input file
:param database: [optional] the database
:param output_path: [optional] the output path
:param output_file: [optional] the output file
:param nproc: [default=1] the number of procs to use
"""
def compose_command(params, check=False, input_file=None, database=None, output_path=None, output_file=None, nproc=1):
    program_name = None
    stdin = None
    stdout = None
    environment = os.environ.copy()
    r_output_path = None
    r_output_file = None
    command_line = params['command_line']

    if 'program_name' in list(params):
        command_line = command_line.replace('#program_name#', params['program_name'])
        program_name = params['program_name']
    else:
        error('Error: something wrong... '+program_name+' not found!', exit=True)

    if check:
        command_line = program_name

        if 'version' in list(params):
            command_line = '{} {}'.format(program_name, params['version'])
    else:
        if 'params' in list(params):
            command_line = command_line.replace('#params#', params['params'])

        if 'threads' in list(params):
            command_line = command_line.replace('#threads#', '{} {}'.format(params['threads'], nproc))

        if output_path:
            r_output_path = output_path

            if 'output_path' in list(params):
                command_line = command_line.replace('#output_path#', '{} {}'.format(params['output_path'], output_path))
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
            command_line = command_line.replace('#database#', '{} {}'.format(params['database'], database))

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
        command_line = command_line.replace(command_line[s + 1:e], command_line[s + 1:e].replace(' ', '#'))

    return {'command_line': [str(a).replace('#', ' ') for a in re.sub(' +', ' ', command_line.replace('"', '')).split(' ') if a],
            'stdin': stdin, 'stdout': stdout, 'env': environment, 'output_path': r_output_path, 'output_file': r_output_file}
