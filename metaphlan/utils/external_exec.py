__author__ = ('Aitor Blanco Miguez (aitor.blancomiguez@unitn.it), '
              'Duy Tin Truong (duytin.truong@unitn.it), '
              'Francesco Asnicar (f.asnicar@unitn.it), '
              'Moreno Zolfo (moreno.zolfo@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)')
__version__ = '3.0.8'
__date__ = '7 May 2021'

import os, sys, re, shutil
import subprocess as sb
try:
    from .util_fun import info, error
except ImportError:
    from util_fun import info, error

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

    exec_res = sb.run(cmd['command_line'], stdin=inp_f, stdout=out_f)
    if exec_res.returncode == 1:
        error("An error was ocurred executing a external tool, exiting...", 
            init_new_line=True, exit=True) 

    if cmd['stdin']:
        inp_f.close()
    if cmd['stdout']:
        out_f.close()


"""
Creates the BLASTn database to align the reference genomes

:param output_dir: the output directory
:param reference: the FASTA with the reference
:returns: the created BLASTn database
"""
def create_blastn_db(output_dir, reference):
    reference_name = os.path.splitext(os.path.basename(reference))[0]    
    params = {
        "program_name" : "makeblastdb",
        "params" : "-parse_seqids -dbtype nucl",
        "input" : "-in",
        "output" : "-out",
        "command_line" : "#program_name# #params# #input# #output#"
    }
    execute(compose_command(params, input_file=reference, output_file=output_dir+reference_name))
    return output_dir+reference_name


"""
Executes BLASTn

:param output_dir: the output directory
:param clade_markers: the FASTA with the markers
:param blastn_db: the BLASTn database
:param nprocs: the number of thread to use
:returns: the BLASTn output file
"""
def execute_blastn(output_dir, clade_markers, blastn_db, nprocs=1):
    db_name = os.path.splitext(os.path.basename(blastn_db))[0]
    params = {
        "program_name" : "blastn",
        "params" : "-outfmt \"6 qseqid sseqid qlen qstart qend sstart send\" -evalue 1e-10 -max_target_seqs 1",
        "input" : "-query",
        "database": "-db",
        "output" : "-out",
        "threads" : "-num_threads",
        "command_line" : "#program_name# #params# #threads# #database# #input# #output#"
    }
    execute(compose_command(params, input_file=clade_markers, database=blastn_db, 
        output_file=output_dir+db_name+".blastn", nproc=nprocs))
    return output_dir+db_name+".blastn"


"""
Creates the PhyloPhlAn database

:param output_dir: the output directory
:param clade: the clade
"""
def create_phylophlan_db(output_dir, clade):
    markers = output_dir+clade
    params = {
        "program_name" :"phylophlan_setup_database",
        "params" : "-d "+clade+" -e fna -t n --overwrite",
        "input" : "-i",
        "command_line" : "#program_name# #input# #params#"
    }
    execute(compose_command(params, input_file=markers))
    #os.rename(output_dir+clade+".fna", markers+"/"+clade+".fna")


"""
Generates the PhyloPhlan configuration file

:param output_dir: the output directory
:returns: the generated configuration file
"""
def generate_phylophlan_config_file(output_dir, configuration):
    conf_file = output_dir+"phylophlan.cfg"
    params = {
        "program_name" : "phylophlan_write_config_file",
        "params" : "-d n --db_dna makeblastdb --map_dna "+configuration['map']+
            " --msa "+configuration['aligner']+" --trim "+configuration['trim']+
            " --tree1 "+configuration['tree1'], #+
            # configuration['tree2'],
        "output" : "-o",
        "command_line" : "#program_name# #output# #params#"
    }
    execute(compose_command(params, output_file=conf_file))
    return conf_file


"""
Executes PhyloPhlAn

:param samples_markers_dir: the temporal samples markers directory
:param conf_file: the PhyloPhlAn configuration file
:param min_entries: the minimun number of entries to consider a good marker 
:param tmp_dir: the temporal output directory
:param output_dir: the output_directory
:param clade: the clade
:param phylogeny_conf: the precision of the phylogenetic analysis
:param mutation_rates: whether get  the mutation rates for the markers
:param nproc: the number of threads to run phylophlan
"""
def execute_phylophlan(samples_markers_dir, conf_file, min_entries, tmp_dir, output_dir, 
    clade, phylogeny_conf, mutation_rates, nprocs):
    accuracy = " --"+phylogeny_conf

    if mutation_rates:
        accuracy = accuracy + " --mutation_rates"
    params = {
        "program_name" : "phylophlan",
        "params" : "-d "+clade[:30]+" --data_folder "+tmp_dir+
            " --databases_folder "+tmp_dir+" -t n -f "+conf_file+
            " --diversity low"+accuracy+" --genome_extension fna"+
            " --force_nucleotides --min_num_entries "+str(min_entries),
        "input" : "-i",
        "output_path" : "--output_folder",
        "output" : "-o",
        "threads" : "--nproc",
        "command_line" : "#program_name# #input# #output# #output_path# #params# #threads#"
    }
    execute(compose_command(params=params, input_file=samples_markers_dir, output_path=output_dir,
        output_file=".", nproc=nprocs))
        

#ToDo: Parametrize this function: default output_dir, remove the compressed file...
"""
Decompressed BZ2 files

:param input: the input BZ2 file to decompress
:param output_dir: the output directory
:returns: the decompressed file
"""
def decompress_bz2(input, output_dir):
    n, _ = os.path.splitext(os.path.basename(input))
    params = {
        "program_name" : "bzip2",
        "input" : "-cdk",
        "command_line" : "#program_name# #input# > #output#"
    }      
    decompressed_file = output_dir + n 
    execute(compose_command(params, input_file=input, output_file=decompressed_file))
    if decompressed_file.endswith('_sam'):
        os.rename(decompressed_file, decompressed_file[:-4] + '.sam')
        decompressed_file = decompressed_file[:-4] + '.sam'
    return decompressed_file


"""
Converts SAM files to sorted BAM files using samtools

:param input: the input SAM file
:param output_dir: the output directory
:returns: the sorted BAM file
"""
def samtools_sam_to_bam(input, output_dir):
    n, _ = os.path.splitext(os.path.basename(input))   
    params = {
        "program_name" : "samtools",
        "params" : "view",
        "input" : "-Sb",
        "command_line" : "#program_name# #params# #input# > #output#"
    }       
    execute(compose_command(params, input_file=input, output_file=output_dir+n+".bam"))
    return output_dir+n+".bam"


"""
Sort BAM files using samtools

:param input: the input BAM file
:param output_dir: the output directory
:returns: the sorted BAM file
"""
def samtools_sort_bam_v0(input, output_dir):
    n, _ = os.path.splitext(os.path.basename(input))
    params = {
        "program_name" : "samtools",
        "params" : "sort",
        "command_line" : "#program_name# #params# #input# #output#"
    }        
    execute(compose_command(params, input_file=input, output_file=output_dir+n+".sorted"))
    return output_dir+n+".sorted.bam"


"""
Sort BAM files using samtools

:param input: the input BAM file
:param output_dir: the output directory
:returns: the sorted BAM file
"""
def samtools_sort_bam_v1(input, output_dir):
    n, _ = os.path.splitext(os.path.basename(input))  
    params = {
        "program_name" : "samtools",
        "params" : "sort",
        "output" : "-o",
        "command_line" : "#program_name# #params# #input# #output#"
    }      
    execute(compose_command(params, input_file=input, output_file=output_dir+n+".sorted.bam"))
    shutil.move(output_dir+n+".sorted.bam", output_dir+n+".bam")
    return output_dir+n+".bam"


"""
Generates a FASTA file with the markers form a MetaPhlAn database

:param database: the MetaPhlan markers database
:param output_dir: the output directory
"""
def generate_markers_fasta(database, output_dir):
    db_markers_faa = output_dir+"db_markers.fna"
    bowtie_database, _ = os.path.splitext(database)
    params = {
        "program_name" : "bowtie2-inspect",
        "command_line" : "#program_name# #input# > #output#"
    }
    execute(compose_command(params, input_file=bowtie_database, output_file=db_markers_faa))
    return db_markers_faa


"""
Compose a command for further executions

:param params: the params of the command
:param check: [default=False] check if program is available
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

