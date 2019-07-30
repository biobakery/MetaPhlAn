#!/usr/bin/env python

__author__ = ('Duy Tin Truong (duytin.truong@unitn.it), '
              'Francesco Asnicar (f.asnicar@unitn.it), '
              'Moreno Zolfo (moreno.zolfo@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it), '
              'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it)')
__version__ = '2.0.0'
__date__ = '29 Jul 2019'


import os, sys, re, pickletools, pickle, time

"""
Prints a message as normal info

:param message: the message
:param init_new_line: inserts a new line before the message
:param exit: exists after print the message
:param exit_value: the exit value if exit=True
"""
def info(message, init_new_line=False, exit=False, exit_value=0):
    if init_new_line:
        sys.stdout.write('\n')
    sys.stdout.write('{}: '.format(time.ctime(int(time.time()))))
    sys.stdout.write('{}'.format(message))
    sys.stdout.flush()

    if exit:
        sys.exit(exit_value)


"""
Prints a message as an error

:param message: the message
:param init_new_line: inserts a new line before the message
:param exit: exists after print the message
:param exit_value: the exit value if exit=True
"""
def error(message, init_new_line=False, exit=False, exit_value=1):
    if init_new_line:
        sys.stderr.write('\n')    
    sys.stderr.write('[e] {}\n'.format(message))
    sys.stderr.flush()

    if exit:
        sys.stderr.write('{}: Stop StrainPhlAn2 execution.\n'.format(time.ctime(int(time.time()))))
        sys.exit(exit_value)


"""
Optimized method for write Pickle files

:param fout: the file to write
:param elem: the element to write on file
"""
def optimized_dump(fout, elem):
    fout.write(pickletools.optimize(pickle.dumps(elem, pickle.HIGHEST_PROTOCOL)))


"""
Creates a folder if the path does not exists,
if not, returns an error

:param path: the path of the new folder
"""
def create_folder(path):
    try:
        os.mkdir(path, 755)
    except Exception as e:
        error('Folder \"'+path+'\" already exists!\n'+str(e), exit=True,
            init_new_line=True)  


"""
Checks if the submitted clade is in the species level

:param clade: the submitted clade
:returns: whether the clade is in the species level
"""
def check_clade(clade):
    species_pattern = re.compile("s__[A-z]*_[A-z]*")
    return species_pattern.match(clade)


"""
Parse marker names to support PhyloPhlAn execution

:param marker_name: the old marker name
:returns: the parsed marker name
"""
def parse_marker_name(marker_name):
    return re.sub('[^a-zA-Z0-9 \n\.]', '-', marker_name)
        

"""
Gets the Breath of Coverage measure for a consensus sequence

:param sequence: the consensus sequence
:returns: the breath of coverage
"""
def get_breath(sequence):
    seq_len = len(sequence)
    return ((seq_len - sequence.count('N')) * 100) / seq_len 