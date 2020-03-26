__author__ = ('Aitor Blanco Miguez (aitor.blancomiguez@unitn.it), '
              'Duy Tin Truong (duytin.truong@unitn.it), '
              'Francesco Asnicar (f.asnicar@unitn.it), '
              'Moreno Zolfo (moreno.zolfo@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)')
__version__ = '3.0'
__date__ = '21 Feb 2020'


import os, sys, re, pickletools, pickle, time, bz2, gzip
from hashlib import sha256

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
        sys.stderr.write('{}: Stop StrainPhlAn 3.0 execution.\n'.format(time.ctime(int(time.time()))))
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
        os.mkdir(path)
    except Exception as e:
        error('Folder \"'+path+'\" already exists!\n'+str(e), exit=True,
            init_new_line=True)  


"""
Parse marker names to support PhyloPhlAn execution

:param marker_name: the old marker name
:returns: the parsed marker name
"""
def parse_marker_name(marker_name):
    return str(int(sha256(marker_name.encode('utf-8')).hexdigest(), 16) % 10**12)
        

"""
Gets the Breath of Coverage measure for a consensus sequence

:param sequence: the consensus sequence
:returns: the breath of coverage
"""
def get_breath(sequence):
    seq_len = len(sequence)
    return ((seq_len - sequence.count('N') - sequence.count('*') - sequence.count('-')) * 100) / seq_len 


def openr(fn, mode="r"):
    if fn is None:
        return sys.stdin

    if fn.endswith(".bz2"):
        if sys.version_info[0] < 3:
            return bz2.BZ2File(fn)
        else:
            return bz2.open(fn, 'rt')
    elif fn.endswith(".gz"):
        if sys.version_info[0] < 3:
            return None  # need to check if gzip is different in Python2
        else:
            return gzip.open(fn, 'rt')
    else:
        return open(fn, mode)


def openw(fn):
    if fn is None:
        return sys.stdout

    if fn.endswith(".bz2"):
        if sys.version_info[0] < 3:
            return bz2.BZ2File(fn, "w")
        else:
            return bz2.open(fn, 'wt')
    elif fn.endswith(".gz"):
        if sys.version_info[0] < 3:
            return None  # need to check if gzip is different in Python2
        else:
            return gzip.open(fn, 'wt')
    else:
        return open(fn, "w")


def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False