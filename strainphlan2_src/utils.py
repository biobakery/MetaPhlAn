#!/usr/bin/env python

__author__ = ('Duy Tin Truong (duytin.truong@unitn.it), '
              'Francesco Asnicar (f.asnicar@unitn.it), '
              'Moreno Zolfo (moreno.zolfo@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it), '
              'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it)')
__version__ = '2.0.0'
__date__ = '17 Jul 2019'


import sys, pickletools, pickle

"""
Prints a message as normal info

:param s: the message
:param init_new_line: inserts a new line before the message
:param exit: exists after print the message
:param exit_value: the exit value if exit=True
"""
def info(s, init_new_line=False, exit=False, exit_value=0):
    if init_new_line:
        sys.stdout.write('\n')

    sys.stdout.write('{}'.format(s))
    sys.stdout.flush()

    if exit:
        sys.exit(exit_value)


"""
Prints a message as an error

:param s: the message
:param init_new_line: inserts a new line before the message
:param exit: exists after print the message
:param exit_value: the exit value if exit=True
"""
def error(s, init_new_line=False, exit=False, exit_value=1):
    if init_new_line:
        sys.stderr.write('\n')

    sys.stderr.write('[e] {}\n'.format(s))
    sys.stderr.flush()

    if exit:
        sys.exit(exit_value)


"""
Optimized method for write Pickle files

:param fout: the file to write
:param element: the element to write on file
"""
def optimized_dump(fout, elem):
    fout.write(pickletools.optimize(pickle.dumps(elem, pickle.HIGHEST_PROTOCOL)))