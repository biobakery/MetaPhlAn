__author__ = ('Aitor Blanco Miguez (aitor.blancomiguez@unitn.it), '
              'Duy Tin Truong (duytin.truong@unitn.it), '
              'Francesco Asnicar (f.asnicar@unitn.it), '
              'Moreno Zolfo (moreno.zolfo@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)')
__version__ = '4.1.2'
__date__ = '16 Oct 2025'

import bz2
import gzip
import os
import sys
import time


def info(message, init_new_line=True, stderr=False, exit=False, exit_value=0):
    """Prints an info message

    Args:
        message (str): The message to print
        init_new_line (bool, optional): Whether to print a new line after the message. Defaults to True.
        exit (bool, optional): Whether to finish the execution after the message. Defaults to False.
        exit_value (int, optional): The exit value. Defaults to 0.
    """
    outw = sys.stdout if not stderr else sys.stderr
    outw.write('{}: '.format(time.ctime(int(time.time()))))
    outw.write('{}'.format(message))
    outw.flush()
    if init_new_line:
        outw.write('\n')
    if exit:
        sys.exit(exit_value)

def warning(message, init_new_line=True, exit=False, exit_value=0):
    """Prints an Warning message

    Args:
        message (str): The message to print
        init_new_line (bool, optional): Whether to print a new line after the message. Defaults to True.
        exit (bool, optional): Whether to finish the execution after the message. Defaults to False.
        exit_value (int, optional): The exit value. Defaults to 0.
    """
    sys.stderr.write('{}: '.format(time.ctime(int(time.time()))))
    sys.stderr.write('[Warning] {}'.format(message))
    sys.stderr.flush()
    if init_new_line:
        sys.stderr.write('\n')
    if exit:
        sys.exit(exit_value)


def error(message, init_new_line=True, exit=False, exit_value=1):
    """Prints an error message

    Args:
        message (str): The message to print
        init_new_line (bool, optional): Whether to print a new line after the message. Defaults to True.
        exit (bool, optional): Whether to finish the execution after the message. Defaults to False.
        exit_value (int, optional): The exit value. Defaults to 1.
    """
    sys.stderr.write('{}: '.format(time.ctime(int(time.time()))))
    sys.stderr.write('[Error] {}'.format(message))
    sys.stderr.flush()
    if init_new_line:
        sys.stdout.write('\n')

    if exit:
        sys.stderr.write('{}: Stop execution.\n'.format(
            time.ctime(int(time.time()))))
        sys.exit(exit_value)


def create_folder(path):
    """Creates a folder and throws errors when not possible

    Args:
        path (str): The path of the folder to create
    """
    try:
        if os.path.exists(path):
            error('Folder \"{}\" already exists!', exit=True)
        else:
            os.mkdir(path)
    except Exception as e:
        error('An error ocurred when creating the \"{}\" folder'.format(e), exit=True)


def openrt(file_path):
    file_path = str(file_path)
    if file_path.endswith('.bz2'):
        return bz2.open(file_path, 'rt')
    elif file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path)
