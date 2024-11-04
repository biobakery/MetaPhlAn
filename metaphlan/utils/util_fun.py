__author__ = ('Aitor Blanco Miguez (aitor.blancomiguez@unitn.it), '
              'Duy Tin Truong (duytin.truong@unitn.it), '
              'Francesco Asnicar (f.asnicar@unitn.it), '
              'Moreno Zolfo (moreno.zolfo@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)')
__version__ = '4.2.4'
__date__    = '21 Oct 2025'

import os
import sys
import time
import bz2
import gzip

def info(message, init_new_line=True, stderr=True, exit=False, exit_value=0):
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
        sys.stderr.write('\n')

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
        
def read_and_split(ofn):
    """Reads and splits a UTF-8 coded input stream

    Args:
        ofn (stream): the coded input stream

    Returns:
        list: the splitted decoded stream as a list
    """
    return (l.decode('utf-8').strip().split('\t') for l in ofn)

def read_and_split_line(line):
    """Reads and splits a UTF-8 coded line

    Args:
        line (str): the coded line

    Returns:
        list: the splitted decoded line
    """
    return line.decode('utf-8').strip().split('\t')

def plain_read_and_split(ofn):
    """Reads and splits an input stream

    Args:
        ofn (stream): the input stream

    Returns:
        list: the splitted stream as a list
    """
    return (l.strip().split('\t') for l in ofn)

def plain_read_and_split_line(l):
    """Reads and splits a line

    Args:
        line (str): the line

    Returns:
        list: the splitted line
    """
    return l.strip().split('\t')

def mybytes(val):
    """Encodes a value to bytes

    Args:
        val (any): the value

    Returns:
        bytes: the encoded value
    """
    return bytes(val, encoding='utf-8')

def byte_to_megabyte(byte):
    """Convert byte value to megabyte

    Args:
        byte (int): the byte value

    Returns:
        int: the megabyte value of the byte
    """
    return byte / (1024.0**2)


def openrt(file_path):
    file_path = str(file_path)
    if file_path.endswith('.bz2'):
        return bz2.open(file_path, 'rt')
    elif file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path)