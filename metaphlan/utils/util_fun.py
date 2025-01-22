__author__ = ('Aitor Blanco Miguez (aitor.blancomiguez@unitn.it), '
              'Duy Tin Truong (duytin.truong@unitn.it), '
              'Francesco Asnicar (f.asnicar@unitn.it), '
              'Moreno Zolfo (moreno.zolfo@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)')
__version__ = '4.1.1'
__date__ = '11 Mar 2024'

import bz2
import gzip
import os
import sys
from datetime import datetime


class GlobalFlags:
    debug = False


global_flags = GlobalFlags()


def message(*args, outw):
    d = datetime.now().strftime("%d/%m/%y %H:%M:%S")
    outw.write(f"[{d}] ")
    outw.write('{}'.format(' '.join(map(str, args))))
    outw.write('\n')
    outw.flush()


def info_debug(*args):
    if not global_flags.debug:
        return
    message('Debug:', *args, outw=sys.stdout)


def info(*args):
    """

    :param args:
    :return:
    """
    message('Info:', *args, outw=sys.stdout)



def warning(*args):
    """

    :param args:
    :return:
    """
    message('Warning:', *args, outw=sys.stderr)


def error(*args, exit=False, exit_value=1):
    """

    :param args:
    :param exit:
    :param exit_value:
    :return:
    """
    message('Error:', *args, outw=sys.stderr)
    if exit:
        message('Exiting.', outw=sys.stderr)
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
