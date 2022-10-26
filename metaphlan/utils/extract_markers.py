#!/usr/bin/env python
__author__ = ('Aitor Blanco Miguez (aitor.blancomiguez@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)')
__version__ = '4.0.3'
__date__ = '24 Oct 2022'


import os
import time
import argparse as ap
try:
    from .util_fun import info, error
    from .database_controller import MetaphlanDatabaseController
except ImportError:
    from util_fun import info, error
    from database_controller import MetaphlanDatabaseController


def read_params():
    """ Reads and parses the command line arguments of the script

    Returns:
        namespace: The populated namespace with the command line arguments
    """
    p = ap.ArgumentParser(
        description="", formatter_class=ap.ArgumentDefaultsHelpFormatter)
    p.add_argument('-d', '--database', type=str, default='latest',
                   help="The input MetaPhlAn database")
    p.add_argument('-c', '--clades', type=str, nargs='+', default=[],
                   help="The clades to investigate")
    p.add_argument('-o', '--output_dir', type=str, default=None,
                   help="The output directory")
    return p.parse_args()


def check_params(args):
    """Checks the mandatory command line arguments of the script

    Args:
        args (namespace): the arguments to check
    """
    if len(args.clades) == 0:
        error('-c (or --clades) must be specified', exit=True)
    elif not args.output_dir:
        error('-o (or --output_dir) must be specified', exit=True)
    elif not os.path.exists(args.output_dir):
        error('The directory {} does not exist'.format(
            args.output_dir), exit=True)
    elif args.database != 'latest' and not os.path.exists(args.database):
        error('The database does not exist', exit=True)


def main():
    t0 = time.time()
    args = read_params()
    info("Start extract markers execution")
    check_params(args)
    database_controller = MetaphlanDatabaseController(args.database)
    database_controller.extract_markers(args.clades, args.output_dir)
    exec_time = time.time() - t0
    info("Finish extract markers execution ({} seconds): Results are stored at \"{}\"".format(
        round(exec_time, 2), args.output_dir))


if __name__ == '__main__':
    main()
