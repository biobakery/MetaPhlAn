#!/usr/bin/env python3


import argparse
import os
import sys
import pandas as pd
from itertools import takewhile


def merge(aaastrIn, ostm):
    """
    Outputs the table join of the given pre-split string collection.

    :param  aaastrIn:   One or more split lines from which data are read.
    :type   aaastrIn:   collection of collections of string collections
    :param  ostm:       Output stream to which matched rows are written.
    :type   ostm:       output stream
    """

    listmpaVersion = set()
    profiles_list = []
    merged_tables = None

    for f in aaastrIn:
        headers = [x.strip() for x in takewhile(lambda x: x.startswith('#'), open(f))]
        listmpaVersion.add(headers[0])

        if len(headers) == 4:
            names = headers[-1].split('#')[1].strip().split('\t')
        else:
            print('merge_metaphlan_tables: wrong header format for "{}", please check your profiles.\n'.format(f))
            return

        if len(listmpaVersion) > 1:
            print('merge_metaphlan_tables: profiles from differrent versions of MetaPhlAn, please profile your '
                  'samples using the same MetaPhlAn version.\n')
            return

        iIn = pd.read_csv(f, sep='\t', skiprows=len(headers), names=names, usecols=range(3), index_col=0)
        profiles_list.append(pd.Series(data=iIn['relative_abundance'], index=iIn.index,
                                       name=os.path.splitext(os.path.basename(f))[0].replace('_profile', '')))

    merged_tables = pd.concat([merged_tables, pd.concat(profiles_list, axis=1).fillna(0)], axis=1).fillna(0)
    ostm.write(list(listmpaVersion)[0]+'\n')
    merged_tables.to_csv(ostm, sep='\t')


argp = argparse.ArgumentParser(prog="merge_metaphlan_tables.py",
                               description="Performs a table join on one or more metaphlan output files.")
argp.add_argument("aistms", metavar="input.txt", nargs="*", help="One or more tab-delimited text tables to join")
argp.add_argument("-l", help="Name of file containing the paths to the files to combine")
argp.add_argument('-o', metavar="output.txt", help="Name of output file in which joined tables are saved")
argp.add_argument('--overwrite', default=False, action='store_true', help="Overwrite output file if exists")

argp.usage = (argp.format_usage() + "\nPlease make sure to supply file paths to the files to combine.\n\n" +
              "If combining 3 files (Table1.txt, Table2.txt, and Table3.txt) the call should be:\n" +
              "   ./merge_metaphlan_tables.py Table1.txt Table2.txt Table3.txt > output.txt\n\n" +
              "A wildcard to indicate all .txt files that start with Table can be used as follows:\n" +
              "    ./merge_metaphlan_tables.py Table*.txt > output.txt")


def main( ):
    args = argp.parse_args()

    if args.l:
        args.aistms = [x.strip().split()[0] for x in open(args.l)]

    if not args.aistms:
        print('merge_metaphlan_tables: no inputs to merge!')
        return

    if args.o and os.path.exists(args.o) and not args.overwrite:
        print('merge_metaphlan_tables: output file "{}" exists, specify the --overwrite param to ovrewrite it!'.format(args.o))
        return

    merge(args.aistms, open(args.o, 'w') if args.o else sys.stdout)


if __name__ == '__main__':
    main()
