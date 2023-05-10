#!/usr/bin/env python3

import argparse
import os
import sys
import re
import glob
import pandas as pd
from itertools import takewhile

def merge( aaastrIn, ostm, options ):
    """
    Outputs the table join of the given pre-split string collection.

    :param  aaastrIn:   One or more split lines from which data are read.
    :type   aaastrIn:   collection of collections of string collections
    :param  ostm:       Output stream to which matched rows are written.
    :type   ostm:       output stream
    :param  options:    Options tuple: ( group_by field, suppress extra_info )
    :type   options:    tuple
    
    """

    listmpaVersion = set()
    merged_tables = pd.DataFrame()

    groupby_field, suppress_info = options


    if len(aaastrIn) == 1 and os.path.isdir(aaastrIn[0]):
        inList = glob.glob(aaastrIn[0]+'/*.tsv')
    else:
        inList = aaastrIn
    for f in inList:
        
        with open(f) as fin:
            headers = [x.strip() for x in takewhile(lambda x: x.startswith('#'), fin)]


        # get the headers. If sample name is not the default, use it.
        # otherwise, use the filename

        smpl_name = None
        if len(headers) >= 3:
            smpl_name = headers[-1].split('#')[1].strip().split('\t')[1]

        if not smpl_name or smpl_name == 'Metaphlan_Analysis': 
            smpl_name = os.path.splitext(os.path.basename(f))[0]


        # take mpaVersion

        mpaVersion = list(filter(re.compile('#mpa_v[0-9]{2,}_CHOCOPhlAn_[0-9]{0,}').match, headers))
        if len(mpaVersion):
            listmpaVersion.add(mpaVersion[0])

        if len(listmpaVersion) > 1:
            print('merge_metaphlan_tables found tables made with different versions of the MetaPhlAn2 database.\nPlease re-run MetaPhlAn2 with the same database.\n')
            return
         
        iIn = pd.read_csv(f,sep='\t',skiprows=len(headers)).assign(sampleID = smpl_name).fillna('')

        if merged_tables.empty:
            merged_tables = iIn
        else:
            merged_tables = pd.concat([iIn, merged_tables])

    indexes_for_pivot = 'M-Group/Cluster' if suppress_info else ['M-Group/Cluster','M-Group-Type [k|u]','First Genome in Cluster','Other Genomes']
    out_table = pd.pivot_table(merged_tables, index = indexes_for_pivot, columns = 'sampleID', values = groupby_field)

    if listmpaVersion:
        ostm.write(list(listmpaVersion)[0]+'\n')
    out_table.fillna('0').reset_index().to_csv(ostm, index=False, sep = '\t')

argp = argparse.ArgumentParser( prog = "merge_metaphlan_tables.py",
    description = """Performs a table join on one or more metaphlan output files.""")
argp.add_argument( "aistms",    metavar = "input.txt", nargs = "+",
    help = "One or more tab-delimited text tables to join" )
argp.add_argument( '-o',    metavar = "output.txt", nargs = 1,
    help = "Name of output file in which joined tables are saved" )
argp.add_argument( '-g',  choices=['breadth_of_coverage','depth_of_coverage_mean','depth_of_coverage_median'],
    help = "Value to group by (default is 'breadth_of_coverage')",
    default = "breadth_of_coverage")
argp.add_argument( '--no_info',  action='store_true',
    help = "Include VSCs info in output file (default is unset, meaning info are included) ")
__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" )

argp.usage = argp.format_usage()[7:]+"\n\n\tPlease make sure to supply file paths to the files to combine. If combining 3 files (Table1.txt, Table2.txt, and Table3.txt) the call should be:\n\n\t\tpython merge_metaphlan_tables.py Table1.txt Table2.txt Table3.txt > output.txt\n\n\tA wildcard to indicate all .txt files that start with Table can be used as follows:\n\n\t\tpython merge_metaphlan_tables.py Table*.txt > output.txt"

def main( ):
    args = argp.parse_args( )

    options = ( args.g, args.no_info )

    if args.o is None:
        merge(args.aistms, sys.stdout, options)
    else:
        with open(args.o[0], 'w') as fout:
            merge(args.aistms, fout, options)

if __name__ == '__main__':
    main()