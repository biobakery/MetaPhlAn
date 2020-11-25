#!/usr/bin/env python

# ============================================================================== 
# Conversion script: from MetaPhlAn output to Krona text input file
# Author: Daniel Brami (daniel.brami@gmail.com)
# ==============================================================================

import sys
import optparse
import re

def main():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-p', '--profile', dest='profile', default='', action='store', help='The input file is the MetaPhlAn standard result file' )
    parser.add_option( '-k', '--krona', dest='krona', default='krona.out', action='store', help='the Krons output file name' )
    ( options, spillover ) = parser.parse_args()

    if not options.profile or not options.krona:
        parser.print_help()
        sys.exit()

    re_candidates = re.compile(r"s__")
    re_replace = re.compile(r"\w__")
    re_bar = re.compile(r"\|")

    metaPhLan = list()
    with open(options.profile,'r') as f:
        metaPhLan = f.readlines()
    f.close()

    krona_tmp = options.krona 
    metaPhLan_FH = open(krona_tmp, 'w')

    for aline in (metaPhLan):
        if(re.search(re_candidates, aline)):
            x=re.sub(re_replace, '\t', aline)
            x=re.sub(re_bar, '', x)

            x_cells = x.split('\t')
            lineage = '\t'.join(x_cells[0:(len(x_cells) -1)])
            abundance = float(x_cells[-1].rstrip('\n')) 

            metaPhLan_FH.write('%s\n'%(str(abundance) + '\t' + lineage))

    metaPhLan_FH.close()

if __name__ == '__main__':
    main()
