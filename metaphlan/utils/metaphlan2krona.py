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
            lin,tax,frac = aline.split()[0:3]
            tax = tax.split('|')[-1]
            x = re.sub(re_replace, '\t', lin)
            lineage = re.sub(re_bar, '', x)
            abundance = float(frac.rstrip('\n'))
            metaPhLan_FH.write(f'{abundance}\t{lineage}\t{tax}\n')

    metaPhLan_FH.close()

if __name__ == '__main__':
    main()
