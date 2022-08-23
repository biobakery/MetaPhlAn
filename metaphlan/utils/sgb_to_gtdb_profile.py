__author__ = 'Aitor Blanco (aitor.blancomiguez@unitn.it'
__version__ = '4.0.0'
__date__ = '22 Aug 2022'

import os, time, sys
import argparse as ap

try:
    from .util_fun import openw, info, error
except ImportError:
    from util_fun import openw, info, error

GTDB_ASSIGNMENT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mpa_vJan21_CHOCOPhlAnSGB_202103_SGB2GTDB.tsv")

def read_params():
    p = ap.ArgumentParser(description="", formatter_class=ap.ArgumentDefaultsHelpFormatter)
    p.add_argument('-i', '--input', type=str, default=None,
                   help="The input profile")
    p.add_argument('-o', '--output', type=str, default=None,
                   help="The output profile")
    
    return p.parse_args()

def check_params(args):
    if not args.input:
        error('-i (or --input) must be specified', exit=True, 
            init_new_line=True)
    if not args.output:
        error('-o (or --output) must be specified', exit=True, 
            init_new_line=True)

def get_gtdb_profile(mpa_profile, gtdb_profile):
    tax_levels = ['d','p','c','o','f','g','s']
    sgb2gtdb = dict()
    with open(GTDB_ASSIGNMENT_FILE, 'r') as read_file:
        for line in read_file:
            line = line.strip().split('\t')
            sgb2gtdb[line[0]] = line[1]

    with open(gtdb_profile, 'w') as wf:
        with open(mpa_profile, 'r') as rf:
            unclassified = 0
            abundances = {x:dict() for x in tax_levels}
            for line in rf:
                if line.startswith('#mpa_'):
                    wf.write(line)
                    wf.write('#clade_name\trelative_abundance\n')
                elif line.startswith('UNCLASSIFIED'):
                    unclassified = float(line.strip().split('\t')[2])
                    wf.write('UNCLASSIFIED\t{}\n'.format(unclassified))
                elif 't__SGB' in line:
                    line = line.strip().split('\t')
                    gtdb_tax = sgb2gtdb[line[0].split('|')[-1][3:]]
                    if gtdb_tax not in abundances:
                        abundances['s'][gtdb_tax] = 0
                    abundances['s'][gtdb_tax] += float(line[2])
        tax_levels.reverse()
        for tax_level in tax_levels[:-1]:
            for tax in abundances[tax_level]:
                new_tax = tax.replace(';{}'.format(tax.split(';')[-1]), '')
                new_level = tax_levels[tax_levels.index(tax_level)+1]
                if new_tax not in abundances[new_level]:
                    abundances[new_level][new_tax] = 0
                abundances[new_level][new_tax] += abundances[tax_level][tax]
        tax_levels.reverse()
        for tax_level in tax_levels:
            for tax in abundances[tax_level]:
                wf.write('{}\t{}\n'.format(tax, abundances[tax_level][tax]))

def main():
    t0 = time.time()
    args = read_params()
    info("Start execution")
    check_params(args)

    get_gtdb_profile(args.input, args.output)

    exec_time = time.time() - t0
    info("Finish execution ("+str(round(exec_time, 2))+" seconds)\n", 
        init_new_line=True)
	
if __name__ == "__main__":
	main()
