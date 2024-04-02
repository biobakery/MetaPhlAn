#!/usr/bin/env python
__author__ = ('Aitor Blanco Miguez (aitor.blancomiguez@unitn.it), '
              'Claudia Mengoni (claudia.mengoni@unitn.it)')
__version__ = '4.1.1'
__date__ = '11 Mar 2024'


import os, time
try:
    from .util_fun import info, error, warning
except ImportError:
    from util_fun import info, error, warning
import argparse as ap

script_install_folder = os.path.dirname(os.path.abspath(__file__))
OCT22_FIXES=os.path.join(script_install_folder,'oct22_fix_tax.tsv')

def read_params():
    """ Reads and parses the command line arguments of the script

    Returns:
        namespace: The populated namespace with the command line arguments
    """
    p = ap.ArgumentParser(formatter_class=ap.RawTextHelpFormatter, add_help=False)   
    requiredNamed = p.add_argument_group('required arguments')
    requiredNamed.add_argument('--input', type=str, default=None, help="The path to the input profile")
    requiredNamed.add_argument('--output', type=str, default=None, help="The path to the output profile")
    return p.parse_args()

def read_oct22_fixes(file):
    """Reads the tab separated file with old and new taxonomies of Oct22
    Args:
        file: file with Oct22 fixes
    """
    oct_fixes=dict()
    with open(file) as inf:
        for l in inf.readlines()[1:]:
            old_tax, new_tax, new_tax_id = l.split('\t') 
            oct_fixes[old_tax]= (new_tax, new_tax_id.strip())
    return oct_fixes


def check_params(args):
    """Checks the mandatory command line arguments of the script

    Args:
        args (namespace): the arguments to check
    """
    if not args.input:
        error('--input must be specified', exit=True)
    elif not os.path.exists(args.input):
        error('The file {} does not exist'.format(
            args.input), exit=True)           
    if not args.output:
        error('--output must be specified', exit=True)

def fix_relab_mpa4(input, output):
    taxa_levs = [{},{},{},{},{},{},{},{}] 
    with open(input, 'r') as rf:
        with open(output, 'w') as wf:
            for line in rf:
                if line.startswith('#mpa_v'):
                    release = line.strip()[1:]
                    wf.write('{}_202403\n'.format(line.strip()))
                elif line.startswith('#') or line.startswith('UNCLASSIFIED'):
                    wf.write(line)
                else:
                    if 't__' in line:
                        if release == 'mpa_vJun23_CHOCOPhlAnSGB_202307':
                            if 'p__Bacillota' in line:
                                line = line.replace('p__Bacillota', 'p__Firmicutes')
                            elif 'f__Saccharomycetales_unclassified' in line:
                                line = line.replace('f__Saccharomycetales_unclassified','f__Debaryomycetaceae')
                            line = line.strip().split('\t') 

                        elif release == 'mpa_vOct22_CHOCOPhlAnSGB_202212':
                                line = line.strip().split('\t')
                                if line[0] in oct_fixes:
                                    line[0],line[1] = oct_fixes[line[0]]
                                    
                        taxa_levs[-1][line[0]] = [line[1], round(float(line[2]),5), line[3] if len(line)==4 else '']

            for i in range(1,8):
                j = i+1
                for ss in taxa_levs[-i]:
                    gg = ss.replace('|{}'.format(ss.split('|')[-1]), '')
                    gg_n = '|'.join(taxa_levs[-i][ss][0].split('|')[:-1])
                    if gg not in taxa_levs[-j]:
                        taxa_levs[-j][gg] = [gg_n, taxa_levs[-i][ss][1], '']
                    else:
                        taxa_levs[-j][gg][1] += taxa_levs[-i][ss][1]
            for level in taxa_levs:
                for tax in level:
                    wf.write(tax + '\t' + '\t'.join([str(x) for x in level[tax]]) + '\n')
    
def main():
    global oct_fixes
    t0 = time.time()
    args = read_params()
    info("Start fixing profile")
    check_params(args)
    oct_fixes = read_oct22_fixes(OCT22_FIXES) 
    fix_relab_mpa4(args.input, args.output)
    exec_time = time.time() - t0
    info("Finish fixing profile ({} seconds)".format(round(exec_time, 2)))


if __name__ == '__main__':
    main()
