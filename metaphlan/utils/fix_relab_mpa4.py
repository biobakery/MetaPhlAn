#!/usr/bin/env python
__author__ = 'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it)'
__version__ = '4.1.0'
__date__ = '23 Aug 2023'


import os, time
try:
    from .util_fun import info, error, warning
except ImportError:
    from util_fun import info, error, warning
import argparse as ap


def read_params():
    """ Reads and parses the command line arguments of the script

    Returns:
        namespace: The populated namespace with the command line arguments
    """
    p = ap.ArgumentParser(formatter_class=ap.RawTextHelpFormatter, add_help=False)   
    requiredNamed = p.add_argument_group('requiered arguments')
    requiredNamed.add_argument('--input', type=str, default=None, help="The path to the input profile")
    requiredNamed.add_argument('--output', type=str, default=None, help="The path to the output profile")
    return p.parse_args()


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
    elif not os.path.exists(args.output):
        error('The file {} does not exist'.format(
            args.output), exit=True)   

def fix_relab_mpa4(input, output):
    taxa_levs = [{},{},{},{},{},{},{},{}] 
    with open(input, 'r') as rf:
        with open(output, 'w') as wf:
            for line in rf:
                if line.startswith('#'):
                    wf.write(line)
                else:
                    if 't__' in line:
                        if 'p__Bacillota' in line:
                            line = line.replace('p__Bacillota', 'p__Firmicutes')
                        line = line.strip().split('\t') 
                        taxa_levs[-1][line[0]] = [line[1], float(line[2]), line[3] if len(line)==4 else '']
                    elif 's__' in line:
                        if 'p__Bacillota' in line:
                            line = line.replace('p__Bacillota', 'p__Firmicutes')
                        line = line.strip().split('\t')
                        taxa_levs[-2][line[0]] = [line[1], float(line[2]), '']
            for i in range(2,8):
                j = i+1
                for ss in taxa_levs[-i]:
                    gg = ss.replace('|{}'.format(ss.split('|')[-1]), '')
                    gg_n = taxa_levs[-i][ss][0].replace('|{}'.format(taxa_levs[-i][ss][0].split('|')[-1]), '')
                    if gg not in taxa_levs[-j]:
                        taxa_levs[-j][gg] = [gg_n, taxa_levs[-i][ss][1], '']
                    else:
                        taxa_levs[-j][gg][1] += taxa_levs[-i][ss][1]
            for level in taxa_levs:
                for tax in level:
                    wf.write(tax + '\t' + '\t'.join([str(x) for x in level[tax]]) + '\n')
    
def main():
    t0 = time.time()
    args = read_params()
    info("Start fixing profile")
    check_params(args)
    fix_relab_mpa4(args.input, args.output)
    exec_time = time.time() - t0
    info("Finish fixing profile ({} seconds)".format(round(exec_time, 2)))


if __name__ == '__main__':
    main()
