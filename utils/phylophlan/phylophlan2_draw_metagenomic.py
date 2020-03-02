#!/usr/bin/env python3


__author__ = ('Francesco Asnicar (f.asnicar@unitn.it), '
              'Claudia Mengoni (claudia.mengoni@studenti.unitn.it)')
__version__ = '0.05'
__date__ = '11 September 2019'


import argparse as ap
import os
import sys
import time
from collections import Counter
import math
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt


OUTPUT_NAME = 'output_heatmap'

if sys.version_info[0] < 3:
    raise Exception("PhyloPhlAn2 requires Python 3, your current Python version is {}.{}.{}"
                    .format(sys.version_info[0], sys.version_info[1], sys.version_info[2]))


def info(s, init_new_line=False, exit=False, exit_value=0):
    if init_new_line:
        sys.stdout.write('\n')

    sys.stdout.write('{}'.format(s))
    sys.stdout.flush()

    if exit:
        sys.exit(exit_value)


def error(s, init_new_line=False, exit=False, exit_value=1):
    if init_new_line:
        sys.stderr.write('\n')

    sys.stderr.write('[e] {}\n'.format(s))
    sys.stderr.flush()

    if exit:
        sys.exit(exit_value)


def read_params():
    p = ap.ArgumentParser(description=("The phylophlan2_draw_metagenomic.py script takes as input the output table generated form the "
                                       "phylophlan2_metagenomic.py script and produces two heatmap figures: (1) presence/absence heatmap "
                                       "of the SGBs and the metagenomic samples of the recontructed input genomes; and (2) heatmap "
                                       "showing the amount of kSGB, uSGB, and unassinged for each metagenome"),
                          formatter_class=ap.ArgumentDefaultsHelpFormatter)

    p.add_argument('-i', '--input', type=str, required=True, help='The input file generated from phylophlan2_metagenomic.py')
    p.add_argument('-m', '--map', required=True, type=str, help='A mapping file that maps each bin to its metagenome')
    p.add_argument('--top', default=20, type=int, help='The number of SGBs to display in the figure')
    p.add_argument('-o', '--output', type=str, default=OUTPUT_NAME, help='Prefix output files')
    p.add_argument('-s', '--separator', type=str, default='\t', help='The separator used in the mapping file')
    p.add_argument('--dpi', type=int, default=200, help='Dpi resolution of the images')
    p.add_argument('-f', type=str, default='svg', help='Images output format')
    p.add_argument('--verbose', action='store_true', default=False, help="Prints more stuff")
    p.add_argument('-v', '--version', action='version',
                   version='phylophlan2_draw_metagenomic.py version {} ({})'.format(__version__, __date__),
                   help="Prints the current phylophlan2_draw_metagenomic.py version and exit")
    return p.parse_args()


def check_params(args, verbose=False):
    if not os.path.isfile(args.input):
        error('input file {} does not exist'.format(args.input), exit=True)

    if not os.path.isfile(args.map):
        error('--map file {} does not exist'.format(args.map), exit=True)

    if args.top < 1:
        error('--top cannot be 0 or negative', exit=True)

    if os.path.dirname(args.output):
        if not os.path.exists(os.path.dirname(args.output)):
            error('output path does not exists: "{}"'.format(args.output), exit=True)

    if not os.path.basename(args.output):
        old = args.output
        args.output = os.path.join(old, OUTPUT_NAME)
        info('No output filename specified "{}", writing output to "{}"\n'.format(old, args.output))

    if verbose:
        info('Arguments: {}\n'.format(vars(args)))


def read_input(inputt, map_dict, verbose=False):
    d = dict([(m, []) for m in map_dict.values()])
    unassigned = dict([(m, []) for m in map_dict.values()])

    with open(inputt) as f:
        for r in f:
            if r.startswith('#'):
                continue

            rc = r.strip().split('\t')

            bin_id = rc[0]
            sgb = rc[1]
            avg_dist = float(sgb.split(':')[3])

            sgb_id = ' '.join(sgb.split(':')[0].split('_'))
            taxa_level = sgb.split(':')[1]
            taxonomy = sgb.split(':')[2]

            label = [sgb_id]

            if avg_dist > 0.05:
                unassigned[map_dict[bin_id]].append(label)
                continue

            if taxa_level == 'Species':
                t = ' '.join([l for l in taxonomy.split('|') if l.startswith('s__')][0].replace('s__', '').split('_'))
                label.append('(' + t + ')')
            elif taxa_level == 'Genus':
                t = [l for l in taxonomy.split('|') if l.startswith('g__')][0].replace('g__', '')
                label.append('(' + t + ' genus)')
            elif taxa_level == 'Family':
                t = [l for l in taxonomy.split('|') if l.startswith('f__')][0].replace('f__', '')
                label.append('(' + t + ' family)')
            elif taxa_level == 'Other':
                t = [l for l in taxonomy.split('|') if l.startswith('p__')][0].replace('p__', '')
                label.append('(' + t + ' phylum)')
            else:
                error('Taxa level {} not valid \n'.format(taxa_level))

            d[map_dict[bin_id]].append(' '.join(label))

    return (d, unassigned)


def bin2met(args, sep):
    return dict([(r.strip().split(sep)[0], r.strip().split(sep)[1]) for r in open(args.map)])


def find_top_SGBs(top, meta_dict, verbose=False):
    c = Counter()

    for x in meta_dict:
        tmpc = Counter(meta_dict[x])
        c = tmpc + c

    if len(c) < top:
        if verbose:
            info('Number of top SGBs was set to {} but there are only {}.\n\n'.format(top, len(c)))

        top = len(c)

    elif len(c) - top <= 5 and len(c) != top:
        if verbose:
            info('Number of top SGBs was set to {} but since there are in total {} SGBs, all were left in.\n\n'.format(top, len(c)))

        top = len(c)

    elif len(c) > top:
        c1 = c.most_common(top + 1)
        i = top

        while(c1[-1][1] == c1[-2][1] and i < len(c)):
            i += 1
            c1 = c.most_common(i + 1)

        if verbose:
            info('Top SGBs is {}, but reporting {} SGBs as they were at the same distance.\n'.format(top, i - top))

        top = i

    c = c.most_common(top)

    return [x[0] for x in c]


def phylophlan2_draw_metagenomic():
    args = read_params()

    if args.verbose:
        info('phylophlan2_draw_metagenomic.py version {} ({})\n'.format(__version__, __date__))
        info('Command line: {}\n\n'.format(' '.join(sys.argv)), init_new_line=True)

    check_params(args, verbose=args.verbose)

    # presence/absence heatmap
    map_dict = bin2met(args, args.separator)
    meta_dict, unass_dict = read_input(args.input, map_dict, args.verbose)
    species_list = find_top_SGBs(args.top, meta_dict, args.verbose)

    df1 = pd.DataFrame(0, index=species_list, columns=meta_dict.keys())

    for x in species_list:
        for y in meta_dict.keys():
            if x in meta_dict[y]:
                df1.at[x, y] = 1

    if args.verbose:
        info('Writing to output file {}_pres_abs.{}\n'.format(args.output, args.f))

    output_file = open(args.output + '_pres_abs.csv', 'w')
    df1.to_csv(args.output + '_pres_abs.csv')
    output_file.close()

    if args.verbose:
        info('Drawing image {}_pres_abs.{}\n'.format(args.output, args.f))

    myColors = ((0., 0.0, 0.0, 0.0), (0.00, 0.00, 0.55))
    cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))

    hm1 = sns.clustermap(df1, cmap=cmap, linewidths=.1, square=True, linecolor='Black', figsize=(20, 8))

    colorbar = hm1.ax_heatmap.collections[0].colorbar
    colorbar.set_ticks([0, 1])
    colorbar.set_ticklabels(['Absent', 'Present'])

    plt.setp(hm1.ax_col_dendrogram, visible=False)
    plt.setp(hm1.ax_row_dendrogram, visible=False)
    plt.setp(hm1.ax_heatmap.yaxis, ticks_position='none')
    plt.setp(hm1.ax_heatmap.xaxis, ticks_position='none')

    hm1.savefig(args.output + '_pres_abs.' + args.f, dpi=args.dpi)

    # counts-heatmap: unassigned, uSGB, and kSGB
    xlabels_clustered = [x.get_text() for x in list(plt.getp(hm1.ax_heatmap.xaxis, 'ticklabels'))]
    df2 = pd.DataFrame(0, index=['unassigned', 'uSGBs', 'kSGBs'], columns=xlabels_clustered)

    for x in xlabels_clustered:
        for md in meta_dict[x]:
            if md.startswith('k'):
                df2.at['kSGBs', x] += 1
            elif md.startswith('u'):
                df2.at['uSGBs', x] += 1

            df2.at['unassigned', x] = len(unass_dict[x])

    if args.verbose:
        info('Writing to output file {}_counts.{}\n'.format(args.output, args.f))

    output_file = open(args.output + '_counts.csv', 'w')
    df2.to_csv(args.output + '_counts.csv')
    output_file.close()

    if args.verbose:
        info(('Drawing image {}_counts.{}').format(args.output, args.f))

    hm2 = sns.clustermap(df2, cmap='Oranges', linewidths=.1, figsize=(40, 10), row_cluster=False, col_cluster=False, square=True,
                         cbar_kws={'orientation': 'horizontal', 'label': 'Number of SGBs'})

    colorbar = hm2.ax_heatmap.collections[0].colorbar
    max_ticks_value = max(df2.max()) + 1
    colorbar.set_ticks(range(0, max_ticks_value, math.ceil(max_ticks_value / 10)))

    plt.setp(hm2.ax_col_dendrogram, visible=False)
    plt.setp(hm2.ax_row_dendrogram, visible=False)
    plt.setp(hm2.ax_heatmap.yaxis, ticks_position='none')
    plt.setp(hm2.ax_heatmap.xaxis, ticks_position='none')

    hm2.savefig(args.output + '_counts.' + args.f, dpi=args.dpi)


if __name__ == '__main__':
    t0 = time.time()
    phylophlan2_draw_metagenomic()
    t1 = time.time()
    info('Total elapsed time {}s\n'.format(int(t1 - t0)), init_new_line=True)
    sys.exit(0)
