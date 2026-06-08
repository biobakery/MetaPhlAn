#!/usr/bin/env python
__author__ = 'Aitor Blanco (aitor.blancomiguez@unitn.it'
__version__ = '4.2.4'
__date__ = '21 Oct 2025'

import os
import glob
import time
import argparse as ap

try:
    from .util_fun import info, error
except ImportError:
    from util_fun import info, error


def read_params():
    """Reads and parses the command line arguments of the script."""
    p = ap.ArgumentParser(description="", formatter_class=ap.ArgumentDefaultsHelpFormatter)
    p.add_argument('-i', '--input', type=str, default=None, help='The input profile')
    p.add_argument('-o', '--output', type=str, default=None, help='The output profile')
    p.add_argument('--gtdb_assignment_file', type=str, default=None,
                   help='Optional GTDB assignment TSV to use instead of auto-detecting it')
    p.add_argument('--merged_profiles', action='store_true', default=False,
                   help='Specify this when the input is a merged MetaPhlAn profile table')
    return p.parse_args()


def check_params(args):
    """Checks the mandatory command line arguments of the script."""
    if not args.input:
        error('-i (or --input) must be specified', exit=True)
    if not args.output:
        error('-o (or --output) must be specified', exit=True)


def get_gtdb_assignment_file(mpa_profile, gtdb_assignment_file=None):
    """Resolve the GTDB assignment TSV for the MetaPhlAn index used by a profile."""
    if gtdb_assignment_file:
        return gtdb_assignment_file

    with open(mpa_profile, 'r') as rf:
        header = rf.readline().strip()

    if not header.startswith('#mpa_'):
        error('Could not infer the MetaPhlAn index from the first header line', exit=True)

    mpa_index = header[1:]
    search_dir = os.path.dirname(os.path.abspath(__file__))
    matches = sorted(glob.glob(os.path.join(search_dir, '{}_SGB2GTDB*.tsv'.format(mpa_index))), reverse=True)

    if not matches:
        error('Could not find an SGB2GTDB assignment file for MetaPhlAn index "{}" in {}'
              .format(mpa_index, search_dir), exit=True)

    return matches[0]


def load_sgb2gtdb(gtdb_assignment_file):
    """Loads the SGB-to-GTDB mapping."""
    sgb2gtdb = {}
    with open(gtdb_assignment_file, 'r') as read_file:
        for line in read_file:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                sgb2gtdb[parts[0]] = parts[1]
    return sgb2gtdb


def get_parent_taxon(taxon):
    """Gets the direct parent taxon in a semicolon-separated taxonomy."""
    parts = taxon.split(';')
    if len(parts) > 1:
        return ';'.join(parts[:-1])
    return parts[0]


def add_abundance(current_value, new_value, merged_profiles):
    """Adds scalar or vector abundances depending on input profile type."""
    if not merged_profiles:
        return current_value + new_value
    return [a + b for a, b in zip(current_value, new_value)]


def format_abundance(value, merged_profiles):
    """Formats abundance values for output writing."""
    if not merged_profiles:
        return '{:.5f}'.format(value)
    return '\t'.join('{:.5f}'.format(x) for x in value)


def parse_abundance(fields, merged_profiles):
    """Parses abundance fields from a profile line."""
    if not merged_profiles:
        return float(fields[2])
    return [float(x) for x in fields[1:]]


def normalize_abundances(abundances, merged_profiles, unclassified_fraction=0):
    """Normalize abundances to sum to 100 using the input profile shape."""
    if not abundances:
        return abundances

    def normalize_vector(values, target):
        total = sum(values)
        rounded = [round(target * value / total, 5) for value in values] if total else list(values)
        correction = round(target - sum(rounded), 5)
        if rounded:
            rounded[-1] = round(rounded[-1] + correction, 5)
        return rounded

    if merged_profiles:
        if not unclassified_fraction:
            unclassified_fraction = [0.0] * len(next(iter(abundances.values())))

        targets = [100 - x for x in unclassified_fraction]
        columns = list(zip(*abundances.values()))
        normalized_columns = [normalize_vector(list(column), target) for column, target in zip(columns, targets)]
        return {
            taxon: [column[i] for column in normalized_columns]
            for i, taxon in enumerate(abundances)
        }

    if isinstance(unclassified_fraction, list):
        unclassified_fraction = unclassified_fraction[0] if unclassified_fraction else 0

    items = list(abundances.items())
    target = 100 - unclassified_fraction
    normalized_values = normalize_vector([value for _, value in items], target)
    return {taxon: value for (taxon, _), value in zip(items, normalized_values)}


def get_gtdb_profile(mpa_profile, gtdb_profile, merged_profiles=False, gtdb_assignment_file=None):
    """Creates the GTDB profile from a MetaPhlAn one."""
    tax_levels = ['d', 'p', 'c', 'o', 'f', 'g', 's']

    gtdb_assignment_file = get_gtdb_assignment_file(mpa_profile, gtdb_assignment_file)
    info('Using GTDB assignment file: {}'.format(gtdb_assignment_file))
    sgb2gtdb = load_sgb2gtdb(gtdb_assignment_file)

    with open(gtdb_profile, 'w') as wf:
        with open(mpa_profile, 'r') as rf:
            abundances = {x: {} for x in tax_levels}
            unclassified_fraction = 0
            for line in rf:
                if line.startswith('#mpa_'):
                    line = line.strip()
                    wf.write(line + '\t' + gtdb_assignment_file.split('_')[-1].split('.')[0] + '\n')
                    if not merged_profiles:
                        wf.write('#clade_name\trelative_abundance\n')
                elif line.startswith('#'):
                    # For merged profiles preserve additional headers
                    if merged_profiles:
                        wf.write(line)
                elif line.startswith('clade_name'):
                    if merged_profiles:
                        wf.write(line)
                elif line.startswith('UNCLASSIFIED'):
                    fields = line.strip().split('\t')
                    if not merged_profiles:
                        unclassified_fraction = float(fields[2])
                        wf.write('UNCLASSIFIED\t{}\n'.format(float(fields[2])))
                    else:
                        unclassified_fraction = [float(x) for x in fields[1:]]
                        wf.write('UNCLASSIFIED\t{}\n'.format('\t'.join(fields[1:])))
                elif 't__SGB' in line:
                    fields = line.strip().split('\t')
                    sgb_id = fields[0].split('|')[-1][3:]
                    gtdb_tax = sgb2gtdb.get(sgb_id, None)
                    if gtdb_tax is None:
                        continue

                    abundance = parse_abundance(fields, merged_profiles)
                    if gtdb_tax not in abundances['s']:
                        if not merged_profiles:
                            abundances['s'][gtdb_tax] = 0.0
                        else:
                            abundances['s'][gtdb_tax] = [0.0] * len(abundance)
                    abundances['s'][gtdb_tax] = add_abundance(
                        abundances['s'][gtdb_tax], abundance, merged_profiles)

        abundances['s'] = normalize_abundances(abundances['s'], merged_profiles, unclassified_fraction)

        tax_levels_rev = list(reversed(tax_levels))
        for i, tax_level in enumerate(tax_levels_rev[:-1]):
            for tax in abundances[tax_level]:
                parent_tax = get_parent_taxon(tax)
                parent_level = tax_levels_rev[i + 1]
                if parent_tax not in abundances[parent_level]:
                    if not merged_profiles:
                        abundances[parent_level][parent_tax] = 0.0
                    else:
                        n_cols = len(abundances[tax_level][tax])
                        abundances[parent_level][parent_tax] = [0.0] * n_cols
                abundances[parent_level][parent_tax] = add_abundance(
                    abundances[parent_level][parent_tax],
                    abundances[tax_level][tax],
                    merged_profiles
                )

        for tax_level in tax_levels:
            for tax in abundances[tax_level]:
                wf.write('{}\t{}\n'.format(
                    tax,
                    format_abundance(abundances[tax_level][tax], merged_profiles)
                ))


def main():
    t0 = time.time()
    args = read_params()
    info("Start execution")
    check_params(args)
    get_gtdb_profile(
        args.input,
        args.output,
        merged_profiles=args.merged_profiles,
        gtdb_assignment_file=args.gtdb_assignment_file
    )
    exec_time = time.time() - t0
    info("Finish execution ({} seconds)".format(round(exec_time, 2)))


if __name__ == "__main__":
    main()
