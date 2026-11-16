#! /usr/bin/env python
__author__ = (
    'Francesco Asnicar (f.asnicar@unitn.it), '
    'Francesco Beghini (francesco.beghini@unitn.it), '
    'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it), '
    'Vitor Heidrich (vitor.heidrich@unitn.it), '
    'Michal Puncochar (michal.puncochar@unitn.it), '
    'Roan Spadazzi (roan.spadazzi@unitn.it), '
    'Duy Tin Truong (duytin.truong@unitn.it), '
    'Moreno Zolfo (moreno.zolfo@unitn.it)'
)
__version__ = '4.3.0'
__date__ = '16 Nov 2026'

import os

# prevent OpenBLAS and co. from creating multiple threads https://pythonspeed.com/articles/concurrency-control/
# we parallelize over samples ourselves
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['GOTO_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['VECLIB_MAXIMUM_THREADS'] = '1'


import argparse as ap
import math
import multiprocessing.pool as mpp
import pathlib
import subprocess as sp
import sys
import traceback
from collections import Counter
from typing import Sequence

from tqdm import tqdm

try:
    import tomllib  # Python 3.11+
except ModuleNotFoundError:
    import tomli as tomllib  # https://github.com/hukkin/tomli

from . import info, fix_hash_seed
from .multistrain.pipeline import run
from .multistrain.prepare_sam import DatabaseMismatch, InputError
from .multistrain.utils import ArgumentType, error, MetaphlanDBInfo
from .database_controller import StrainphlanDatabaseController

DEFAULT_CONFIG = pathlib.Path(__file__).parent / 'multistrain' / 'config-default.toml'


class ArgTypes:
    input: Sequence[pathlib.Path] | None
    input_list: pathlib.Path | None
    output_dir: pathlib.Path
    database: str
    allele_counts_only: bool
    threads: int
    config: pathlib.Path | None
    config_params: list[str] | None
    reuse: str
    save_bam_file: bool
    output_suffix: str
    tmp: pathlib.Path | None


def existing_input(path):
    """An existing file, as given: a symbolic link keeps its own name (the sample name comes from it)."""
    p = pathlib.Path(path).absolute()
    if not p.is_file():
        raise ap.ArgumentTypeError(f'The file does not exist ({path}).')
    return p


def read_params():
    p = ap.ArgumentParser(
        description="StrainPhlAn sample2markers: per-sample allele counts and strain-resolved consensus markers "
                    "(major and, where a second strain is detected, minor) from MetaPhlAn's SAM or BAM output. The "
                    "*.major.json.bz2 / *.minor.json.bz2 outputs are the input of strainphlan.",
        formatter_class=ap.ArgumentDefaultsHelpFormatter)

    input_g = p.add_mutually_exclusive_group(required=True)
    input_g.add_argument('-i', '--input', type=existing_input, nargs='+',
                         help="MetaPhlAn mapping outputs (.sam, .sam.bz2, .bam, .bam.bz2), separated by spaces")
    input_g.add_argument('--input_list', type=ArgumentType.existing_file,
                         help='A text file with one input path per line (as in --input)')
    input_g.add_argument('--print_config', action='store_true',
                         help='Print the configuration (defaults overridden by --config and --config_params) as a commented '
                              'TOML file, and exit')
    p.add_argument('-o', '--output_dir', type=ArgumentType.creatable_dir,
                   help="Output directory (one subfolder per sample)")
    p.add_argument('-d', '--database', type=str, default='latest',
                   help="MetaPhlAn database used for the mapping: a .pkl path, a database name (looked up in the "
                        "MetaPhlAn database folder, $METAPHLAN_DB_DIR), or 'latest'")
    p.add_argument('-t', '-n', '--threads', '--nprocs', dest='threads', type=ArgumentType.positive_int, default=1,
                   help="Number of samples processed in parallel (one process each)")
    p.add_argument('--config', type=ArgumentType.existing_file, default=None,
                   help=f"TOML file with the parameters to change; the others keep their defaults ({DEFAULT_CONFIG})")
    p.add_argument('--config_params', type=str, action='append', default=None, metavar='"KEY=VALUE,..."',
                   help="Parameters to change on the command line, e.g. \"min_base_coverage=3,min_read_pident=95\" "
                        "(the option can be repeated); applied after the defaults and --config. Values as in the "
                        "TOML file (3, 0.5, true, nan; text without quotes)")
    p.add_argument('--tmp', type=ArgumentType.existing_dir, default=None,
                   help="Folder for the temporary uncompressed copy of the kept alignments (one per running sample); "
                        "default: $TMPDIR, else the sample's output folder")
    p.add_argument('--reuse', type=str, default='bam', choices=['none', 'bam'],
                   help="'bam' reuses the preprocessed BAM of an earlier run (kept with --save_bam_file) if present; "
                        "'none' redoes everything")
    p.add_argument('--save_bam_file', action='store_true', default=False,
                   help="Keep the preprocessed BAM file (and its read lengths)")
    p.add_argument('--allele_counts_only', action='store_true',
                   help="Stop before the genotype reconstruction: after the pileup and SNP calling, write the allele "
                        "counts and results.tsv (coverage, polyallelic sites, multi-strain call), skip linkage, model "
                        "fit and consensus markers (no JSON files)")
    p.add_argument('--output_suffix', type=str, default='', help="Suffix appended to the JSON file names")
    p.add_argument('-v', '--version', action='version',
                   version=f"StrainPhlAn sample2markers version {__version__} ({__date__})")
    return p


# ----------------------------------------------------------------------------------------------
# configuration: defaults from config-default.toml, overridden by --config
# ----------------------------------------------------------------------------------------------
def _type_name(v):
    return 'boolean' if isinstance(v, bool) else 'number' if isinstance(v, (int, float)) else type(v).__name__


def _layer_problems(defaults, user):
    """Unknown parameters and values of the wrong type in one configuration layer (file or --config_params)."""
    problems = []
    unknown = sorted(set(user) - set(defaults))
    if unknown:
        problems.append('unknown parameters: ' + ', '.join(unknown))
    for k, v in user.items():
        if k not in defaults:
            continue
        d = defaults[k]
        if k == 'error_rate':
            ok = v in ('phred', 'estimate') or (isinstance(v, (int, float)) and not isinstance(v, bool))
        elif isinstance(d, bool) or isinstance(v, bool):
            ok = isinstance(v, bool) and isinstance(d, bool)
        elif isinstance(d, (int, float)):
            ok = isinstance(v, (int, float))
        else:
            ok = type(v) is type(d)
        if not ok and k == 'error_rate':
            problems.append(f"error_rate = {v!r} should be 'phred', 'estimate' or a number")
        elif not ok:
            problems.append(f'{k} = {v!r} should be a {_type_name(d)} like the default ({d!r})')
    for k in ('min_base_quality', 'min_mapping_quality'):
        v = user.get(k)
        if v is not None and not (isinstance(v, int) and not isinstance(v, bool) and 0 <= v < 2 ** 31):
            problems.append(f'{k} = {v!r} should be a non-negative integer')
    return problems


def parse_config_params(items):
    """--config_params values ("key=value,key=value", the option can be repeated) as a dict. Each value is read as
    a TOML value (3, 0.5, true, nan, 'text'); a bare word is taken as text (error_rate=estimate)."""
    params, problems = {}, []
    for item in items or []:
        for part in item.split(','):
            part = part.strip()
            if not part:
                continue
            if '=' not in part:
                problems.append(f'"{part}" is not key=value')
                continue
            k, v = (x.strip() for x in part.split('=', 1))
            try:
                params[k] = tomllib.loads(f'v = {v}')['v']
            except tomllib.TOMLDecodeError:
                params[k] = v
    return params, problems


def load_config(path, config_params=None):
    """The configuration in effect: config-default.toml, then the `path` TOML file (any subset of the parameters),
    then --config_params. Unknown parameters and values of the wrong type stop the run.
    Returns (config, changed): changed maps every parameter that differs from its default to the source of its
    value (the file's path, or 'the command line')."""
    with open(DEFAULT_CONFIG, 'rb') as f:
        defaults = tomllib.load(f)
    layers = []
    if path is not None:
        with open(path, 'rb') as f:
            layers.append((f'Configuration file {path}', str(path), tomllib.load(f)))
    if config_params:
        params, problems = parse_config_params(config_params)
        if problems:
            error('--config_params: ' + '; '.join(problems), exit=True)
        layers.append(('--config_params', 'the command line', params))

    config, changed = dict(defaults), {}
    for label, source, user in layers:
        problems = _layer_problems(defaults, user)
        if problems:
            error(f'{label}: ' + '; '.join(problems), exit=True)
        for k, v in user.items():
            config[k] = v
            if v == defaults[k] or (_isnan(v) and _isnan(defaults[k])):
                changed.pop(k, None)
            else:
                changed[k] = source
    return config, changed


def describe_changes(config, changed):
    """'N changed by FILE: k = v, ...; M set on the command line: k = v, ...' for the log."""
    by_source = {}
    for k, source in changed.items():
        by_source.setdefault(source, []).append(f'{k} = {_toml_value(config[k])}')
    parts = [f'{len(v)} {"set on" if s == "the command line" else "changed by"} {s}: ' + ', '.join(v)
             for s, v in by_source.items()]
    return 'defaults' + ('; ' + '; '.join(parts) if parts else ' (nothing changed)')


def _isnan(v):
    return isinstance(v, float) and math.isnan(v)


# What each parameter does, by pipeline step (printed as comments by --print_config)
PARAMETER_DOCS = [
    ('Reads (SAM preprocessing)', [
        ('fix_pairend', 'Whether to fix paired-end information from MetaPhlAn output; requires SAM input.'),
        ('dedup_stats', 'Whether to run samtools markdup to mark and count duplicate reads; requires fix_pairend = true.'),
        ('remove_viral_hits', 'Strain reconstruction and tracking of viral SGBs have not been tested.'),
        ('min_mapping_quality', 'Minimum MAPQ.'),
        ('filter_indels', 'Remove indel and clipped alignments.'),
        ('filter_duplicates', 'Works only if dedup_stats is set to true.'),
        ('min_read_pident', 'Minimum percentage identity required to keep a read.'),
    ]),
    ('Bases and sequencing errors', [
        ('min_base_coverage', 'Minimum base coverage.'),
        ('min_base_quality', 'Minimum Phred base quality.'),
        ('seq_error_stats_min_cov', 'Minimum coverage of a marker position for its bases to contribute to sequencing-error statistics.'),
        ('null_err_rate_pseudocount', 'Pseudocount added to the sample-wide null error estimate.'),
        ('null_err_rate_excess_pseudocount', 'Pseudocount added to the error-excess estimate for each read position.'),
        ('max_null_err_rate', 'Set to 0 to turn off the null error estimate and use the Phred sequencing qualities.'),
        ('discard_read_ends', 'Whether to check for an excess error rate at the ends of reads.'),
        ('discard_read_ends_support', 'Threshold for the number of counted bases at a read position to enter the error-excess outlier test.'),
        ('error_rate', 'Expected sequencing error rate of bases with qualities >= min_base_quality; options: phred, estimate, or a float.'),
    ]),
    ('Species presence', [
        ('min_marker_breadth', 'Minimum breadth of coverage (at coverage depth >= 1x) for the marker to be kept.'),
        ('quasi_marker_frac', 'Exclude a quasi-marker if any other SGB associated with it has at least this fraction of its database markers represented in the pileup.'),
        ('min_markers_abs', 'Minimum absolute number of markers required to retain an SGB; combined with min_markers_rel using AND.'),
        ('min_markers_rel', 'Minimum fraction of available database markers required to retain an SGB; combined with min_markers_abs using AND.'),
    ]),
    ('Positions kept for SNP calling', [
        ('outlier_coverage_p', "Tail probability for the negative-binomial coverage filter; exclude positions outside the model's lower/upper coverage bounds."),
        ('trim_marker_ends', 'Number of positions to trim from each end of a marker.'),
        ('outlier_coverage_marker', 'Threshold for the fraction of positions with outlier coverage above which the whole marker is discarded.'),
    ]),
    ('SNPs and strains', [
        ('polymorphism_q_alpha', 'The maximum FDR-corrected p-value (q-value) in the polymorphism test.'),
        ('multistrain_min_markers_abs', 'Minimum number of markers containing a significantly polyallelic position required to classify an SGB as multi-strain.'),
        ('multistrain_min_markers_rel', 'Minimum fraction of all database markers for the SGB that must contain a significantly polyallelic position for multi-strain classification.'),
        ('linkage_graph_merging', 'Whether to force physically linked polyallelic alleles (linked by reads) to stay together in the genotype reconstruction model; recommended off by default.'),
        ('max_w_across', 'Maximum summed read/read-pair support contradicting a proposed merge.'),
        ('min_w_along', 'Minimum summed read/read-pair support for a proposed merge.'),
    ]),
    ('Model fit (SNP rate and strain ratio r)', [
        ('max_iter', 'Maximum number of optimizer iterations.'),
        ('init_r_min', 'Minimum r to use as an initial guess (to avoid boundary problems).'),
        ('init_r_max', 'Maximum r to use as an initial guess (to avoid boundary problems).'),
        ('acceptable_tol', 'Acceptable convergence error tolerance (gradient norm).'),
        ('negative_variance_tol', 'If the calculated variance is negative but close to zero within this tolerance, accept it.'),
        ('init_snp_rate_min', 'Minimum snp_rate to use as an initial guess for the optimization (to avoid boundary problems).'),
        ('snp_rate_fallback', 'SNP rate to use when fitting fails.'),
        ('r_fallback', 'The major strain fraction to use when fitting fails.'),
    ]),
    ('Output markers and the minor-strain release', [
        ('min_output_quality', 'Minimum Phred base quality required to call a base; otherwise, insert a gap (e.g. 10 ~ 90%, 20 ~ 99%, 30 ~ 99.9% confidence).'),
        ('min_output_base_coverage', 'Minimum counted coverage for a position to enter model fitting and genotype reconstruction.'),
        ('min_output_breadth', 'Minimum breadth of coverage required to output a marker [%].'),
        ('min_output_snp_rate', 'Minimum SNP rate required to report the second strain [%].'),
        ('min_recall', 'Minimum recall required to report the second strain [%].'),
        ('min_precision', 'Minimum precision required to report the second strain [%].'),
        ('min_precision_support', 'Minimum precision support required to report the second strain.'),
        ('max_r_var', 'Maximum variance of the strain ratio for the fit to be considered reliable and the second strain to be reported; nan turns off this filter.'),
        ('max_snp_rate_var', 'Maximum variance of the SNP rate for the fit to be considered reliable and the second strain to be reported; nan turns off this filter.'),
    ]),
]


def _toml_value(v):
    if isinstance(v, bool):
        return 'true' if v else 'false'
    if isinstance(v, float) and math.isnan(v):
        return 'nan'
    if isinstance(v, str):
        return f"'{v}'"
    return repr(v)


def config_to_toml(config, changed=None):
    """The configuration as a commented TOML file (usable as --config), grouped by pipeline step; values that differ
    from the default are marked with where they come from and the default."""
    with open(DEFAULT_CONFIG, 'rb') as f:
        defaults = tomllib.load(f)
    changed = changed or {}
    lines = ['# sample2markers configuration',
             f'# defaults: {DEFAULT_CONFIG}',
             f'# in effect: {describe_changes(config, changed)}',
             '# a --config file needs only the parameters to change', '']
    sections = PARAMETER_DOCS + [('Other', [(k, '') for k in config if k not in
                                            {k for _, ps in PARAMETER_DOCS for k, _ in ps}])]
    for title, params in sections:
        params = [(k, d) for k, d in params if k in config]
        if not params:
            continue
        lines += [f'# --- {title} ---']
        for k, doc in params:
            if doc:
                lines.append(f'# {doc}')
            line = f'{k} = {_toml_value(config[k])}'
            if k in changed:
                line += f'  # {"set on" if changed[k] == "the command line" else "changed by"} {changed[k]} ' \
                        f'(default {_toml_value(defaults[k])})'
            lines.append(line)
        lines.append('')
    return '\n'.join(lines)


# ----------------------------------------------------------------------------------------------
def resolve_database(database):
    """A .pkl path, a database name in the MetaPhlAn database folder, or 'latest'."""
    if database != 'latest' and not os.path.isfile(database):
        candidate = os.path.join(StrainphlanDatabaseController(database).default_db_folder, database + '.pkl')
        if os.path.isfile(candidate):
            return candidate
        error(f'The database {database} is neither a file nor {candidate}', exit=True)
    return database


def check_params(argp: ap.ArgumentParser):
    args = argp.parse_args(namespace=ArgTypes())

    if args.print_config:
        return args
    if args.output_dir is None:
        argp.error('the following arguments are required: -o/--output_dir')

    if args.input_list is not None:
        with open(args.input_list) as f:
            paths = [pathlib.Path(line.strip()) for line in f if line.strip()]
        info(f'Found {len(paths)} samples in the specified input file')
        args.input = paths

    names = Counter(sample_name(s) for s in args.input)
    dup = sorted(n for n, c in names.items() if c > 1)
    if dup:
        error('Several inputs have the same sample name, their outputs would overwrite each other: '
              + ', '.join(dup), exit=True)

    if args.tmp is not None:
        os.environ['MS_TMPDIR'] = str(args.tmp)

    args.database = resolve_database(args.database)
    return args


def sample_name(path):
    path = pathlib.Path(path)
    return path.with_suffix('').stem if path.suffix == '.bz2' else path.stem


def check_samtools():
    try:
        samtools_version_output = sp.run("samtools", capture_output=True).stderr.decode().split('\n')
    except FileNotFoundError:
        error('samtools (>= 1.17) is required but was not found in PATH', exit=True)
    for line in samtools_version_output:
        if line.startswith('Version:'):
            samtools_version_info = line
            break
    else:
        error("Couldn't get information about the samtools version", exit=True)

    samtools_version = samtools_version_info.rstrip('\n').split(' ')[1]
    samtools_version_parts = [int(x) for x in samtools_version.split('.')[:2]]
    if samtools_version_parts < [1, 17]:
        error(f'The samtools version required is >= 1.17, installed is {samtools_version}', exit=True)

    info(f'Using samtools version {samtools_version}')


def try_run(args):
    """One sample; any failure (including a SystemExit raised deep inside) is reported and returns False, so a
    worker never takes the pool down."""
    try:
        return run(*args, mp_db_info=try_run.mp_db_info)
    except KeyboardInterrupt:
        raise
    except BaseException as e:
        error(f'Error running sample {args[0]}: {type(e).__name__}: {e}')
        if not isinstance(e, (SystemExit, DatabaseMismatch, InputError)):
            error(traceback.format_exc())
        return False


def main():
    fix_hash_seed()  # reproducible JSON marker order
    argp = read_params()
    args = check_params(argp)

    if args.print_config:
        config, changed = load_config(args.config, args.config_params)
        sys.stdout.write(config_to_toml(config, changed))
        return

    info('Starting sample2markers')
    check_samtools()

    config, changed = load_config(args.config, args.config_params)
    if args.config is not None or args.config_params:
        info(f'Configuration: {describe_changes(config, changed)}')

    missing = [s for s in args.input if not s.is_file()]
    for s in missing:
        error(f'Input file not found, skipping: {s}')
    samples = [s for s in args.input if s.is_file()]

    mp_db_controller = StrainphlanDatabaseController(args.database)
    mp_db_controller.load_database()
    info('Getting information from the MetaPhlAn database')
    mp_db_info = MetaphlanDBInfo.from_mp_controller(mp_db_controller)

    ss_args = [(sample_path, args.output_dir, config, args.save_bam_file, args.reuse, args.output_suffix,
                args.allele_counts_only)
               for sample_path in samples]

    n_threads = max(1, min(args.threads, len(ss_args)))

    info(f'Running on {len(ss_args)} samples using {n_threads} processes')

    # bind constant data to the function so that they are shared on fork (copy-on-write)
    try_run.mp_db_info = mp_db_info

    if n_threads == 1:
        successes = [try_run(ss_arg) for ss_arg in ss_args]
    else:
        with mpp.Pool(n_threads) as pool:
            successes = list(tqdm(pool.imap(try_run, ss_args), total=len(ss_args)))

    failed = [sample_name(s) for s, ok in zip(samples, successes) if not ok] + [sample_name(s) for s in missing]
    if not failed:
        info(f'Successfully finished sample2markers, results are stored at {args.output_dir}')
    else:
        error(f'Finished sample2markers; {len(failed)} of {len(args.input)} samples failed: ' + ', '.join(failed))
        sys.exit(1)


if __name__ == '__main__':
    main()
