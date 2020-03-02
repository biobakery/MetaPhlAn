#!/usr/bin/env python3


__author__ = ('Francesco Asnicar (f.asnicar@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it), '
              'Mattia Bolzan (mattia.bolzan@unitn.it), '
              'Nicola Segata (nicola.segata@unitn.it)')
__version__ = '0.14'
__date__ = '6 November 2019'


import sys
import bz2
import os
import argparse as ap
from urllib.request import urlretrieve
import time
import ftplib


if sys.version_info[0] < 3:
    raise Exception("PhyloPhlAn2 requires Python 3, your current Python version is {}.{}.{}"
                    .format(sys.version_info[0], sys.version_info[1], sys.version_info[2]))

DB_TYPE_CHOICES = ['n', 'a']
DOWNLOAD_URL = "https://bitbucket.org/nsegata/phylophlan/downloads/"
GB_ASSEMBLY_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"
TAXA2GENOMES_FILE = "taxa2genomes.txt"
GB_ASSEMBLY_FILE = "assembly_summary_genbank.txt"


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
    p = ap.ArgumentParser(description=("The phylophlan2_get_reference.py script allows to download a specified number (-n/--how_many) of "
                                       "reference genomes from the Genbank repository. Special case \"all\" allows to download a specified "
                                       "number of reference genomes for all available taxonomic species. With the -l/--list_clades params "
                                       "the phylophlan2_get_reference.py scripts returns the list of all species in the database"),
                          formatter_class=ap.ArgumentDefaultsHelpFormatter)

    group = p.add_mutually_exclusive_group()
    group.add_argument('-g', '--get', type=str,
                       help=('Specify the taxonomic label for which download the set of reference genomes. '
                             'The label must represent a valid taxonomic level or the special case "all"'))
    group.add_argument('-l', '--list_clades', action='store_true', default=False,
                       help='Print for all taxa the total number of species and reference genomes available')

    p.add_argument('--database_update', action='store_true', default=False, help="Update the databases file")
    p.add_argument('-e', '--output_file_extension', type=str, default='.fna.gz',
                   help="Specify path to the extension of the output files")
    p.add_argument('-o', '--output', type=str,
                   help="Specify path to the output folder where to save the files, required when -g/--get is specified")
    p.add_argument('-n', '--how_many', type=int, default=4,
                   help='Specify how many reference genomes to download, where -1 stands for "all available"')
    p.add_argument('-m', '--genbank_mapping', type=str, default=GB_ASSEMBLY_FILE,
                   help='The local GenBank mapping file, if not found it will be automatically downloaded')
    p.add_argument('--verbose', action='store_true', default=False, help="Prints more stuff")
    p.add_argument('-v', '--version', action='version',
                   version='phylophlan2_get_reference.py version {} ({})'.format(__version__, __date__),
                   help="Prints the current phylophlan2_get_reference.py version and exit")
    return p.parse_args()


def check_params(args, verbose=False):
    if args.database_update:
        database_update(update=args.database_update, verbose=args.verbose)
        args.database_update = False

    if not args.get and not args.list_clades:
        error('either -g/--get or -l/--list must be specified', init_new_line=True, exit=True)

    if args.list_clades:
        return

    if args.get[:3] not in ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__', 'all']:
        error('Taxonomic label provided "{}" is not in the correct format, it should starts with either: '
              '"k__", "p__", "c__", "o__", "f__", "g__", "s__", or "all"'.format(args.get), exit=True)

    if not args.output:
        error('-o/--output is required', exit=True)

    if os.path.exists(args.output):
        if not os.path.isdir(args.output):
            error('output param is not a directory', exit=True)

    if args.how_many < 0:
        args.how_many = None

    if not args.output_file_extension.startswith('.'):
        args.output_file_extension = '.' + args.output_file_extension

    if args.output_file_extension.endswith('.'):
        args.output_file_extension = args.output_file_extension[:-1]

    if verbose:
        info('Arguments: {}\n'.format(vars(args)))


def create_folder(output, verbose=False):
    if not os.path.exists(output):
        if verbose:
            info('Creating output folder "{}"\n'.format(output))

        os.mkdir(output, mode=0o775)
    elif verbose:
        info('Output folder "{}" present\n'.format(output))


def database_update(update=False, verbose=False):
    taxa2genomes_file_latest = None
    download(os.path.join(DOWNLOAD_URL, TAXA2GENOMES_FILE), TAXA2GENOMES_FILE, overwrite=update, verbose=verbose)

    with open(TAXA2GENOMES_FILE) as f:
        for r in f:
            if not r.startswith('#'):
                taxa2genomes_file_latest = r.strip()
                break  # file should contains only one line, i.e., the name of the latest taxa2genomes file

    download(os.path.join(DOWNLOAD_URL, taxa2genomes_file_latest), taxa2genomes_file_latest, overwrite=update, verbose=verbose)

    return taxa2genomes_file_latest


def byte_to_megabyte(byte):
    """
    Convert byte value to megabyte
    """

    return (byte / 1048576)


class ReportHook():

    def __init__(self):
        self.start_time = time.time()

    def report(self, blocknum, block_size, total_size):
        """
        Print download progress message
        """

        if blocknum == 0:
            self.start_time = time.time()

            if total_size > 0:
                info("Downloading file of size: {:.2f} MB\n".format(byte_to_megabyte(total_size)))
        else:
            total_downloaded = blocknum * block_size
            status = "{:3.2f} MB ".format(byte_to_megabyte(total_downloaded))

            if total_size > 0:
                percent_downloaded = total_downloaded * 100.0 / total_size
                # use carriage return plus sys.stderr to overwrite stderr
                download_rate = total_downloaded / (time.time() - self.start_time)
                estimated_time = (total_size - total_downloaded) / download_rate
                estimated_minutes = int(estimated_time / 60.0)
                estimated_seconds = estimated_time - estimated_minutes * 60.0
                status += ("{:3.2f} %  {:5.2f} MB/sec {:2.0f} min {:2.0f} sec "
                           .format(percent_downloaded, byte_to_megabyte(download_rate),
                                   estimated_minutes, estimated_seconds))

            status += "        \r"
            info(status)


def download(url, download_file, overwrite=False, verbose=False):
    """
    Download a file from a url
    """

    if (not os.path.isfile(download_file)) or overwrite:
        try:
            if verbose:
                info('Downloading "{}" to "{}"\n'.format(url, download_file))

            urlretrieve(url, download_file, reporthook=ReportHook().report)
            info('\n')
        except EnvironmentError:
            error('unable to download "{}"'.format(url), exit=True)
    elif verbose:
        info('File "{}" present\n'.format(download_file))


def retrieve_refseq_url(gcx_id):
    refseq_base_ftp_url = 'ftp.ncbi.nlm.nih.gov'
    refseq_genomes_url = 'genomes/all'

    gcx, number = gcx_id.split('.')[0].split('_')
    gcx_url = '/'.join([gcx] + [number[i:i + 3] for i in range(0, len(number), 3)])

    ftp = ftplib.FTP(refseq_base_ftp_url)
    _ = ftp.login()
    _ = ftp.cwd(refseq_genomes_url + '/' + gcx_url)
    folder = ftp.nlst()[0]
    _ = ftp.cwd(folder)
    files = ftp.nlst()
    _ = ftp.quit()

    url = None

    for ff in files:
        if (folder + '_genomic.fna.gz') == ff:
            url = 'https://' + '/'.join([refseq_base_ftp_url, refseq_genomes_url, gcx_url, folder, ff])
            break

    return url


def list_available_clades(taxa2proteomes_file, verbose=False):
    clades = {}
    metadata = None

    for r in bz2.open(taxa2proteomes_file, 'rt'):
        if r.startswith('#'):
            metadata = r.strip()
            continue

        taxa = r.strip().split('\t')[1].split('|')
        num_ref_gen = len(r.strip().split('\t')[-1].split(';'))

        for i, c in enumerate(taxa):
            cl = '|'.join(taxa[:i + 1])

            if cl in clades:
                num_spp, num_ref = clades[cl]
                clades[cl] = (num_spp + 1, num_ref + num_ref_gen)
            else:
                clades[cl] = (1, num_ref_gen)

    info('#taxa\tspecies\treference_genomes\n')
    for k, v in sorted(clades.items(), key=lambda x: x[0]):
        info('\t'.join([k, str(v[0]), str(v[1])]) + '\n')


def get_reference_genomes(gb_assembly_file, taxa2genomes_file, taxa_label, num_ref, out_file_ext, output, update, verbose=False):
    core_genomes = {}
    metadata = None

    download(GB_ASSEMBLY_URL, gb_assembly_file, overwrite=update, verbose=verbose)

    # load GenBank assembly summary
    gb_assembly_summary = dict([(r.strip().split('\t')[0],
                                 (r.strip().split('\t')[19].replace('ftp://', 'https://') + '/' +
                                  r.strip().split('\t')[19].split('/')[-1] + '_genomic.fna.gz'))
                                for r in open(gb_assembly_file) if not r.startswith('#')])

    # taxa2genomes format
    #
    #   # taxid   taxonomy      genomes_list
    #   12345     k__|...|s__   UP123;UP456;...

    for r in bz2.open(taxa2genomes_file, 'rt'):
        if r.startswith('#'):
            metadata = r.strip()
            continue

        r_clean = r.strip().split('\t')

        if (taxa_label in r_clean[1].split('|')) or (taxa_label == 'all'):
            core_genomes[r_clean[1]] = [(g.split('.')[0],
                                         gb_assembly_summary[g] if g in gb_assembly_summary else retrieve_refseq_url(g))
                                        for g in r_clean[2].split(';')[:num_ref]]

    if not len(core_genomes):
        error('no reference genomes found for "{}", please check the taxonomic label provided'.format(taxa_label),
              exit=True)

    for lbl, genomes in core_genomes.items():
        if verbose:
            info('Downloading {} reference genomes for {}\n'.format(len(genomes), lbl))

        for t in genomes:
            genome, url = t

            if url:
                download(url, os.path.join(output, genome + out_file_ext), verbose=verbose)
            else:
                error('no URL found for "{}"'.format(genome))


def phylophlan2_get_reference():
    taxa2genomes_file_latest = None
    args = read_params()

    if args.verbose:
        info('phylophlan2_get_reference.py version {} ({})\n'.format(__version__, __date__))
        info('Command line: {}\n\n'.format(' '.join(sys.argv)), init_new_line=True)

    check_params(args, verbose=args.verbose)
    taxa2genomes_file_latest = database_update(update=args.database_update, verbose=args.verbose)

    if args.list_clades:
        list_available_clades(taxa2genomes_file_latest, verbose=args.verbose)
        sys.exit(0)

    create_folder(os.path.join(args.output), verbose=args.verbose)
    get_reference_genomes(args.genbank_mapping, taxa2genomes_file_latest, args.get, args.how_many,
                          args.output_file_extension, args.output, args.database_update, verbose=args.verbose)


if __name__ == '__main__':
    t0 = time.time()
    phylophlan2_get_reference()
    t1 = time.time()
    info('Total elapsed time {}s\n'.format(int(t1 - t0)), init_new_line=True)
    sys.exit(0)
