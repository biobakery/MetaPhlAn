#!/usr/bin/env python3


__author__ = ('Francesco Asnicar (f.asnicar@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it), '
              'Claudia Mengoni (claudia.mengoni@studenti.unitn.it), '
              'Mattia Bolzan (mattia.bolzan@unitn.it), '
              'Nicola Segata (nicola.segata@unitn.it)')
__version__ = '0.17'
__date__ = '6 November 2019'


import sys
import bz2
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import argparse as ap
from urllib.request import urlretrieve
from urllib.parse import urlencode
from urllib.request import Request
from urllib.request import urlopen
import time


if sys.version_info[0] < 3:
    raise Exception("PhyloPhlAn2 requires Python 3, your current Python version is {}.{}.{}"
                    .format(sys.version_info[0], sys.version_info[1], sys.version_info[2]))

DB_TYPE_CHOICES = ['n', 'a']
GENOME_EXTENSION = '.fna'
PROTEOME_EXTENSION = '.faa'
DOWNLOAD_URL = "https://bitbucket.org/nsegata/phylophlan/downloads/"
TAXA2CORE_FILE = "taxa2core.txt"


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
    p = ap.ArgumentParser(description=("The phylophlan2_setup_database.py script can be used to either format an input folder or "
                                       "multi-fasta file to be used as database in phylophlan2.py, or automatically download a "
                                       "pre-identified set of core UniRef90 proteins for the taxonomic label of a given species"),
                          formatter_class=ap.ArgumentDefaultsHelpFormatter)

    group = p.add_mutually_exclusive_group()
    group.add_argument('-i', '--input', type=str,
                       help=("Specify the path to either the folder containing the marker files or the file of markers, in "
                             "(multi-)fasta format"))
    group.add_argument('-g', '--get_core_proteins', type=str, default=None,
                       help=('Specify the taxonomic label for which download the set of core proteins. The label must represent a '
                             'species: "--get_core_proteins s__Escherichia_coli"'))

    p.add_argument('--database_update', action='store_true', default=False, help="Update the databases file")
    p.add_argument('-o', '--output', type=str, default=None, help="Specify path to the output folder where to save the database")
    p.add_argument('-d', '--db_name', type=str, help="Specify the name of the output database")
    p.add_argument('-e', '--input_extension', type=str, default=None,
                   help="Specify the extension of the input file(s) specified via -i/--input")
    p.add_argument('-t', '--db_type', default=None, choices=DB_TYPE_CHOICES,
                   help='Specify the type of the database, where "n" stands for nucleotides and "a" for amino acids')
    p.add_argument('-x', '--output_extension', type=str, default=None, help="Set the database output extension")
    p.add_argument('--overwrite', action='store_true', default=False, help="If specified and the output file exists it will be overwritten")
    p.add_argument('--verbose', action='store_true', default=False, help="Prints more stuff")
    p.add_argument('-v', '--version', action='version',
                   version='phylophlan2_setup_database.py version {} ({})'.format(__version__, __date__),
                   help="Prints the current phylophlan2_setup_database.py version and exit")
    return p.parse_args()


def check_params(args, verbose=False):
    if args.database_update:
        database_update(update=args.database_update, verbose=args.verbose)
        args.database_update = False

    if not args.input and not args.get_core_proteins:
        error('either -i/--input or -g/--get_core_proteins must be specified', init_new_line=True, exit=True)

    if args.input:
        if (not os.path.isdir(args.input)) and (not os.path.isfile(args.input)):
            error('input must be either a folder or a file', exit=True)

        if not args.db_name:
            error('-d/--db_name must be specified', exit=True)

        if os.path.isdir(args.input):
            if not args.input_extension:
                error('input is a folder, hence --input_extension must be specified', exit=True)

    if args.get_core_proteins:
        # not handling keys different than species (Tue 19 Jun 2018)
        if args.get_core_proteins[:3] != 's__':
            error('The taxonomic label provided "{}" does not starts with "s__"'.format(args.get_core_proteins),
                  exit=True)

        if not args.output:
            args.output = args.get_core_proteins

            if verbose:
                info('Setting output folder "{}"\n'.format(args.output))
        elif os.path.isdir(args.output) and (not args.output.endswith(args.get_core_proteins)):
            args.output = os.path.join(args.output, args.get_core_proteins)

            if verbose:
                info('Setting output folder "{}"\n'.format(args.output))

        args.input = args.output
        args.input_extension = PROTEOME_EXTENSION
        args.output_extension = PROTEOME_EXTENSION

        if not args.db_name:
            args.db_name = args.get_core_proteins

            if verbose:
                info('Setting db_name to "{}"\n'.format(args.db_name))

    if not args.db_type:
        if not args.output_extension:
            error('either -t/--db_type or -x/--output_extension must be specified', exit=True)
        else:
            if not args.output_extension.startswith('.'):
                args.output_extension = '.' + args.output_extension

            if args.output_extension.endswith('.'):
                args.output_extension = args.output_extension[:-1]
    else:
        if not args.output_extension:
            if 'n' == args.db_type:
                args.output_extension = GENOME_EXTENSION
            if 'a' == args.db_type:
                args.output_extension = PROTEOME_EXTENSION

            if verbose:
                info('Setting output extension to "{}"\n'.format(args.output_extension))
        else:
            error("both -t/--db_type and -x/--output_extension were specified, don't know which one to use!",
                  exit=True)

    if not args.output:
        args.output = args.input

        if os.path.isdir(args.input):
            args.output = os.path.dirname(args.input)

        if verbose:
            info('Output folder not specified, setting to "{}"\n'.format(args.output))

    if verbose:
        info('Arguments: {}\n'.format(vars(args)))


def database_update(update=False, verbose=False):
    taxa2core_file_latest = None
    download(os.path.join(DOWNLOAD_URL, TAXA2CORE_FILE), TAXA2CORE_FILE, overwrite=update, verbose=verbose)

    with open(TAXA2CORE_FILE) as f:
        for r in f:
            if not r.startswith('#'):
                taxa2core_file_latest = r.strip()
                break  # file should contains only one line, i.e., the name of the latest taxa2core file

    download(os.path.join(DOWNLOAD_URL, taxa2core_file_latest), taxa2core_file_latest, overwrite=update, verbose=verbose)

    return taxa2core_file_latest


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

    if overwrite or not os.path.isfile(download_file):
        try:
            if verbose:
                info('Downloading "{}" to "{}"\n'.format(url, download_file))

            urlretrieve(url, download_file, reporthook=ReportHook().report)
            info('\n')
        except EnvironmentError:
            error('unable to download "{}"'.format(url))
    elif verbose:
        info('File "{}" present\n'.format(download_file))


def create_folder(output, verbose=False):
    if not os.path.exists(output):
        if verbose:
            info('Creating output folder "{}"\n'.format(output))

        os.mkdir(output, mode=0o775)
    elif verbose:
        info('Output "{}" folder present\n'.format(output))


def get_core_proteins(taxa2core_file, taxa_label, output, output_extension, verbose=False):
    core_proteins = {}
    url = None
    metadata = None
    retry2download = []
    not_mapped = []
    not_mapped_again = []

    # taxa2core format
    #
    #   # taxid   taxonomy      UniRefURL       core_list
    #   12345     k__|...|s__   http://..{}..   UP123;UP456;...

    for r in bz2.open(taxa2core_file, 'rt'):
        if r.startswith('#'):
            metadata = r.strip()
            continue

        r_clean = r.strip().split('\t')

        if taxa_label in r_clean[1]:
            if taxa_label == r_clean[1].split('|')[-1]:
                url = r_clean[2]
                core_proteins[r_clean[1]] = r_clean[3].split(';')
                break
            else:
                core_proteins[r_clean[1]] = None

    if not len(core_proteins):
        error('no entry found for "{}", please check the taxonomic label provided'.format(taxa_label),
              exit=True)
    elif len([k for k, v in core_proteins.items() if v is not None]) > 1:
        error('{} entries found for "{}":\n{}    please check the taxonomic label provided'
              .format(len(core_proteins), taxa_label, '    - {}\n'.join(core_proteins.keys())), exit=True)

    for lbl, core_prots in core_proteins.items():
        if core_prots is None:
            continue

        if verbose:
            info('Downloading {} core proteins for {}\n'.format(len(core_prots), lbl))

        for core_prot in core_prots:
            local_prot = os.path.join(output, core_prot + output_extension)
            download(url.format(core_prot), local_prot, verbose=verbose)

            if not os.path.exists(local_prot):
                retry2download.append(core_prot)

    # try to re-map the ids in case "not mapped" store in not_mapped
    if retry2download:
        if verbose:
            info("Re-trying to download {} core proteins that just failed, please wait as it might take some time\n"
                 .format(len(retry2download)))
        idmapping_url = 'https://www.uniprot.org/uploadlists/'
        contact = "phylophlan@cibiocm.com"
        params = {'from': 'ACC+ID',
                  'to': 'NF90',
                  'format': 'tab',  # or 'list' for only converted clusters
                  'query': ' '.join(retry2download)}
        data = urlencode(params).encode('utf-8')
        request = Request(idmapping_url, data, headers={'User-Agent': 'Python {}'.format(contact)})

        try:
            response = urlopen(request, data)
            uniprotkb2uniref90 = [line.decode().strip().split('\t')[:2] for line in response.readlines()][1:]  # skip ['From', 'To']
        except Exception:
            error('unable convert UniProtKB ID to UniRef90 ID')

        for uniref90_id in (x[1].split('_')[-1] for x in uniprotkb2uniref90[1:]):
            local_prot = os.path.join(output, uniref90_id + output_extension)
            download(url.format(uniref90_id), local_prot, verbose=verbose)

            if not os.path.exists(local_prot):
                not_mapped.append(uniref90_id)

        if len(uniprotkb2uniref90) != len(retry2download):
            # probably deleted proteins in the Uniprot versions, try to download their latest version in any case
            for ur90 in set([ur90 for ur90, _ in uniprotkb2uniref90]) - set(retry2download):
                local_prot = os.path.join(output, ur90 + output_extension)
                download('https://www.uniprot.org/uniprot/{}.fasta?version=*'.format(ur90), local_prot, verbose=verbose)

                if not os.path.exists(local_prot):
                    not_mapped_again.append(core_prot)

        # really don't know what else to try... I'm sorry!
        if not_mapped_again:
            nd_out = os.path.join(output, taxa_label + '_core_proteins_not_mapped.txt')

            if verbose:
                info('There are {} core proteins that could not be downloaded, writing thier IDs to "{}"\n'
                     .format(len(not_mapped_again), nd_out))

            with open(nd_out, 'w') as f:
                f.write('\n'.join(not_mapped_again) + '\n')


def create_database(db_name, inputt, input_ext, output, overwrite, verbose=False):
    seqs = []

    if os.path.exists(output) and (not overwrite):
        error('output file exists and --overwrite not specified', exit=True)

    if os.path.isdir(inputt):
        for marker in glob.iglob(os.path.join(inputt, '*' + input_ext + '*')):
            seqs += [SeqRecord(record.seq,
                               id='_'.join([db_name.replace('_', '-').replace(':', ''),
                                            record.id.replace('_', '-').replace(',', '-').replace(':', ''),
                                            str(count)]),
                               description='')
                     for count, record in enumerate(SeqIO.parse(bz2.open(marker, 'rt') if marker.endswith('.bz2')
                                                                else open(marker), "fasta"))]
    else:
        seqs = [SeqRecord(record.seq,
                          id='_'.join([db_name.replace('_', '-').replace(':', ''),
                                       record.id.replace('_', '-').replace(',', '-').replace(':', ''),
                                       str(count)]),
                          description='')
                for count, record in enumerate(SeqIO.parse(bz2.open(inputt, 'rt') if inputt.endswith('.bz2')
                                                           else open(inputt), "fasta"))]

    if not seqs:
        error('no sequences found, make sure the input folder/file provided is not empty', exit=True)

    if verbose:
        info('Writing output database "{}"\n'.format(output))

    with open(output, 'w') as f:
        SeqIO.write(seqs, f, "fasta")


def phylophlan2_setup_database():
    args = read_params()

    if args.verbose:
        info('phylophlan2_setup_database.py version {} ({})\n'.format(__version__, __date__))
        info('Command line: {}\n\n'.format(' '.join(sys.argv)), init_new_line=True)

    check_params(args, verbose=args.verbose)
    create_folder(args.output, verbose=args.verbose)

    if args.get_core_proteins:
        taxa2core_file_latest = database_update(update=args.database_update, verbose=args.verbose)
        get_core_proteins(taxa2core_file_latest, args.get_core_proteins, args.output, args.output_extension, verbose=args.verbose)

    create_database(args.db_name, args.input, args.input_extension, os.path.join(args.output, args.db_name + args.output_extension),
                    args.overwrite, verbose=args.verbose)


if __name__ == '__main__':
    t0 = time.time()
    phylophlan2_setup_database()
    t1 = time.time()
    info('Total elapsed time {}s\n'.format(int(t1 - t0)), init_new_line=True)
    sys.exit(0)
