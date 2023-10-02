#!/usr/bin/env python
__author__ = ('Aitor Blanco Miguez (aitor.blancomiguez@unitn.it), '
              'Duy Tin Truong (duytin.truong@unitn.it), '
              'Francesco Asnicar (f.asnicar@unitn.it), '
              'Moreno Zolfo (moreno.zolfo@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it),'
              'Michal Puncochar (michal.puncochar@unitn.it')
__version__ = '4.1.0'
__date__ = '23 Aug 2023'


import argparse as ap
import collections
import io
import os
import re
import tempfile
import time
from shutil import copyfile, rmtree
from typing import Iterable

import numpy as np
import pandas as pd
from Bio import SeqIO, Seq

try:
    from .utils import *
    from .utils.external_exec import run_command
except ImportError:
    from utils import *
    from utils.external_exec import run_command


class Strainphlan:

    @staticmethod
    def sample_path_to_name(sample_path):
        if sample_path.endswith('.bz2'):
            sample_path = sample_path[:-len('.bz2')]

        sample_file = os.path.basename(sample_path)
        name, ext = os.path.splitext(sample_file)  # ext can be pkl, json, fna, ...
        return name


    def get_markers_matrix_from_samples(self):
        """Gets a binary matrix representing the presence/absence of the clade markers in the uploaded samples

        Returns:
            list: the list containing the samples-to-markers information
        """

        markers_matrix = execute_pool(((Strainphlan.get_matrix_for_sample, sample, self.clade_markers_names,
                                        self.breadth_thres) for sample in self.samples), self.nprocs)
        return markers_matrix


    @staticmethod
    def get_matrix_for_sample(sample_path, clade_markers, breadth_thres):
        """Returns the matrix with the presence / absence of the clade markers in a samples

        Args:
            sample_path (str): the path to the sample
            clade_markers (Iterable): the list with the clade markers names
            breadth_thres:

        Returns:
            dict: dictionary containing the sample-to-markers information as a binary matrix
        """
        sample = ConsensusMarkers.from_file(sample_path)

        markers = {"sample": sample_path}
        markers.update({m: 0 for m in clade_markers})
        markers.update({marker.name: 1 for marker in sample.consensus_markers
                        if marker.name in clade_markers and marker.breadth >= breadth_thres})
        return markers


    def filter_markers_matrix(self, markers_matrix, messages=True):
        """Filters the primary samples, references and markers based on the user defined thresholds

        Args:
            markers_matrix (list[dict]): The list with the sample-to-markers information
            messages (bool): Whether to be verbose and halt execution when less than 4 samples are left.
             Defaults to True.
        """
        mm = pd.DataFrame.from_records(markers_matrix, index='sample')  # df with index samples and columns markers

        # Step I: samples with not enough markers are treated as secondary
        #   here the percentage is calculated from the number of markers in at least one sample, it can be less than
        #   the total number of available markers in the database for the clade
        n_samples_0, n_markers_0 = mm.shape
        min_markers = max(self.sample_with_n_markers, n_markers_0 * self.sample_with_n_markers_perc / 100)
        mm_primary = mm.loc[mm.sum(axis=1) >= min_markers]

        # Step II: filter markers in not enough primary samples
        n_samples_primary, _ = mm_primary.shape
        min_samples = n_samples_primary * self.marker_in_n_samples_perc / 100
        mm = mm.loc[:, mm_primary.sum(axis=0) >= min_samples]

        # Step III: filter samples with not enough markers
        #   here the percentage is calculated from the remaining markers
        _, n_markers_2 = mm.shape
        min_markers = max(self.sample_with_n_markers_after_filt,
                          n_markers_2 * self.sample_with_n_markers_after_filt_perc / 100)
        mm = mm.loc[mm.sum(axis=1) >= min_markers]
        n_samples_3 = len(mm)
        if n_samples_3 < 4:
            if messages:
                error(f"Phylogeny can not be inferred. Less than 4 samples remained after filtering.\n"
                      f"{n_samples_3} / {n_samples_0} samples ({n_samples_primary} primary) "
                      f"and {n_markers_2} / {n_markers_0 } markers remained.",
                      exit=True)
            return

        self.cleaned_markers_matrix = mm


    def copy_filtered_references(self, markers_tmp_dir, filtered_samples):
        """Copies the FASTA files of the filtered references to be processed by PhyloPhlAn

        Args:
            markers_tmp_dir (str): the temporary folder where to copy the reference genomes
            filtered_samples (set): the set of samples after filtering
        """
        for reference in self.references:
            reference_name = self.sample_path_to_name(reference)

            if reference_name in filtered_samples:
                reference_marker = os.path.join(self.tmp_dir, "reference_markers", f'{reference_name}.fna.bz2')
                copyfile(reference_marker, os.path.join(markers_tmp_dir, f"{reference_name}.fna.bz2"))


    def matrix_markers_to_fasta(self):
        """For each sample, writes the FASTA files with the sequences of the filtered markers

        Returns:
            str: the temporary folder where the FASTA files were written
        """
        markers_tmp_dir = os.path.join(self.tmp_dir, "{}.StrainPhlAn4".format(self.clade))
        create_folder(markers_tmp_dir)
        filtered_samples = set(self.cleaned_markers_matrix.index)
        filtered_markers = set(self.cleaned_markers_matrix.columns)
        self.copy_filtered_references(markers_tmp_dir, filtered_samples)
        execute_pool(((Strainphlan.sample_markers_to_fasta, sample, filtered_samples, filtered_markers,
                       self.trim_sequences, markers_tmp_dir) for sample in self.samples), self.nprocs)
        return markers_tmp_dir


    @classmethod
    def sample_markers_to_fasta(cls, sample_path, filtered_samples, filtered_markers, trim_sequences, markers_tmp_dir):
        """Writes a FASTA file with the filtered clade markers of a sample

        Args:
            sample_path (str): the path to the sample
            filtered_markers:
            filtered_samples:
            trim_sequences:
            markers_tmp_dir (str): the temporary folder were the FASTA file is written
        """
        if sample_path in filtered_samples:
            sample_name = cls.sample_path_to_name(sample_path)
            marker_output_file = os.path.join(markers_tmp_dir, f'{sample_name}.fna.bz2')
            sample = ConsensusMarkers.from_file(sample_path)
            sample.consensus_markers = [m for m in sample.consensus_markers if m.name in filtered_markers]
            sample.to_fasta(marker_output_file, trim_ends=trim_sequences)


    def get_markers_from_references(self):
        """Gets markers from reference files and returns the marker matrix with the reference markers

        Args:

        Returns:
            list: the list with the samples-to-markers information of the main samples and references
        """
        if not self.clade_markers_file:
            self.database_controller.extract_markers([self.clade], self.tmp_dir)
            self.clade_markers_file = os.path.join(self.tmp_dir, "{}.fna".format(self.clade))
        elif self.clade_markers_file.endswith(".bz2"):
            self.clade_markers_file = decompress_bz2(self.clade_markers_file, self.tmp_dir)

        return execute_pool(((Strainphlan.process_reference, reference, self.tmp_dir, self.clade_markers_file,
                              self.clade_markers_names, self.trim_sequences)
                             for reference in self.references), self.nprocs)


    @classmethod
    def process_reference(cls, reference_file, tmp_dir, clade_markers_file, clade_markers, trim_sequences):
        """Processes each reference file and get a markers dictionary to add to the markers matrix

        Args:
            reference_file (str): path to the reference file
            tmp_dir (str): the temporary folder where the BLASTn results where saved
            clade_markers_file (str):
            clade_markers (Iterable): the list with the clade markers names
            trim_sequences:

        Returns:
            dict: the dictionary with the reference-to-markers information
        """
        if reference_file.endswith(".bz2"):
            uncompressed_refernces_dir = os.path.join(tmp_dir, "uncompressed_references")
            os.makedirs(uncompressed_refernces_dir, exist_ok=True)
            reference_file = decompress_bz2(reference_file, uncompressed_refernces_dir)

        ext_markers = cls.extract_markers_from_genome(reference_file, clade_markers_file)

        reference_markers_dir = os.path.join(tmp_dir, "reference_markers")
        os.makedirs(reference_markers_dir, exist_ok=True)

        consensus_markers = ConsensusMarkers([ConsensusMarker(m, s) for m, s in ext_markers.items()])
        reference_name = cls.sample_path_to_name(reference_file)
        consensus_markers.to_fasta(os.path.join(reference_markers_dir, f'{reference_name}.fna.bz2'), trim_ends=trim_sequences)

        markers_matrix = {'sample': reference_file}
        markers_matrix.update({m: int(m in ext_markers) for m in clade_markers})

        return markers_matrix


    @staticmethod
    def extract_markers_from_genome(reference_file, clade_markers_file):
        """

        Args:
            reference_file (str):
            clade_markers_file (str):

        Returns:

        """
        # load the raw fasta data
        with util_fun.openrt(reference_file) as f:
            input_file_data = f.read()

        # parse the fasta
        input_seqs = {seq.id: seq for seq in SeqIO.parse(io.StringIO(input_file_data), 'fasta')}

        # we need the additional btop column
        columns = 'qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send evalue bitscore btop'
        cmd = f'blastn -query {clade_markers_file} -subject - -num_threads {1} -outfmt "6 {columns}"'

        # run blastn and pass the raw data to stdin
        r = run_command(cmd, input=input_file_data, text=True)

        # load the blast output
        df = pd.read_csv(io.StringIO(r.stdout), sep='\t', names=columns.split(' '))

        ext_markers = {}
        for _, row in df.iterrows():
            sseq = str(input_seqs[row['sseqid']].seq)
            btop = re.split(r'(\d+|\D{2})', row['btop'])[1::2]  # blast trace-back operations

            assert row['qend'] >= row['qstart']
            strand = 1 if row['send'] >= row['sstart'] else -1  # whether reverse-complemented

            qi = row['qstart'] - 1
            si = row['sstart'] - 1
            ext_s = '-' * qi
            for op in btop:
                if op.isnumeric():
                    op = int(op)
                    ext_s += sseq[si: si + strand * op: strand]
                    qi += op
                    si += strand * op
                else:
                    if strand == -1:
                        op = str(Seq.Seq(op).complement())

                    if op[0] == '-':  # query gap
                        si += strand
                    elif op[1] == '-':  # subject gap
                        ext_s += '-'
                        qi += 1
                    else:
                        qi += 1
                        si += strand
                        ext_s += op[1]

            ext_s += '-' * (row['qlen'] - row['qend'])

            assert qi == row['qend']
            assert len(ext_s) == row['qlen']

            if strand == -1:
                ext_s = str(Seq.Seq(ext_s).complement())

            ext_markers[row['qseqid']] = ext_s

        return ext_markers


    def calculate_polymorphic_rates(self):
        """Generates a file with the polymorphic rates of the species for each sample"""
        rows = []
        for sample_path in self.samples:
            sample = ConsensusMarkers.from_file(sample_path)
            p_stats, p_count, m_len = [], 0, 0
            for marker in sample.consensus_markers:
                if marker.name in self.clade_markers_names:
                    p_count += marker.get_polymorphisms()
                    m_len += marker.get_sequence_length()
                    p_stats.append(marker.get_polymorphism_perc())
            if m_len > 0:
                rows.append({
                    'sample':  self.sample_path_to_name(sample_path),
                    'percentage_of_polymorphic_sites': p_count * 100 / m_len,
                    'avg_by_marker': np.mean(p_stats),
                    'median_by_marker': np.median(p_stats),
                    'std_by_marker': np.std(p_stats),
                    'min_by_marker': np.min(p_stats),
                    'max_by_marker': np.max(p_stats),
                    'q25_by_marker': np.percentile(p_stats, 25),
                    'q75_by_marker': np.percentile(p_stats, 75),
                })

        df = pd.DataFrame(rows)
        df.to_csv(os.path.join(self.output_dir, f"{self.clade}.polymorphic"), sep='\t', index=False)


    def write_info(self):
        """Writes the information file for the execution"""
        filtered_names = [self.sample_path_to_name(sample) for sample in self.cleaned_markers_matrix.index]
        with open(os.path.join(self.output_dir, "{}.info".format(self.clade)), 'w') as info_file:
            info_file.write("Clade: {}\n".format(self.clade))
            info_file.write(
                "Number of samples: {}\n".format(len(self.samples)))
            info_file.write(
                "Number of references: {}\n".format(len(self.references)))
            info_file.write("Number of available markers for the clade: {}\n".format(len(self.clade_markers_names)))
            info_file.write("Filtering parameters:\n")
            info_file.write("\tNumber of bases to remove when trimming markers: {}\n".format(
                self.trim_sequences))
            info_file.write(f"\tMinimum number of markers to make a sample primary: {self.sample_with_n_markers}\n")
            info_file.write(f"\tMinimum percentage of markers to make a sample primary: {self.sample_with_n_markers_perc}\n")
            info_file.write(f"\tMinimum number of markers to keep a sample after filtering: {self.sample_with_n_markers_after_filt}\n")
            info_file.write(f"\tMinimum percentage of markers to keep a sample after filtering: {self.sample_with_n_markers_after_filt_perc}\n")
            info_file.write(f"\tMinimum percentage of samples to keep a marker: {self.marker_in_n_samples_perc}\n")
            info_file.write("Number of markers selected after filtering: {}\n".format(
                len(self.cleaned_markers_matrix.columns)))
            info_file.write("Number of samples after filtering: {}\n".format(len(
                [sample for sample in self.samples if sample in self.cleaned_markers_matrix.index])))
            info_file.write("Number of references after filtering: {}\n".format(len([reference for reference in self.references if self.sample_path_to_name(reference) in filtered_names])))
            info_file.write(
                "PhyloPhlan phylogenetic precision mode: {}\n".format(self.phylophlan_mode))
            info_file.write(
                "Number of processes used: {}\n".format(self.nprocs))


    def detect_clades(self):
        """Checks the clades that can be reconstructed from the pkl files

        Returns:
            dict: dictionary containing the number of samples a clade can be reconstructed from
        """
        markers2species = self.database_controller.get_markers2species()
        species2samples = {}
        species_to_check = set()
        info('Detecting clades...')
        for sample_path in self.samples:
            sample = ConsensusMarkers.from_file(sample_path)
            species_to_check.update((markers2species[marker.name] for marker in sample.consensus_markers if (
                marker.name in markers2species and marker.breadth >= self.breadth_thres)))
        info(f'  Will check {len(species_to_check)} species')
        for species in species_to_check:
            self.cleaned_markers_matrix = pd.DataFrame()
            self.clade = species
            self.clade_markers_names = self.database_controller.get_markers_for_clade(species)
            self.filter_markers_samples(print_clades=True)
            if len(self.cleaned_markers_matrix) >= 4:
                species2samples[species] = len(self.cleaned_markers_matrix)
        info('Done.')
        return species2samples


    def print_clades(self):
        """Prints the clades detected in the reconstructed markers"""
        species2samples = self.detect_clades()
        info('Detected clades: ')
        sorted_species2samples = collections.OrderedDict(sorted(species2samples.items(), key=lambda kv: kv[1],
                                                                reverse=True))
        with open(os.path.join(self.output_dir, 'print_clades_only.tsv'), 'w') as wf:
            wf.write('Clade\tNumber_of_samples\n')
            for species in sorted_species2samples:
                info('\t{}: in {} samples.'.format(species, sorted_species2samples[species]))
                wf.write('{}\t{}\n'.format(species, sorted_species2samples[species]))
        info('Done.')


    def interactive_clade_selection(self):
        """Allows the user to interactively select the SGB-level clade when specifying the clade at the species level"""
        if not self.non_interactive:
            info("The clade has been specified at the species level, starting interactive clade selection...")
        species2sgbs = self.database_controller.get_species2sgbs()
        if self.clade not in species2sgbs:
            error('The specified species "{}" is not present in the database. Exiting...'.format(
                self.clade), exit=True)
        sgbs_in_species = dict(sorted(
            species2sgbs[self.clade].items(), key=lambda item: item[1], reverse=True))
        if self.non_interactive or len(sgbs_in_species) == 1:
            selected_clade = next(iter(sgbs_in_species))
            if len(sgbs_in_species) == 1:
                info('Only one SGB ("{}") is available for the species "{}".'.format(
                    selected_clade, self.clade))
            info('Clade "{}" has been automatically selected.'.format(selected_clade))
            info("Done.")
            self.clade = selected_clade
            return
        info('Available SGBs for species "{}":'.format(self.clade))
        option_counter = 1
        for sgb in sgbs_in_species:
            info('[{}] {} ({} genomes)'.format(
                option_counter, sgb, sgbs_in_species[sgb]))
            option_counter += 1
        info('[{}] Exit'.format(option_counter))
        selected_option = input('Select option: ')
        tries = 1
        while not selected_option.isdigit() or int(selected_option) > option_counter or int(selected_option) == 0:
            info('"{}" is not a valid option'.format(selected_option))
            if tries == 3:
                error('Too many tries (3). Exiting...', exit=True)
            tries += 1
            selected_option = input('Select option: ')
        selected_option = int(selected_option)
        if selected_option == option_counter:
            error('Exiting...', exit=True)
        else:
            selected_clade = list(sgbs_in_species.keys())[selected_option-1]
            info('Clade "{}" has been selected.'.format(selected_clade))
            info("Done.")
            self.clade = selected_clade


    def filter_markers_samples(self, print_clades=False):
        """Retrieves the filtered markers matrix with the filtered samples and references

        Args:
            print_clades (bool, optional): Whether it was run in the print_clade_only mode. Defaults to False.
        """
        if not print_clades:
            info("Getting markers from samples...")
        markers_matrix = self.get_markers_matrix_from_samples()
        if not print_clades:
            info("Done.")
            if len(self.references) > 0:
                info("Getting markers from references...")
                markers_matrix += self.get_markers_from_references()
                info("Done.")
            info("Removing bad markers / samples...")
        self.filter_markers_matrix(markers_matrix, messages=not print_clades)
        if not print_clades:
            info("Done.")


    def run_strainphlan(self):
        """Runs the full StrainPhlAn pipeline"""
        if self.print_clades_only:
            self.print_clades()
            return

        if self.clade.startswith('s__') and 'CHOCOPhlAnSGB' in self.database_controller.get_database_name():
            self.interactive_clade_selection()
        self.clade_markers_names = self.database_controller.get_markers_for_clade(self.clade)
        info("Creating temporary directory...")
        self.tmp_dir = tempfile.mkdtemp(dir=self.tmp_dir)
        info("Done.")
        info("Filtering markers and samples...")
        self.filter_markers_samples()
        info("Done.")
        info("Writing samples as markers' FASTA files...")
        samples_as_markers_dir = self.matrix_markers_to_fasta()
        info("Done.")
        info("Calculating polymorphic rates...")
        self.calculate_polymorphic_rates()
        info("Done.")
        info("Computing phylogeny...")
        self.phylophlan_controller.compute_phylogeny(samples_as_markers_dir, len(self.cleaned_markers_matrix),
                                                     self.tmp_dir)
        info("Done.")
        info("Writing information file...")
        self.write_info()
        info("Done.")
        if not self.debug:
            info("Removing temporary files...")
            rmtree(self.tmp_dir, ignore_errors=False, onerror=None)
            info("Done.")


    def __init__(self, args):
        self.clade_markers_names = None
        self.cleaned_markers_matrix = None
        self.database_controller = MetaphlanDatabaseController(args.database)
        self.clade_markers_file = args.clade_markers
        self.samples = args.samples
        self.references = args.references if not args.print_clades_only else []
        self.clade = args.clade
        self.output_dir = args.output_dir
        self.trim_sequences = args.trim_sequences
        self.sample_with_n_markers = args.sample_with_n_markers
        self.sample_with_n_markers_perc = args.sample_with_n_markers_perc
        self.marker_in_n_samples_perc = args.marker_in_n_samples_perc
        self.sample_with_n_markers_after_filt = args.sample_with_n_markers_after_filt
        self.sample_with_n_markers_after_filt_perc = args.sample_with_n_markers_after_filt_perc
        self.breadth_thres = args.breadth_thres
        self.print_clades_only = args.print_clades_only
        self.non_interactive = args.non_interactive
        self.phylophlan_mode = args.phylophlan_mode
        self.phylophlan_params = args.phylophlan_params
        self.debug = args.debug
        self.nprocs = args.nprocs
        self.tmp_dir = args.output_dir if args.tmp is None else args.tmp
        self.phylophlan_controller = Phylophlan3Controller(args)


def read_params():
    """
    Reads and parses the command line arguments of the script
    
    Returns:
        namespace: The populated namespace with the command line arguments
    """
    p = ap.ArgumentParser(description="", formatter_class=ap.ArgumentDefaultsHelpFormatter)
    p.add_argument('-d', '--database', type=str, default='latest',
                   help="The input MetaPhlAn {} database".format(__version__))
    p.add_argument('-m', '--clade_markers', type=str, default=None,
                   help="The clade markers as FASTA file")
    p.add_argument('-s', '--samples', type=str, nargs='+',
                   help="The reconstructed markers for each sample")
    p.add_argument('-r', '--references', type=str, nargs='+', default=[],
                   help="The reference genomes")
    p.add_argument('-c', '--clade', type=str, default=None,
                   help="The clade to investigate")
    p.add_argument('-o', '--output_dir', type=str, required=True,
                   help="The output directory")
    p.add_argument('-n', '--nprocs', type=int, default=1,
                   help="The number of threads to use")
    p.add_argument('--trim_sequences', type=int, default=50,
                   help="The number of bases to remove from both ends when trimming markers")
    p.add_argument('--sample_with_n_markers', type=int, default=20,
                   help="Threshold defining the minimum absolute number of markers for a sample to be primary. "
                        "This rule is combined with AND logic with --sample_with_n_markers_perc")
    p.add_argument('--sample_with_n_markers_perc', type=float, default=25,
                   help="Threshold defining the minimum percentage of markers to for a sample to be primary. "
                        "This rule is combined with AND logic with --sample_with_n_markers")
    p.add_argument('--marker_in_n_samples_perc', type=float, default=50,
                   help="Threshold defining the minimum percentage of primary samples to keep a marker")
    p.add_argument('--sample_with_n_markers_after_filt', type=int, default=20,
                   help="Threshold defining the minimum absolute number of markers after filtering to keep a sample. "
                        "This rule is combined with AND logic with --sample_with_n_markers_after_filt_perc")
    p.add_argument('--sample_with_n_markers_after_filt_perc', type=float, default=25,
                   help="Threshold defining the minimum percentage of markers kept after filtering to keep a sample. "
                        "This rule is combined with AND logic with --sample_with_n_markers_after_filt")
    p.add_argument('--breadth_thres', type=int, default=80,
                   help="Threshold defining the minimum breadth of coverage for the markers")
    p.add_argument('--phylophlan_mode', choices=['accurate', 'fast'], default='fast',
                   help="The presets for fast or accurate phylogenetic analysis")
    p.add_argument('--phylophlan_configuration', type=str, default=None,
                   help="The PhyloPhlAn configuration file")
    p.add_argument('--tmp', type=str, default=None,
                   help="If specified, the directory where to store the temporary files.")
    p.add_argument('--mutation_rates', action='store_true', default=False,
                   help="If specified, StrainPhlAn will produce a mutation rates table for each of the aligned markers"
                        " and a summary table for the concatenated MSA. This operation can take long time to finish")
    p.add_argument('--phylophlan_params', type=str, default=None, help="Additional phylophlan parameters")
    p.add_argument('--print_clades_only', action='store_true', default=False,
                   help="If specified, StrainPhlAn will only print the potential clades and stop the execution")
    p.add_argument('--non_interactive', action='store_true', default=False,
                   help="If specified, StrainPhlAn will select the first SGB available when the clade is specified at"
                        " the species level")
    p.add_argument('--treeshrink', action='store_true', default=False,
                   help="If specified, StrainPhlAn will execute TreeShrink after building the tree")
    p.add_argument('--debug', action='store_true', default=False,
                   help="If specified, StrainPhlAn will not remove the temporary folders")
    p.add_argument('-v', '--version', action='store_true',
                   help="Shows this help message and exit")

    return p.parse_args()


def check_params(args):
    """Checks the mandatory command line arguments of the script

    Args:
        args (namespace): the arguments to check
    """
    if args.print_clades_only and args.clade_markers:
        error('-m (or --clade_markers) cannot be specified together with --print_clades_only', exit=True)
    elif not args.print_clades_only and not args.clade:
        error('-c (or --clade) must be specified', exit=True)
    elif not os.path.exists(args.output_dir):
        error('The directory {} does not exist'.format(args.output_dir), exit=True)
    elif args.tmp is not None and not os.path.exists(args.tmp):
        error('The directory {} does not exist'.format(args.tmp), exit=True)
    elif args.database != 'latest' and not os.path.exists(args.database):
        error('The database does not exist', exit=True)
    elif args.clade_markers and not os.path.exists(args.clade_markers):
        error('The clade markers file does not exist', exit=True)
    elif args.phylophlan_configuration and not os.path.exists(args.phylophlan_configuration):
        error('The phylophlan configuration file does not exist', exit=True)
    elif len(args.samples)+len(args.references) < 4:
        error('The main inputs samples + references are less than 4', exit=True)
    for s in args.samples:
        if not os.path.exists(s):
            error('The input sample file \"{}\" does not exist'.format(s), exit=True)
    for r in args.references:
        if not os.path.exists(r):
            error('The reference file \"{}\" does not exist'.format(r), exit=True)


def main():
    t0 = time.time()
    args = read_params()
    if args.version:
        info('StrainPhlAn version {} ({})'.format(__version__, __date__))
        exit(0)
    check_params(args)
    info("Start StrainPhlAn {} execution".format(__version__))
    strainphlan_runner = Strainphlan(args)
    strainphlan_runner.run_strainphlan()
    exec_time = time.time() - t0
    info("Finish StrainPhlAn {} execution ({} seconds): Results are stored at \"{}\"".format(
        __version__, round(exec_time, 2), args.output_dir))


if __name__ == '__main__':
    main()
    