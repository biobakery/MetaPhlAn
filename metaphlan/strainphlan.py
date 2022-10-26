#!/usr/bin/env python
__author__ = ('Aitor Blanco Miguez (aitor.blancomiguez@unitn.it), '
              'Duy Tin Truong (duytin.truong@unitn.it), '
              'Francesco Asnicar (f.asnicar@unitn.it), '
              'Moreno Zolfo (moreno.zolfo@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)')
__version__ = '4.0.3'
__date__ = '24 Oct 2022'


import argparse as ap
import collections
import os
import tempfile
import time
from shutil import copyfile, rmtree
import numpy
import pandas as pd
from Bio import SeqIO
try:
    from .utils import *
except ImportError:
    from utils import *


class Strainphlan:
    """StrainPhlAn class"""

    def add_secondary_samples_and_references(self, messages=True):
        """Adds secondary samples and references to the marker matrix

        Args:
            messages (bool, optional): Whether the function is verbose. Defaults to True.
        """
        if len(self.secondary_samples) > 0:
            if messages:
                info("Getting markers from secondary sample files...")
            self.add_secondary_samples()
            if messages:
                info("Done.")
        if len(self.secondary_references) > 0:
            info("Getting markers from secondary reference files...")
            self.add_secondary_references()
            info("Done.")

    def add_secondary_samples(self):
        """Adds secondary samples to the marker matrix"""
        filtered_clade_markers = self.cleaned_markers_matrix.columns.tolist()
        markers_matrix = execute_pool(((self.get_matrix_for_sample, sample, filtered_clade_markers)
                                      for sample in self.secondary_samples), self.nprocs)
        self.include_secondary_matrix(markers_matrix)

    def add_secondary_references(self):
        """Adds secondary references to the marker matrix"""
        filtered_clade_markers = self.cleaned_markers_matrix.columns.tolist()
        markers_matrix = execute_pool(((self.process_reference, reference, os.path.join(
            self.tmp_dir, "blastn"), filtered_clade_markers) for reference in self.secondary_references), self.nprocs)
        self.include_secondary_matrix(markers_matrix)

    def include_secondary_matrix(self, markers_matrix):
        """Returns the cleaned matrix adding secondary samples/references

        Args:
            markers_matrix (list): List containing the samples-to-markers information for the secondary entries
        """
        mm = pd.DataFrame.from_records(markers_matrix, index='sample')
        mm = mm.loc[mm.sum(axis=1) >= (self.secondary_sample_with_n_markers if self.abs_n_markers_thres else len(
            self.db_clade_markers) * self.secondary_sample_with_n_markers // 100)]
        self.cleaned_markers_matrix = pd.concat(
            [self.cleaned_markers_matrix, mm])

    def get_markers_matrix(self, print_clades=False):
        """Gets a binary matrix representing the presence/ausence of the clade markers in the uploaded samples

        Args:
            print_clades (bool, optional): Whether the function has been called from the print_clades_only function. Defaults to False.

        Returns:
            list: the list containing the samples-to-markers information
        """
        if not print_clades:
            if not self.clade_markers_file:
                self.database_controller.extract_markers([self.clade], self.tmp_dir)
                self.clade_markers_file = os.path.join(self.tmp_dir, "{}.fna".format(self.clade))
            else:
                base_name = os.path.basename(self.clade_markers_file)
                name, ext = os.path.splitext(base_name)
                if ext == '.bz2':
                    decompress_bz2(self.clade_markers_file, self.tmp_dir)
                    self.clade_markers_file = os.path.join(self.tmp_dir, name)
                else:
                    copyfile(self.clade_markers_file,
                             os.path.join(self.tmp_dir, base_name))
                    self.clade_markers_file = os.path.join(
                        self.tmp_dir, base_name)
            self.db_clade_markers = {rec.id: rec.seq for rec in SeqIO.parse(
                open(self.clade_markers_file, 'r'), 'fasta')}
        markers_matrix = execute_pool(((self.get_matrix_for_sample, sample, list(
            self.db_clade_markers.keys())) for sample in self.samples), self.nprocs)
        return markers_matrix

    def get_matrix_for_sample(self, sample_path, clade_markers):
        """Returns the matrix with the presence / absence of the clade markers in a samples

        Args:
            sample_path (str): the path to the sample
            clade_markers (list): the list with the clade markers names

        Returns:
            dict: dictionary containing the sample-to-markers information as a binary matrix
        """
        sample = ConsensusMarkers(pkl_file=sample_path)
        markers = {"sample": sample_path}
        markers.update({m: 0 for m in clade_markers})
        markers.update(
            {marker.name: 1 for marker in sample.consensus_markers if marker.name in clade_markers and marker.breadth >= self.breadth_thres})
        return markers

    def clean_markers_matrix(self, markers_matrix, messages=True):
        """Filters the primary samples, references and markers based on the user defined thresholds

        Args:
            markers_matrix (list): The list with the sample-to-markers information
            messages (bool, optional): Whether to be verbose. Defaults to True.
        """
        mm = pd.DataFrame.from_records(markers_matrix, index='sample')
        self.min_markers = self.sample_with_n_markers_after_filt if self.abs_n_markers_thres else len(
            self.db_clade_markers) * self.sample_with_n_markers_after_filt // 100
        # Checks if the percentage of markers of a sample sample reaches sample_with_n_markers, if not, removes the sample
        mm = mm.loc[mm.sum(axis=1) >= (self.sample_with_n_markers if self.abs_n_markers_thres else len(
            self.db_clade_markers) * self.sample_with_n_markers // 100)]
        if not self.check_matrix_length(mm.index, 4, "Phylogeny can not be inferred. Too many samples were discarded.", messages):
            return
        # Checks if the percentage of samples that contain a marker reaches marker_in_n_samples, if not, removes the marker
        mm = mm.loc[:, mm.sum(axis=0) >= (
            self.marker_in_n_samples if self.abs_n_samples_thres else self.marker_in_n_samples * len(mm) // 100)]
        if not self.check_matrix_length(mm.columns, self.min_markers, "Phylogeny can not be inferred. Too many markers were discarded.", messages):
            return
        # Checks again if the percentage of markers of a sample reaches sample_with_n_markers_after_filt, if not, removes the sample
        mm = mm.loc[mm.sum(axis=1) >= self.min_markers]
        if not self.check_matrix_length(mm.index, 4, "Phylogeny can not be inferred. No enough markers were kept for the samples.", messages):
            return
        self.cleaned_markers_matrix = mm

    def check_matrix_length(self, matrix, threshold, message, verbose):
        """Checks the length of the markers' matrix

        Args:
            matrix (list): the matrix to check the length
            threshold (int): size threshold
            message (str): message to print if verbose
            verbose (bool): whether to be verbose

        Returns:
            bool: whether the matrix size reaches the length threshold
        """
        if len(matrix) < threshold:
            if verbose:
                error(message, exit=True)
            return False
        return True

    def copy_filtered_references(self, markers_tmp_dir):
        """Copies the FASTA files of the filtered references to be processed by PhyloPhlAn

        Args:
            markers_tmp_dir (str): the temporal folder where to copy the reference genomes
        """
        for reference in self.references+self.secondary_references:
            sample_name = os.path.splitext(os.path.basename(reference[:-4]))[0] if reference.endswith(
                '.bz2') else os.path.splitext(os.path.basename(reference))[0]
            if sample_name in [os.path.splitext(os.path.basename(sample))[0] for sample in self.cleaned_markers_matrix.index.tolist()]:
                if reference.endswith('.bz2'):
                    decompress_bz2(reference, markers_tmp_dir)
                else:
                    copyfile(reference, os.path.join(
                        markers_tmp_dir, "{}.fna".format(sample_name)))

    def matrix_markers_to_fasta(self):
        """For each sample, writes the FASTA files with the sequences of the filtered markers

        Returns:
            str: the temporal folder where the FASTA files were written
        """
        markers_tmp_dir = os.path.join(
            self.tmp_dir, "{}.StrainPhlAn4".format(self.clade))
        create_folder(markers_tmp_dir)
        self.copy_filtered_references(markers_tmp_dir)
        execute_pool(((self.sample_markers_to_fasta, sample, markers_tmp_dir)
                     for sample in self.samples+self.secondary_samples), self.nprocs)
        return markers_tmp_dir

    def sample_markers_to_fasta(self, sample_path, markers_tmp_dir):
        """Writes a FASTA file with the filtered clade markers of a sample

        Args:
            sample_path (str): the path to the sample
            markers_tmp_dir (str): the temporal folder were the FASTA file is written
        """
        if sample_path in self.cleaned_markers_matrix.index.tolist():
            sample_name = os.path.splitext(os.path.basename(sample_path))[
                0].replace(".pkl", "")
            with open(os.path.join(markers_tmp_dir, '{}.fna'.format(sample_name)), 'w') as marker_fna:
                sample = ConsensusMarkers(pkl_file=sample_path)
                for marker in sample.consensus_markers:
                    if marker.name in self.cleaned_markers_matrix.columns.tolist():
                        SeqIO.write(marker.get_sequence(
                            self.trim_sequences), marker_fna, 'fasta')

    def cleaned_clade_markers_to_fasta(self):
        """Writes a FASTA file with the sequences of the filtered clade markers"""
        clade_tmp_dir = os.path.join(self.tmp_dir, self.clade[:30])
        create_folder(clade_tmp_dir)
        for m in self.cleaned_markers_matrix.columns.tolist():
            marker = ConsensusMarker(m, self.db_clade_markers[m])
            with open(os.path.join(clade_tmp_dir, '{}.fna'.format(marker.parse_marker_name())), 'w') as marker_fna:
                SeqIO.write(marker.get_sequence(
                    self.trim_sequences), marker_fna, 'fasta')

    def get_markers_from_references(self, markers_matrix):
        """Gets markers from reference files and returns the marker matrix with the reference markers

        Args:
            markers_matrix (list): the list with the samples-to-markers information of the main samples

        Returns:
            list: the list with the samples-to-markers information of the main samples and references
        """
        create_folder(os.path.join(self.tmp_dir, "blastn"))
        results = execute_pool(((self.process_reference, reference, os.path.join(self.tmp_dir, "blastn"), list(
            self.db_clade_markers.keys())) for reference in self.references), self.nprocs)
        return markers_matrix + results


    def process_reference(self, reference, blastn_dir, clade_markers):
        """Processes each reference file and get a markers dictionary to add to the markers matrix

        Args:
            reference (str): path to the reference file
            blastn_dir (str): the temporal folder where the BLASTn results where saved
            clade_markers (list): the list with the clade markers names

        Returns:
            dict: the dictionary with the reference-to-markers information
        """
        name = os.path.splitext(os.path.basename(reference))[0]
        if reference.endswith(".bz2"):
            name = os.path.splitext(name)[0]
            reference = decompress_bz2(reference, blastn_dir)
        blastn_db = create_blastn_db(blastn_dir, reference)
        blastn_file = execute_blastn(
            blastn_dir, self.clade_markers_file, blastn_db)
        return self.parse_blastn_results(os.path.join(blastn_dir, '{}.pkl'.format(name)), blastn_file, clade_markers)


    def parse_blastn_results(self, sample, blastn_file, clade_markers):
        """Parses BLASTn results and gets the presence of the clade markers in the reference file

        Args:
            sample (str): path to the reference file
            blastn_file (str): path to the BLASTn result
            clade_markers (list): the list with the clade markers names

        Returns:
            dict: the dictionary with the reference-to-markers information
        """
        markers = {"sample": sample}
        markers.update({m: 0 for m in clade_markers})
        with open(blastn_file, "r") as rf:
            for line in rf:
                if line == '':
                    break
                markers[line.split("\t")[0]] = 1
        return markers

    def calculate_polymorphic_rates(self):
        """Generates a file with the polymorphic rates of the species for each sample"""
        with open(os.path.join(self.output_dir, "{}.polymorphic".format(self.clade)), 'w') as polymorphic_file:
            polymorphic_file.write("sample\tpercentage_of_polymorphic_sites\tavg_by_marker\tmedian_by_marker" +
                                   "\tstd_by_marker\tmin_by_marker\tmax_by_marker\tq25_by_marker\tq75_by_marker")
            for sample_path in self.samples+self.secondary_samples:
                sample = ConsensusMarkers(pkl_file=sample_path)
                p_stats, p_count, m_len = [], 0, 0
                for marker in sample.consensus_markers:
                    if marker.name in self.db_clade_markers:
                        p_count += marker.get_polymorphisms()
                        m_len += marker.get_sequence_length()
                        p_stats.append(marker.get_polymorphism_perc())
                if m_len > 0:
                    polymorphic_file.write("\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(os.path.splitext(os.path.basename(sample_path))[0],
                                                                                         p_count * 100 / m_len, numpy.average(p_stats), numpy.percentile(
                                                                                             p_stats, 50), numpy.std(p_stats), numpy.min(p_stats),
                                                                                         numpy.max(p_stats), numpy.percentile(p_stats, 25), numpy.percentile(p_stats, 75)))

    def write_info(self):
        """Writes the information file for the execution"""
        filtered_names = [os.path.splitext(os.path.basename(sample))[
            0] for sample in self.cleaned_markers_matrix.index.tolist()]
        with open(os.path.join(self.output_dir, "{}.info".format(self.clade)), 'w') as info_file:
            info_file.write("Clade: {}\n".format(self.clade))
            info_file.write(
                "Number of main samples: {}\n".format(len(self.samples)))
            info_file.write("Number of secondary samples: {}\n".format(
                len(self.secondary_samples)))
            info_file.write(
                "Number of main references: {}\n".format(len(self.references)))
            info_file.write("Number of secondary references: {}\n".format(
                len(self.secondary_references)))
            info_file.write("Number of available markers for the clade: {}\n".format(
                len(self.db_clade_markers)))
            info_file.write("Filtering parameters:\n")
            info_file.write("\tNumber of bases to remove when trimming markers: {}\n".format(
                self.trim_sequences))
            info_file.write("\tMinimum {} of markers to keep a main sample: {}\n".format(
                'number' if self.abs_n_markers_thres else 'percentage', self.sample_with_n_markers))
            info_file.write("\tMinimum {} of markers to keep a main sample after first filtering: {}\n".format(
                'number' if self.abs_n_markers_thres else 'percentage', self.sample_with_n_markers_after_filt))
            info_file.write("\tMinimum {} of markers to keep a secondary sample: {}\n".format(
                'number' if self.abs_n_markers_thres else 'percentage', self.secondary_sample_with_n_markers))
            info_file.write("\tMinimum {} of samples to keep a marker: {}\n".format(
                'number' if self.abs_n_samples_thres else 'percentage', self.marker_in_n_samples))
            info_file.write("Number of markers selected after filtering: {}\n".format(
                len(self.cleaned_markers_matrix.columns)))
            info_file.write("Number of main samples after filtering: {}\n".format(len(
                [sample for sample in self.samples if sample in self.cleaned_markers_matrix.index.tolist()])))
            info_file.write("Number of secondary samples after filtering: {}\n".format(len(
                [sample for sample in self.secondary_samples if sample in self.cleaned_markers_matrix.index.tolist()])))
            info_file.write("Number of main references after filtering: {}\n".format(len([reference for reference in self.references if (os.path.splitext(
                os.path.basename(reference[:-4]))[0] if reference.endswith('.bz2') else os.path.splitext(os.path.basename(reference))[0]) in filtered_names])))
            info_file.write("Number of secondary references after filtering: {}\n".format(len([reference for reference in self.secondary_references if (os.path.splitext(
                os.path.basename(reference[:-4]))[0] if reference.endswith('.bz2') else os.path.splitext(os.path.basename(reference))[0]) in filtered_names])))
            info_file.write(
                "PhyloPhlan phylogenetic precision mode: {}\n".format(self.phylophlan_mode))
            info_file.write(
                "Number of processes used: {}\n".format(self.nprocs))

    def detect_clades(self, markers2species):
        """Checks the clades that can be reconstructed from the pkl files

        Args:
            markers2species (dict): dictionary containing the clade each marker belong to

        Returns:
            dict: dictionary containing the number of samples a clade can be reconstructed from
        """
        species2samples = {}
        species_to_check = set()
        info('Detecting clades...')
        for sample_path in self.samples:
            sample = ConsensusMarkers(pkl_file=sample_path)
            species_to_check.update({markers2species[marker.name] for marker in sample.consensus_markers if (
                marker.name in markers2species and marker.breadth >= self.breadth_thres)})
        for species in species_to_check:
            self.cleaned_markers_matrix = pd.DataFrame()
            self.clade = species
            self.db_clade_markers = {
                marker: marker for marker in markers2species if markers2species[marker] == self.clade}
            self.filter_markers_samples(print_clades=True)
            if len(self.cleaned_markers_matrix) >= 4:
                species2samples[species] = len(self.cleaned_markers_matrix)
        info('Done.')
        return species2samples

    def print_clades(self):
        """Prints the clades detected in the reconstructed markers"""
        markers2species = self.database_controller.get_markers2species()
        species2samples = self.detect_clades(markers2species)
        info('Detected clades: ')
        sorted_species2samples = collections.OrderedDict(
            sorted(species2samples.items(), key=lambda kv: kv[1], reverse=True))
        with open(os.path.join(self.output_dir, 'print_clades_only.tsv'), 'w') as wf:
            wf.write('Clade\tNumber_of_samples\n')
            for species in sorted_species2samples:
                info('\t{}: in {} samples.'.format(
                    species, sorted_species2samples[species]))
                wf.write('{}\t{}\n'.format(
                    species, sorted_species2samples[species]))
        info('Done.')

    def interactive_clade_selection(self):
        """Allows the user to interactively select the SGB-level clade when specifing the clade at the species level"""
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
            print_clades (bool, optional): Whether it was runned in the print_clade_only mode. Defaults to False.
        """
        if not print_clades:
            info("Getting markers from main samples...")
        markers_matrix = self.get_markers_matrix(print_clades)
        if not print_clades:
            info("Done.")
            info("Getting markers from main references...")
            markers_matrix = self.get_markers_from_references(markers_matrix)
            info("Done.")
            info("Removing bad markers / samples...")
        self.clean_markers_matrix(markers_matrix, messages=not print_clades)
        if not print_clades:
            info("Done.")
            info("Getting markers from secondary samples and references...")
        self.add_secondary_samples_and_references(messages=not print_clades)
        if not print_clades:
            info("Done.")

    def run_strainphlan(self):
        """Runs the full StrainPhlAn pipeline"""
        if self.print_clades_only:
            self.print_clades()
        else:
            if self.clade.startswith('s__') and 'CHOCOPhlAnSGB' in self.database_controller.get_database_name():
                self.interactive_clade_selection()
            info("Creating temporary directory...")
            self.tmp_dir = tempfile.mkdtemp(dir=self.tmp_dir)
            info("Done.")
            info("Filtering markers and samples...")
            self.filter_markers_samples()
            info("Done.")
            info("Writing samples as markers' FASTA files...")
            samples_as_markers_dir = self.matrix_markers_to_fasta()
            info("Done.")
            info("Writing filtered clade markers as FASTA file...")
            self.cleaned_clade_markers_to_fasta()
            info("Done.")
            info("Calculating polymorphic rates...")
            self.calculate_polymorphic_rates()
            info("Done.")
            info("Executing PhyloPhlAn...")
            self.phylophlan_controller.compute_phylogeny(samples_as_markers_dir, len(
                self.cleaned_markers_matrix), self.min_markers, self.tmp_dir)
            info("Done.")
            info("Writing information file...")
            self.write_info()
            info("Done.")
            if not self.debug:
                info("Removing temporary files...")
                rmtree(self.tmp_dir, ignore_errors=False, onerror=None)
                info("Done.")
                
    def __init__(self, args):
        self.database_controller = MetaphlanDatabaseController(args.database)
        self.clade_markers_file = args.clade_markers
        self.samples = args.samples
        self.references = args.references if not args.print_clades_only else []
        self.secondary_samples = args.secondary_samples
        self.secondary_references = args.secondary_references if not args.print_clades_only else []
        self.clade = args.clade
        self.output_dir = args.output_dir
        self.trim_sequences = args.trim_sequences
        self.sample_with_n_markers = args.sample_with_n_markers
        self.marker_in_n_samples = args.marker_in_n_samples
        self.secondary_sample_with_n_markers = args.secondary_sample_with_n_markers
        self.sample_with_n_markers_after_filt = args.sample_with_n_markers * \
            2 // 3 if (args.sample_with_n_markers_after_filt is None or args.sample_with_n_markers <
                       args.sample_with_n_markers_after_filt) else args.sample_with_n_markers_after_filt
        self.abs_n_markers_thres = args.abs_n_markers_thres
        self.abs_n_samples_thres = args.abs_n_samples_thres
        self.breadth_thres = args.breadth_thres
        self.print_clades_only = args.print_clades_only
        self.non_interactive = args.non_interactive
        self.phylophlan_mode = args.phylophlan_mode
        self.debug = args.debug
        self.nprocs = args.nprocs
        self.tmp_dir = args.output_dir if args.tmp is None else args.tmp
        self.phylophlan_controller = Phylophlan3Controller(args)
        self.cleaned_markers_matrix = pd.DataFrame()

def read_params():
    """ Reads and parses the command line arguments of the script
    
    Returns:
        namespace: The populated namespace with the command line arguments
    """
    p = ap.ArgumentParser(
        description="", formatter_class=ap.ArgumentDefaultsHelpFormatter)
    p.add_argument('-d', '--database', type=str, default='latest',
                   help="The input MetaPhlAn {} database".format(__version__))
    p.add_argument('-m', '--clade_markers', type=str, default=None,
                   help="The clade markers as FASTA file")
    p.add_argument('-s', '--samples', type=str, nargs='+', default=[],
                   help="The reconstructed markers for each sample")
    p.add_argument('-r', '--references', type=str, nargs='+', default=[],
                   help="The reference genomes")
    p.add_argument('-c', '--clade', type=str, default=None,
                   help="The clade to investigate")
    p.add_argument('-o', '--output_dir', type=str, default=None,
                   help="The output directory")
    p.add_argument('-n', '--nprocs', type=int, default=1,
                   help="The number of threads to use")
    p.add_argument('--secondary_samples', type=str, nargs='+', default=[],
                   help="The reconstructed markers for each secondary sample")
    p.add_argument('--secondary_references', type=str, nargs='+', default=[],
                   help="The secondary reference genomes")
    p.add_argument('--trim_sequences', type=int, default=50,
                   help="The number of bases to remove from both ends when trimming markers")
    p.add_argument('--marker_in_n_samples', type=int, default=80,
                   help="Theshold defining the minimum percentage of samples to keep a marker")
    p.add_argument('--sample_with_n_markers', type=int, default=80,
                   help="Threshold defining the minimun percentage of markers to keep a sample")
    p.add_argument('--secondary_sample_with_n_markers', type=int, default=80,
                   help="Threshold defining the minimun percentage of markers to keep a secondary sample")
    p.add_argument('--sample_with_n_markers_after_filt', type=int, default=None,
                   help="Threshold defining the minimun percentage of markers to keep a sample after filtering the markers [only for dev]")
    p.add_argument('--abs_n_markers_thres', action='store_true', default=False,
                   help="If specified, the *sample_with_n_markers* thresholds will be specified as absolute numbers")
    p.add_argument('--abs_n_samples_thres', action='store_true', default=False,
                   help="If specified, the marker_in_n_samples threshold will be specified as absolute numbers")
    p.add_argument('--breadth_thres', type=int, default=80,
                   help="Threshold defining the minimum breadth of coverage for the markers")
    p.add_argument('--phylophlan_mode', choices=['accurate', 'fast'], default='fast',
                   help="The presets for fast or accurate phylogenetic analysis")
    p.add_argument('--phylophlan_configuration', type=str, default=None,
                   help="The PhyloPhlAn configuration file")
    p.add_argument('--tmp', type=str, default=None,
                   help="If specified, the directory where to store the temporal files.")
    p.add_argument('--mutation_rates', action='store_true', default=False,
                   help="If specified, StrainPhlAn will produce a mutation rates table for each of the aligned markers and a summary table "
                   "for the concatenated MSA. This operation can take long time to finish")
    p.add_argument('--print_clades_only', action='store_true', default=False,
                   help="If specified, StrainPhlAn will only print the potential clades and stop the execution")
    p.add_argument('--non_interactive', action='store_true', default=False,
                   help="If specified, StrainPhlAn will select the first SGB available when the clade is specified at the species level")
    p.add_argument('--treeshrink', action='store_true', default=False,
                   help="If specified, StrainPhlAn will execute TreeShrink after building the tree")
    p.add_argument('--debug', action='store_true', default=False,
                   help="If specified, StrainPhlAn will not remove the temporal folders")
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
    elif not args.samples:
        error('-s (or --samples) must be specified', exit=True)
    elif not args.print_clades_only and not args.clade:
        error('-c (or --clade) must be specified', exit=True)
    elif not args.output_dir:
        error('-o (or --output_dir) must be specified', exit=True)
    elif not os.path.exists(args.output_dir):
        error('The directory {} does not exist'.format(args.output_dir), exit=True)
    elif not (args.tmp is None) and not os.path.exists(args.tmp):
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
    for s in args.secondary_samples:
        if not os.path.exists(s):
            error('The secondary input sample file \"{}\" does not exist'.format(
                s), exit=True)
    for r in args.secondary_references:
        if not os.path.exists(r):
            error('The secondary reference file \"{}\" does not exist'.format(
                r), exit=True)

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
    