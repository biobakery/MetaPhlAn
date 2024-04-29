#!/usr/bin/env python
__author__ = ('Aitor Blanco Miguez (aitor.blancomiguez@unitn.it), '
              'Duy Tin Truong (duytin.truong@unitn.it), '
              'Francesco Asnicar (f.asnicar@unitn.it), '
              'Moreno Zolfo (moreno.zolfo@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it),'
              'Michal Puncochar (michal.puncochar@unitn.it')
__version__ = '4.1.1'
__date__ = '11 Mar 2024'


import argparse as ap
import io
import os
import re
import tempfile
import time
from collections import OrderedDict
from shutil import copyfile, rmtree
from typing import Iterable

import numpy as np
import pandas as pd
from Bio import SeqIO, Seq

try:
    from .utils import *
except ImportError:
    from utils import *


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

        return execute_pool(((Strainphlan.get_matrix_for_sample, sample_path, self.clade_markers_names,
                              self.breadth_thres) for sample_path in self.samples), self.nprocs)


    @classmethod
    def get_matrix_for_sample(cls, sample_path, clade_markers, breadth_thres):
        """Returns the matrix with the presence / absence of the clade markers in a samples

        Args:
            sample_path (str): the path to the sample
            clade_markers (Iterable): the list with the clade markers names
            breadth_thres:

        Returns:
            dict: dictionary containing the sample-to-markers information as a binary matrix
        """
        sample = ConsensusMarkers.from_file(sample_path)
        sample_name = cls.sample_path_to_name(sample_path)

        markers = {"sample_name": sample_name}
        markers.update({m: 0 for m in clade_markers})
        markers.update({marker.name: 1 for marker in sample.consensus_markers
                        if marker.name in clade_markers and marker.breadth >= breadth_thres})
        return markers


    def filter_markers_matrix(self, markers_matrix, messages=True):
        """Filters the primary samples, references and markers based on the user defined thresholds

        Args:
            markers_matrix (pd.DataFrame): The presence-absence data frame with index as samples and columns as markers
            messages (bool): Whether to be verbose and halt execution when less than 4 samples are left.
             Defaults to True.
        """

        # Step I: samples with not enough markers are treated as secondary
        #   here the percentage is calculated from the number of markers in at least one sample, it can be less than
        #   the total number of available markers in the database for the clade
        n_samples_0, n_markers_0 = markers_matrix.shape
        min_markers = max(self.sample_with_n_markers, n_markers_0 * self.sample_with_n_markers_perc / 100)
        mm_primary = markers_matrix.loc[markers_matrix.sum(axis=1) >= min_markers]

        # Step II: filter markers in not enough primary samples
        n_samples_primary, _ = mm_primary.shape
        min_samples = n_samples_primary * self.marker_in_n_samples_perc / 100
        markers_matrix = markers_matrix.loc[:, mm_primary.sum(axis=0) >= min_samples]

        # Step III: filter samples with not enough markers
        #   here the percentage is calculated from the remaining markers
        _, n_markers_2 = markers_matrix.shape
        min_markers = max(self.sample_with_n_markers_after_filt,
                          n_markers_2 * self.sample_with_n_markers_after_filt_perc / 100)
        markers_matrix = markers_matrix.loc[markers_matrix.sum(axis=1) >= min_markers]

        n_samples_3 = len(markers_matrix)
        if n_samples_3 < 4:
            if messages:
                error(f"Phylogeny can not be inferred. Less than 4 samples remained after filtering.\n"
                      f"{n_samples_3} / {n_samples_0} samples ({n_samples_primary} primary) "
                      f"and {n_markers_2} / {n_markers_0 } markers remained.",
                      exit=True)
            return markers_matrix

        return markers_matrix


    def copy_filtered_references(self, markers_tmp_dir, filtered_sample_names):
        """Copies the FASTA files of the filtered references to be processed by PhyloPhlAn

        Args:
            markers_tmp_dir (str): the temporary folder where to copy the reference genomes
            filtered_sample_names (set): the set of samples after filtering
        """
        for reference in self.references:
            reference_name = self.sample_path_to_name(reference)

            if reference_name in filtered_sample_names:
                reference_marker = os.path.join(self.tmp_dir, "reference_markers", f'{reference_name}.fna.bz2')
                copyfile(reference_marker, os.path.join(markers_tmp_dir, f"{reference_name}.fna.bz2"))


    def matrix_markers_to_fasta(self, markers_matrix):
        """For each sample, writes the FASTA files with the sequences of the filtered markers

        Args:
            markers_matrix (pd.DataFrame):

        Returns:
            str: the temporary folder where the FASTA files were written
        """
        markers_tmp_dir = os.path.join(self.tmp_dir, "{}.StrainPhlAn4".format(self.clade))
        create_folder(markers_tmp_dir)
        filtered_sample_names = set(markers_matrix.index)
        filtered_markers = set(markers_matrix.columns)
        self.copy_filtered_references(markers_tmp_dir, filtered_sample_names)
        filtered_sample_paths = [sample_path for sample_path in self.samples
                                 if Strainphlan.sample_path_to_name(sample_path) in filtered_sample_names]
        execute_pool(((Strainphlan.sample_markers_to_fasta, sample_path, filtered_markers, self.trim_sequences,
                       markers_tmp_dir) for sample_path in filtered_sample_paths), self.nprocs)
        return markers_tmp_dir


    @classmethod
    def sample_markers_to_fasta(cls, sample_path, filtered_markers, trim_sequences, markers_tmp_dir):
        """Writes a FASTA file with the filtered clade markers of a sample

        Args:
            sample_path (str): the path to the sample
            filtered_markers:
            trim_sequences:
            markers_tmp_dir (str): the temporary folder were the FASTA file is written
        """
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
            clade_markers_file = os.path.join(self.tmp_dir, "{}.fna".format(self.clade))
        elif self.clade_markers_file.endswith(".bz2"):
            clade_markers_file = decompress_bz2(self.clade_markers_file, self.tmp_dir)
        else:
            clade_markers_file = self.clade_markers_file

        return execute_pool(((Strainphlan.process_reference, reference, self.tmp_dir, clade_markers_file,
                              self.clade_markers_names, self.trim_sequences)
                             for reference in self.references), self.nprocs)


    @classmethod
    def process_reference(cls, reference_path, tmp_dir, clade_markers_file, clade_markers, trim_sequences):
        """Processes each reference file and get a markers dictionary to add to the markers matrix

        Args:
            reference_path (str): path to the reference file
            tmp_dir (str): the temporary folder where the BLASTn results where saved
            clade_markers_file (str):
            clade_markers (Iterable): the list with the clade markers names
            trim_sequences:

        Returns:
            dict: the dictionary with the reference-to-markers information
        """
        if reference_path.endswith(".bz2"):
            uncompressed_refernces_dir = os.path.join(tmp_dir, "uncompressed_references")
            os.makedirs(uncompressed_refernces_dir, exist_ok=True)
            reference_path = decompress_bz2(reference_path, uncompressed_refernces_dir)

        ext_markers = cls.extract_markers_from_genome(reference_path, clade_markers_file)

        reference_markers_dir = os.path.join(tmp_dir, "reference_markers")
        os.makedirs(reference_markers_dir, exist_ok=True)

        consensus_markers = ConsensusMarkers([ConsensusMarker(m, s) for m, s in ext_markers.items()])
        reference_name = cls.sample_path_to_name(reference_path)
        consensus_markers.to_fasta(os.path.join(reference_markers_dir, f'{reference_name}.fna.bz2'),
                                   trim_ends=trim_sequences)

        markers_matrix = {'sample_name': reference_name}
        markers_matrix.update({m: int(m in ext_markers) for m in clade_markers})

        return markers_matrix


    @staticmethod
    def extract_with_btop(sseq, btop, qstart, qend, qlen, sstart, send):
        btop = re.split(r'(\d+|\D{2})', btop)[1::2]  # blast trace-back operations

        assert qend >= qstart  # query should be forward
        strand = 1 if send >= sstart else -1  # whether reverse-complemented

        qi = qstart - 1
        si = sstart - 1
        ext_s = '-' * qi
        for op in btop:
            if op.isnumeric():  # match
                op = int(op)
                b = si + strand * op
                if b == -1 and strand == -1:  # we should go to the beginning, but -1 gets interpreted as the end
                    b = 0
                    ext_s += sseq[si: b: strand] + sseq[0]
                else:
                    ext_s += sseq[si: b: strand]

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

        ext_s += '-' * (qlen - qend)

        # Check we parsed everything correctly
        assert qi == qend
        assert si + 1 - strand == send
        assert len(ext_s) == qlen

        if strand == -1:
            ext_s = str(Seq.Seq(ext_s).complement())

        return ext_s


    @classmethod
    def extract_markers_from_genome(cls, reference_file, clade_markers_file):
        """

        Args:
            reference_file (str):
            clade_markers_file (str):

        Returns:

        """
        # load the raw fasta data
        with openrt(reference_file) as f:
            input_file_data = f.read()

        # parse the fasta
        input_seqs = {seq.id: seq for seq in SeqIO.parse(io.StringIO(input_file_data), 'fasta')}

        # we need the additional btop column
        columns = 'qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send sstrand evalue bitscore btop'
        blast_params = '-task megablast ' \
                       '-word_size 28 ' \
                       '-reward 1 -penalty -3 ' \
                       '-gapopen 5 -gapextend 2 ' \
                       '-perc_identity 90 -qcov_hsp_perc 10'
        cmd = f'blastn -query {clade_markers_file} -subject - {blast_params} -num_threads {1} -outfmt "6 {columns}"'

        # run blastn and pass the raw data to stdin
        r = run_command(cmd, input=input_file_data, text=True)

        # load the blast output
        df = pd.read_csv(io.StringIO(r.stdout), sep='\t', names=columns.split(' '), dtype={'btop': str})
        # df['qcov'] = (df['qend'] - df['qstart'] + 1) / df['qlen']
        # df = df.query('qcov >= .1 and pident >= 90').copy()
        df = df.sort_values('bitscore', ascending=False)

        ext_ss = []
        for idx, row in df.iterrows():
            sseq = str(input_seqs[row['sseqid']].seq)
            ext_s = cls.extract_with_btop(sseq, **row[['btop', 'qstart', 'qend', 'qlen', 'sstart', 'send']])
            ext_ss.append(ext_s)
        df['ext_s'] = ext_ss

        def segments_overlap(sa, sb):
            """Intervals are closed, i.e. s = [s0, s1]"""
            if sa[0] > sb[1] or sb[0] > sa[1]:
                return False
            return True

        def segments_intersection(sa, sb):
            return max(sa[0], sb[0]), min(sa[1], sb[1])

        def segment_difference(sa, sb):
            si = segments_intersection(sa, sb)
            if si[0] > si[1]:
                return sa
            return (sa[0], sb[0] - 1) if sa[0] < sb[0] else (sb[1] + 1, sa[1])

        ext_markers = {}
        for query, idx in df.groupby('qseqid').groups.items():  # for each query/marker
            df_query = df.loc[idx]

            # take the best hit and consider only that contig/reference and that strand
            best_ref = df_query.iloc[0]['sseqid']
            best_strand = df_query.iloc[0]['sstrand']
            df_query = df_query.query(f'sseqid=="{best_ref}" and sstrand=="{best_strand}"')

            covered_regions = []
            ext_s = ['-'] * df_query.iloc[0]['qlen']
            for _, row in df_query.iterrows():
                reg = (row['qstart'], row['qend'])

                # trim by already covered regions
                for cr in covered_regions:
                    reg = segment_difference(reg, cr)

                if reg[0] > reg[1]:  # nothing left of the region
                    continue

                covered_regions.append(reg)

                for i, b in enumerate(row['ext_s']):
                    if i + 1 < reg[0] or i + 1 > reg[1]:
                        continue
                    assert ext_s[i] == '-'  # make sure the trimming works
                    ext_s[i] = b

                # if any(segments_overlap(reg, cr) for cr in covered_regions):
                #     continue
                #
                # # non-overlaping hit => expand
                # covered_regions.append(reg)
                # assert all('-' in [a, b] for a, b in zip(ext_s, row['ext_s']))  # make sure they really don't overlap
                # ext_s = [a if b == '-' else b for a, b in zip(ext_s, row['ext_s'])]

            ext_markers[query] = ''.join(ext_s)

        return ext_markers


    def calculate_polymorphic_rates(self):
        """Generates a file with the polymorphic rates of the species for each sample"""
        rows = []
        consensus_markers = execute_pool(((ConsensusMarkers.from_file, sample_path) for sample_path in self.samples),
                                         nprocs=self.nprocs, return_generator=True, ordered=True)
        for sample_path, cm in zip(self.samples, consensus_markers):
            p_stats, p_count, m_len = [], 0, 0
            for marker in cm.consensus_markers:
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


    def write_info(self, markers_matrix):
        """Writes the information file for the execution"""
        with open(os.path.join(self.output_dir, "{}.info".format(self.clade)), 'w') as info_file:
            info_file.write("Clade: {}\n".format(self.clade))
            info_file.write("Number of samples: {}\n".format(len(self.samples)))
            info_file.write("Number of references: {}\n".format(len(self.references)))
            info_file.write("Number of available markers for the clade: {}\n".format(len(self.clade_markers_names)))
            info_file.write("Filtering parameters:\n")
            info_file.write("\tNumber of bases to remove when trimming markers: {}\n".format(self.trim_sequences))
            info_file.write(f"\tMinimum number of markers to make a sample primary: "
                            f"{self.sample_with_n_markers}\n")
            info_file.write(f"\tMinimum percentage of markers to make a sample primary: "
                            f"{self.sample_with_n_markers_perc}\n")
            info_file.write(f"\tMinimum number of markers to keep a sample after filtering: "
                            f"{self.sample_with_n_markers_after_filt}\n")
            info_file.write(f"\tMinimum percentage of markers to keep a sample after filtering: "
                            f"{self.sample_with_n_markers_after_filt_perc}\n")
            info_file.write(f"\tMinimum percentage of samples to keep a marker: {self.marker_in_n_samples_perc}\n")
            info_file.write("Number of markers selected after filtering: {}\n".format(len(markers_matrix.columns)))
            n_samples = len([sample for sample in self.samples
                             if Strainphlan.sample_path_to_name(sample) in markers_matrix.index])
            info_file.write("Number of samples after filtering: {}\n".format(n_samples))
            n_refs = len([reference for reference in self.references
                          if Strainphlan.sample_path_to_name(reference) in markers_matrix.index])
            info_file.write("Number of references after filtering: {}\n".format(n_refs))
            info_file.write("PhyloPhlan phylogenetic precision mode: {}\n".format(self.phylophlan_mode))
            info_file.write("Number of processes used: {}\n".format(self.nprocs))


    def detect_clades(self):
        """Checks the clades that can be reconstructed from the pkl files

        Returns:
            dict: dictionary containing the number of samples a clade can be reconstructed from
        """
        markers2clade = self.database_controller.get_markers2clade()
        clade2markers = self.database_controller.get_clade2markers()
        sample2markers = {}
        clades_to_check = set()
        all_markers = set()
        info('Processing samples...')
        consensus_markers = execute_pool(((ConsensusMarkers.from_file, sample_path) for sample_path in self.samples),
                                         nprocs=self.nprocs, return_generator=True, ordered=True)
        for sample_path, cm in zip(self.samples, consensus_markers):
            sample_markers = [marker.name for marker in cm.consensus_markers
                              if marker.name in markers2clade and marker.breadth >= self.breadth_thres]
            sample2markers[sample_path] = sample_markers
            all_markers.update(sample_markers)
            clades_to_check.update((markers2clade[m] for m in sample_markers))

        info('Constructing the big marker matrix')
        all_samples = list(sample2markers.keys())
        all_markers = list(all_markers)
        markers_matrix_big = np.zeros(shape=(len(all_samples), len(all_markers)), dtype=int)
        for i, sample in enumerate(all_samples):
            sample_markers = set(sample2markers[sample])
            for j, marker in enumerate(all_markers):
                if marker in sample_markers:
                    markers_matrix_big[i, j] = 1

        markers_matrix_big = pd.DataFrame(markers_matrix_big, index=all_samples, columns=all_markers)

        info(f'Checking {len(clades_to_check)} species')
        species2samples = {}
        for clade in clades_to_check:
            markers_for_clade = markers_matrix_big.columns.intersection(clade2markers[clade])
            markers_matrix = markers_matrix_big.loc[:, markers_for_clade]
            markers_matrix_filtered = self.filter_markers_matrix(markers_matrix, messages=False)
            if len(markers_matrix_filtered) >= 4:
                species2samples[clade] = markers_matrix_filtered.index
        info('Done.')
        return species2samples


    def print_clades(self):
        """Prints the clades detected in the reconstructed markers"""
        species2samples = self.detect_clades()
        info(f'Detected {len(species2samples)} clades: ')
        with open(os.path.join(self.output_dir, 'print_clades_only.tsv'), 'w') as f:
            f.write('Clade\tNumber_of_samples\tSamples\n')
            for species, samples in sorted(species2samples.items(), key=lambda kv: len(kv[1]), reverse=True):
                info('\t{}: in {} samples.'.format(species, len(samples)))
                f.write('{}\t{}\t{}\n'.format(species, len(samples),
                                              ','.join(map(Strainphlan.sample_path_to_name, samples))))
        info('Done.')


    def interactive_clade_selection(self):
        """Allows the user to interactively select the SGB-level clade when specifying the clade at the species level"""
        if not self.non_interactive:
            info("The clade has been specified at the species level, starting interactive clade selection...")
        species2sgbs = self.database_controller.get_species2sgbs()
        if self.clade not in species2sgbs:
            error('The specified species "{}" is not present in the database. Exiting...'.format(self.clade), exit=True)
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


    def filter_markers_samples(self):
        """Retrieves the filtered markers matrix with the filtered samples and references

        """
        info("Getting markers from samples...")
        markers_matrix = self.get_markers_matrix_from_samples()
        info("Done.")
        if len(self.references) > 0:
            info("Getting markers from references...")
            markers_matrix += self.get_markers_from_references()
            info("Done.")
        info("Removing markers / samples...")
        # df with index samples and columns markers
        markers_matrix = pd.DataFrame.from_records(markers_matrix, index='sample_name')
        markers_matrix_filtered = self.filter_markers_matrix(markers_matrix, messages=True)
        info("Done.")

        return markers_matrix_filtered


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
        markers_matrix = self.filter_markers_samples()
        info("Done.")
        info("Writing samples as markers' FASTA files...")
        samples_as_markers_dir = self.matrix_markers_to_fasta(markers_matrix)
        info("Done.")
        info("Calculating polymorphic rates...")
        self.calculate_polymorphic_rates()
        info("Done.")
        info("Computing phylogeny...")
        self.phylophlan_controller.compute_phylogeny(samples_as_markers_dir, len(markers_matrix), self.tmp_dir)
        info("Done.")
        info("Writing information file...")
        self.write_info(markers_matrix)
        info("Done.")
        if not self.debug:
            info("Removing temporary files...")
            rmtree(self.tmp_dir, ignore_errors=False, onerror=None)
            info("Done.")


    @staticmethod
    def get_input_samples(args_samples):
        samples = []
        for s in args_samples:
            if os.path.isfile(s):
                samples.append(s)
            elif os.path.isdir(s):
                dir_files = [os.path.join(s, f) for f in os.listdir(s)]
                samples.extend([f for f in dir_files if os.path.isfile(f)])
            elif not os.path.exists(s):
                error(f'Sample file/folder {s} does not exist', exit=True)
            else:
                error(f'Neither file nor directory: {s}', exit=True)

        return samples


    def __init__(self, args):
        self.clade_markers_names = None
        self.database_controller = MetaphlanDatabaseController(args.database)
        self.clade_markers_file = args.clade_markers
        self.samples = Strainphlan.get_input_samples(args.samples)
        self.references = Strainphlan.get_input_samples(args.references)
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

    sample_names = [Strainphlan.sample_path_to_name(s) for s in args.samples]
    ref_names = [Strainphlan.sample_path_to_name(s) for s in args.references]
    if len(sample_names) + len(ref_names) != len(set(sample_names + ref_names)):
        error('Some sample or reference names are duplicated', exit=True)


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
    info("Finish StrainPhlAn {} execution ({} seconds): Results are stored at "
         "\"{}\"".format(__version__, round(exec_time, 2), args.output_dir))


if __name__ == '__main__':
    main()
    
