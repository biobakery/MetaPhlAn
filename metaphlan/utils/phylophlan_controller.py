__author__ = 'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it'
__version__ = '4.0.3'
__date__ = '24 Oct 2022'


import os
import bz2
from shutil import move, copy
try:
    from .util_fun import info, error
    from .external_exec import generate_phylophlan_config_file, create_phylophlan_db, execute_phylophlan, execute_treeshrink
except ImportError:
    from util_fun import info, error
    from external_exec import generate_phylophlan_config_file, create_phylophlan_db, execute_phylophlan, execute_treeshrink


class PhylophlanController:
    """PhyloPhlAnController interface class"""

    def compute_phylogeny(self):
        """Executes PhyloPhlAn to compute phylogeny"""
        pass

    def get_phylophlan_configuration(self):
        """Gets PhyloPhlAn configuration

        Returns:
            dict: The dictionary with the PhyloPhlAn configuration
        """
        configuration = {}
        configuration['map'] = 'blastn'
        configuration['aligner'] = 'mafft'
        configuration['trim'] = 'trimal'
        configuration['tree1'] = 'raxml'
        return configuration

class Phylophlan3Controller(PhylophlanController):
    """PhyloPhlAn3Controller class"""

    def fake_phylophlan_inputs(self, samples_markers_dir):
        """Fakes the PhyloPhlAn inputs for the reconstructed markers

        Args:
            samples_markers_dir (str): Directory containing the samples as FASTA
        """
        samples = [s.split('/')[-1].replace('.pkl', '')
                   for s in self.samples+self.secondary_samples]
        os.mkdir(os.path.join(self.tmp_dir, 'markers_dna'))
        os.mkdir(os.path.join(self.tmp_dir, 'map_dna'))
        os.mkdir(os.path.join(self.tmp_dir, 'clean_dna'))
        for sample in os.listdir(samples_markers_dir):
            if sample.replace('.fna', '') in samples:
                open(os.path.join(self.tmp_dir, 'map_dna',
                     sample.replace('.fna', '.b6o.bkp')), 'a').close()
                open(os.path.join(self.tmp_dir, 'map_dna',
                     sample.replace('.fna', '.b6o.bz2')), 'a').close()
                copy(os.path.join(samples_markers_dir, sample),
                     os.path.join(self.tmp_dir, 'clean_dna', sample))
                with bz2.open(os.path.join(self.tmp_dir, 'markers_dna', '{}.bz2'.format(sample)), 'wt') as write_file:
                    with open(os.path.join(samples_markers_dir, sample), 'r') as read_file:
                        for line in read_file:
                            write_file.write(line)

    def compute_phylogeny(self, samples_markers_dir, num_samples, min_markers, temporal_dir):
        """Executes PhyloPhlAn to compute phylogeny

        Args:
            samples_markers_dir (str): Directory containing the samples as FASTA
            num_samples (int): Minimum number or samples containing a marker
            min_markers (int): Mininum number of markers per sample
            temporal_dir (str): Temporal directory
        """
        info("\tCreating PhyloPhlAn database...")
        self.tmp_dir = temporal_dir
        create_phylophlan_db(self.tmp_dir, self.clade[:30])
        info("\tDone.")
        if not self.phylophlan_configuration:
            info("\tGenerating PhyloPhlAn configuration file...")
            self.phylophlan_configuration = generate_phylophlan_config_file(
                self.tmp_dir, self.get_phylophlan_configuration())
            info("\tDone.")
        self.fake_phylophlan_inputs(samples_markers_dir)
        info("\tProcessing samples...")
        min_entries = self.marker_in_n_samples if self.abs_n_markers_thres else self.marker_in_n_samples * num_samples // 100
        execute_phylophlan(samples_markers_dir, self.phylophlan_configuration, min_entries, int(round(min_markers, 0)),
                           self.tmp_dir, self.output_dir, self.clade, self.phylophlan_mode, self.mutation_rates, self.nprocs)
        if self.mutation_rates:
            move(os.path.join(self.output_dir, "mutation_rates.tsv"), os.path.join(
                self.output_dir, "{}.mutation".format(self.clade)))
            move(os.path.join(self.output_dir, "mutation_rates"), os.path.join(
                self.output_dir, "{}_mutation_rates".format(self.clade)))
        info("\tDone.")
        if self.treeshrink:
            info("Executing TreeShrink...")
            execute_treeshrink(os.path.join(self.output_dir, 'RAxML_bestTree.{}.StrainPhlAn4.tre'.format(
                self.clade)), self.output_dir, tmp=self.tmp_dir, centroid=True if num_samples >= 100 else False)
            info("Done.")
            info('This StrainPhlAn analysis ran TreeShrink, do not forget to cite:')
            info('Mai, Uyen, and Siavash Mirarab. 2018. “TreeShrink: Fast and Accurate Detection of Outlier Long Branches in Collections of Phylogenetic Trees.” BMC Genomics 19 (S5): 272. https://doi.org/10.1186/s12864-018-4620-2.')

    def __init__(self, args):
        self.samples = args.samples
        self.secondary_samples = args.secondary_samples
        self.marker_in_n_samples = args.marker_in_n_samples
        self.abs_n_markers_thres = args.abs_n_markers_thres
        self.phylophlan_mode = args.phylophlan_mode
        self.phylophlan_configuration = args.phylophlan_configuration
        self.mutation_rates = args.mutation_rates
        self.output_dir = args.output_dir
        self.clade = args.clade
        self.treeshrink = args.treeshrink
        self.nprocs = args.nprocs
