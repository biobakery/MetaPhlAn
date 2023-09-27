__author__ = 'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it'
__version__ = '4.1.0'
__date__ = '23 Aug 2023'

import abc
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

    @abc.abstractmethod
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


    def compute_phylogeny(self, samples_markers_dir, num_samples, tmp_dir):
        """Executes PhyloPhlAn to compute phylogeny

        Args:
            samples_markers_dir (str): Directory containing the samples as FASTA
            num_samples (int): Minimum number of samples containing a marker
            tmp_dir (str): Temporal directory
        """
        if not self.phylophlan_configuration:
            info("\tGenerating PhyloPhlAn configuration file...")
            self.phylophlan_configuration = generate_phylophlan_config_file(
                tmp_dir, self.get_phylophlan_configuration())
            info("\tDone.")
        info("\tProcessing samples...")
        execute_phylophlan(samples_markers_dir, self.phylophlan_configuration, tmp_dir, self.output_dir,
                           self.phylophlan_mode, self.phylophlan_params, self.mutation_rates, self.nprocs)
        if self.mutation_rates:
            move(os.path.join(self.output_dir, "mutation_rates.tsv"), os.path.join(
                self.output_dir, "{}.mutation".format(self.clade)))
            move(os.path.join(self.output_dir, "mutation_rates"), os.path.join(
                self.output_dir, "{}_mutation_rates".format(self.clade)))
        info("\tDone.")
        if self.treeshrink:
            info("Executing TreeShrink...")
            execute_treeshrink(os.path.join(self.output_dir, 'RAxML_bestTree.{}.StrainPhlAn4.tre'.format(
                self.clade)), self.output_dir, tmp=tmp_dir, centroid=True if num_samples >= 100 else False)
            info("Done.")
            info('This StrainPhlAn analysis ran TreeShrink, do not forget to cite:')
            info('Mai, Uyen, and Siavash Mirarab. 2018. “TreeShrink: Fast and Accurate Detection of Outlier Long Branches in Collections of Phylogenetic Trees.” BMC Genomics 19 (S5): 272. https://doi.org/10.1186/s12864-018-4620-2.')

    def __init__(self, args):
        self.samples = args.samples
        self.secondary_samples = args.secondary_samples
        self.marker_in_n_samples = args.marker_in_n_samples
        self.abs_n_samples_thres = args.abs_n_samples_thres
        self.phylophlan_mode = args.phylophlan_mode
        self.phylophlan_params = args.phylophlan_params
        self.phylophlan_configuration = args.phylophlan_configuration
        self.mutation_rates = args.mutation_rates
        self.output_dir = args.output_dir
        self.clade = args.clade
        self.treeshrink = args.treeshrink
        self.nprocs = args.nprocs
