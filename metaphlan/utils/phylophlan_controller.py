__author__ = 'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it),' \
             'Michal Puncochar (michal.puncochar@unitn.it)'
__version__ = '4.1.1'
__date__ = '11 Mar 2024'

import abc
import os
from shutil import move

try:
    from .util_fun import info, error
    from .external_exec import generate_phylophlan_config_file, execute_treeshrink, run_command
except ImportError:
    from util_fun import info, error
    from external_exec import generate_phylophlan_config_file, execute_treeshrink, run_command


class PhylophlanController(abc.ABC):
    """PhyloPhlAnController interface class"""

    @abc.abstractmethod
    def compute_phylogeny(self, *args, **kwargs):
        """Executes PhyloPhlAn to compute phylogeny"""
        pass


    @staticmethod
    def get_phylophlan_configuration():
        """Gets PhyloPhlAn configuration

        Returns:
            dict: The dictionary with the PhyloPhlAn configuration
        """
        configuration = {
            'map': 'blastn',
            'aligner': 'mafft',
            'tree1': 'raxml'
        }
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
            self.phylophlan_configuration = generate_phylophlan_config_file(tmp_dir, self.get_phylophlan_configuration())
            info("\tDone.")
        info("\tExecuting PhyloPhlAn...")
        self.execute_phylophlan(samples_markers_dir, tmp_dir)
        if self.mutation_rates:
            move(os.path.join(self.output_dir, "mutation_rates.tsv"), os.path.join(
                self.output_dir, "{}.mutation".format(self.clade)))
            move(os.path.join(self.output_dir, "mutation_rates"), os.path.join(
                self.output_dir, "{}_mutation_rates".format(self.clade)))
        info("\tDone.")
        if self.treeshrink:
            info("\tExecuting TreeShrink...")
            execute_treeshrink(os.path.join(self.output_dir, 'RAxML_bestTree.{}.StrainPhlAn4.tre'.format(
                self.clade)), self.output_dir, tmp=tmp_dir, centroid=True if num_samples >= 100 else False)
            info("\tDone.")
            info('\tThis StrainPhlAn analysis ran TreeShrink, do not forget to cite:')
            info('\tMai, Uyen, and Siavash Mirarab. 2018. “TreeShrink: Fast and Accurate Detection of Outlier Long'
                 ' Branches in Collections of Phylogenetic Trees.” BMC Genomics 19 (S5): 272.'
                 ' https://doi.org/10.1186/s12864-018-4620-2.')


    def execute_phylophlan(self, samples_markers_dir, tmp_dir):
        """
        Executes PhyloPhlAn

        Args:
            samples_markers_dir:
            tmp_dir:

        Returns:

        """
        cmd = f'phylophlan' \
              f' -i {samples_markers_dir} -o . --output_folder {self.output_dir} --nproc {self.nprocs}' \
              f' --strainphlan --{self.phylophlan_mode} --data_folder {tmp_dir} -f {self.phylophlan_configuration}' \
              f' -t n --diversity low --genome_extension fna --min_num_entries 1 --min_num_markers 1' \
              f' --fragmentary_threshold 1.0 --not_variant_threshold 1.0 --gap_perc_threshold 0.8'

        if self.phylophlan_params is not None:
            cmd += " " + self.phylophlan_params
        if self.mutation_rates:
            cmd += " --mutation_rates"

        r = run_command(cmd)
        with open(os.path.join(tmp_dir, "phylophlan_log.stdout"), 'wb') as f:
            f.write(r.stdout)
        with open(os.path.join(tmp_dir, "phylophlan_log.stderr"), 'wb') as f:
            f.write(r.stderr)


    def __init__(self, args):
        self.samples = args.samples
        self.phylophlan_mode = args.phylophlan_mode
        self.phylophlan_params = args.phylophlan_params
        self.phylophlan_configuration = args.phylophlan_configuration
        self.mutation_rates = args.mutation_rates
        self.output_dir = args.output_dir
        self.clade = args.clade
        self.treeshrink = args.treeshrink
        self.nprocs = args.nprocs
