__author__ = 'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it)'
__version__ = '4.0.3'
__date__ = '24 Oct 2022'


import os
import pickle
import bz2
from Bio import SeqIO
try:
    from .external_exec import generate_markers_fasta
    from .util_fun import info, error
except ImportError:
    from external_exec import generate_markers_fasta
    from util_fun import info, error


class MetaphlanDatabaseController:
    """MetaphlanDatabaseController class"""

    def load_database(self, verbose=True):
        """Loads the MetaPhlAn PKL database"""
        if self.database_pkl is None:
            if verbose:
                info('Loading MetaPhlAn {} database...'.format(self.get_database_name()))
            self.database_pkl = pickle.load(bz2.BZ2File(self.database))
            if verbose:
                info('Done.')

    def get_database_name(self):
        """Gets database name

        Returns:
            str: the database name
        """
        return self.database.split('/')[-1][:-4]

    def get_markers2species(self):
        """Retrieve information from the MetaPhlAn database

        Returns:
            dict: the dictionary assigning markers to clades
        """
        self.load_database()
        return {marker: self.database_pkl['markers'][marker]['clade'] for marker in self.database_pkl['markers']}

    def get_filtered_markers(self, clades):
        """Retrieve the markers belonging to a list of clades

        Args:
            clades (list): the list of clades

        Returns:
            list: the list of markers from the input clades
        """
        self.load_database()
        return [marker for marker in self.database_pkl['markers'] if self.database_pkl['markers'][marker]['clade'] in clades]

    def get_species2sgbs(self):
        """Retrieve information from the MetaPhlAn database

        Returns:
            dict: the dictionary with the SGBs spanning each species
        """
        self.load_database()
        species2sgbs = {}
        sgb2size = self.get_sgbs_size()
        for taxa in self.database_pkl['taxonomy']:
            if taxa.split('|')[-2] not in species2sgbs:
                species2sgbs[taxa.split('|')[-2]] = {}
            species2sgbs[taxa.split('|')[-2]][taxa.split('|')
                                              [-1]] = sgb2size[taxa.split('|')[-1]]
        return species2sgbs

    def extract_markers(self, clades, output_dir):
        """Extracts markers as a FASTA file for the specified clades

        Args:
            clades (list): the list with the clades to extract markers from
            output_dir (str): the output folder
        """
        info('\tExtracting markers from the Bowtie2 database...')
        fasta_markers = generate_markers_fasta(self.database, output_dir)
        info('\tDone.')
        for clade in clades:
            markers = self.get_filtered_markers([clade])
            if len(markers) == 0:
                exit_value = True if len(clades) == 1 else False
                error('No markers were found for the clade "{}".'.format(
                    clade), exit=exit_value)
            else:
                info('\tNumber of markers for the clade "{}": {}'.format(
                    clade, len(markers)))
                info('\tExporting markers for clade {}...'.format(clade))
                with open(os.path.join(output_dir, "{}.fna".format(clade)), 'w') as ofile:
                    for rec in SeqIO.parse(open(fasta_markers, 'r'), 'fasta'):
                        if rec.name in markers:
                            SeqIO.write(rec, ofile, 'fasta')
                info('\tDone.')
        info('\tRemoving temporal FASTA files...')
        os.remove(fasta_markers)
        info('\tDone.')

    def resolve_database(self, database):
        """Resolves the path to the MPA database

        Args:
            database (str): the name or path of the database

        Returns:
            str: the resolved path to the database
        """
        if database == 'latest':
            if os.path.exists(os.path.join(self.default_db_folder, 'mpa_latest')):
                with open(os.path.join(self.default_db_folder, 'mpa_latest'), 'r') as mpa_latest:
                    return '{}/{}.pkl'.format(self.default_db_folder, [line.strip() for line in mpa_latest if not line.startswith('#')][0])
            else:
                error('The default MetaPhlAn database cannot be found at: {}'.format(
                    os.path.join(self.default_db_folder, 'mpa_latest')), exit=True)
        else:
            return database

    def resolve_index(self):
        """Resolves the name to the MPA database

        Args:
            database (str): the name or path of the database

        Returns:
            str: the resolved name to the database
        """
        return self.database.split('/')[-1][:-4]

    def get_sgbs_size(self):
        """Returns the size of the SGBs in the database

        Returns:
            dict: the dictionary with the size of each SGB
        """
        with bz2.open(os.path.join(self.mpa_script_folder, '{}_size.txt.bz2'.format(self.get_database_name())), 'rt') as rf:
            sgb2size = {line.strip().split('\t')[0]: int(
                line.strip().split('\t')[1]) for line in rf}
        return sgb2size

    def __init__(self, database, bowtie2db=None):
        self.mpa_script_folder = os.path.dirname(os.path.abspath(__file__))
        if bowtie2db == None:
            self.default_db_folder = os.path.join(
                self.mpa_script_folder, "..", "metaphlan_databases")
            self.default_db_folder = os.environ.get(
                'METAPHLAN_DB_DIR', self.default_db_folder)
        else:
            self.default_db_folder = bowtie2db
        self.database = self.resolve_database(database)
        self.database_pkl = None
