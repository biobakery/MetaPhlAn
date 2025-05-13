__author__ = 'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it)'
__version__ = '4.1.1'
__date__ = '11 Mar 2024'


import os
import pickle
import bz2
import pathlib
from collections.abc import Iterable

import pandas as pd
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
            with bz2.open(self.database, 'rb') as f:
                self.database_pkl = pickle.load(f)
            if verbose:
                info('Done.')


    def get_database_name(self):
        """Gets database name

        Returns:
            str: the database name
        """
        return self.database.split('/')[-1][:-4]


    def get_markers2clade(self):
        """
        Get a dictionary mapping marker names to clade

        Returns:
            dict[str, str]: the dictionary assigning markers to clades
        """
        self.load_database()
        return {marker_name: marker_info['clade'] for marker_name, marker_info in self.database_pkl['markers'].items()}


    def get_clade2markers(self):
        """
        Get a dictionary mapping clades to lists of markers

        Returns:
            dict[str, Iterable]:
        """
        markers2clade = self.get_markers2clade()
        markers2clade = pd.Series(markers2clade)
        return markers2clade.groupby(markers2clade).groups


    def get_markers_for_clade(self, clade):
        """

        Args:
            clade (str):

        Returns:
            set: marker names for the given clade

        """
        db_sgbs = self.get_all_sgbs()
        if clade not in db_sgbs:
            error(f'Clade {clade} not found in the database', exit=True)

        self.load_database()
        return set(marker_name for marker_name, marker_info in self.database_pkl['markers'].items()
                   if marker_info['clade'] == clade)


    def get_all_markers(self):
        self.load_database()
        return list(self.database_pkl['markers'].keys())

    def get_all_sgbs(self):
        self.load_database()
        return set(tax.split('|')[-1] for tax in self.database_pkl['taxonomy'].keys())

    def get_markers2ext(self):
        self.load_database()
        return {marker_name: ['t__' + sgb for sgb in marker_info['ext']]
                for marker_name, marker_info in self.database_pkl['markers'].items()}


    def get_filtered_markers(self, clades):
        """Retrieve the markers belonging to a list of clades

        Args:
            clades (Iterable): the list of clades

        Returns:
            set: the list of markers from the input clades
        """
        self.load_database()

        db_sgbs = self.get_all_sgbs()
        for clade in clades:
            if clade not in db_sgbs:
                error(f'Clade {clade} not found in the database', exit=True)

        clades_set = set(clades)
        return set((marker for marker in self.database_pkl['markers']
                    if self.database_pkl['markers'][marker]['clade'] in clades_set))

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
            species2sgbs[taxa.split('|')[-2]][taxa.split('|')[-1]] = sgb2size[taxa.split('|')[-1]]
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
        db_sgbs = self.get_all_sgbs()
        for clade in clades:
            if clade not in db_sgbs:
                error(f'Clade {clade} not found in the database', exit=True)
            markers = self.get_filtered_markers([clade])
            info(f'\tExporting {len(markers)} markers for clade {clade}...')
            with open(os.path.join(output_dir, f"{clade}.fna"), 'w') as ofile:
                for rec in SeqIO.parse(open(fasta_markers, 'r'), 'fasta'):
                    if rec.name in markers:
                        SeqIO.write(rec, ofile, 'fasta')
            info('\tDone.')
        info('\tRemoving temporary FASTA files...')
        os.remove(fasta_markers)
        info('\tDone.')

    def resolve_database(self, database):
        """Resolves the path to the MPA database

        Args:
            database (pathlib.Path|str): the name or path of the database

        Returns:
            str: the resolved path to the database
        """
        database = str(database)
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
        """

        Args:
            database (pathlib.Path|str):
            bowtie2db:
        """
        self.mpa_script_folder = os.path.dirname(os.path.abspath(__file__))
        if bowtie2db is None:
            self.default_db_folder = os.path.join(
                self.mpa_script_folder, "..", "metaphlan_databases")
            self.default_db_folder = os.environ.get(
                'METAPHLAN_DB_DIR', self.default_db_folder)
        else:
            self.default_db_folder = bowtie2db
        self.database = self.resolve_database(database)
        self.database_pkl = None
