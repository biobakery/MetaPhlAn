__author__ = 'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it)'
__version__ = '4.0.3'
__date__ = '24 Oct 2022'


import os
import pickle
import bz2
from hashlib import sha256
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


class ConsensusMarker:
    """ConsensusMarker class"""

    def parse_marker_name(self):
        """Parses the marker name to avoid long marker names issue

        Returns:
            str: the parsed marker name
        """
        return str(int(sha256(self.name.encode('utf-8')).hexdigest(), 16) % 10**12)

    def get_sequence(self, trim_sequences=0):
        """Gets FASTA sequence as a biopython object

        Args:
            trim_sequences (int, optional): The number of nt to trim from both ends of the sequence. Defaults to 0.

        Returns:
            SeqRecord: the parsed and trimmed sequence
        """
        marker_name = self.parse_marker_name()
        return SeqRecord(Seq(self.sequence[trim_sequences:-trim_sequences].replace("*", "-").replace('-', 'N')), id=marker_name, description=marker_name)

    def get_polymorphisms(self):
        """ Gets the number of polymorphic positions in the marker

        Returns:
            int: the number of polymorphisms in the marker sequence
        """
        return self.sequence.count('*')

    def get_sequence_length(self):
        """Gets the length of the marker

        Returns:
            int: the length of the marker sequence
        """
        return len(self.sequence)

    def get_polymorphism_perc(self):
        """Gets the percentage of the polymorphic positions in the marker

        Returns:
            float: the percentage of the marker sequence containing polymorphisms
        """
        return self.sequence.count('*') * 100 / self.get_sequence_length()

    def get_breadth(self):
        """Returns the breadth of coverage of the marker

        Returns:
            float: the breadth of coverage of the marker sequence
        """
        seq_len = len(self.sequence)
        return ((seq_len - self.sequence.count('N') - self.sequence.count('*') - self.sequence.count('-')) * 100) / seq_len

    def __init__(self, name, sequence, breadth=None):
        self.name = name
        self.sequence = sequence
        if breadth == None:
            self.breadth = self.get_breadth()
        else:
            self.breadth = breadth


class ConsensusMarkers:
    """ConsensusMarkers class"""

    def to_fasta(self, output_file):
        """Writes the consensus markers to FASTA"""
        pass

    def from_pkl(self, pkl_file):
        """Init from PKL file

        Args:
            pkl_file (str): the path to the PKL file
        """
        self.consensus_markers = []
        sample_as_pkl = pickle.load(bz2.BZ2File(pkl_file)) if os.path.splitext(
            pkl_file)[1] == ".bz2" else pickle.load(open(pkl_file, "rb"))
        for marker in sample_as_pkl:
            self.consensus_markers.append(ConsensusMarker(
                marker['marker'], marker['sequence'], marker['breath']))

    def __init__(self, pkl_file=None):
        if pkl_file != None:
            self.from_pkl(pkl_file)
        else:
            self.consensus_markers = []
