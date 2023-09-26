__author__ = 'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it)'
__version__ = '4.1.0'
__date__ = '23 Aug 2023'

import os
import pickle
import bz2
import pickletools
from hashlib import sha256

from Bio import SeqIO
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

    def to_seq_record(self, trim_sequences=0):
        """Gets FASTA sequence as a biopython object

        Args:
            trim_sequences (int, optional): The number of nt to trim from both ends of the sequence. Defaults to 0.

        Returns:
            SeqRecord: the parsed and trimmed sequence
        """
        marker_name = self.parse_marker_name()
        return SeqRecord(Seq(self.sequence[trim_sequences:-trim_sequences].replace("*", "-")), id=marker_name, description=marker_name)

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


    def to_dict(self):
        return {"marker": self.name, "breath": self.breadth, "avg_depth": self.avg_depth, "sequence": self.sequence}


    @classmethod
    def from_dict(cls, d):
        cls(d['marker'], d['sequence'], breadth=d['breath'], avg_depth=d['avg_depth'] if 'avg_depth' in d else None)


    def __init__(self, name, sequence, breadth=None, avg_depth=None):
        self.name = name
        self.sequence = sequence
        if breadth is None:
            self.breadth = self.get_breadth()
        else:
            self.breadth = breadth
        self.avg_depth = avg_depth


class ConsensusMarkers:
    """ConsensusMarkers class"""

    @classmethod
    def from_pkl(cls, pkl_file):
        """Init from PKL file

        Args:
            pkl_file (str): the path to the PKL file
        """
        sample_as_pkl = pickle.load(bz2.BZ2File(pkl_file)) if os.path.splitext(
            pkl_file)[1] == ".bz2" else pickle.load(open(pkl_file, "rb"))
        return cls([ConsensusMarker.from_dict(marker) for marker in sample_as_pkl])


    def to_pkl(self, output_path):
        marker_dicts = [marker.to_dict() for marker in self.consensus_markers]
        with open(output_path, 'wb') as markers_pkl:
            markers_pkl.write(pickletools.optimize(pickle.dumps(marker_dicts, pickle.HIGHEST_PROTOCOL)))


    def to_fasta(self, output_file, trim_ends=0):
        """Writes the consensus markers to FASTA"""
        seq_records = [m.to_seq_record(trim_ends) for m in self.consensus_markers]
        with open(output_file) as f:
            SeqIO.write(seq_records, f, 'fasta')


    def __init__(self, consensus_markers):
        self.consensus_markers = consensus_markers
