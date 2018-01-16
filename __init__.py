# from . import metaphlan2
import metaphlan2
from ._metaphlan2 import profile_single_fastq
from ._metaphlan2 import profile_paired_fastq


__author__ = metaphlan2.__author__
__version__ = metaphlan2.__version__
__date__ = metaphlan2.__date__

__all__ = ['profile_single_fastq', 'profile_paired_fastq']
