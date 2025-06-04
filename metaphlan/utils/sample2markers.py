#!/usr/bin/env python
__author__ = (
    'Michal Puncochar (michal.puncochar@unitn.it), '
    'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it), '
    'Duy Tin Truong (duytin.truong@unitn.it), '
    'Francesco Asnicar (f.asnicar@unitn.it), '
    'Moreno Zolfo (moreno.zolfo@unitn.it), '
    'Francesco Beghini (francesco.beghini@unitn.it)'
)
__version__ = '4.2.2'
__date__ = '4 Jun 2025'



import argparse as ap
import bz2
import os
import subprocess as sp
import tempfile
import time
import itertools as it
from collections import Counter
from dataclasses import dataclass
from shutil import rmtree
from typing import TextIO

import numpy as np
import pysam
import scipy.stats as sps

try:
    from .external_exec import samtools_sam_to_bam, samtools_sort_bam_v1, decompress_bz2
    from .util_fun import info, error, warning
    from .parallelisation import execute_pool
    from .database_controller import StrainphlanDatabaseController
    from .consensus_markers import ConsensusMarker, ConsensusMarkers
except ImportError:
    from external_exec import samtools_sam_to_bam, samtools_sort_bam_v1, decompress_bz2
    from util_fun import info, error, warning
    from parallelisation import execute_pool
    from database_controller import StrainphlanDatabaseController
    from consensus_markers import ConsensusMarker, ConsensusMarkers


@dataclass
class MapperSpecificDefaults:
    min_reads_aligning: int
    min_mapping_quality: int


class GenericSamFilter:
    def __init__(self, min_mapping_quality, markers_subset):
        self.min_mapping_quality = min_mapping_quality
        self.markers_subset = markers_subset


    def filter_mapping_line(self, aln):
        marker = aln.reference_name

        if marker.startswith('VDB'):
            return False

        if self.markers_subset is not None and marker not in self.markers_subset:
            return False

        if aln.is_secondary or aln.is_qcfail or aln.is_unmapped:
            return False

        if aln.mapping_quality < self.min_mapping_quality:
            return False

        return True


class Bowtie2SamFilter(GenericSamFilter):
    pass


class MinimapSamFilter(GenericSamFilter):
    def __init__(self, min_mapping_quality, markers_subset, max_gcsd):
        super().__init__(min_mapping_quality, markers_subset)
        self.max_gcsd = max_gcsd

    def filter_mapping_line(self, aln):
        if not super().filter_mapping_line(aln):
            return False

        tags = dict(aln.tags)
        gcsd = tags['de']

        if gcsd > self.max_gcsd:
            return False

        return True



class SampleToMarkers:
    MAPPER_SPECIFIC_DEFAULTS = {
        'bowtie2': MapperSpecificDefaults(min_reads_aligning=8, min_mapping_quality=10),
        'minimap2': MapperSpecificDefaults(min_reads_aligning=1, min_mapping_quality=50),
        None: MapperSpecificDefaults(min_reads_aligning=1, min_mapping_quality=0),
    }

    class DEFAULTS:
        depth_avg_q = 0.2
        quasi_marker_frac = 0.33
        min_breadth = 80
        min_base_coverage = 1
        min_base_quality = 30
        max_gcsd = 0.1
        poly_dominant_frq_thrsh = 0.8

    class CONSTANTS:
        """Settings not modifiable through arguments"""
        allowed_bases = set(list('ACTGactg'))
        pileup_stepper = 'nofilter'
        poly_error_rate = 0.001
        poly_pvalue_threshold = 0.05


    def decompress_from_bz2(self):
        """ Decompressed SAM.BZ2 files

        Returns:
            (list, str): tuple with the list of decompressed files and their format
        """
        results = execute_pool(((SampleToMarkers.decompress_bz2_file, i, self.tmp_dir)
                                for i in self.input), self.nprocs)
        decompressed = [r[0] for r in results]
        decompressed_format = [r[1] for r in results]
        if decompressed_format[1:] == decompressed_format[:-1]:
            if decompressed_format[0][1:].lower() == "sam":
                return decompressed, "sam"
            elif decompressed_format[0][1:].lower() == "bam":
                return decompressed, "bam"
            else:
                error("Decompressed files are not in SAM or BAM format", exit=True)
        else:
            error("Decompressed files have different formats", exit=True)

    @staticmethod
    def decompress_bz2_file(input_file, output_dir):
        """Decompress a BZ2 file and returns the decompressed file and the decompressed file format

        Args:
            input_file (str): the input file to decompress
            output_dir (str): the output directory

        Returns:
            (str, str): tuple with the decompressed file and its format
        """
        decompressed_file = decompress_bz2(input_file, output_dir)
        return decompressed_file, os.path.splitext(decompressed_file)[1]

    def convert_inputs(self):
        """Convert input sample files to sorted BAMs"""
        if self.input_format.lower() == "bz2":
            info("\tDecompressing samples...")
            self.input, self.input_format = self.decompress_from_bz2()
            info("\tDone.")
        if self.input_format.lower() == "sam":
            info("\tConverting samples to BAM format...")
            self.input = execute_pool(((samtools_sam_to_bam, i, self.tmp_dir) for i in self.input), self.nprocs)
            info("\tDone.")
            info("\tSorting BAM samples...")
            self.input = execute_pool(((samtools_sort_bam_v1, i, self.tmp_dir) for i in self.input), self.nprocs)
            info("\tDone.")
        elif not self.sorted:
            info("\tSorting BAM samples...")
            self.input = execute_pool(((samtools_sort_bam_v1, i, self.tmp_dir) for i in self.input), self.nprocs)
            info("\tDone.")

        info('\tIndexing BAM samples...')
        execute_pool(((pysam.index, i) for i in self.input), self.nprocs)
        info('\tDone.')

    def build_consensus_markers(self, filtered):
        """Gets the markers for each sample and writes the Pickle files

        Args:
            filtered (str): string to append when filtered for clades
        """
        for i in self.input:
            info("\tProcessing sample: {}".format(i))
            consensuses, coverages = self.get_consensuses_for_sample(i)
            consensuses_filtered = self.filter_consensuses(consensuses, coverages)
            if len(consensuses_filtered) == 0:
                warning(f'\t\tSkipping sample as it contains no markers after filtering')

            consensus_markers = ConsensusMarkers(consensuses_filtered, self.database_controller.get_database_name())
            output_filename = f'{os.path.splitext(os.path.basename(i))[0]}{filtered}.json.bz2'
            output_path = os.path.join(self.output_dir, output_filename)
            consensus_markers.to_json(output_path)
            info("\tDone.")


    def get_consensuses_for_sample(self, input_bam):
        """Pileup on the specified BAM file

        Args:
            input_bam (str): the path to the input BAM file

        Returns:
            tuple[dict[str, str], dict[str, np.ndarray[int]]]: the list with the marker names and consensus sequences
        """

        stepper = SampleToMarkers.CONSTANTS.pileup_stepper
        error_rate = SampleToMarkers.CONSTANTS.poly_error_rate
        poly_pvalue_threshold = SampleToMarkers.CONSTANTS.poly_pvalue_threshold

        info('\t\tLoading the bam file and extracting information...')
        sam_file = pysam.AlignmentFile(input_bam, check_sq=False)
        if not sam_file.references:
            warning(f'\t\t\tSample {input_bam} has empty BAM file, skipping')
            return {}, {}


        all_markers = sam_file.references
        marker_lengths = sam_file.lengths
        marker_to_length = dict(zip(all_markers, marker_lengths))
        info('\t\tDone.')

        info('\t\tRunning the pileup...')
        consensuses = {m: bytearray(b'-' * marker_to_length[m]) for m in all_markers}
        coverages = {m: np.zeros(marker_to_length[m], dtype=int) for m in all_markers}
        for base_pileup in sam_file.pileup(contig=None, stepper=stepper, min_base_quality=self.min_base_quality):
            marker = base_pileup.reference_name
            pos = base_pileup.pos
            bases = [b.upper() for b in base_pileup.get_query_sequences()
                     if b in SampleToMarkers.CONSTANTS.allowed_bases]
            base_coverage = len(bases)

            coverages[marker][pos] = base_coverage

            if base_coverage < self.min_base_coverage:
                continue

            base_frequencies = Counter(bases)
            max_frequency = max(base_frequencies.values())
            max_ratio = max_frequency / base_coverage
            base_consensus = base_frequencies.most_common(1)[0][0]

            if max_ratio < self.dominant_frq_threshold:
                p_value = sps.binom.cdf(max_frequency, base_coverage, 1 - error_rate)
                if p_value <= poly_pvalue_threshold:
                    base_consensus = '*'  # mask out likely polymorphic positions

            consensuses[marker][pos] = ord(base_consensus)

        info('\t\tDone.')

        consensuses = {m: c.decode() for m, c in consensuses.items()}  # convert bytearrays to strings
        return consensuses, coverages


    def filter_consensuses(self, consensuses, coverages):
        """
        Filters the markers for colliding quasi-markers and the breadth of coverage.
        Also computes the depth of coverage

        Args:
            consensuses (dict[str, str]): dictionary marker name => sequence
            coverages: (dict[str, np.ndarray[int]]): dictionary marker name => per-base coverages

        Returns:
            list[ConsensusMarker]:
        """
        markers2ext = self.database_controller.get_markers2ext()
        markers2clade = self.database_controller.get_markers2clade()
        clade2nmarkers = Counter(markers2clade.values())

        markers_all = list(consensuses.keys())
        clades_all = [markers2clade[m] for m in markers_all]
        clades_counter = Counter(clades_all)  # number of marker with non-zero coverage for each clade

        markers_non_quasi = [m for m in markers_all if len(markers2ext[m]) == 0]
        markers_quasi = [m for m in markers_all if len(markers2ext[m]) > 0]
        # if any external SGB has >= quasi_marker_frac non-zero markers ==> discard
        markers_quasi_filtered = [m for m in markers_quasi
                                  if not any(clades_counter[ext_sgb] / clade2nmarkers[ext_sgb] >= self.quasi_marker_frac
                                             for ext_sgb in markers2ext[m])]

        consensus_markers_filtered = []
        for m in markers_non_quasi + markers_quasi_filtered:
            seq = consensuses[m]
            cov = coverages[m]
            assert len(seq) == len(cov)

            # calculate robust average by trimming upper and lower quantiles at depth_avg_q
            cov = sorted(cov)
            t = int(self.depth_avg_q * len(cov))
            cov_trimmed = cov[t:-t]
            avg_depth = np.mean(cov_trimmed)

            c = ConsensusMarker(m, seq, avg_depth=avg_depth)
            if c.breadth < self.breadth_threshold:
                continue
            consensus_markers_filtered.append(c)

        return consensus_markers_filtered


    @classmethod
    def parallel_filter_sam(cls, input_file, tmp_dir, input_format, min_mapping_quality, min_reads_aligning,
                            max_gcsd, filtered_markers, all_markers, db_name):
        """
        Filters an input SAM file
            * filters out viral markers (VDB)
            * filters out hits with low mapQ (argument --min_mapping_quality)
            * filters out markers with not enough mapped reads (argument --min_reads_aligning)

        Args:
            input_file (str): the input SAM file
            tmp_dir:
            input_format:
            min_reads_aligning:
            min_mapping_quality:
            max_gcsd:
            filtered_markers (set): the set with the markers of the filtered clades, None if to use all markers
            all_markers:
            db_name:

        Returns:
            str: the path to the output file
        """
        output_file = os.path.join(tmp_dir, input_file.split('/')[-1])
        ifn: TextIO
        if input_format.lower() == "bz2":
            ifn = bz2.open(input_file, 'rt')
            output_file = output_file.replace('.bz2', '')
        else:
            assert input_format.lower() == "sam"
            ifn = open(input_file, 'rt')

        # Pass through header
        mapper_name = None
        first_aln_line = None
        header_lines = []
        for line in ifn:
            if line.startswith('@'):
                if line.startswith('@CO'):
                    for item in line[len('@CO'):].strip().split('\t'):
                        if item.startswith('index:'):
                            db_name_sam = item[len('index:'):].strip().strip('"')
                            if db_name_sam != db_name:
                                error(f'The database of the sample {db_name_sam} does not match {db_name}', exit=True)
                elif line.startswith('@PG'):
                    for item in line[len('@PG'):].strip().split('\t'):
                        if item.startswith('PN:'):
                            mapper_name = item[len('PN:'):].strip().strip('"')
                header_lines.append(line)
            else:  # already in the alignment part
                first_aln_line = line
                break


        # Resolve mapper specific default arguments and sam filtering
        if min_reads_aligning is None or min_mapping_quality is None:
            if mapper_name is None:
                warning('Could not infer the mapper from the SAM file header, using default parameters')
            elif mapper_name not in cls.MAPPER_SPECIFIC_DEFAULTS:
                warning(f'The mapper name {mapper_name} is not recognized, using default parameters')
                mapper_name = None
            else:
                info(f'Setting default parameters for mapper {mapper_name}')

            if min_reads_aligning is None:
                min_reads_aligning = cls.MAPPER_SPECIFIC_DEFAULTS[mapper_name].min_reads_aligning
            if min_mapping_quality is None:
                min_mapping_quality = cls.MAPPER_SPECIFIC_DEFAULTS[mapper_name].min_mapping_quality

        if mapper_name == 'bowtie2':
            sam_filter = Bowtie2SamFilter(min_mapping_quality, filtered_markers)
        elif mapper_name == 'minimap2':
            sam_filter = MinimapSamFilter(min_mapping_quality, filtered_markers, max_gcsd)
        else:
            sam_filter = GenericSamFilter(min_mapping_quality, filtered_markers)


        # Pass through alignment lines just to count hits
        header = pysam.AlignmentHeader.from_text(''.join(header_lines))
        marker_to_reads = Counter()
        if first_aln_line is not None:
            for line in it.chain([first_aln_line], ifn):
                aln = pysam.AlignedSegment.fromstring(line, header)
                if sam_filter.filter_mapping_line(aln):
                    marker = aln.reference_name
                    marker_to_reads[marker] += 1
                    if marker not in all_markers:
                        error(f'Marker {marker} not in the metaphlan database', exit=True)

        selected_markers = set((m for m, c in marker_to_reads.items() if c >= min_reads_aligning))
        sam_filter.markers_subset = selected_markers


        # Second pass of the file to write the filtered output SAM file
        ifn.seek(0)
        with open(output_file, 'wt') as ofn:
            for line in ifn:
                line_fields = line.rstrip('\n').split('\t')
                if line.startswith('@SQ'):
                    assert line_fields[1].startswith('SN:')  # in bowtie2 output SN is always the first tag
                    marker = line_fields[1][3:]
                    if marker in selected_markers:
                        ofn.write(line)
                elif line.startswith('@'):  # other header lines like @HD and @PG
                    ofn.write(line)
                else:  # mapping lines
                    aln = pysam.AlignedSegment.fromstring(line, header)
                    if sam_filter.filter_mapping_line(aln):
                        ofn.write(line)

        ifn.close()
        return output_file


    def filter_sam_files(self):
        """Filters the input SAM files with the hits against markers of specific clades and low quality reads"""
        filtered_markers = self.database_controller.get_filtered_markers(self.clades) if len(self.clades) > 0 else None
        all_markers = set(self.database_controller.get_all_markers())
        db_name = self.database_controller.get_database_name()
        self.input = execute_pool(((SampleToMarkers.parallel_filter_sam, i, self.tmp_dir, self.input_format,
                                    self.min_mapping_quality, self.min_reads_aligning, self.max_gcsd, filtered_markers,
                                    all_markers, db_name)
                                   for i in self.input), self.nprocs)
        self.input_format = 'sam'
        self.sorted = False


    def run_sample2markers(self):
        """Runs the full sample2markes pipeline"""
        info("Creating temporary directory...")
        self.tmp_dir = tempfile.mkdtemp(dir=self.tmp_dir)
        info("Done.")
        if self.input_format in ['sam', 'bz2']:
            info("Filtering SAM files...")
            self.filter_sam_files()
            info("Done.")
        info("Converting input files...")
        self.convert_inputs()
        info("Done.")
        info("Getting consensus markers from samples...")
        self.build_consensus_markers(filtered='_filtered' if len(self.clades) > 0 else '')
        info("Done.")
        if not self.debug:
            info("Removing temporary files...")
            rmtree(self.tmp_dir)
            info("Done.")


    def __init__(self, args):
        self.input = args.input
        self.output_dir = args.output_dir
        self.database_controller = StrainphlanDatabaseController(args.database)
        self.tmp_dir = args.output_dir if args.tmp is None else args.tmp
        self.breadth_threshold = args.breadth_threshold
        self.input_format = args.input_format
        self.clades = args.clades
        self.sorted = args.sorted
        self.min_reads_aligning = args.min_reads_aligning
        self.max_gcsd = args.max_gcsd
        self.min_base_coverage = args.min_base_coverage
        self.min_base_quality = args.min_base_quality
        self.min_mapping_quality = args.min_mapping_quality
        self.dominant_frq_threshold = args.dominant_frq_threshold
        self.quasi_marker_frac = args.quasi_marker_frac
        self.depth_avg_q = args.depth_avg_q
        self.debug = args.debug
        self.nprocs = args.nprocs


def read_params():
    """ Reads and parses the command line arguments of the script

    Returns:
        namespace: The populated namespace with the command line arguments
    """
    p = ap.ArgumentParser(description="", formatter_class=ap.ArgumentDefaultsHelpFormatter, add_help=False)
    # required
    p.add_argument('-i', '--input', type=str, nargs='+', required=True, help="The input samples as SAM or BAM files")
    p.add_argument('-o', '--output_dir', type=str, required=True, help="The output directory")

    # optional
    p.add_argument('-d', '--database', type=str, default='latest',
                   help="The input MetaPhlAn " + __version__ + " database (path to the pkl file)")
    p.add_argument('--clades', type=str, nargs='+', default=[],
                   help="Restricts the reconstruction of the markers to the specified clades")
    p.add_argument('-f', '--input_format', type=str, default="bz2", help="The input samples format {bam, sam, bz2}")
    p.add_argument('--sorted', action='store_true', default=False, help="Whether the BAM input files are sorted")
    p.add_argument('--min_reads_aligning', type=int, default=None,
                   help="The minimum number of reads to cover a marker."
                        "Default 8 for bowtie2, 1 for minimap2 nad 1 otherwise.")
    p.add_argument('--min_base_coverage', type=int, default=SampleToMarkers.DEFAULTS.min_base_coverage,
                   help="The minimum depth of coverage for a base to be considered")
    p.add_argument('--min_base_quality', type=int, default=SampleToMarkers.DEFAULTS.min_base_quality,
                   help="The minimum quality for a base to be considered. This is performed BEFORE --min_base_coverage")
    p.add_argument('--min_mapping_quality', type=int, default=None,
                   help="The minimum quality for a mapping of the read to be considered. "
                        "Default 10 for bowtie2, 50 for minimap2 and 0 otherwise.")
    p.add_argument('--max_gcsd', type=float, default=SampleToMarkers.DEFAULTS.max_gcsd,
                   help="The maximum gap-compressed sequence divergence threshold to use in case of Minimap2 mapper.")
    p.add_argument('--dominant_frq_threshold', type=float, default=SampleToMarkers.DEFAULTS.poly_dominant_frq_thrsh,
                   help="The cutoff for degree of 'allele dominance' for a position to be considered polymorphic")
    p.add_argument('--quasi_marker_frac', type=float, default=SampleToMarkers.DEFAULTS.quasi_marker_frac,
                   help="Fraction [0-1] of markers with a hit of an external SGB to disqualify a quasi-marker.")
    p.add_argument('--depth_avg_q', type=float, default=SampleToMarkers.DEFAULTS.depth_avg_q,
                   help="A quantile to cut from both ends of the coverage distributions to calculate robust average.")
    p.add_argument('--tmp', type=str, default=None,
                   help="If specified, the directory where to store the temporary files. "
                        "Otherwise the output directory will be used.")
    p.add_argument('-b', '--breadth_threshold', type=int, default=SampleToMarkers.DEFAULTS.min_breadth,
                   help="The breadth of coverage threshold for the consensus markers")
    p.add_argument('--debug', action='store_true', default=False,
                   help="If specified, StrainPhlAn will not remove the temporary folder. "
                        "Not available with inputs in BAM format")
    p.add_argument('-n', '--nprocs', type=int, default=1, help="The number of threads to execute the script")
    p.add_argument('-v', '--version', action='version',
                   version=f"StrainPhlAn sample2markers version {__version__} ({__date__})")
    p.add_argument('-h', '--help', action='help', help="Show help.")

    return p.parse_args()


def check_params(args):
    """Checks the mandatory command line arguments of the script

    Args:
        args (namespace): the arguments to check
    """
    if not os.path.exists(args.output_dir):
        error('The directory {} does not exist'.format(args.output_dir), exit=True)
    if args.tmp is not None and not os.path.exists(args.tmp):
        error('The directory {} does not exist'.format(args.tmp), exit=True)
    if args.database != 'latest' and not os.path.exists(args.database):
        error('The database does not exist', exit=True)
    if args.input_format.lower() not in ['bam', 'sam', 'bz2']:
        error('The input format must be SAM, BAM, or compressed in BZ2 format', exit=True)
    if args.input_format.lower() == "bam" and len(args.clades) > 0:
        error('The --clades option cannot be used with inputs in BAM format', exit=True)

    check_input_files(args.input, args.input_format)

    return args


def check_input_files(input_files, input_format):
    """Checks the format input sample files

    Args:
        input_files (list): list containing all the input files
        input_format (str): the format to check
    """
    for s in input_files:
        _, extension = os.path.splitext(s)
        if not os.path.exists(s):
            error('The input file \"{}\" does not exist'.format(s), exit=True)
        elif not input_format.lower() == extension[1:].lower():
            error('The the input file \"{}\" must be in \"{}\" format'.format(
                s, input_format.upper()), exit=True)


def check_samtools():
    samtools_version_output = sp.run("samtools", capture_output=True).stderr.decode().split('\n')
    for line in samtools_version_output:
        if line.startswith('Version:'):
            samtools_version_info = line
            break
    else:
        error("Couldn't get information about the samtools version", exit=True)
        return

    samtools_version = samtools_version_info.rstrip('\n').split(' ')[1]
    samtools_version_parts = [int(x) for x in samtools_version.split('.')]
    if samtools_version_parts < [1]:
        error(f'The samtools version required is >= 1.x.x, installed is {samtools_version}', exit=True)
    info(f'Using samtools version {samtools_version}')


def main():
    t0 = time.time()
    args = read_params()
    info("Start samples to markers execution")
    check_samtools()
    check_params(args)
    sampletomarkers = SampleToMarkers(args)
    sampletomarkers.run_sample2markers()
    exec_time = time.time() - t0
    info("Finish samples to markers execution ({} seconds): Results are stored at "
         "\"{}\"".format(round(exec_time, 2), args.output_dir))


if __name__ == '__main__':
    main()
