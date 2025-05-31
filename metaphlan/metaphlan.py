#!/usr/bin/env python
__author__ = ('Aitor Blanco-Miguez (aitor.blancomiguez@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it), '
              'Nicola Segata (nicola.segata@unitn.it), '
              'Duy Tin Truong, '
              'Francesco Asnicar (f.asnicar@unitn.it), '
              'Claudia Mengoni (claudia.mengoni@unitn.it), '
              'Linda Cova (linda.cova@unitn.it)')
__version__ = '4.2.0'
__date__ = '14 May 2025'


import time
import bz2, gzip
import re
import os
from glob import glob
import sys
import stat
import random
import tempfile
import argparse as ap
import subprocess as subp
import numpy as np 

import pickle as pkl
from collections import Counter
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pysam
import pandas as pd
import shutil
from packaging import version
import biom
import biom.table
import json

from collections import defaultdict as defdict
try:
    from .utils import *
except ImportError:
    from utils import *

metaphlan_script_install_folder = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DB_FOLDER = os.path.join(metaphlan_script_install_folder, "metaphlan_databases")
DEFAULT_DB_FOLDER= os.environ.get('METAPHLAN_DB_DIR', DEFAULT_DB_FOLDER)

class TaxClade:
    """TaxClade class"""
    
    def add_child(self, name, tax_id ):
        """Adds a new child to a taxonomy

        Args:
            name (str): The taxonomy label of the child
            tax_id (str): The taxID of the child
        """
        child = TaxClade( name, tax_id, self)
        self.children[name] = child
    
    def get_full_taxids( self ):
        """Returns the full taxID

        Returns:
            str: the full taxID
        """
        fullname = ['']
        if self.tax_id:
            fullname = [self.tax_id]
        clade = self.father
        while clade:
            fullname = [clade.tax_id] + fullname
            clade = clade.father
        return "|".join(fullname[1:])

    def get_full_name( self ):
        """Returns the full taxonomy label

        Returns:
            str: the full taxonomy label
        """
        fullname = [self.name]
        clade = self.father
        while clade:
            fullname = [clade.name] + fullname
            clade = clade.father
        return "|".join(fullname[1:])
    
    
    def normalize_coverage(self, rat_nreads):
        """Normalizes the clade coverage based on the selected stat

        Args:
            rat_nreads (list): the list of the filtered markers as a tuple marker length and number of reads mapping
        """
        rat_nreads = sorted(rat_nreads, key = lambda x: x[1])  
        rat_v,nreads_v = zip(*rat_nreads)                
        quant = int(self.stat_q * len(rat_nreads))
        ql,qr,qn = (quant, -quant, quant)       
        if self.stat == 'avg_g' or (not qn and self.stat in ['wavg_g','tavg_g']):
            self.coverage = sum(nreads_v) / float(sum(rat_v))
        elif self.stat == 'avg_l' or (not qn and self.stat in ['wavg_l','tavg_l']):
            self.coverage = np.mean([float(n)/(np.absolute(r - self.avg_read_length) + 1) for r,n in rat_nreads])
        elif self.stat == 'tavg_g':
            wnreads = sorted([(float(n)/(np.absolute(r-self.avg_read_length)+1),(np.absolute(r - self.avg_read_length)+1) ,n) for r,n in rat_nreads], key=lambda x:x[0])
            den,num = zip(*[v[1:] for v in wnreads[ql:qr]]) #if wnreads[ql:qr] else ([],[])
            self.coverage = float(sum(num))/float(sum(den)) if any(den) else 0.0
        elif self.stat == 'tavg_l':
            self.coverage = np.mean(sorted([float(n)/(np.absolute(r - self.avg_read_length) + 1) for r,n in rat_nreads])[ql:qr])
        elif self.stat == 'wavg_g':
            vmin, vmax = nreads_v[ql], nreads_v[qr]
            wnreads = [vmin]*qn+list(nreads_v[ql:qr])+[vmax]*qn
            self.coverage = float(sum(wnreads)) / float(sum(rat_v))
        elif self.stat == 'wavg_l':
            wnreads = sorted([float(n)/(np.absolute(r - self.avg_read_length) + 1) for r,n in rat_nreads])
            vmin, vmax = wnreads[ql], wnreads[qr]
            wnreads = [vmin]*qn+list(wnreads[ql:qr])+[vmax]*qn
            self.coverage = np.mean(wnreads)
        elif self.stat == 'med':
            self.coverage = np.median(sorted([float(n)/(np.absolute(r - self.avg_read_length) +1) for r,n in rat_nreads])[ql:qr])

    def estimate_number_reads(self, rat_nreads):
        """Estimates the number of reads mapping to the clade

        Args:
            rat_nreads (list): the list of the filtered makers as a tuple marker length and number of reads mapping
        """
        self.nreads = int(np.mean([float(n)/(np.absolute(r - self.avg_read_length) + 1) for r,n in rat_nreads if n > 0]) * self.glen)
    
    def filter_markers_for_ext(self):
        """Filters the makers for external hits

        Returns:
            list: the list of the filtered makers as a tuple marker length and number of reads mapping
        """
        rat_nreads = []
        for marker, nreads in sorted(self.markers2nreads.items(),key=lambda x:x[0]):
            misidentified = False
            if not self.avoid_disqm:
                for ext in self.markers2exts[marker]:
                    ext_markers2nreads = self.taxa2clades[ext].markers2nreads
                    if len(ext_markers2nreads) > 0 and float(sum([nr > 0 for nr in ext_markers2nreads.values()])) / len(ext_markers2nreads) > self.perc_nonzero:
                        misidentified = True
                        break
            if not misidentified:
                rat_nreads.append( (self.markers2lens[marker], nreads) ) 
        return rat_nreads
    
    def get_estimated_number_reads(self):  
        """Retrieves the estimated number of reads mapped to the clade   
        
        Returns:
            int: the estimated number of reads
        """      
        if self.nreads is None and len(self.children) > 0:
            self.nreads = sum([child.get_estimated_number_reads() for child in self.children.values()])
        return self.nreads
    
    def get_normalized_counts( self ):
        """Normalize markers counts as RPKs
        
        Returns:
            list: the list of the normalized marker counts as tuple (marker, normalized counts)
        """
        return [(marker, float(nreads)*1000.0/(np.absolute(self.markers2lens[marker] - self.avg_read_length) +1) )
                    for marker,nreads in self.markers2nreads.items()]
            

    def compute_coverage(self):
        """Computes the coverage of the clade

        Returns:
            float: the coverage of the clade
        """
        if self.coverage is not None: 
            pass
        elif len(self.children) > 0:
            self.coverage = sum([child.compute_coverage() for child in self.children.values()])
            self.nreads = self.get_estimated_number_reads()
        else:
            self.coverage = 0
            self.nreads = 0
            rat_nreads = self.filter_markers_for_ext()
            if sum([item[1] for item in rat_nreads]) > 0:                                    
                self.normalize_coverage(rat_nreads)
                if self.coverage > 0:
                    self.estimate_number_reads(rat_nreads)
        return self.coverage

    def __init__(self, name, tax_id, father=None):
        self.name = name
        self.tax_id = tax_id
        self.father = father
        self.children = {}
        self.markers2nreads = {}
        self.nreads = None
        self.coverage = None
        self.rel_abundance = None
        self.glen = None
        
    stat = None
    stat_q = None
    perc_nonzero = None
    avoid_disqm = None        
    markers2lens = None
    taxa2clades = None
    markers2exts = None
    avg_read_length = None
    

class TaxTree:
    """TaxTree class"""

    def add_tree_lens(self, node):
        """Adds the genome length to a node in the tree and its children
        Args:
            node (TaxClade): the clade to add the length to

        Returns:
            float: The genome length of the taxa
        """
        if not node.children:
            return node.glen
        lens = []
        for clade in node.children.values():
            lens.append(self.add_tree_lens(clade))
        node.glen = min(np.mean(lens), np.median(lens))
        return node.glen
    
    def add_tree_nodes(self):
        """Adds the nodes to the tree"""
        for clade, c_info in self.mpa['taxonomy'].items():
            clade = clade.strip().split("|")
            taxids, clade_len = c_info
            taxids = taxids.strip().split("|")
            father = self.root
            for i in range(len(clade)):
                clade_lev = clade[i]
                clade_taxid = taxids[i] if i < 8 and taxids is not None else None
                if not clade_lev in father.children:
                    father.add_child(clade_lev, tax_id=clade_taxid)
                    self.all_clades[clade_lev] = father.children[clade_lev]
                father = father.children[clade_lev]
                if clade_lev[0] == "t":
                    self.taxa2clades[clade_lev[3:]] = father
                if clade_lev[0] == "t":
                    father.glen = clade_len
                    
    def add_tree_markers(self):
        """Adds the information of the marker sequences to the tree"""
        for marker, m_info in self.mpa['markers'].items():
            if marker in self.ignore_markers:
                continue
            self.markers2lens[marker] = m_info['len']
            self.markers2clades[marker] = m_info['clade']
            self.markers2exts[marker] = m_info['ext']
            self.all_clades[m_info['clade']].markers2nreads[marker] = 0

    def build_tree(self):
        """Builds the taxonomy tree"""
        self.add_tree_nodes()
        self.add_tree_lens(self.root)
        self.add_tree_markers()
            
    def add_reads(self, marker, nreads,
                  ignore_eukaryotes = False,
                  ignore_bacteria = False, ignore_archaea = False, 
                  ignore_ksgbs = False, ignore_usgbs = False):
        """Assigns a number of reads to a marker in the database

        Args:
            marker (str): name of the marker
            nreads (int): number of read
            ignore_eukaryotes (bool, optional): whether to add hits to eukaryotic markers. Defaults to False.
            ignore_bacteria (bool, optional): whether to add hits to bacterial markers. Defaults to False.
            ignore_archaea (bool, optional): whether to add hits to archaeal markers. Defaults to False.
            ignore_ksgbs (bool, optional): whether to add hits to kSGBs markers. Defaults to False.
            ignore_usgbs (bool, optional): whether to add hits to uSGBs markers. Defaults to False.
        """
        clade = self.all_clades[self.markers2clades[marker]]
        if ignore_bacteria or ignore_archaea or ignore_eukaryotes:
            clade_name = clade.get_full_name()
            if ignore_archaea and clade_name.startswith("k__Archaea"):
                return
            if ignore_bacteria and clade_name.startswith("k__Bacteria"):
                return
            if ignore_eukaryotes and clade_name.startswith("k__Eukaryota"):
                return
        if ignore_ksgbs or ignore_usgbs:
            clade_name = clade.get_full_name()
            if ignore_ksgbs and not '_SGB' in clade_name.split('|')[-2]:
                return
            if ignore_usgbs and '_SGB' in clade_name.split('|')[-2]:
                return
        
        clade.markers2nreads[marker] = nreads
        
    def relative_abundances(self):
        """Compute the relative abundances for the taxa present in the sample

        Returns:
            int: the number of mapped reads (or bases for long reads)
        """
        total_ab = 0
        total = 0 # total instead of total_reads because it is bases for long reads
        for clade in self.all_clades.values():
            clade.compute_coverage()
            if len(clade.children) == 0 and clade.coverage > 0:
                total_ab += clade.coverage
                total += clade.nreads        
        for clade in self.all_clades.values():
            if clade.coverage > 0:
                clade.rel_abundance = round(100 * clade.coverage / total_ab, 5) if total_ab > 0 else 0
        
        if total_ab == 0:
            warning('Warning: No species were detected.', init_new_line = True)
            
        return total
    
    def clade_profiles(self):
        """Retrieves the clade profiles

        Returns:
            dict: the clade profiles dictionary
        """
        clade2profiles = {}
        for clade in self.all_clades.values():
            normalized_counts = clade.get_normalized_counts()
            if len(normalized_counts) == 0 or sum([marker[1] for marker in normalized_counts]) <= 0.0:
                continue
            clade2profiles[clade.get_full_name()] = normalized_counts
        return clade2profiles

    def __init__(self, database_pkl, ignore_markers, stat, stat_q, perc_nonzero, avoid_disqm, avg_read_length):
        self.mpa = database_pkl
        self.ignore_markers = ignore_markers
        self.root = TaxClade('root', 0) #if not long_reads else TaxClade_longreads('root', 0)
        self.all_clades = {}
        self.markers2clades = {}
        self.markers2lens = {}
        self.taxa2clades = {}
        self.markers2exts = {}
        self.build_tree()
        TaxClade.stat = stat
        TaxClade.stat_q = stat_q
        TaxClade.perc_nonzero = perc_nonzero
        TaxClade.avoid_disqm = avoid_disqm
        TaxClade.markers2lens = self.markers2lens
        TaxClade.taxa2clades = self.taxa2clades
        TaxClade.markers2exts = self.markers2exts
        TaxClade.avg_read_length = avg_read_length


class VSCController():
    """ Class for controlling the viral profiling"""

    def vsc_samtools(self):
        """Convert SAM to BAM and sort it"""
        try:
            # stv_command = ['samtools','view','-bS', os.path.basename(self.vscBamFile).split('.')[0]+'.sam','-@',str(self.nproc)]
            stv_command = ['samtools','view','-bS', self.vscSamFile,'-@',str(self.nproc)]
            sts_command = ['samtools','sort','-','-@',str(self.nproc),'-o',self.vscBamFile]
            sti_command = ['samtools','index',self.vscBamFile]

            p3 = subp.Popen(sts_command, stdin=subp.PIPE)
            p2 = subp.Popen(stv_command, stdin=subp.PIPE, stdout=p3.stdin)

            p2.communicate()
            p3.communicate()

            subp.check_call(sti_command)

        except Exception as e:
            error('Error: "{}"\nFatal error creating BAM file for Viruses\n'.format(e), init_new_line = True, exit = True) 

        if not os.path.exists(self.vscBamFile):
            error('Error:\nUnable to create BAM FILE for viruses.', init_new_line = True, exit = True) 

    def vsc_parsing(self, rpkm=True, total_metagenome=None):
        """Parsing of the output sam file from mapping viral markers to reads"""
        try:
            bamHandle = pysam.AlignmentFile(self.vscBamFile, "rb")
        except Exception as e:
            error('Error: "{}"\nCheck PySam is correctly working\n'.format(e), exit = True)
        
        VSC_report=list()
        coverage_positions = defdict(dict) 
        ref_to_len = dict(zip(bamHandle.references,bamHandle.lengths))
        for pileupcolumn in bamHandle.pileup():
            c = pileupcolumn.reference_name
            
            tCoverage = 0
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip \
                        and (self.input_type == 'fasta' or (self.input_type == 'fastq' and pileupread.alignment.query_qualities[pileupread.query_position] >= 20)) \
                        and pileupread.alignment.query_sequence[pileupread.query_position].upper() in ('A','T','C','G'):
                        tCoverage +=1
            if tCoverage >= 1:
                coverage_positions[c][pileupcolumn.pos] = tCoverage
            
        for c in coverage_positions.keys():
            length = ref_to_len[c]
            breadth = float(len(coverage_positions[c].keys()))/float(length)
            if breadth > 0:
                cvals=list(coverage_positions[c].values())
                if rpkm:
                    n_reads_map = int(bamHandle.count(c))
                    calc_rpkm = n_reads_map /  ((length/(10**3)) * (total_metagenome / (10**6)))
                    VSC_report.append({'M-Group/Cluster':c.split('|')[2].split('-')[0], 'genomeName':c, 'len':length, 'breadth_of_coverage':breadth, 'mapping_reads_count':n_reads_map, 'RPKM':calc_rpkm, 'depth_of_coverage_mean': np.mean(cvals), 'depth_of_coverage_median': np.median(cvals)})
                else:
                    VSC_report.append({'M-Group/Cluster':c.split('|')[2].split('-')[0], 'genomeName':c, 'len':length, 'breadth_of_coverage':breadth, 'depth_of_coverage_mean': np.mean(cvals), 'depth_of_coverage_median': np.median(cvals)})
        
        if bamHandle:
            bamHandle.close()
        
        self.create_vsc_report(VSC_report)

    def create_vsc_report(self, VSC_report):
        """Create a report of the VSC that mapped"""
        with open(self.vsc_out,'w') as outf:
            outf.write('#{}\n'.format(self.index))
            outf.write('#{}\n'.format(' '.join(sys.argv)))
            outf.write('#SampleID\t{}\n'.format(self.sample_id))
            vsc_out_df = pd.DataFrame.from_dict(VSC_report).query('breadth_of_coverage >= {}'.format(self.vsc_breadth))
            vsc_out_df = self.filter_vsc_report(vsc_out_df)
            vsc_info_df = pd.read_table(self.vsc_vinfo, sep='\t')
            vsc_out_df = vsc_out_df.merge(vsc_info_df, on='M-Group/Cluster').sort_values(by='breadth_of_coverage', ascending=False).set_index('M-Group/Cluster')
            vsc_out_df.to_csv(outf,sep='\t',na_rep='-')
            if vsc_out_df.shape[0] == 0:
                warning('No viral clusters detected, the output report is empty', init_new_line = True)

    def __init__(self, args, mapping_controller, database_controller):
        # parsing previous mapout run
        self.nproc = args.nproc
        self.sample_id = args.sample_id
        self.tmp_dir = tempfile.mkdtemp(dir=args.tmp_dir)
        self.input_type = args.input_type 
        
        # mapping 
        self.vsc_out = args.vsc_out
        self.vsc_breadth = args.vsc_breadth
        self.vscBamFile = os.path.join(self.tmp_dir,'v_align.bam')
        self.vscSamFile = os.path.join(self.tmp_dir,'v_align.sam')
        self.mapping_controller = mapping_controller
        self.index = self.mapping_controller.get_index()
        self.database_controller = database_controller

        #database  
        self.db_dir = args.db_dir
        self.vsc_fna = os.path.join(self.db_dir,"{}_VSG.fna".format(self.database_controller.get_index()))
        self.vsc_vinfo = os.path.join(self.db_dir,"{}_VINFO.csv".format(self.database_controller.get_index()))

class VSC_bt2_controller(VSCController):
    """VSC controller class for short reads"""

    def filter_vsc_report(self, vsc_df):
        """Returns the VSC report (no need to filter by depth the bowtie2 mapping)"""
        return vsc_df

    def extract_viral_mappings_line(self, o):
        """"Extraction of reads mapping to viral markers
        
        Args:   
            o (str): the line of the SAM file
        
        Returns:
            str: the marker group
            str: the marker cluster
            SeqRecord: the read record
        """
        mCluster = o[2]
        mGroup = o[2].split('|')[2].split('-')[0]
                            
        if (hex(int(o[1]) & 0x10) == '0x0'): #front read
            rr=SeqRecord(Seq(o[9]),letter_annotations={'phred_quality':[ord(_)-33 for _ in o[10][::-1]]}, id=o[0])
        else:
            rr=SeqRecord(Seq(o[9]).reverse_complement(),letter_annotations={'phred_quality':[ord(_)-33 for _ in o[10][::-1]]}, id=o[0])
        if o[10] == "*":
            rr.letter_annotations['phred_quality'] = None # Remove phred quality if not present in SAM

        return mGroup, mCluster, rr
    
    def get_viral_mapping(self, sam_line):
        """When reading a SAM file, save viral mappings to the VSCController variables
        
        Args:   
            sam_line (str): the line of the SAM file
        """
        mGroup, mCluster, rr = self.extract_viral_mappings_line(sam_line)
        self.viral_markers[mGroup].append(mCluster)
        self.viral_reads.append(rr)

    def infer_sam_input_type(self):
        """Infer the input type for the viral profiling"""
        if not self.input_type in ['fastq','fasta']:
            if any([bool(rr.letter_annotations['phred_quality']) for rr in self.viral_reads]):
                self.input_type = 'fastq'
            else:
                self.input_type = 'fasta'

    def process_SGB_mapping(self):
        """From reads detected in the first bowtie2 run, extracts the reads/markers for a second bowtie2 run"""
        self.infer_sam_input_type()
        SeqIO.write(self.viral_reads,self.reads_file, self.input_type) # write viral reads to remap

        VSCs_markers = SeqIO.index(self.vsc_fna, "fasta")
        selectedMarkers=[]
        for _,v in self.viral_markers.items():

            cv=Counter(v)
            if (len(cv) > 1):
                normalized_occurrencies=sorted([(mk,mk_occurrencies,len(VSCs_markers[mk].seq)) for (mk,mk_occurrencies) in cv.items()],key=lambda x: x[1]/x[2],reverse=True)
                topMarker= VSCs_markers[normalized_occurrencies[0][0]] #name of the viralGenome
            else:
                topMarker=VSCs_markers[v[0]] # the only hitted viralGenome is the winning one

            selectedMarkers.append(topMarker)
        SeqIO.write(selectedMarkers, self.top_marker_file, 'fasta') # write top viral markers to use for remapping

    def check_vsc_files(self):
        """Check if files are present to map to viral database"""
        if not os.path.exists(self.top_marker_file) or not os.path.exists(self.reads_file):
            error('There was an error in the VSCs file lookup.\n \
                It may be that there are not enough reads to profile (not enough depth).\n \
                Passing without reporting any virus.', init_new_line = True, exit=True)

        if os.stat(self.top_marker_file).st_size == 0:
            warning('No reads aligning to VSC markers in this file.', init_new_line = True, exit=True)        
    
    def initialize_mapping(self):
        """Initialize parameters for mapping to viral database"""
        self.database_controller.set_db_dir(self.tmp_dir)
        self.database_controller.set_index(os.path.basename(self.top_marker_file).split('.')[0])
        self.database_controller.prepare_indexes()

        # set parameters for bowtie mapping
        self.mapping_controller.set_inp(self.reads_file)
        self.mapping_controller.set_db_dir(self.tmp_dir)
        self.mapping_controller.set_index(os.path.basename(self.top_marker_file).split('.')[0])        
        self.mapping_controller.set_input_type(self.input_type)
        self.mapping_controller.set_samout(self.vscSamFile)
        self.mapping_controller.set_mapout(os.path.join(self.tmp_dir, 'v_mapout.txt'))
        self.mapping_controller.set_mapping_parameters(None)
        self.mapping_controller.run_mapping()

    def vsc_bowtie2(self):
        """Run bowtie2 only on viral markers and interested reads"""
        try:
            self.check_vsc_files()        
            self.initialize_mapping()

        except Exception as e:
            error('Error: "{}"\nFatal error running BowTie2 for Viruses\n'.format(e), init_new_line = True, exit = True) 

        ## check if the output file is empty
        if not os.path.exists(self.mapping_controller.samout) or os.stat(self.mapping_controller.samout).st_size == 0:
            warning('Problems when running bowtie2 on VSC reads ({}):\nSam file empty or not existing ({}).'.format(self.reads_file, self.mapping_controller.samout), init_new_line = True, exit=True)

        self.vsc_samtools()

    def run_analysis(self, total_metagenome):
        """All the steps to run analysis on VSC"""
        self.process_SGB_mapping()
        self.vsc_bowtie2()
        self.vsc_parsing(rpkm=True, total_metagenome=total_metagenome)
        shutil.rmtree(self.tmp_dir)


    def __init__(self, args, mapping_controller, database_controller):
        super().__init__(args, mapping_controller, database_controller)
        self.viral_reads = []
        self.viral_markers = defdict(list)
        self.top_marker_file = os.path.join(self.tmp_dir,'v_top_markers.fna')
        self.reads_file= os.path.join(self.tmp_dir,'v_reads.fx')

class VSC_mm2_controller(VSCController):
    """VSC controller class for long reads"""

    def infer_sam_input_type(self):
        """Infers the input type for the viral profiling"""
        if not self.input_type in ['fastq','fasta']:
            with open(self.vscSamFile) as samfile:
                for line in samfile:
                    if line.startswith('@'):
                        continue
                    line = line.strip().split('\t')
                    self.input_type = 'fasta' if line[10] == "*" else 'fastq'
                    break

    def check_vsc_files(self):
        """Checks if viral sam is present"""
        if not os.path.exists(self.vscSamFile):
            error('There was an error in the VSCs file lookup.\n \
                It may be that there are not enough reads to profile (not enough depth).\n \
                Passing without reporting any virus.', init_new_line = True, exit=True)

        if os.stat(self.vscSamFile).st_size == 0:
            warning('No viral mappings in the viral sam file: {}'.format(self.vscSamFile), init_new_line = True, exit=True)        

    def filter_vsc_report(self, vsc_df):
        """Filter the VSC report to keep the top marker by depth of coverage for each cluster
        
        Args:
            vsc_df (DataFrame): the VSC report
            
        Returns:
            DataFrame: the filtered VSC report
        """
        return vsc_df.loc[vsc_df.groupby('M-Group/Cluster')['depth_of_coverage_mean'].idxmax()].reset_index(drop=True)
    
    def run_analysis(self, total_metagenome):
        """All the steps to run analysis on VSC"""
        self.infer_sam_input_type()
        self.check_vsc_files()
        self.vsc_samtools()
        self.vsc_parsing(rpkm=False)

    def __init__(self, args, mapping_controller, database_controller):
        super().__init__(args, mapping_controller, database_controller)


class MappingController:
    """MappingController interface"""

    def set_inp(self, value):
        """Sets the input file
            
        Args:   
            value (str): the input file
        """
        self.inp = value

    def set_db_dir(self, value):
        """Sets the path to the MetaPhlAn database

        Args:   
            value (str): path to the MetaPhlAn database
        """    
        self.db_dir = value

    def set_input_type(self, value):
        """Sets the input type  

        Args:   
            value (str): the input type
        """
        self.input_type = value    
    
    def set_samout(self, value):
        """Sets the SAM output file

        Args:   
            value (str): the SAM output file
        """
        self.samout = value

    def set_mapout(self, value):
        """Sets the mapping output file

        Args:   
            value (str): the mapping output file
        """        
        self.mapout = value    

    def get_index(self):
        """Get the index of the database""" 
        return self.index
    
    def set_index(self, value):
        """Set the index of the database
            
        Args:
            value (str): the index of the database
        """
        self.index = value
    
    #def get_sample_id(self):
    #    """Get the sample id"""
    #   return self.sample_id
    
    def set_vsc_controller(self, value):
        """Set the VSC controller

        Args:
            value (VSCController): the VSC controller
        """
        self.vsc_controller = value

    def init_mapout(self):
        """Inits the mapping output file"""
        if self.no_map:
            self.mapout = tempfile.NamedTemporaryFile(dir=self.tmp_dir).name
        elif self.mapout is None:
            if self.inp is None:
                self.mapout = 'stdin_map.mapout.txt'
            elif stat.S_ISFIFO(os.stat(self.inp).st_mode):
                self.mapout = 'fifo_map.mapout.txt'
            else:
                self.mapout = 'fifo_map.mapout.txt'.format(self.inp)
    
    def __init__(self, args, index):
        self.inp = args.inp
        self.input_type = args.input_type
        self.mapout = args.mapout
        self.db_dir = args.db_dir
        self.samout = args.samout
        self.min_alignment_len = args.min_alignment_len
        self.min_mapq_val = args.min_mapq_val
        self.index = index
        self.nproc = args.nproc
        self.no_map = args.no_map
        self.read_min_len = args.read_min_len
        self.tmp_dir = args.tmp_dir
        self.vsc_controller = None
        self.verbose = args.verbose
    
class Bowtie2Controller(MappingController):
    """Bowtie2Controller class"""
    
    def set_mapping_parameters(self, value):
        self.min_alignment_len = value
        self.nreads = value
        self.min_mapq_val = value
        self.vsc_controller = value

    def check_bowtie2_database(self):
        """Checks the presence and consistency of the Bowtie2 database"""
        if glob(os.path.join(self.db_dir, "{}*.{}".format(self.index, 'bt2*'))):
            
            if glob(os.path.join(self.db_dir, "{}*.{}".format(self.index, 'bt2*')))[0].endswith('l'):
                bt2_ext = 'bt2l'
            else:
                bt2_ext = 'bt2'

            self.db_dir = os.path.join(self.db_dir, "{}".format(self.index))

        else:
            error('No MetaPhlAn BowTie2 database was found at: {}'.format(
                self.db_dir), exit=True)
        #bt2_ext = 'bt2l'
        if not all([os.path.exists(".".join([str(self.db_dir), p]))
                    for p in ["1." + bt2_ext, "2." + bt2_ext, "3." + bt2_ext, "4." + bt2_ext, "rev.1." + bt2_ext, "rev.2." + bt2_ext]]):
            error('No MetaPhlAn BowTie2 database found (--index option)!\nExpecting location {}'.format(self.db_dir), exit=True)
        if not (abs(os.path.getsize(".".join([str(self.db_dir), "1." + bt2_ext])) - os.path.getsize(".".join([str(self.db_dir), "rev.1." + bt2_ext]))) <= 1000):
            error('Partial MetaPhlAn BowTie2 database found at {}. Please remove and rebuild the database'.format(
                self.db_dir), exit=True)

    def init_mapping_arguments(self):
        """Automatically set the mapping arguments"""
        if (self.bt2_ps in ["sensitive-local", "very-sensitive-local"]) and (self.min_alignment_len is None):
            self.min_alignment_len = 100
            warning('bt2_ps is set to local mode, and min_alignment_len is None, min_alignment_len has automatically been set to 100')
        self.init_mapout()

    def get_bowtie2cmd(self):
        """Gets the command for Bowtie2 execution

        Returns:
            list: the command for the bowtie2 execution
        """
        bowtie2_cmd = [self.bowtie2_exe if self.bowtie2_exe else 'bowtie2', "--seed", "1992", "--quiet", "--no-unal", "--{}".format(self.bt2_ps),
                       "-S", "-", "-U", "-", "-x", self.db_dir]
        if int(self.nproc) > 1:
            bowtie2_cmd += ["-p", str(self.nproc)]
        if self.input_type == "fasta":
            bowtie2_cmd += ["-f"]
        return bowtie2_cmd

    def run_bowtie2(self):
        """Runs Bowtie2"""
        try:
            if self.verbose:
                info("Running BowTie2", init_new_line=True)
            read_fastx = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'utils', 'read_fastx.py')

            if self.inp:
                readin = subp.Popen([read_fastx, '-l', str(self.read_min_len), '--split_reads', str(self.split_reads), self.inp], stdout=subp.PIPE, stderr=subp.PIPE)
            else:
                readin = subp.Popen([read_fastx, '-l', str(self.read_min_len), '--split_reads', str(self.split_reads)], stdin=sys.stdin, stdout=subp.PIPE, stderr=subp.PIPE)
            p = subp.Popen(self.get_bowtie2cmd(), stdout=subp.PIPE, stdin=readin.stdout)
            readin.stdout.close()
            outf = bz2.open(self.mapout, "wt") if self.mapout.endswith(".bz2") else open(self.mapout, "wt")
            
            try:
                if self.samout:
                    if self.samout[-4:] == '.bz2':
                        sam_file = bz2.BZ2File(self.samout, 'wb')
                    else:
                        sam_file = open(self.samout, 'wb')
            except IOError as e:
                error('IOError: "{}"\nUnable to open sam output file.\n'.format(e), exit=True)

            if self.samout:
                sam_file.write(mybytes('@CO\tindex:{}\n'.format(self.index)))
            for line in p.stdout:
                if self.samout:
                    sam_file.write(line)
                o = read_and_split_line(line)
                if self.check_hq_mapping(o):
                    outf.write("\t".join([o[0], o[2].split('/')[0]]) + "\n")
                    ## check for VSC markers
                    if self.vsc_controller and o[2].startswith('VDB|'):
                        self.vsc_controller.get_viral_mapping(o)

            if self.samout:
                sam_file.close()
            p.communicate()

            nreads, avg_read_len = self.get_nreads_and_avg_rlen(readin.stderr.readlines())
            outf.write('#nreads\t{}\n'.format(nreads))
            outf.write('#avg_read_length\t{}'.format(avg_read_len))
            outf.close()
            self.input_type = 'mapout'
            self.inp = self.mapout
        except OSError as e:
            error('OSError: "{}"\nFatal error running BowTie2.'.format(e), exit=True)
        except IOError as e:
            error('IOError: "{}"\nFatal error running BowTie2.'.format(e), exit=True)
        if p.returncode == 13:
            error("Permission Denied Error: fatal error running BowTie2. Is the BowTie2 file in the path with execution and read permissions?", exit=True)
        elif p.returncode != 0:
            error("Error while running bowtie2.\n", exit=True)
            
    def check_hq_mapping(self, sam_line):
        """Checks whether a hit in the SAM file is of high quality

        Args:
            sam_line (list): SAM file line as a list

        Returns:
            bool: Whether the hit was of HQ
        """
        if not sam_line[0].startswith('@'):
            if not sam_line[2].endswith('*'):
                if (hex(int(sam_line[1]) & 0x100) == '0x0'):  # no secondary
                    if self.mapq_filter(sam_line[2], int(sam_line[4]), self.min_mapq_val):  # filter low mapq reads
                        if ((self.min_alignment_len is None) or
                                (max([int(x.strip('M')) for x in re.findall(r'(\d*M)', sam_line[5]) if x]) >= self.min_alignment_len)):
                            return True
        return False

    def mapq_filter(self, marker_name, mapq_value, min_mapq_val):
        """Checks whether a mapping hit pass a mapq quality filter

        Args:
            marker_name (str): the marker name of the hit
            mapq_value (int): the mapq value of the hit
            min_mapq_val (int): the mapq threshold

        Returns:
            bool: whether the hit passed the mapq filter
        """
        if 'GeneID:' in marker_name  or 'VDB' in marker_name:
            return True
        else:
            if mapq_value > min_mapq_val:
                return True
        return False  
    
    def get_nreads_and_avg_rlen(self, read_fastx_stderr):
        """Gets the number of reads and the average read lenght from the readfastx execution

        Args:
            read_fastx_stderr (list): standard error from the readfastx execution
        """
        try:
            nreads, avg_read_len, _ = list(
                map(float, read_fastx_stderr[0].decode().split()))
            if not nreads or int(nreads)==0:
                error('Fatal error running MetaPhlAn. Total metagenome size was not estimated or is zero.\nPlease check your input files.', exit=True)
            if not avg_read_len:
                error('Fatal error running MetaPhlAn. The average read length was not estimated.\nPlease check your input files.', exit=True)
            return int(nreads), avg_read_len
        except ValueError:
            os.unlink(self.mapout)
            error(b''.join(read_fastx_stderr).decode(), exit=True)
    
    def get_reads2markers(self):
        """Retrieves the reads to markers information from the mapping results

        Returns:
            dict: dictionary containing the reads to mapping results (dict of lists)
        """
        if not self.inp:
            ras, ras_line, inpf = plain_read_and_split, plain_read_and_split_line, sys.stdin
        elif self.inp.endswith(".bz2"):
            ras, ras_line, inpf = read_and_split, read_and_split_line, bz2.BZ2File(self.inp, "r")
        else:
            ras, ras_line, inpf = plain_read_and_split, plain_read_and_split_line, open(self.inp)
        reads2markers = defdict(list) ## dict of lists to allow multiple markers per read (for long reads)
        if self.input_type == 'mapout':
            for r, c in ras(inpf):
                if r.startswith('#') and 'nreads' in r:
                    nreads = int(c)
                elif r.startswith('#') and 'avg_read_length' in r:
                    avg_read_len = float(c)
                else:
                    reads2markers[r].append(c)
        elif self.input_type == 'sam':
            nreads = int(self.nreads)
            read_lengths = []
            for line in inpf:
                sam_line = ras_line(line)
                if self.check_hq_mapping(sam_line):
                    reads2markers[sam_line[0]].append(sam_line[2].split('/')[0])
                    read_lengths.append(len(sam_line[9]))
                    ## check for VSC markers
                    if self.vsc_controller and sam_line[2].startswith('VDB|'):
                        self.vsc_controller.get_viral_mapping(sam_line)
            avg_read_len = sum(read_lengths) / len(read_lengths)
        inpf.close()
        return nreads, avg_read_len, reads2markers, None

    def run_mapping(self):
        """Runs all the steps of the Bowtie2 mapping"""
        self.check_bowtie2_database()
        self.init_mapping_arguments()
        self.run_bowtie2()

    def __init__(self, args, index):
        super().__init__(args, index)
        self.bt2_ps = args.bt2_ps
        self.bowtie2_exe = args.bowtie2_exe
        self.bowtie2_build = args.bowtie2_build
        self.nreads = args.nreads
        self.split_reads = args.split_readlen if args.split_reads else 0

class Minimap2Controller(MappingController):
    """Minimap2Controller class"""

    def set_nbases(self, value):
        self.nbases = value

    def build_mm2_index(self): 
        """Build mm2 index with the selected parameters"""
        # Check fasta file
        fna_file = os.path.join(self.db_dir, self.index + ".fna")
        if not glob(fna_file):
            error('No Metaphlan database found at: {} to build the minimap index'.format(fna_file), exit=True)
                            
        # Build mm2 index
        if self.verbose:
            info('Building minimap2 index for parameters:{}'.format(self.mm2_ps_str.replace("_"," -")), init_new_line = True)
        mmi_index = os.path.join(self.db_dir, self.index + self.mm2_ps_str + ".mmi")
        # mm2_cmd = [self.minimap2_exe]+self.mm2_ps_list+["-d", mmi_index, fna_file]
        mm2_cmd = [self.minimap2_exe]+self.mm2_ps_list+["-v","0","-I","20G","-d", mmi_index, fna_file]

        try:
            subp.run(mm2_cmd)
        except Exception as e:
            error("Fatal error running '{}'\nError message: '{}'\n\n".format(' '.join(mm2_cmd), e), exit = True)

        # check mm2 index        
        try:
            os.chmod(mmi_index, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)  # change permissions to 664
        except PermissionError as e:
            error('PermissionError: "{}"\nCannot change permission for {}. Make sure the files are readable.'.format(e, mmi_index), exit=True)

        return mmi_index

    def check_mm2_database(self): 
        """Checks the presence and consistency of the mpa database"""
        mmi_index = os.path.join(self.db_dir, self.index + self.mm2_ps_str + ".mmi")

        if not glob(mmi_index):
            ## check minimap2 
            try:
                subp.check_call([self.minimap2_exe, "-h"], stdout=subp.DEVNULL)
            except Exception as e:
                error('OSError: "{}"\nFatal error running minimap2 at {}. Please check minimap2 installation and path'.format(e,self.minimap2_exe), exit=True)
            ## build the index
            self.db_dir = self.build_mm2_index()

        elif os.path.getsize(mmi_index) < (30*10**9):
            warning('Small minimap index found at: {}, please check and remove if needed'.format(mmi_index), init_new_line = True)
        else:
            info("Minimap2 index found")
            self.db_dir = mmi_index

    def init_mapping_arguments(self):
        """Automatically set the mapping arguments"""
        ## add here default changes in mapping arguments if needed for mm2
        self.init_mapout()

    def get_minimap2cmd(self):
        """Gets the command for Minimap2 execution

        Returns:
            list: the command for the Minimap2 execution
        """
        ## old version for input in file
        # sp_tmp = ".".join(self.mapout.split(".")[:-2]) if self.mapout.endswith("bz2") else ".".join(self.mapout.split(".")[-1])
        # mm2_cmd = [self.minimap2_exe, "-x", "asm20", "-B", "3", "-O", "3,12", "--sam-hit-only", "--split-prefix", sp_tmp, "-a", self.db_dir, self.inp]
        # mm2_cmd = [self.minimap2_exe]+self.mm2_ps_list+["-v","0","--sam-hit-only", "--split-prefix", sp_tmp, "-a", self.db_dir, self.inp]
        ## version for input in stdin
        mm2_cmd = [self.minimap2_exe]+self.mm2_ps_list+["-v","0","--sam-hit-only", "-a", self.db_dir]
        if int(self.nproc) > 3:
            mm2_cmd += ["-t", str(self.nproc)]

        return mm2_cmd + ["-"] # "-" needed for stdin input
    
    def get_gscd(self, sam_line):
        """Gets the GCSD value from the sam line
        
        Returns:
            float: the GCSD value
        """
        gcsd = [float(x.split(":")[-1]) for x in sam_line[11:] if x.startswith("de")]
        if len(gcsd):
            return gcsd[0]
        else:
            error('Problem with the SAM file: --long_reads option is only tested for SAM produced by Minimap2', exit=True)


    def run_minimap2(self):
        """Runs Minimap2"""
        try:
            if self.verbose:
                info("Running Minimap2", init_new_line=True)
            read_fastx = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'utils', 'read_fastx.py')

            if self.inp:
                readin = subp.Popen([read_fastx, '-l', str(self.read_min_len), self.inp], stdout=subp.PIPE, stderr=subp.PIPE)
            else:
                readin = subp.Popen([read_fastx, '-l', str(self.read_min_len)], stdin=sys.stdin, stdout=subp.PIPE, stderr=subp.PIPE)
            p = subp.Popen(self.get_minimap2cmd(), stdout=subp.PIPE, stdin=readin.stdout) #, stderr=subp.DEVNULL)
            readin.stdout.close()

            ## run without read_fastx
            # p = subp.Popen(self.get_minimap2cmd(), stdout=subp.PIPE, stderr=subp.DEVNULL)
            outf = bz2.open(self.mapout, "wt") if self.mapout.endswith(".bz2") else open(self.mapout, "wt")
            
            try:
                if self.samout:
                    if self.samout[-4:] == '.bz2':
                        sam_file = bz2.BZ2File(self.samout, 'wb')
                    else:
                        sam_file = open(self.samout, 'wb')
                if self.vsc_controller:
                    viral_sam = open(self.vsc_controller.vscSamFile, 'wb')
            except IOError as e:
                error('IOError: "{}"\nUnable to open sam output file.\n'.format(e), exit=True)
            
            if self.samout:
                sam_file.write(mybytes('@CO\tindex:{}\n'.format(self.index))) 
            for line in p.stdout:
                if self.samout:
                    sam_file.write(line)
                o = read_and_split_line(line)
                # check for VSC markers
                if self.vsc_controller and ((o[2].startswith('VDB|')) or (o[0].startswith('@') and 'VDB|' in o[1])):
                    viral_sam.write(line)
                if self.check_hq_mapping(o):
                    gcsd = self.get_gscd(o)
                    nb = sum([int(x[:-1]) for x in re.findall(r"\d*[DM]", o[5])])
                    outf.write("\t".join([o[0], o[2].split('/')[0], str(nb), str(gcsd), str(len(o[9]))]) + "\n")
                    ## outfile line: readname, marker, bases_covered, GCSD, read_length

            if self.samout:
                sam_file.close()
            if self.vsc_controller:
                viral_sam.close()
            p.communicate()

            self.nbases = self.get_nbases(readin.stderr.readlines())
            outf.write('#\tnbases\t{}\t#\t#\n'.format(self.nbases))
            outf.close()
            self.input_type = 'mapout'
            self.inp = self.mapout
        except OSError as e:
            error('OSError: "{}"\nFatal error running Minimap2.'.format(e), exit=True)
        except IOError as e:
            error('IOError: "{}"\nFatal error running Minimap2.'.format(e), exit=True)
        if p.returncode == 13:
            error("Permission Denied Error: fatal error running Minimap2. Is the Minimap2 file in the path with execution and read permissions?", exit=True)
        elif p.returncode != 0:
            error("Error while running Minimap2.\n", exit=True)
            
    def check_hq_mapping(self, sam_line):
        """Checks whether a hit in the SAM file is of high quality

        Args:
            sam_line (list): SAM file line as a list

        Returns:
            bool: Whether the hit was of HQ
        """
        if not sam_line[0].startswith('@'):
            if not sam_line[2].endswith('*'):
                if len(sam_line[9]) >= self.read_min_len: # check read length
                    if (hex(int(sam_line[1]) & 0x100) == '0x0'):  # no secondary
                        if self.mapq_filter(sam_line[2], int(sam_line[4]), self.min_mapq_val):  # filter low mapq reads
                            if ((self.min_alignment_len is None) or (max([int(x.strip('M')) for x in re.findall(r'(\d*M)', sam_line[5]) if x]) >= self.min_alignment_len)): # check alignment length
                                if self.get_gscd(sam_line) <= self.max_gcsd: # check GCSD
                                    return True
        return False

    def mapq_filter(self, marker_name, mapq_value, min_mapq_val):
        """Checks whether a mapping hit pass a mapq quality filter

        Args:
            marker_name (str): the marker name of the hit
            mapq_value (int): the mapq value of the hit
            min_mapq_val (int): the mapq threshold

        Returns:
            bool: whether the hit passed the mapq filter
        """
        if 'GeneID:' in marker_name  or 'VDB' in marker_name:
            return True
        else:
            if mapq_value > min_mapq_val:
                return True
        return False
    
    def discard_multiple_clades(self, reads2markers, reads2gcsd, reads2nbases, m2c):
        """Assigns reads to clades based on markers by majority vote, GCSD is used as tie-breaker.
        Only markers of the same clade are kept for each read.

        Args:
            reads2markers (dict): dictionary containing the reads to markers
            reads2gcsd (dict): dictionary containing the reads to GCSD values
            reads2nbases (dict): dictionary containing the reads to number of bases
            m2c (dict): dictionary containing the markers to clades information

        Returns:
            dict, dict: the reads to markers and reads to number of bases dictionaries, after removing markers from different clades
        """

        for r, m in reads2markers.items():
            clades = []
            for marker in m:
                ## if marker in mpa_pkl append the clade, otherwise append "viral"
                if marker in m2c.keys():
                    clades.append(m2c[marker]['clade'])
                else: clades.append("viral")
                
            top = pd.Series(clades).value_counts()
            if len(top)>1 and top.iloc[0] == top.iloc[1]:
                ## if there is a tie, use GCSD as tie-breaker
                top = top[top == top.max()].index.tolist()
                top = min([(reads2gcsd[r][i],clades[i]) for i in range(len(clades)) if clades[i] in top], key=lambda x: x[0])[1]
                # top = clades[reads2gcsd[r].index(min(reads2gcsd[r]))]
            else:
                ## keep the most represented clade
                top = top.index[0]
            
            if top != "viral": ## Do not keep reads assigned to viral clusters
                reads2markers[r] = [m[i] for i,c in enumerate(clades) if c == top]
                reads2nbases[r] = [reads2nbases[r][i] for i,c in enumerate(clades) if c == top]

        return reads2markers, reads2nbases
    
    def get_nbases(self, read_fastx_stderr):
        """Gets the total number of bases from the readfastx execution

        Args:
            read_fastx_stderr (list): standard error from the readfastx execution
        """
        try:
            _, _, nbases = list(map(float, read_fastx_stderr[0].decode().split()))
            if not nbases or int(nbases)==0:
                error('Fatal error running MetaPhlAn. Total metagenome size was not estimated or is zero.\nPlease check your input files.', exit=True)
            return int(nbases)
        except ValueError:
            os.unlink(self.mapout)
            error(b''.join(read_fastx_stderr).decode(), exit=True)
    
    def get_reads2markers(self):
        """Retrieves the reads to markers information from the mapping results

        Returns:
            dict: dictionary containing the reads to mapping results
        """
        if not self.inp:
            ras, ras_line, inpf = plain_read_and_split, plain_read_and_split_line, sys.stdin
        elif self.inp.endswith(".bz2"):
            ras, ras_line, inpf = read_and_split, read_and_split_line, bz2.BZ2File(self.inp, "r")
        else:
            ras, ras_line, inpf = plain_read_and_split, plain_read_and_split_line, open(self.inp)
        reads2markers = defdict(list) ## dict of lists to allow multiple markers per read (for long reads)
        reads2gcsd = defdict(list) ## LONG READS: store GCSD values for long reads
        reads2nbases = defdict(list) ## LONG READS: number of marker positions covered by the read
        reads2readlen = defdict(int) if self.mapping_subsampling else 1 ## will be needed for mapping subsampling

        if self.input_type == 'mapout':
            for r, m, b, g, l in ras(inpf):
                if r.startswith('#') and 'nbases' in m:
                    nbases = int(b)
                else:
                    reads2markers[r].append(m)
                    reads2nbases[r].append(int(b))
                    reads2gcsd[r].append(float(g))
                    if self.mapping_subsampling and int(l) > reads2readlen[r]:
                        reads2readlen[r] = int(l)
        elif self.input_type == 'sam':
            nbases = int(self.nbases)
            if self.vsc_controller:
                viral_sam = open(self.vsc_controller.vscSamFile, 'wb')
            for line in inpf:
                sam_line = ras_line(line)
                # check for VSC markers
                if self.vsc_controller and ((sam_line[2].startswith('VDB|')) or (sam_line[0].startswith('@') and 'VDB|' in sam_line[1])):
                    viral_sam.write(line)
                if self.check_hq_mapping(sam_line):
                    gcsd = self.get_gscd(sam_line)
                    reads2markers[sam_line[0]].append(sam_line[2].split('/')[0])
                    reads2gcsd[sam_line[0]].append(gcsd)
                    # Number of positions covered by the read (sum of the M and D of the CIGAR string)
                    reads2nbases[sam_line[0]].append(sum([int(x[:-1]) for x in re.findall(r"\d*[DM]", sam_line[5])]))
                    if self.mapping_subsampling and len(sam_line[9]) > reads2readlen[sam_line[0]]:
                        reads2readlen[sam_line[0]] = len(sam_line[9])
            if self.vsc_controller:
                viral_sam.close()
        inpf.close()
        reads2markers, reads2nbases = self.discard_multiple_clades(reads2markers, reads2gcsd, reads2nbases, self.database_controller.database_pkl['markers'])
        return nbases, reads2readlen, reads2markers, reads2nbases

    def run_mapping(self):
        """Runs all the steps of the Minimap2 mapping"""
        self.check_mm2_database()
        self.init_mapping_arguments()
        self.run_minimap2()

    def __init__(self, args, index, database_controller):
        super().__init__(args, index)
        self.minimap2_exe = args.minimap2_exe if args.minimap2_exe else 'minimap2'
        self.mm2_ps_list = args.minimap2_ps.strip().split(" ")
        self.mm2_ps_str = args.minimap2_ps.strip().replace(" ","").replace(",","_").replace("-","_")
        self.max_gcsd = args.max_gcsd
        self.nbases = args.nbases ## nbases instead of nreads for long reads
        self.database_controller = database_controller ## needed to discard multiple clades
        self.mapping_subsampling = args.mapping_subsampling ## needed to know if readlength is needed from sam

class MetaphlanAnalysis:
    def get_mapped_fraction(self):
        """Gets the estimated fraction of mapped reads (bases for long reads)""" 
        self.mapped = self.tree.relative_abundances()
        if self.unclassified_estimation:
            self.fraction_mapped = min(self.mapped/float(self.total_metagenome), 1.0)
        else:
            self.fraction_mapped = 1.0
            
    def write_common_headers(self, outf):
        """Writes the common headers of the MetaPhlAn output

        Args:
            outf (stream): output stream
        """
        outf.write('#{}\n'.format(self.index))
        outf.write('#{}\n'.format(' '.join(sys.argv)))
        outf.write('#{} {} processed\n'.format(self.total_metagenome, 'bases' if self.avg_read_length == 1 else 'reads'))   
        outf.write('#{}\n'.format('\t'.join([self.sample_id_key, self.sample_id])))
    
    def report_results(self, tree, total_metagenome, avg_read_length):
        """Reports the MetaPhlAn results
        
        Args:
            tree (Tree): the MetaPhlAn tree
            total_metagenome (int): the total metagenome size
            avg_read_length (int): the average read length
        """
        self.tree = tree
        self.total_metagenome = total_metagenome
        self.avg_read_length = avg_read_length
        
    def __init__(self, args, database_controller, index):
        ## n_metagenome_reads substituted with total_metagenome (n of reads for short reads and n of bases for long reads)
        self.output = args.output_file
        self.sample_id_key = args.sample_id_key
        self.sample_id = args.sample_id
        self.index = index
        self.database_controller = database_controller
        self.total_metagenome = None
        self.fraction_mapped = None
        self.mapped = None
        self.avg_read_length = None        
        self.tree = None
        

class RelativeAbundanceAnalysis(MetaphlanAnalysis):
    def report_cami_output(self):
        """Reports the MetaPhlAn results in the CAMI output format"""
        ranks2code = { 'k' : 'superkingdom', 'p' : 'phylum', 'c':'class', 'o' : 'order', 'f' : 'family', 'g' : 'genus', 's' : 'species'}
        out_stream = open(self.output,"w") if self.output else sys.stdout
        with out_stream as outf:
            outf.write('''@SampleID:{}\n@Version:0.10.0\n@Ranks:superkingdom|phylum|class|order|family|genus|species|strain\n@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n'''.format(self.sample_id))                 
            clade2abundance = self.get_clade2abundance()
            for clade, values in clade2abundance.items():
                taxid, relab = values
                if taxid and clade.split('|')[-1][0] != 't': 
                    rank = ranks2code[clade.split('|')[-1][0]]       
                    leaf_taxid = taxid.split('|')[-1]
                    taxpathsh = '|'.join([re.sub(r'^[a-z]__', '', name) for name in clade.split('|')])
                    outf.write( '\t'.join( [ leaf_taxid, rank, taxid, taxpathsh, str(relab*self.fraction_mapped) ] ) + '\n' )
    


    def to_biomformat(self, clade_name):
        """Converts the clade name to a BIOM format 

        Args:   
            clade_name (str): the clade name
        """
        return {'taxonomy': clade_name.split(self.biom_mdelim)}
    
    def report_biom_output(self):
        """Reports the MetaPhlAn results in the BIOM output format"""
        json_key = "MetaPhlAn"
        clade2abundance = self.get_clade2abundance()
        out_stream = open(self.output,"w") if self.output else sys.stdout
        if len(clade2abundance) == 0:
            biom_table = biom.Table([], [], [])
            biom_table.to_json(json_key, direct_io=out_stream)
        else:
            packed=list()
            for clade, values in clade2abundance.items():
                taxid, relab = values
                if clade.split(self.biom_mdelim)[-1].startswith('s__'): 
                    packed.append([[relab], clade, taxid])
            data, clade_names, _ = zip(*packed)
            data = np.array(data)
            sample_ids = [self.sample_id]
            table_id = 'MetaPhlAn_Analysis'
            if version.parse(biom.__version__) < version.parse('2.0.0'):
                biom_table = biom.table.table_factory(data, sample_ids, clade_names, sample_metadata=None, 
                                                      observation_metadata=list(map(self.to_biomformat, clade_names)),
                                                      table_id=table_id, constructor= biom.table.DenseOTUTable)
                json.dump(biom_table.getBiomFormatObject(json_key), out_stream)

            else:  
                biom_table = biom.table.Table(data, clade_names, sample_ids, sample_metadata=None, 
                                              observation_metadata=list(map(self.to_biomformat, clade_names)),
                                              table_id=table_id, input_is_dense=True)
                biom_table.to_json(json_key, direct_io=out_stream)


    def get_clade2abundance(self):
        """Gets the filtered clade to relative abundance dictionary

        Returns:
            dict: The clade to relative abundance dictionary
        """
        clade2abundance = {}
        for clade in self.tree.all_clades.values():
            if clade.coverage > 0 and (self.tax_lev == 'a' or clade.name.startswith(self.tax_lev + '__')):
                clade2abundance[clade.get_full_name()] = (clade.get_full_taxids(), clade.rel_abundance)
        if len(clade2abundance) > 0:
            clade2abundance = dict(sorted(clade2abundance.items(), reverse=True, key=lambda x:x[1][1]+(100.0*(8-(x[0].count("|"))))))
        return clade2abundance
    
    def report_results(self, tree, total_metagenome, avg_read_length):
        """Reports the MetaPhlAn results"""
        super().report_results(tree, total_metagenome, avg_read_length)
        self.get_mapped_fraction()
        if self.cami_output:
            self.report_cami_output()
        elif self.biom_output:
            self.report_biom_output()
        else:
            out_stream = open(self.output,"w") if self.output else sys.stdout
            with out_stream as outf:
                add_repr = None
                self.write_common_headers(outf)
                if not self.use_group_representative:
                    outf.write('#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species\n')
                else:
                    outf.write('#clade_name\tNCBI_tax_id\trelative_abundance\n')
                if self.unclassified_estimation:
                    outf.write( "\t".join(["UNCLASSIFIED", "-1", str(round((1-self.fraction_mapped)*100,5)),""]) + "\n" )                   
                clade2abundance = self.get_clade2abundance()
                for clade, values in clade2abundance.items():
                    taxid, relab = values
                    if (clade, taxid) in self.database_controller.database_pkl['merged_taxon'] and not self.use_group_representative:
                        add_repr = '{}'.format(','.join( [ n[0] for n in self.database_controller.database_pkl['merged_taxon'][(clade, taxid)]] ))
                        outf.write( "\t".join( [clade,  taxid, str(relab*self.fraction_mapped), add_repr ] ) + "\n" )
                    else:                        
                        outf.write( "\t".join( [clade,  taxid, str(relab*self.fraction_mapped)] ) + "\n" )                            
            if add_repr is not None:
                warning("The metagenome profile contains clades that represent multiple species merged into a single representant. "
                        "An additional column listing the merged species is added to the MetaPhlAn output.")
            
    def __init__(self, args, database_controller, index):
        super().__init__(args, database_controller, index)
        self.tax_lev = args.tax_lev
        self.cami_output = args.CAMI_format_output
        self.biom_output = args.biom_format_output
        self.biom_mdelim = args.biom_mdelim
        self.use_group_representative = args.use_group_representative
        self.unclassified_estimation = not args.skip_unclassified_estimation
        
        
class RelativeAbundanceReadStatsAnalysis(MetaphlanAnalysis):
    def get_clade2abundance(self):
        """Gets the filtered clade to relative abundance dictionary

        Returns:
            dict: The clade to relative abundance dictionary
        """
        clade2abundance = {}
        for clade in self.tree.all_clades.values():
            if clade.coverage > 0 and (self.tax_lev == 'a' or clade.name.startswith(self.tax_lev + '__')):
                clade2abundance[clade.get_full_name()] = (clade.get_full_taxids(), clade.rel_abundance, clade.coverage, clade.nreads)
        if len(clade2abundance) > 0:
            clade2abundance = dict(sorted(clade2abundance.items(), reverse=True, key=lambda x:x[1][1]+(100.0*(8-(x[0].count("|"))))))
        return clade2abundance
    
    def report_results(self, tree, total_metagenome, avg_read_length):
        """Reports the MetaPhlAn results"""
        super().report_results(tree, total_metagenome, avg_read_length)
        self.get_mapped_fraction()
        unit = 'bases' if self.avg_read_length == 1 else 'reads'
        out_stream = open(self.output,"w") if self.output else sys.stdout
        with out_stream as outf:
            self.write_common_headers(outf)
            outf.write( "#Estimated {} mapped to known clades: {}\n".format(unit, int(self.mapped)))
            outf.write( "\t".join( ["#clade_name","clade_taxid","relative_abundance","coverage","estimated_number_of_{}_from_the_clade".format(unit) ]) +"\n" )
            if self.unclassified_estimation:
                outf.write( "\t".join(["UNCLASSIFIED", "-1", str(round((1-self.fraction_mapped)*100,5)),"-", str(self.total_metagenome - self.mapped)]) + "\n" )                   
            clade2abundance = self.get_clade2abundance()
            for clade, values in clade2abundance.items():
                taxid, relab, coverage, nreads = values
                outf.write( "\t".join( [clade,  taxid, str(relab*self.fraction_mapped), str(coverage), str(nreads)] ) + "\n" )         
            
    def __init__(self, args, database_controller, index):
        super().__init__(args, database_controller, index)
        self.tax_lev = args.tax_lev
        self.unclassified_estimation = args.unclassified_estimation
    
class CladeProfilesAnalysis(MetaphlanAnalysis):
    def report_results(self, tree, total_metagenome, avg_read_length):
        """Reports the MetaPhlAn results"""
        super().report_results(tree, total_metagenome, avg_read_length)
        clade2profiles = self.tree.clade_profiles()
        out_stream = open(self.output,"w") if self.output else sys.stdout
        with out_stream as outf:
            self.write_common_headers(outf)
            for clade, profile in clade2profiles.items():
                mn,n = zip(*profile)
                outf.write( "\t".join( [""]+[str(s) for s in mn] ) + "\n" )
                outf.write( "\t".join( [clade]+[str(s) for s in n] ) + "\n" )
            
    def __init__(self, args, database_controller, index):
        super().__init__(args, database_controller, index)
    
class MarkerAbundanceTableAnalysis(MetaphlanAnalysis):
    def report_results(self, tree, total_metagenome, avg_read_length):
        """Reports the MetaPhlAn results"""
        super().report_results(tree, total_metagenome, avg_read_length)
        clade2profiles = self.tree.clade_profiles()
        out_stream = open(self.output,"w") if self.output else sys.stdout
        with out_stream as outf:
            self.write_common_headers(outf)     
            for v in clade2profiles.values():
                outf.write( "\n".join(["\t".join([str(a),str(b/float(self.total_metagenome)) if self.total_metagenome else str(b)]) for a,b in v if b > 0.0]) + "\n" )
            
    def __init__(self, args, database_controller, index):
        super().__init__(args, database_controller, index)
    
class MarkerPresenceTableAnalysis(MetaphlanAnalysis):
    def report_results(self, tree, total_metagenome, avg_read_length):
        """Reports the MetaPhlAn results"""
        super().report_results(tree, total_metagenome, avg_read_length)
        clade2profiles = self.tree.clade_profiles()
        out_stream = open(self.output,"w") if self.output else sys.stdout
        with out_stream as outf:
            self.write_common_headers(outf)       
            for v in clade2profiles.values():
                strout = ["\t".join([str(a),"1"]) for a,b in v if b > self.pres_th]
                if strout:
                    outf.write( "\n".join(strout) + "\n" )
            
    def __init__(self, args, database_controller, index):
        super().__init__(args, database_controller, index)
        self.pres_th = args.pres_th
        

class Metaphlan:
    """MetaPhlAn class"""

    def get_markers_to_ignore(self, ignore_markers):
        """Retrieves the markers to ignore from a file

        Args:
            ignore_markers (str): path to the file with the markers to ignore

        Returns:
            set: the set of markers to ignore
        """
        if ignore_markers:
            with open(ignore_markers) as ignv:
                return set([l.strip() for l in ignv])
        else:
            return set()


    def build_taxonomy_tree(self):
        """Build the MetaPhlAn taxonomy tree"""
        self.tree = TaxTree(self.database_controller.database_pkl, self.ignore_markers,
                       self.stat, self.stat_q, self.perc_nonzero, self.avoid_disqm, self.avg_read_length) #, self.long_reads)
            
    def separate_reads2markers(self, reads2markers):
        """Separates the viral hits from the reads to markers dictionary

        Args:
            reads2markers (dict): full reads to markers dictionary

        Returns:
            (dict, dict): tuple of dictionaries for the SGBs/EUKs and viral contigs
        """
        ## changed to work with dict
        return {r: m for r, m in reads2markers.items() if ('SGB' in m[0] or 'EUK' in m[0]) and not 'VDB' in m[0]}, {r: m for r, m in reads2markers.items() if 'VDB' in m[0] and not ('SGB' in m[0] or 'EUK' in m[0])}

    def make_gen_fastq(self, reader):
        """Read fastq file in chunks
        
        Args:
            reader (function): the reader function
        """
        b = reader(1024 * 1024) 
        while (b):
            yield b
            b = reader(1024 * 1024)

    def rawpycount(self, filename, intype):
        """Counts lines in a fastq file
        
        Args:
            filename (str): the file name
            intype (str): the input type
        
        Returns:
            int: the number of bases
        """
        f = bz2.BZ2File(filename, 'rb') if filename.endswith(".bz2") else gzip.open(filename, 'rb') if filename.endswith('.gz') else open(filename, 'rb')
        f_gen = self.make_gen_fastq(f.read)
        total_metagenome = sum( buf.count(b'\n') for buf in f_gen) if intype == 'fastq' else sum( buf.count(b'>') for buf in f_gen)
        f.close()
        return int(total_metagenome/4) if intype == 'fastq' else int(total_metagenome)
    
    def rawpycount_bases(self, filename, intype, tot=True):
        """Counts bases in a fastq file
        
        Args:
            filename (str): the filename
            intype (str): the input type
            tot (bool): whether to count the total number of bases or the number of bases per read
        
        Returns:
            int or dict: the total number of bases or the number of bases per read"""
        f = bz2.BZ2File(filename, 'rb') if filename.endswith(".bz2") else gzip.open(filename, 'rb') if filename.endswith('.gz') else open(filename, 'rb')
        f_gen = self.make_gen_fastq(f.read)
        
        # Determine which lines to count based on the input type
        i = 4 if intype == 'fastq' else 2
        j = 2 if intype == 'fastq' else 0

        nbases = 0 if tot else dict()
        read_num = -1
        line_num = 0
        buffer = b''
        
        for chunk in f_gen:
            lines = (buffer + chunk).split(b'\n')
            buffer = lines[-1]  # Save the last incomplete line to be processed in the next chunk   
            for line in lines[:-1]:  # Process all but the last line
                line_num += 1
                if line_num % i == j:
                    if tot:
                        nbases += len(line)
                    else:
                        read_num += 1
                        nbases[read_num] = len(line)
        
        # Check if there's remaining content in the buffer
        if buffer:
            line_num += 1
            if line_num % i == j:
                    if tot:
                        nbases += len(line)
                    else:
                        read_num += 1
                        nbases[read_num] = len(line)
        f.close()
        return nbases
    
    def subsample_file(self, sample, out):
        """Subsamples the reads of the file

        Args:
            sample (set): the set of reads to sample
            out (file): the output file
        """
        length, read_num = 0, -1
        counter=0
        out_l = out
        x = 3 if self.input_type == 'fastq' else 1
        for inp_f in self.inp.split(','):
            if self.subsampling_paired:
                length, read_num = 0, -1
                out=out_l[counter]
                counter+=1
            line = 1
            reader = bz2.open(inp_f, 'rt') if inp_f.endswith(".bz2") else gzip.open(inp_f, 'rt') if inp_f.endswith('.gz') else open(inp_f, 'r')
            while (line and length < self.subsampling):
                line = reader.readline()
                read_num += 1
                if read_num in sample:
                    length += 1
                    out.write(line)       
                    for _ in range(x):
                        out.write(reader.readline())
                else:
                    for _ in range(x):
                        reader.readline()
            reader.close()
        out.close()

    def prepare_subsample_output(self):
        """Prepares the output file for subsampling
        
        Returns:
            file: the output file
        """
        if self.subsampling_output is None:
            if self.subsampling_paired:
                self.subsampling_output = list()
                out = list()
                for _ in self.inp.split(','):
                    out_f = tempfile.NamedTemporaryFile(dir=self.tmp_dir, mode='w', delete=False)
                    out.append(out_f)
                    self.subsampling_output.append(out_f.name)
            else:
                out= tempfile.NamedTemporaryFile(dir=self.tmp_dir, mode='w', delete=False)
                self.subsampling_output = out.name
        else:
            if self.subsampling_paired:
                r, ext = os.path.splitext(self.subsampling_output)

                self.subsampling_output = list()
                out = list()

                self.subsampling_output.append('.'.join([r, 'R1'+ext]))
                self.subsampling_output.append('.'.join([r, 'R2'+ext]))

                for s in self.subsampling_output:
                    out.append(bz2.open(s, 'wt') if s.endswith(".bz2") else gzip.open(s, 'wt') if s.endswith('.gz') else open(s, 'w'))
            else:
                out = bz2.open(self.subsampling_output, 'wt') if self.subsampling_output.endswith(".bz2") else gzip.open(self.subsampling_output, 'wt') if self.subsampling_output.endswith('.gz') else open(self.subsampling_output, 'w')
        return out


    def subsample_reads(self):
        """Subsamples the reads of the input file"""
        self.total_metagenome = execute_pool(((self.rawpycount, inp_f, self.input_type) for inp_f in self.inp.split(',')), 2)

        if self.subsampling >= sum(self.total_metagenome):
            warning("The specified subsampling ({}) is equal or higher than the original number of reads ({}). Subsampling will be skipped.\n".format(self.subsampling, self.total_metagenome)) 
            self.total_metagenome = sum(self.total_metagenome)
            self.subsampling_output = self.inp
            return

        if self.subsampling_paired:
            self.subsampling //= 2
            if self.total_metagenome[0] != self.total_metagenome[1]:
                error("The specified reads file are not the same length! Make sure the forward and reverse reads are files are not damaged and reads are in the same order.", init_new_line = True, exit = True) 
            self.total_metagenome = self.total_metagenome[0]
        else:
            self.total_metagenome = sum(self.total_metagenome)


        if self.subsampling_seed.lower() != 'random':
            random.seed(int(self.subsampling_seed))

        sample = set(random.sample(range(self.total_metagenome), self.subsampling))

        out = self.prepare_subsample_output()
        self.subsample_file(sample, out)

        if isinstance(self.subsampling_output, list):
            self.subsampling_output = ','.join(self.subsampling_output)
            self.subsampling = self.subsampling*2

        # update number of reads and input file
        self.total_metagenome = self.subsampling
        self.inp = self.subsampling_output


    def subsample_bases(self):
        """Subsamples the bases of the input file"""
        r2rl = self.rawpycount_bases(self.inp, self.input_type, tot=False)
        self.total_metagenome = sum(r2rl.values())

        if self.subsampling >= self.total_metagenome:
            warning("The specified subsampling ({}) is equal or higher than the original number of bases ({}). Subsampling will be skipped.\n".format(self.subsampling, self.total_metagenome)) 
            self.subsampling_output = self.inp
            return

        if self.subsampling_seed.lower() != 'random':
            random.seed(int(self.subsampling_seed))

        cur_bases = 0
        r = list(r2rl.keys())
        random.shuffle(r)
        i = 0
        while cur_bases < self.subsampling:
            cur_bases += r2rl[r[i]]
            i += 1

        sample = set(r[:i])
        self.subsampling = len(sample)

        out = self.prepare_subsample_output()
        self.subsample_file(sample, out)

        # update number of reads and input file
        self.total_metagenome = cur_bases
        self.inp = self.subsampling_output
        

    def mapping_subsample_reads(self, reads2markers):
        """Subsamples the reads using the mapping results 

        Args:
            reads2markers (dict): full reads to markers dictionary

        Returns:
            dict: the subsampled reads to markers dictionary
        """
        if self.subsampling >= self.total_metagenome:
            warning("WARNING: The specified subsampling ({}) is equal or higher than the original number of reads ({}). Subsampling will be skipped.".format(self.subsampling, self.total_metagenome), init_new_line=True)
        else:
            reads2markers =  dict(sorted(reads2markers.items()))
            if self.subsampling_seed.lower() != 'random':
                random.seed(int(self.subsampling_seed))
            reads2filtmarkers = {}
            sgb_reads2markers, viral_reads2markers = self.separate_reads2markers(reads2markers)            
            n_sgb_mapped_reads = int((len(sgb_reads2markers) * self.subsampling) / self.total_metagenome)
            reads2filtmarkers = { r:sgb_reads2markers[r] for r in random.sample(list(sgb_reads2markers.keys()), n_sgb_mapped_reads) }         
            n_viral_mapped_reads = int((len(viral_reads2markers) * self.subsampling) / self.total_metagenome)
            reads2filtmarkers.update({ r:viral_reads2markers[r] for r in random.sample(list(viral_reads2markers.keys()), n_viral_mapped_reads) })            
            reads2markers = reads2filtmarkers
            sgb_reads2markers.clear()
            viral_reads2markers.clear()
            self.total_metagenome = self.subsampling

        return reads2markers
    
    def mapping_subsample_longreads(self, reads2readlen,reads2markers,reads2nbases):
        """Subsamples the reads using the mapping results 

        Args:
            reads2readlen (dict): dictionary containing the reads to read length information
            reads2markers (dict): full reads to markers dictionary
            reads2nbases (dict): full reads to number of bases covered per marker

        Returns:
            dict: the subsampled reads to markers dictionary
        """
        if self.subsampling >= self.total_metagenome:
            warning("WARNING: The specified subsampling ({}) is equal or higher than the original number of reads ({}). Subsampling will be skipped.".format(self.subsampling, self.total_metagenome), init_new_line=True)
        else:
            reads2markers =  dict(sorted(reads2markers.items()))
            if self.subsampling_seed.lower() != 'random':
                random.seed(int(self.subsampling_seed))
            reads2keep = []
            sgb_reads2markers, viral_reads2markers = self.separate_reads2markers(reads2markers)

            n_sgb_mapped_bases = int((sum([reads2readlen[r] for r in  sgb_reads2markers.keys()]) * self.subsampling) / self.total_metagenome)
            n_viral_mapped_bases = int((sum([reads2readlen[r] for r in viral_reads2markers.keys()]) * self.subsampling) / self.total_metagenome)
            cur_bases = 0
            while cur_bases < n_sgb_mapped_bases:
                r = random.choice(list(sgb_reads2markers.keys()))
                if r not in reads2keep:
                    reads2keep.append(r)
                    cur_bases += reads2readlen[r]
            while cur_bases < (n_sgb_mapped_bases + n_viral_mapped_bases):
                r = random.choice(list(viral_reads2markers.keys()))
                if r not in reads2keep:
                    reads2keep.append(r)
                    cur_bases += reads2readlen[r]

            reads2markers = {r:reads2markers[r] for r in reads2keep}
            reads2nbases = {r:reads2nbases[r] for r in reads2keep}
            sgb_reads2markers.clear()
            viral_reads2markers.clear()
            self.total_metagenome = cur_bases

        return reads2markers, reads2nbases

    def parse_mapping(self):
        """Parses the mapping results into a marker to reads dictionary

        Returns:
            dict: marker to reads dictionary
        """
        self.total_metagenome, self.avg_read_length, reads2markers, reads2nbases = self.mapping_controller.get_reads2markers()
        self.build_taxonomy_tree()     
        if not self.long_reads and self.subsampling is not None and self.mapping_subsampling:
            reads2markers = self.mapping_subsample_reads(reads2markers)
        elif self.long_reads and self.subsampling is not None and self.mapping_subsampling:
            reads2markers, reads2nbases = self.mapping_subsample_longreads(self.avg_read_length, reads2markers, reads2nbases)
            self.avg_read_length = 1
        elif not self.long_reads and self.subsampling is None and self.total_metagenome < 10000:
            warning("The number of reads in the sample ({}) is below the recommended minimum of 10,000 reads.".format(self.total_metagenome))
        elif self.long_reads and self.subsampling is None and self.total_metagenome < (10000*150):
            warning("The number of bases in the sample ({}) is below the recommended minimum of 1,500,000 bases.".format(self.total_metagenome))  
        markers2reads = defdict(list) ## list instead of set for long reads
        for r, m in reads2markers.items(): ## changed to work with dict of lists
            for i in range(len(m)):
                if not reads2nbases:
                    markers2reads[m[i]].append(r) ## short reads --> append read names
                else:
                    # markers2nbases[m[i]] += reads2nbases[r][i]
                    markers2reads[m[i]].append(reads2nbases[r][i]) ## long reads --> append number of bases covered by the read
        self.add_reads_to_tree(markers2reads)

        if self.no_map:
            os.remove(self.mapping_controller.inp)

    def add_reads_to_tree(self, markers2reads):
        """Adds reads mapping to the markers to the tree structure

        Args:
            markers2reads (dict): dictionary with the markers to reads information
        """
        for marker,basesxread in sorted(markers2reads.items(), key=lambda pars: pars[0]):
            if marker not in self.tree.markers2lens:
                continue
            self.tree.add_reads( marker, len(basesxread) if not self.long_reads else sum(basesxread), ## add number of reads if short reads , number of bases if long reads
                ignore_eukaryotes = self.ignore_eukaryotes,
                ignore_bacteria = self.ignore_bacteria,
                ignore_archaea = self.ignore_archaea,
                ignore_ksgbs = self.ignore_ksgbs,
                ignore_usgbs = self.ignore_usgbs
            )
                
    def init_metaphlan_analysis(self, args):
        """Initializes the MetaPhlAn analysis to the specified analysis type

        Args:
            args (namespace): the arguments of the metaphlan execution

        Returns:
            MetaphlanAnalysis: the selected MetaPhlAn analysis object
        """
        if args.t == 'rel_ab':
            return RelativeAbundanceAnalysis(args, self.database_controller, self.index)
        elif args.t == 'rel_ab_w_read_stats':
            return RelativeAbundanceReadStatsAnalysis(args, self.database_controller, self.index)
        elif args.t == 'clade_profiles':        
            return CladeProfilesAnalysis(args, self.database_controller, self.index)
        elif args.t == 'marker_ab_table':   
            return MarkerAbundanceTableAnalysis(args, self.database_controller, self.index)
        elif args.t == 'marker_pres_table':  
            return MarkerPresenceTableAnalysis(args, self.database_controller, self.index)
        else:
            error('The specified analysis type is not available', exit=True)


    def run_metaphlan(self):
        """Runs the MetaPhlAn pipeline"""    
        self.index=self.database_controller.check_and_install_database()
        if self.verbose or self.install:
            info('The database is installed ({})'.format(self.index), stderr=True, exit=self.install)
        if (self.subsampling_paired or self.subsampling) and not self.mapping_subsampling:
            if not (self.long_reads or self.split_reads):
                self.subsample_reads()
            else:
                self.subsample_bases()
            self.mapping_controller.set_inp(self.inp)
            if self.long_reads:
                self.mapping_controller.set_nbases(self.total_metagenome)

        if self.input_type in ['fastq', 'fasta']:
            ## This is not needed anymore as the file is now read by read_fastx, that counts the number of bases
            # if self.long_reads and not self.total_metagenome:
                ## calculate nbases for long reads if needed because minimap skips read_fastx
                # self.mapping_controller.set_nbases(self.rawpycount_bases(self.inp, self.input_type))
            self.mapping_controller.run_mapping()       
        self.parse_mapping()
        self.metaphlan_analysis.report_results(self.tree, self.total_metagenome, self.avg_read_length)
        if self.profile_vsc:
            self.vsc_controller.run_analysis(self.total_metagenome)

    def __init__(self, args):
        self.inp = args.inp
        self.verbose = args.verbose
        self.long_reads = args.long_reads
        self.split_reads = args.split_reads
        self.database_controller = MetaphlanDatabaseController(args)
        self.index = self.database_controller.resolve_index()
        self.mapping_controller = Bowtie2Controller(args, self.index) if not self.long_reads else Minimap2Controller(args, self.index, self.database_controller)
        self.metaphlan_analysis = self.init_metaphlan_analysis(args)
        self.input_type = args.input_type
        self.ignore_eukaryotes = args.ignore_eukaryotes
        self.ignore_bacteria = args.ignore_bacteria
        self.ignore_archaea = args.ignore_archaea
        self.ignore_ksgbs = args.ignore_ksgbs
        self.ignore_usgbs = args.ignore_usgbs
        self.perc_nonzero = args.perc_nonzero
        self.no_map = args.no_map
        self.stat = args.stat
        self.stat_q = args.stat_q
        self.ignore_markers = self.get_markers_to_ignore(args.ignore_markers)
        self.avoid_disqm = args.avoid_disqm
        self.subsampling = args.subsampling
        self.subsampling_paired = args.subsampling_paired
        self.subsampling_output = args.subsampling_output
        self.mapping_subsampling = args.mapping_subsampling
        self.subsampling_seed = args.subsampling_seed
        self.install = args.install
        self.offline = args.offline
        self.force_download = args.force_download
        self.tree = None
        self.total_metagenome = None
        self.avg_read_length = None
        self.profile_vsc = args.profile_vsc
        self.tmp_dir = args.tmp_dir
        if args.profile_vsc:
            self.vsc_controller = VSC_bt2_controller(args, self.mapping_controller, self.database_controller) if not self.long_reads else VSC_mm2_controller(args, self.mapping_controller, self.database_controller)
            self.mapping_controller.set_vsc_controller(self.vsc_controller)


def read_params(args):
    """ Reads and parses the command line arguments of the script

    Returns:
        namespace: The populated namespace with the command line arguments
    """
    p = ap.ArgumentParser(description=" MetaPhlAn version "+__version__+" ("+__date__+"): \n"
                          " METAgenomic PHyLogenetic ANalysis for metagenomic taxonomic profiling.\n\n"
                          " AUTHORS: "+__author__+"\n\n", formatter_class=ap.RawTextHelpFormatter, add_help=False)
    arg = p.add_argument
    arg('inp', metavar='INPUT_FILE', type=str, nargs='?', default=None, help="the input file can be:\n"
        "* a fastq file containing metagenomic reads\n"
        "OR\n"
        "* a BowTie2 produced SAM file. \n"
        "OR\n"
        "* a minimap2 produced SAM file (for long reads). \n"
        "OR\n"
        "* an intermediary mapping file of the metagenome generated by a previous MetaPhlAn run (mapout)\n"
        "If the input file is missing, the script assumes that the input is provided using the standard \n"
        "input, or named pipes.\n"
        "IMPORTANT: the type of input needs to be specified with --input_type")
    g = p.add_argument_group('Required arguments')
    arg = g.add_argument
    input_type_choices = ['fastq', 'fasta', 'mapout', 'sam']
    arg('--input_type', choices=input_type_choices, required = '--install' not in args, help="set whether the input is the FASTA file of metagenomic reads or \n"
        "the SAM file of the mapping of the reads against the MetaPhlAn db.\n")
    g = p.add_argument_group('Mapping arguments')
    arg = g.add_argument
    arg('--force', action='store_true',
        help="Force profiling of the input file by removing the mapout file")
    arg('--db_dir', metavar="METAPHLAN_DB", type=str, default=DEFAULT_DB_FOLDER,
        help=("Folder containing the MetaPhlAn database. You can specify the location by exporting the DEFAULT_DB_FOLDER variable in the shell (old --bowtie2db option)."
              "[default "+DEFAULT_DB_FOLDER+"]\n"))
    arg('-x', '--index', type=str, default='latest',
        help=("Specify the id of the database version to use. "
              "If \"latest\", MetaPhlAn will get the latest version.\n"
              "If an index name is provided, MetaPhlAn will try to use it, if available, and skip the online check.\n"
              "If the database files are not found on the local MetaPhlAn installation they\n"
              "will be automatically downloaded [default latest]\n"))
    arg('--mapout', metavar="FILE_NAME", type=str, default=None,
        help="The file for saving the mapping output (old --bowtie2out option)")
    arg('--min_mapq_val', type=int, default=None,
        help="Minimum mapping quality value (MAPQ) [default 5 for short reads, 50 for long reads]")
    arg('--no_map', action='store_true',
        help="Avoid storing the --mapout map file")
    arg('--tmp_dir', metavar="", default=None, type=str,
        help="The folder used to store temporary files [default is the OS "
             "dependent tmp dir]")
    g = p.add_argument_group('Bowtie2 arguments')
    arg = g.add_argument
    bt2ps = ['sensitive', 'very-sensitive',
             'sensitive-local', 'very-sensitive-local']
    arg('--bt2_ps', metavar="BowTie2 presets", default='very-sensitive',
        choices=bt2ps, help="Presets options for BowTie2 (applied only when a "
                            "FASTA file is provided)\n"
                            "The choices enabled in MetaPhlAn are:\n"
                            " * sensitive\n"
                            " * very-sensitive\n"
                            " * sensitive-local\n"
                            " * very-sensitive-local\n"
                            "[default very-sensitive]\n")
    arg('--bowtie2_exe', type=str, default=None,
        help='Full path and name of the BowTie2 executable. This option allows'
             'MetaPhlAn to reach the executable even when it is not in the '
             'system PATH or the system PATH is unreachable')
    arg('--bowtie2_build', type=str, default='bowtie2-build',
        help="Full path to the bowtie2-build command to use, deafult assumes "
             "that 'bowtie2-build is present in the system path")
    g = p.add_argument_group('Post-mapping arguments')
    arg = g.add_argument
    stat_choices = ['avg_g', 'avg_l', 'tavg_g',
                    'tavg_l', 'wavg_g', 'wavg_l', 'med',
                    'npos_lr', 'nreads_lr']
    arg('--tax_lev', metavar='TAXONOMIC_LEVEL', type=str,
        choices='akpcofgst', default='a', help="The taxonomic level for the relative abundance output:\n"
        "'a' : all taxonomic levels\n"
        "'k' : kingdoms\n"
        "'p' : phyla only\n"
        "'c' : classes only\n"
        "'o' : orders only\n"
        "'f' : families only\n"
        "'g' : genera only\n"
        "'s' : species only\n"
        "'t' : SGBs only\n"
        "[default 'a']")
    arg('--min_alignment_len', metavar="", default=None, type=int, help="The sam records for aligned reads with the longest subalignment\n"
        "length smaller than this threshold will be discarded.\n"
        "[default None]\n")
    arg('--ignore_eukaryotes', action='store_true',
        help="Do not profile eukaryotic organisms")
    arg('--ignore_bacteria', action='store_true',
        help="Do not profile bacterial organisms")
    arg('--ignore_archaea', action='store_true',
        help="Do not profile archeal organisms")
    arg('--ignore_ksgbs', action='store_true',
        help="Do not profile known SGBs (together with --sgb option)")
    arg('--ignore_usgbs', action='store_true',
        help="Do not profile unknown SGBs (together with --sgb option)")
    arg('--stat_q', metavar="", type=float, default=0.2, help="Quantile value for the robust average\n"
        "[default 0.2]")
    arg('--perc_nonzero', metavar="", type=float, default=0.33, help="Percentage of markers with a non zero relative abundance for misidentify a species\n"
        "[default 0.33]")
    arg('--ignore_markers', type=str, default=None,
        help="File containing a list of markers to ignore. \n")
    arg('--avoid_disqm', action="store_true", help="Deactivate the procedure of disambiguating the quasi-markers based on the \n"
        "marker abundance pattern found in the sample. It is generally recommended \n"
        "to keep the disambiguation procedure in order to minimize false positives\n")
    arg('--stat', metavar="", choices=stat_choices, default="tavg_g", type=str, help="Statistical approach for converting marker abundances into clade abundances\n"
        "'avg_g'  : clade global (i.e. normalizing all markers together) average\n"
        "'avg_l'  : average of length-normalized marker counts\n"
        "'tavg_g' : truncated clade global average at --stat_q quantile\n"
        "'tavg_l' : truncated average of length-normalized marker counts (at --stat_q)\n"
        "'wavg_g' : winsorized clade global average (at --stat_q)\n"
        "'wavg_l' : winsorized average of length-normalized marker counts (at --stat_q)\n"
        "'med'    : median of length-normalized marker counts\n"
        "[default tavg_g]")
    arg = p.add_argument
    g = p.add_argument_group('Additional analysis types and arguments')
    arg = g.add_argument
    analysis_types = ['rel_ab', 'rel_ab_w_read_stats', 'clade_profiles',
                      'marker_ab_table', 'marker_pres_table']
    arg('-t', metavar='ANALYSIS TYPE', type=str, choices=analysis_types,
        default='rel_ab', help="Type of analysis to perform: \n"
        " * rel_ab: profiling a metagenomes in terms of relative abundances\n"
        " * rel_ab_w_read_stats: profiling a metagenomes in terms of relative abundances and estimate the number of reads coming from each clade.\n"
        " * clade_profiles: normalized marker counts for clades with at least a non-null marker\n"
        " * marker_ab_table: normalized marker counts (only when > 0.0 and normalized by metagenome size if --nreads is specified)\n"
        " * marker_pres_table: list of markers present in the sample (threshold at 1.0 if not differently specified with --pres_th\n"
        "[default 'rel_ab']")
    arg('--nreads', metavar="NUMBER_OF_READS", type=int, default=None, help="The total number of reads in the original metagenome. \n"
        "It is mandatory when the --input_type is a SAM file.")
    arg('--pres_th', metavar="PRESENCE_THRESHOLD", type=int, default=1.0,
        help='Threshold for calling a marker present by the -t marker_pres_table option')
    g = p.add_argument_group('Output arguments')
    arg = g.add_argument
    arg('-o', '--output_file',  metavar="output file", type=str, default=None,
        help="The output file (if not specified stdout)\n")
    arg('--sample_id_key',  metavar="name", type=str, default="SampleID",
        help=("Specify the sample ID key for this analysis."
              " Defaults to 'SampleID'."))
    arg('--use_group_representative', action='store_true',
        help=("Use a species as representative for species groups."))
    arg('--sample_id',  metavar="value", type=str,
        default="Metaphlan_Analysis",
        help=("Specify the sample ID for this analysis."
              " Defaults to 'Metaphlan_Analysis'."))
    arg('-s', '--samout', metavar="sam_output_file",
        type=str, default=None, help="The sam output file\n")
    arg('--CAMI_format_output', action='store_true',
        help="Report the profiling using the CAMI output format\n")
    arg('--skip_unclassified_estimation', action='store_true',
        help="Do not scale relative abundances to the estimate unclassified taxa\n")
    arg('--biom_format_output',action='store_true',
        help="Report the profiling using the biom output format\n")
    arg('--biom_mdelim',  metavar="mdelim", type=str, default="|",
        help="Delimiter for metadata in the biom output format [default '|'] \n")
    g = p.add_argument_group('Viral Sequence Clusters Analisys')
    arg = g.add_argument
    arg("--profile_vsc", action="store_true",help="Add this parameter to profile Viruses with VSCs approach.")
    arg("--vsc_out", help="Path to the VSCs breadth-of-coverage output file", default="mp3_viruses.csv")
    arg("--vsc_breadth", help="Minimum Breadth of Coverage for a Viral Group to be reported.\n"
    "Default is 0.75 (at least 75 percent breadth to report), 0.5 for long reads", default=None,type=float)
    g = p.add_argument_group('Long reads arguments')
    arg = g.add_argument
    arg('--long_reads', action="store_true",help="Add this parameter to profile long reads.")
    arg('--max_gcsd', type=float, default=0.10,
        help="Specify the max gap-compressed sequence divergence (minimap2) allowed between marker and read, default value is 0.10")
    arg('--minimap2_exe', type=str, default=None,
        help='Full path and name of the Minimap2 executable. This option allows'
             'MetaPhlAn to reach the executable even when it is not in the '
             'system PATH or the system PATH is unreachable')    
    arg('--minimap2_ps', metavar="minimap2 presets", type=str, default="-x asm20 -B 3 -O 3,12",
        help="Presets options for Minimap2 (applied only when a FASTA/Q file is provided)\n"
        "Please refer to the minimap2 documentation and provide the options in a string such as: '-x asm20 -B 3 -O 3,12'\n"
        "[default '-x asm20 -B 3 -O 3,12']\n")
    arg('--nbases', metavar="NUMBER_OF_BASES", type=int, default = None, help =
         "The total number of bases in the original metagenome. It is mandatory when the --input_type is a SAM file from long reads mapping (--long_reads is selected)" )
    arg('--split_reads', action="store_true",help="Add this parameter to profile long reads by splitting into short reads.")
    arg('--split_readlen', type=int, default=150,
        help="Specify length of the reads to be split into when using --split_reads, default value is 150")
    g = p.add_argument_group('Other arguments')
    arg = g.add_argument
    arg('--nproc', metavar="N", type=int, default=4,
        help="The number of CPUs to use for parallelizing the mapping [default 4]")
    arg('--subsampling', type=int, default=None,
        help="Specify the number of reads to be considered from the input metagenomes (number of bases if --long_reads) [default None]")
    arg('--mapping_subsampling', action='store_true',
        help="If used, the subsamping will be done on the mapping results instead of on the reads.")
    arg('--subsampling_seed', type=str, default='1992',
        help="Random seed to use in the selection of the subsampled reads. Choose \"random\r for a random behaviour")
    arg('--subsampling_output', type=str, default=None,
        help="The output file for the subsampled reads. If not specified the subsampled reads will not be saved.")
    arg('--subsampling_paired',  type=int, default=None,
        help="Specify the number of paired reads to be considered from the input metagenomes [default None]")
    arg('-1', type=str, default=None, metavar='FORWARD_READS', dest='forward',
        help="Specify the fastq file with forward reads of the input metagenomes. Reads are assumed to be in the same order in the forward and reverse files! [default None]")
    arg('-2', type=str, default=None, metavar='REVERSE_READS', dest='reverse',
        help="Specify the fastq file with reverse reads of the input metagenomes. Reads are assumed to be in the same order in the forward and reverse files! [default None]")
    arg('--install', action='store_true',
        help="Only checks if the MetaPhlAn DB is installed and installs it if not. All other parameters are ignored.")
    arg('--offline', action='store_true',
        help="If used, MetaPhlAn will not check for new database updates.")
    arg('--force_download', action='store_true',
        help="Force the re-download of the latest MetaPhlAn database.")
    arg('--read_min_len', type=int, default=70,
        help="Specify the minimum length of the reads to be considered when parsing the input file with "
             "'read_fastx.py' script, default value is 70")
    arg('--verbose', action='store_true',
        help="Makes MetaPhlAn verbose")
    arg('-v', '--version', action='version',
        version="MetaPhlAn version {} ({})".format(__version__, __date__),
        help="Prints the current MetaPhlAn version and exit")
    arg("-h", "--help", action="help", help="show this help message and exit")
    return p.parse_args()


def check_params(args):
    """Checks the mandatory command line arguments of the script

    Args:
        args (namespace): the arguments to check
    """
    if args.inp:
        for file in args.inp.split(','):
            if not os.path.exists(file):
                error("Error: Input file not found: {}".format(file), exit=True)
    if (args.subsampling_paired is None) and (args.forward or args.reverse):
        error("You specified forward and reverse reads, but did not specify --subsampling_paired. Provide the input as standard input or with -inp if you don't want to subsample the reads with --subsampling_paired", init_new_line = True, exit = True)
    if not (args.subsampling_seed.lower() == 'random' or args.subsampling_seed.isdigit()):
        error("The --subsampling_seed parameter is not accepted. It should contain an integer number or \"random\".", exit=True)
    if args.inp and ',' in args.inp and not args.mapout:
        error("--mapout needs to be specified when multiple FASTQ or FASTA files (comma separated) are provided", exit=True)
    if args.mapout and os.path.exists(args.mapout) and not args.force and not args.profile_vsc:
        error("Mapping output file detected: {}\n. Please use it as input or remove it if you want to re-perform the BowTie2 run".format(args.mapout), exit=True)
    if args.input_type == 'sam' and not args.nreads and not args.long_reads:
        error('The --nreads parameter must be specified when using input files in SAM format', exit=True)
    if args.input_type not in ['fasta', 'fastq'] and args.no_map:
        error('The --no_map parameter can only be used with FASTA or FASTQ input formats', exit=True)
    if args.CAMI_format_output and args.t != 'rel_ab':
        error('The --CAMI_format_output parameter can only be used with the default analysis type (rel_ab)', exit=True)    
    if args.biom_format_output and args.t != 'rel_ab':
        error('The --biom_format_output parameter can only be used with the default analysis type (rel_ab)', exit=True)
    if args.force and args.input_type not in ['fastq', 'fasta']:
        warning("The --force parameter can only be used with FASTA or FASTQ input formats, --force option will be ignored", init_new_line = True)
        args.force = False       
    if args.force and os.path.exists(args.mapout):
        os.remove(args.mapout)
        warning("Previous mapping output file has been removed from: {}".format(args.mapout)) 
    if args.profile_vsc and args.input_type == 'mapout':
        error("The Viral Sequence Clusters mode requires fasta/q or sam input!", init_new_line = True, exit = True)
    if args.profile_vsc and (args.input_type == 'fasta' or  args.input_type == 'fastq') and not args.samout:
        error("The Viral Sequence Clusters mode with fasta files requires to specify a SAM output file with the -s parameter", init_new_line = True, exit = True)
    if args.profile_vsc and not args.vsc_out:
        error("The Viral Sequence Clusters mode requires to specify a profiling output file with the --vsc_out parameter", init_new_line = True, exit = True)    
    if args.mapping_subsampling and args.subsampling is None:
        error("The --mapping_subsampling parameter should be used together with the --subsampling parameter.", init_new_line = True, exit = True)
    if args.subsampling and args.subsampling_paired:
        error("You specified both --subsampling and --subsampling_paired. Choose only one of the two options.", init_new_line = True, exit = True) 
    if not args.long_reads and not args.split_reads and not args.mapping_subsampling and (args.subsampling or args.subsampling_paired) and args.subsampling_paired < 10000:
        warning("The specified subsampling is below the recommended minimum of 10,000 reads.", init_new_line = True) 
    if not args.mapping_subsampling and ((args.subsampling is not None) or (args.subsampling_paired is not None)):
        if args.input_type not in ['fastq', 'fasta']:
            error("The --subsampling/--subsampling_paired parameter requires FASTQ or FASTA input or --mapping_subsampling", init_new_line = True, exit = True)
        if args.subsampling is not None:
            if not args.long_reads:
                warning("If you use --subsampling the reads will be subsampled NOT taking into account paired information. Use --subsampling_paired if you have paired ends reads.", init_new_line = True)
            args.subsampling_paired=False
            if not args.inp:
                error("Input reads for the subsampling must be provided as parameter. Stdin input is not allowed.", init_new_line = True, exit = True)
        elif args.subsampling_paired is not None and not args.long_reads:
            args.subsampling = args.subsampling_paired
            args.subsampling_paired = True

            if args.forward is None or args.reverse is None:
                error("If you specify --subsampling_paired you have to provide forward and reverse reads as -1 and -2 respectively. \n Reads are assumed to be in the same order in the two files.", init_new_line = True, exit = True)
            if not os.path.exists(args.forward) or not os.path.exists(args.reverse):
                error("Error: Files passed with -1 ({}) or -2 ({}) not found.".format(args.forward, args.reverse), init_new_line = True, exit = True)
            if args.inp is not None:
                warning("Since --subsampling_paired has been specified, reads are taken from -1 ({}) and -2 ({}), not from -inp.".format(args.forward, args.reverse), init_new_line = True)
            args.inp = args.forward + ',' + args.reverse
    if args.min_mapq_val is None:
        args.min_mapq_val = 50 if args.long_reads else 5
    if args.vsc_breadth is None:
        args.vsc_breadth = 0.5 if args.long_reads else 0.75
    ## checks for long reads
    if args.long_reads:
        if args.inp is None and args.input_type in ['fastq', 'fasta']:
            error("The input cannot be stdin when using long reads, please provide the path to the fasta or fastq file", exit=True)
        if args.split_reads:
            error("The --split_reads parameter is not accepted with --long_reads, please remove --long_reads if you want to split your long reads", exit=True)
        if args.input_type == 'sam' and not args.nbases:
            error('The --nbases parameter must be specified when using input files in SAM format with long reads', exit=True)
        ## This is not needed anymore as the file is now read by read_fastx, that accepts multiploe comma-separated files
        # if args.inp and "," in args.inp:
        #     warning("Multiple input files are not supported with long reads, only the first file will be used", init_new_line = True)
        #     args.inp = args.inp.split(",")[0]
        if args.subsampling_paired:
            error("The --subsampling_paired parameter is not accepted with --long_reads", exit=True)
        if not args.mapping_subsampling and args.subsampling is not None and args.subsampling < 10000*150: # TODO check threshold
            warning("The specified subsampling is below the recommended minimum of 1,500,000 bases.", init_new_line = True)
    if args.profile_vsc and (args.long_reads or args.split_reads):
        warning("Long-read sequencing techniques may not be suited for viral profiling", init_new_line = True)
    return args

def main():
    t0 = time.time()
    args = read_params(sys.argv)
    if args.verbose:
        info("Start MetaPhlAn execution", stderr=True)
    args = check_params(args)
    metaphlan_runner = Metaphlan(args)
    metaphlan_runner.run_metaphlan()
    exec_time = time.time() - t0    
    if args.verbose:
        info("Finish MetaPhlAn execution ({} seconds)".format(round(exec_time, 2)), stderr=True)


if __name__ == '__main__':
    main()
