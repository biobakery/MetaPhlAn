#!/usr/bin/env python
__author__ = ('Aitor Blanco-Miguez (aitor.blancomiguez@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it), '
              'Nicola Segata (nicola.segata@unitn.it), '
              'Duy Tin Truong, '
              'Francesco Asnicar (f.asnicar@unitn.it)')
__version__ = '4.0.3'
__date__ = '24 Oct 2022'


import time
import bz2
import re
import os
from glob import glob, iglob
import sys
import stat
import random
import tempfile
import argparse as ap
import subprocess as subp
import numpy as np 
import hashlib
import urllib.request
import tarfile
import pickle as pkl
from collections import Counter
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pysam
import pandas as pd
import shutil

from collections import defaultdict as defdict
try:
    from .utils import *
except ImportError:
    from utils import *

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
        if self.stat == 'avg_g':
            self.coverage = sum(nreads_v) / float(sum(rat_v))
        elif self.stat == 'avg_l':
            self.coverage = np.mean([float(n)/(np.absolute(r - self.avg_read_length) + 1) for r,n in rat_nreads])
        elif self.stat == 'tavg_g':
            wnreads = sorted([(float(n)/(np.absolute(r-self.avg_read_length)+1),(np.absolute(r - self.avg_read_length)+1) ,n) for r,n in rat_nreads], key=lambda x:x[0])
            den,num = zip(*[v[1:] for v in wnreads[ql:qr]])
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
        for marker, nreads in self.markers2nreads.items():
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
            int: the number of mapped reads
        """
        total_ab = 0
        total_reads = 0
        for clade in self.all_clades.values():
            clade.compute_coverage()
            if len(clade.children) == 0 and clade.coverage > 0:
                total_ab += clade.coverage
                total_reads += clade.nreads        
        for clade in self.all_clades.values():
            if clade.coverage > 0:
                clade.rel_abundance = round(100 * clade.coverage / total_ab, 5)
        return total_reads
    
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
        self.root = TaxClade('root', 0)
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


class MappingController:
    """MappingController interface"""
    def run_mapping(self):
        pass
    
    def get_reads2markers(self):
        pass
    
    def __init__(self):
        pass
    
class MetaphlanDatabaseController(): 
    """MetaphlanDatabaseController class"""
    def prepare_bwt_indexes(self):
        """Prepare for building BowTie indexes"""

        if not glob(os.path.join(self.bowtie2db, self.index + "*.bt2l")):
            self.build_bwt_indexes()
        else:
            try:
                subp.check_call(['bowtie2-inspect', '-n', os.path.join(self.bowtie2db, self.index)], stdout=subp.DEVNULL, stderr=subp.DEVNULL)
            except Exception as e:
                warning('Downloaded indexes are not compatible with the installed version of Bowtie2', init_new_line = True)
                info('Building indexes from the FASTA files', init_new_line = True)
                for btw_file in iglob(os.path.join(self.bowtie2db, self.index + "*.bt2l")):
                    os.remove(btw_file)
                self.build_bwt_indexes()
        try:
            for bt2 in glob(os.path.join(self.bowtie2db, self.index + "*.bt2l")):
                os.chmod(bt2, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)  # change permissions to 664
        except PermissionError as e:
            error('PermissionError: "{}"\nCannot change permission for {}. Make sure the files are readable.'.format(e, os.path.join(self.bowtie2db, self.self.index + "*.bt2l")))
        
        #remove all the individual FASTA files but ViralDB
        info('Removing uncompressed databases', init_new_line = True)
        for fna_file in iglob(os.path.join(self.bowtie2db, self.index + "_*.fna")):
            if not fna_file.endswith('_VSG.fna'):
                os.remove(fna_file)

    def set_bowtie2db(self, value):
        self.bowtie2db = value

    def set_index(self, value):
        self.index = value

    def build_bwt_indexes(self):
        """Build BowTie indexes"""
        info('Joining FASTA databases', init_new_line = True )
        if len(glob(os.path.join(self.bowtie2db, self.index + "*.fna"))) > 1:
            with open(os.path.join(self.bowtie2db, self.index + ".fna"), 'w') as fna_h:
                for fna_file in iglob(os.path.join(self.bowtie2db, self.index + "_*.fna")):
                    with open(fna_file, 'r') as fna_r:
                        for line in fna_r:
                            fna_h.write(line)
        fna_file = os.path.join(self.bowtie2db, self.index + ".fna")
        
        bt2_base = os.path.join(self.bowtie2db, self.index)
        bt2_cmd = [self.bowtie2_build, '--quiet']

        if self.nproc > 1:
            bt2_build_output = subp.check_output([self.bowtie2_build, '--usage'], stderr=subp.STDOUT)

            if 'threads' in str(bt2_build_output):
                bt2_cmd += ['--threads', str(self.nproc)]

        bt2_cmd += ['-f', fna_file, bt2_base]

        info('Building Bowtie2 indexes', init_new_line = True)

        try:
            subp.check_call(bt2_cmd)
        except Exception as e:
            error("Fatal error running '{}'\nError message: '{}'\n\n".format(' '.join(bt2_cmd), e), exit = True)


    def report(self, blocknum, block_size, total_size):
        """Print download progress message"""
        if blocknum == 0:
            self.start_time = time.time()
            if total_size > 0:
                info("Downloading file of size: {:.4f} MB".format(byte_to_megabyte(total_size)), init_new_line = True)
        else:
            total_downloaded = blocknum * block_size
            status = "{:3.2f} MB ".format(byte_to_megabyte(total_downloaded))

            if total_size > 0.001: # don't do it if very small file
                percent_downloaded = total_downloaded * 100.0 / total_size
                # use carriage return plus sys.stderr to overwrite stderr
                download_rate = total_downloaded / (time.time() - self.start_time)
                estimated_time = (total_size - total_downloaded) / download_rate
                estimated_minutes = int(estimated_time / 60.0)
                estimated_seconds = estimated_time - estimated_minutes * 60.0
                status += ("{:3.2f} %  {:5.2f} MB/sec {:2.0f} min {:2.0f} sec "
                           .format(min(percent_downloaded,100),
                                   byte_to_megabyte(download_rate),
                                   estimated_minutes, estimated_seconds))

            status += "        \r"
            sys.stderr.write(status)


    def calculate_md5(self, file_path, md5):
        """Calculate md5 of .tar.bz2 and read md5"""
        if os.path.isfile(file_path):
            # read md5
            if md5:
                with open(file_path) as f:
                    for row in f:
                        md5_md5 = row.strip().split(' ')[0]
                        return md5_md5
            # calculate md5 of .tar.bz2
            else:
                hash_md5 = hashlib.md5()
                with open(file_path, "rb") as f:
                    for chunk in iter(lambda: f.read(4096), b""):
                        hash_md5.update(chunk)
                return hash_md5.hexdigest()[:32]
        else:
            error('File "{}" not found!'.format(file_path), init_new_line = True)


    def download(self, url, download_file, force=False):
        """Download a file from a url"""
        if not os.path.isfile(download_file) or force:
            try:
                info("Downloading " + url, init_new_line = True)
                urllib.request.urlretrieve(url, download_file, reporthook=self.report)
            except EnvironmentError as e:
                error('EnvironmentError "{}"\n Unable to download {}'.format(e, url), init_new_line = True)
        else:
            warning("File {} already present!".format(download_file), init_new_line = True)


    def download_and_untar(self, download_file_name, folder, origin):
        """Download a file and untar it"""
        # local path of the tarfile and md5file
        tar_file = os.path.join(folder, download_file_name + ".tar")
        md5_file = os.path.join(folder, download_file_name + ".md5")
        # download the list of all the files in the FPT    
        url_tar_file = "{}/{}.tar".format(origin, download_file_name)
        url_md5_file = "{}/{}.md5".format(origin, download_file_name)
        # download tar and MD5 checksum
        self.download(url_tar_file, tar_file)
        self.download(url_md5_file, md5_file)

        # compute MD5 of .tar.bz2
        md5_md5 = self.calculate_md5(md5_file, md5 = True)
        md5_tar = self.calculate_md5(tar_file, md5 = False)

        if (md5_tar is None) or (md5_md5 is None):
            error("MD5 checksums not found, something went wrong!", init_new_line = True, exit = True)

        # compare checksums
        if md5_tar != md5_md5:
            error("MD5 checksums do not correspond! If this happens again, "
                  "you should remove the database files and rerun MetaPhlAn "
                  "so they are re-downloaded", init_new_line = True, exit = True)

        # untar
        try:
            tarfile_handle = tarfile.open(tar_file)
            tarfile_handle.extractall(path=folder)
            tarfile_handle.close()
            os.remove(tar_file)
            os.remove(md5_file)
        except EnvironmentError as e:
            error('EnvironmentError: "{}"\n Unable to extract {}'.format(e, tar_file), exit = True)


    def download_unpack_tar(self):
        """Download the url to the file and decompress into the folder"""
        # Create the folder if it does not already exist
        if not os.path.isdir(self.bowtie2db):
            try:
                os.makedirs(self.bowtie2db)
            except EnvironmentError as e:
                error('EnvironmentError "{}"\n Unable to create folder for database install: {}'.format(e, self.bowtie2db), exit = True)

        # Check the directory permissions
        if not os.access(self.bowtie2db, os.W_OK):
            error("The directory is not writable: {}\n Please modify the permissions.".format(self.bowtie2db), exit = True)
            
        info('Downloading and uncompressing indexes', init_new_line = True)
        self.download_and_untar("{}_bt2".format(self.index), self.bowtie2db, os.path.join(DB_URL,"bowtie2_indexes"))
        info('Downloading and uncompressing additional files', init_new_line = True)
        self.download_and_untar(self.index, self.bowtie2db, DB_URL)

        # uncompress sequences
        for bz2_file in iglob(os.path.join(self.bowtie2db, self.index + "_*.fna.bz2")):
            fna_file = bz2_file[:-4]

            if not os.path.isfile(fna_file):
                info('Decompressing {} into {}'.format(bz2_file, fna_file), init_new_line = True)

                with open(fna_file, 'wb') as fna_h, \
                    bz2.BZ2File(bz2_file, 'rb') as bz2_h:
                    for data in iter(lambda: bz2_h.read(100 * 1024), b''):
                        fna_h.write(data)
            os.remove(bz2_file)  


    def check_and_install_database(self):
        # Create the folder if it does not already exist
        if not os.path.isdir(self.bowtie2db):
            try:
                os.makedirs(self.bowtie2db)
            except EnvironmentError as e:
                error('EnvironmentError "{}"\n Unable to create folder for database install: '.format(e, self.bowtie2db), exit = True)
        
        # database present locally and not force download, return 
        if len(glob(os.path.join(self.bowtie2db, "*{}*".format(self.index)))) >= 7 and not self.force_download:
            return self.index
        
        # not enough database files present locally and offline option is on
        if self.offline:
            error("Database cannot be downloaded with the --offline option activated and database files for {} were not detected in {}".format(self.index, self.bowtie2db), init_new_line = True, exit = True)

        # database not present, download and install
        info("Downloading MetaPhlAn database\n Please note due to the size this might take a few minutes.", init_new_line = True)
        self.download_unpack_tar()
        self.prepare_bwt_indexes()
      
        info("Download complete.", init_new_line = True)
        
        return self.index

    def set_pkl(self):
        #mpa_pkl = 'mpa_pkl'
        #bowtie2db = 'bowtie2db'
        bt2_ext = 'bt2l'  

        if os.path.isfile(os.path.join(self.bowtie2db, "{}.pkl".format(self.index))):
            mpa_pkl = os.path.join(self.bowtie2db, "{}.pkl".format(self.index))
        else:
            error('Unable to find the mpa_pkl file at {}'.format(os.path.isfile(os.path.join(self.bowtie2db, "{}.pkl".format(self.index)))), exit=True)

        with bz2.BZ2File(mpa_pkl, 'r') as handle:
            database_pkl = pkl.load(handle)

        #if glob(os.path.join(self.bowtie2db, "{}*.{}".format(self.index, bt2_ext))):
        #    bowtie2db = os.path.join(self.bowtie2db, "{}".format(self.index))

        return database_pkl 
    
    def get_index(self):
        return self.index
    

    def resolve_index(self):
        '''Find out what is the index of the latest mpa DB available online or locally''' 

        if self.index == 'latest':
            if not self.offline:
                # check internet connection
                try:
                    if urllib.request.urlopen(os.path.join(DB_URL,'mpa_latest')).getcode() == 200:
                        pass
                except EnvironmentError as e:
                    warning('It seems that you do not have Internet access.', init_new_line = True)
                    # if you do not have internet access
                    if os.path.exists(os.path.join(self.bowtie2db,'mpa_latest')):
                        with open(os.path.join(self.bowtie2db,'mpa_latest')) as mpa_latest:
                            latest_db_version = ''.join([line.strip() for line in mpa_latest if not line.startswith('#')])
                            self.index = latest_db_version
                            warning('Cannot connect to the database server. The latest available local database will be used.'.format(self.index), init_new_line = True)
                            return self.index
                    else:
                        error('Cannot find a local database. Please run MetaPhlAn using option "-x <database_name>". '
                              'You can download the MetaPhlAn database from: \n {}'.format(DB_URL), init_new_line = True, exit = True)

                # download latest if not available locally
                if not os.path.exists(os.path.join(self.bowtie2db, 'mpa_latest')):
                    self.download(os.path.join(DB_URL, 'mpa_latest'), os.path.join(self.bowtie2db, 'mpa_latest'), force=True)

                else:
                    # if available, check how old it is. If too old, download a new mpa_latest
                    ctime_latest_db = int(os.path.getctime(os.path.join(self.bowtie2db, 'mpa_latest')))
                    if int(time.time()) - ctime_latest_db > 31536000:         #1 year in epoch
                        os.rename(os.path.join(self.bowtie2db, 'mpa_latest'),os.path.join(self.bowtie2db, 'mpa_previous'))
                        self.download(os.path.join(DB_URL, 'mpa_latest'), os.path.join(self.bowtie2db, 'mpa_latest'), force=True)

                        # if mpa_previous present, make the user choose to proceed with it or newer version
                        if not self.force_download:         
                            with open(os.path.join(self.bowtie2db,'mpa_previous')) as mpa_previous:
                                previous_db_version = ''.join([line.strip() for line in mpa_previous if not line.startswith('#')])
                            with open(os.path.join(self.bowtie2db, 'mpa_latest')) as mpa_latest:
                                latest_db_version = ''.join([line.strip() for line in mpa_latest if not line.startswith('#')])
                                
                            choice = ''
                            while choice.upper() not in ['Y','N']:
                                choice = input('A newer version of the database ({}) is available. Do you want to download it and replace the current one ({})?\t[Y/N]'.format(self.index, previous_db_version))

                            # if not, rename mpa_previous to mpa_latest and use it
                            if choice.upper() == 'N':
                                os.rename(os.path.join(self.bowtie2db,'mpa_previous'),os.path.join(self.bowtie2db,'mpa_latest')) 
            else:
                if not os.path.exists(os.path.join(self.bowtie2db, 'mpa_latest')):
                    error("Database cannot be downloaded with the --offline option activated and no existing database was detected in {}".format(self.bowtie2db), init_new_line = True, exit = True) 
        
        with open(os.path.join(self.bowtie2db, 'mpa_latest')) as mpa_latest:
            latest_db_version = ''.join([line.strip() for line in mpa_latest if not line.startswith('#')])
            self.index = latest_db_version

        self.database_pkl = self.set_pkl()

        return self.index
    

    def __init__(self, args):
        self.index = args.index
        self.bowtie2db = args.bowtie2db
        self.nproc = args.nproc
        self.force_download = args.force_download
        self.offline = args.offline
        self.database_pkl = None
        self.bowtie2_build = args.bowtie2_build


class VSCController():

    def extract_viral_mappings_line(self, o):
        """"Extraction of reads mapping to viral markers"""
        mCluster = o[2]
        mGroup = o[2].split('|')[2].split('-')[0]
                            
        if (hex(int(o[1]) & 0x10) == '0x0'): #front read
            rr=SeqRecord(Seq(o[9]),letter_annotations={'phred_quality':[ord(_)-33 for _ in o[10][::-1]]}, id=o[0])
        else:
            rr=SeqRecord(Seq(o[9]).reverse_complement(),letter_annotations={'phred_quality':[ord(_)-33 for _ in o[10][::-1]]}, id=o[0])

        return mGroup, mCluster, rr
    

    def extract_viral_mappings(self):
        """"Line by line extraction of reads mapping to viral markers"""
        CREAD=[] 
        if self.inp.startswith(".bz2"):
            ras, ras_line, inpf = read_and_split, read_and_split_line, bz2.BZ2File(self.inp, "r")
        else:
            ras, ras_line, inpf = plain_read_and_split, plain_read_and_split_line, open(self.inp)
        
        with open(self.marker_file, 'w') as outf:
            for line in inpf:
                sam_line = ras_line(line)
                #if self.check_hq_mapping(sam_line): #redo this in case sam file as input? it has never been filtered?to check
                if sam_line[2].startswith('VDB|'):
                    mGroup, mCluster, rr = self.extract_viral_mappings_line(sam_line)
                    outf.write(mGroup+'\t'+mCluster+'\n')
                    CREAD.append(rr) 
                        
            inpf.close()                
            SeqIO.write(CREAD,self.reads_file,'fastq')

    def process_SGB_mapping(self):
        """From reads detected in the first bowtie2 run, extracts the reads/markers for a second bowtie2 run"""
        self.extract_viral_mappings()
        VSCs_markers = SeqIO.index(self.vsc_fna, "fasta")

        with open(self.marker_file) as marker_file_to_remap:

            allViralMarkers = {}
            for vmarker in marker_file_to_remap:
                group,marker = vmarker.strip().split()
                if group not in allViralMarkers:
                    allViralMarkers[group] = [marker]
                else:
                    allViralMarkers[group].append(marker)

            selectedMarkers=[]
            for grp,v in allViralMarkers.items():

                cv=Counter(v)
                if (len(cv) > 1):
                    normalized_occurrencies=sorted([(mk,mk_occurrencies,len(VSCs_markers[mk].seq)) for (mk,mk_occurrencies) in cv.items()],key=lambda x: x[1]/x[2],reverse=True)
                    topMarker= VSCs_markers[normalized_occurrencies[0][0]] #name of the viralGenome
                else:
                    topMarker=VSCs_markers[v[0]] # the only hitted viralGenome is the winning one

                selectedMarkers.append(topMarker)

            SeqIO.write(selectedMarkers, self.top_marker_file, 'fasta')

    def check_vsc_files(self):
        """Check if files are present to map to viral database"""
        if not os.path.exists(self.marker_file) or not os.path.exists(self.reads_file):
            error('There was an error in the VSCs file lookup.\n \
                It may be that there are not enough reads to profile (not enough depth).\n \
                Passing without reporting any virus.', init_new_line = True, exit=True)

        if os.stat(self.marker_file).st_size == 0:
            warning('No reads aligning to VSC markers in this file.', init_new_line = True, exit=True)        

    
    def initialize_mapping(self):
        """Initialize parameters for mapping to viral database"""
        self.database_controller.set_bowtie2db(self.tmp_dir)
        self.database_controller.set_index(os.path.basename(self.top_marker_file).split('.')[0])
        self.database_controller.prepare_bwt_indexes()

        # set parameters for bowtie mapping
        self.mapping_controller.set_inp(self.reads_file)
        self.mapping_controller.set_bowtie2db(self.tmp_dir)
        self.mapping_controller.set_index(os.path.basename(self.top_marker_file).split('.')[0])        
        self.mapping_controller.set_input_type('fastq')
        self.mapping_controller.set_samout(os.path.basename(self.vscBamFile).split('.')[0]+'.sam')
        self.mapping_controller.set_bowtie2out(None)
        self.mapping_controller.set_mapping_parameters(None)
        self.mapping_controller.run_mapping()

    def vsc_bowtie2(self):
        """Run bowtie2 only on viral markers and interested reads"""

        self.check_vsc_files()        
        self.initialize_mapping()

        try:
            stv_command = ['samtools','view','-bS', os.path.basename(self.vscBamFile).split('.')[0]+'.sam','-@',str(self.nproc)]
            sts_command = ['samtools','sort','-','-@',str(self.nproc),'-o',self.vscBamFile]
            sti_command = ['samtools','index',self.vscBamFile]

            p3 = subp.Popen(sts_command, stdin=subp.PIPE)
            p2 = subp.Popen(stv_command, stdin=subp.PIPE, stdout=p3.stdin)

            p2.communicate()
            p3.communicate()

            subp.check_call(sti_command)

        except Exception as e:
            error('Error: "{}"\nFatal error running BowTie2 for Viruses\n'.format(e), init_new_line = True, exit = True) 

        if not os.path.exists(self.vscBamFile):
            error('Error:\nUnable to create BAM FILE for viruses.', init_new_line = True, exit = True) 

    def vsc_parsing(self):
        """Parsing of theoutput sam file from mapping viral markers to reads"""
        try:
            bamHandle = pysam.AlignmentFile(self.vscBamFile, "rb")
        except Exception as e:
            error('Error: "{}"\nCheck PySam is correctly working\n'.format(e), exit = True)

        VSC_report=[]

        for c, length in zip(bamHandle.references,bamHandle.lengths):
            coverage_positions = {}
            for pileupcolumn in bamHandle.pileup(c):
                tCoverage = 0
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip \
                            and pileupread.alignment.query_qualities[pileupread.query_position] >= 20 \
                            and pileupread.alignment.query_sequence[pileupread.query_position].upper() in ('A','T','C','G'):
                            tCoverage +=1
                if tCoverage >= 1:
                    coverage_positions[pileupcolumn.pos] = tCoverage
            breadth = float(len(coverage_positions.keys()))/float(length)
            if breadth > 0:
                cvals=list(coverage_positions.values())
                VSC_report.append({'M-Group/Cluster':c.split('|')[2].split('-')[0], 'genomeName':c, 'len':length, 'breadth_of_coverage':breadth, 'depth_of_coverage_mean': np.mean(cvals), 'depth_of_coverage_median': np.median(cvals)})
        
        if bamHandle:
            bamHandle.close()

        self.create_vsc_report(VSC_report)
        shutil.rmtree(self.tmp_dir)



    def create_vsc_report(self, VSC_report):
        """Create a report of the VSC that mapped"""
        with open(self.vsc_out,'w') as outf:
            outf.write('#{}\n'.format(self.index))
            outf.write('#{}\n'.format(' '.join(sys.argv)))
            outf.write('#SampleID\t{}\n'.format(self.sample_id))
            vsc_out_df = pd.DataFrame.from_dict(VSC_report).query('breadth_of_coverage >= {}'.format(self.vsc_breadth))
            vsc_info_df = pd.read_table(self.vsc_vinfo, sep='\t')
            vsc_out_df = vsc_out_df.merge(vsc_info_df, on='M-Group/Cluster').sort_values(by='breadth_of_coverage', ascending=False).set_index('M-Group/Cluster')
            vsc_out_df.to_csv(outf,sep='\t',na_rep='-')

    def run_analysis(self):
        """All the steps to run analysis on VSC"""
        self.process_SGB_mapping()
        self.vsc_bowtie2()
        self.vsc_parsing()


    def __init__(self, args, mapping_controller, database_controller):
        # parsing previous bowtie2out run
        self.nproc = args.nproc
        self.sample_id = args.sample_id
        self.tmp_dir = tempfile.mkdtemp(dir=args.tmp_dir) 
        self.marker_file = os.path.join(self.tmp_dir,'v_markers.fa')
        self.top_marker_file = os.path.join(self.tmp_dir,'v_top_markers.fna')
        self.reads_file= os.path.join(self.tmp_dir,'v_reads.fq')
        self.inp = args.input_type if args.input_type == 'sam' else args.samout
        
        # mapping 
        self.vsc_out = args.vsc_out
        self.vsc_breadth = args.vsc_breadth     
        self.vscBamFile = os.path.join(self.tmp_dir,'v_align.bam')
        self.mapping_controller = mapping_controller
        self.index = self.mapping_controller.get_index()
        self.database_controller = database_controller

        #database  
        self.bowtie2db = args.bowtie2db
        self.vsc_fna = os.path.join(self.bowtie2db,"{}_VSG.fna".format(self.database_controller.get_index()))
        self.vsc_vinfo = os.path.join(self.bowtie2db,"{}_VINFO.csv".format(self.database_controller.get_index()))



class Bowtie2Controller(MappingController):
    """Bowtie2Controller class"""

    def set_inp(self, value):
        self.inp = value

    def set_bowtie2db(self, value):
        self.bowtie2db = value

    def set_input_type(self, value):
        self.input_type = value    
    
    def set_samout(self, value):
        self.samout = value

    def set_bowtie2out(self, value):
        self.bowtie2out = value 
    
    def set_mapping_parameters(self, value):
        self.min_alignment_len = value 
        self.nreads = value 
        self.min_mapq_val = value 
        self.no_map = value 

    def get_index(self):
        return self.index
    
    def set_index(self, value):
        self.index = value
    
    def get_sample_id(self):
        return self.sample_id

    def check_bowtie2_database(self):
        """Checks the presence and consistency of the Bowtie2 database"""
        if glob(os.path.join(self.bowtie2db, "{}*.{}".format(self.index, 'bt2*'))):
            
            if glob(os.path.join(self.bowtie2db, "{}*.{}".format(self.index, 'bt2*')))[0].endswith('l'):
                bt2_ext = 'bt2l'
            else:
                bt2_ext = 'bt2'

            self.bowtie2db = os.path.join(self.bowtie2db, "{}".format(self.index))

        else:
            error('No MetaPhlAn BowTie2 database was found at: {}'.format(
                self.bowtie2db), exit=True)
        #bt2_ext = 'bt2l'
        if not all([os.path.exists(".".join([str(self.bowtie2db), p]))
                    for p in ["1." + bt2_ext, "2." + bt2_ext, "3." + bt2_ext, "4." + bt2_ext, "rev.1." + bt2_ext, "rev.2." + bt2_ext]]):
            error('No MetaPhlAn BowTie2 database found (--index option)!\nExpecting location {}'.format(self.bowtie2db), exit=True)
        if not (abs(os.path.getsize(".".join([str(self.bowtie2db), "1." + bt2_ext])) - os.path.getsize(".".join([str(self.bowtie2db), "rev.1." + bt2_ext]))) <= 1000):
            error('Partial MetaPhlAn BowTie2 database found at {}. Please remove and rebuild the database'.format(
                self.bowtie2db), exit=True)

    def init_mapping_arguments(self):
        """Automatically set the mapping arguments"""
        if (self.bt2_ps in ["sensitive-local", "very-sensitive-local"]) and (self.min_alignment_len is None):
            self.min_alignment_len = 100
            warning(
                'bt2_ps is set to local mode, and min_alignment_len is None, min_alignment_len has automatically been set to 100')
        self.init_bowtie2out()

    def init_bowtie2out(self):
        """Inits the Bowtie2 output file"""
        if self.no_map:
            self.bowtie2out = tempfile.NamedTemporaryFile(
                dir=self.tmp_dir).name
        elif self.bowtie2out is None:
            if self.inp is None:
                self.bowtie2out = 'stdin_map.bowtie2out.txt'
            elif stat.S_ISFIFO(os.stat(self.inp).st_mode):
                self.bowtie2out = 'fifo_map.bowtie2out.txt'
            else:
                self.bowtie2out = 'fifo_map.bowtie2out.txt'.format(self.inp)

    def get_bowtie2cmd(self):
        """Gets the command for Bowtie2 execution

        Returns:
            list: the command for the bowtie2 execution
        """
        bowtie2_cmd = [self.bowtie2_exe if self.bowtie2_exe else 'bowtie2', "--seed", "1992", "--quiet", "--no-unal", "--{}".format(self.bt2_ps),
                       "-S", "-", "-U", "-", "-x", self.bowtie2db]
        if int(self.nproc) > 1:
            bowtie2_cmd += ["-p", str(self.nproc)]
        if self.input_type == "fasta":
            bowtie2_cmd += ["-f"]
        return bowtie2_cmd

    # To refactor?
    def run_bowtie2(self):
        """Runs Bowtie2"""
        try:
            if self.inp:
                readin = subp.Popen(
                    ['/shares/CIBIO-Storage/CM/scratch/users/claudia.mengoni/tools/MetaPhlAn/metaphlan/utils/read_fastx.py', '-l', str(self.read_min_len), self.inp], stdout=subp.PIPE, stderr=subp.PIPE)
            else:
                readin = subp.Popen(['/shares/CIBIO-Storage/CM/scratch/users/claudia.mengoni/tools/MetaPhlAn/metaphlan/utils/read_fastx.py', '-l', str(self.read_min_len)],
                                    stdin=sys.stdin, stdout=subp.PIPE, stderr=subp.PIPE)
            p = subp.Popen(self.get_bowtie2cmd(),
                           stdout=subp.PIPE, stdin=readin.stdout)
            readin.stdout.close()
            lmybytes, outf = (mybytes, bz2.BZ2File(self.bowtie2out, "w")) if self.bowtie2out.endswith(
                ".bz2") else (str, open(self.bowtie2out, "w"))
            
            try:
                if self.samout:
                    if self.samout[-4:] == '.bz2':
                        sam_file = bz2.BZ2File(self.samout, 'w')
                    else:
                        sam_file = open(self.samout, 'wb')
            except IOError as e:
                error('IOError: "{}"\nUnable to open sam output file.\n'.format(
                    e), exit=True)
            for line in p.stdout:
                if self.samout:
                    sam_file.write(line)
                o = read_and_split_line(line)
                if self.check_hq_mapping(o):
                    outf.write(
                        lmybytes("\t".join([o[0], o[2].split('/')[0]]) + "\n"))

            if self.samout:
                sam_file.close()
            p.communicate()
            nreads, avg_read_len = self.get_nreads_and_avg_rlen(readin.stderr.readlines())
            outf.write(lmybytes('#nreads\t{}\n'.format(nreads)))
            outf.write(
                lmybytes('#avg_read_length\t{}'.format(avg_read_len)))
            outf.close()
            self.input_type = 'bowtie2out'
            self.inp = self.bowtie2out
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
            nreads, avg_read_len = list(
                map(float, read_fastx_stderr[0].decode().split()))
            if not nreads:
                error('Fatal error running MetaPhlAn. Total metagenome size was not estimated.\nPlease check your input files.', exit=True)
            if not avg_read_len:
                error('Fatal error running MetaPhlAn. The average read length was not estimated.\nPlease check your input files.', exit=True)
            return int(nreads), avg_read_len
        except ValueError:
            os.unlink(self.bowtie2out)
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
        reads2markers = {}
        if self.input_type == 'bowtie2out':
            for r, c in ras(inpf):
                if r.startswith('#') and 'nreads' in r:
                    nreads = int(c)
                elif r.startswith('#') and 'avg_read_length' in r:
                    avg_read_len = float(c)
                else:
                    reads2markers[r] = c
        elif self.input_type == 'sam':
            nreads = int(self.nreads)
            read_lengths = []
            for line in inpf:
                sam_line = ras_line(line)
                if self.check_hq_mapping(sam_line):
                    reads2markers[sam_line[0]] = sam_line[2].split('/')[0]
                    read_lengths.append(len(sam_line[9]))
            avg_read_len = sum(read_lengths) / len(read_lengths)
        inpf.close()
        return nreads, avg_read_len, reads2markers

    def run_mapping(self):
        """Runs all the steps of the Bowtie2 mapping"""
        self.check_bowtie2_database()
        self.init_mapping_arguments()
        self.run_bowtie2()

    def __init__(self, args, index):
        self.inp = args.inp
        self.input_type = args.input_type
        self.bt2_ps = args.bt2_ps
        self.bowtie2_exe = args.bowtie2_exe
        self.bowtie2_build = args.bowtie2_build
        self.bowtie2out = args.bowtie2out
        self.bowtie2db = args.bowtie2db
        self.samout = args.samout
        self.min_alignment_len = args.min_alignment_len
        self.min_mapq_val = args.min_mapq_val
        self.index = index
        self.nreads = args.nreads
        self.nproc = args.nproc
        self.no_map = args.no_map
        self.read_min_len = args.read_min_len
        self.tmp_dir = args.tmp_dir

class MetaphlanAnalysis:
    def get_mapped_fraction(self):
        """Gets the estimated fraction of mapped reads"""
        self.mapped_reads = self.tree.relative_abundances()
        if self.unclassified_estimation:
            self.fraction_mapped_reads = min(self.mapped_reads/float(self.n_metagenome_reads), 1.0)
        else:
            self.fraction_mapped_reads = 1.0
            
    def write_common_headers(self, outf):
        """Writes the common headers of the MetaPhlAn output

        Args:
            outf (stream): output stream
        """
        outf.write('#{}\n'.format(self.index))
        outf.write('#{}\n'.format(' '.join(sys.argv)))
        outf.write('#{} reads processed\n'.format(self.n_metagenome_reads))   
        outf.write('#{}\n'.format('\t'.join([self.sample_id_key, self.sample_id])))
    
    def report_results(self, tree, n_metagenome_reads, avg_read_length):
        """Reports the MetaPhlAn results"""
        self.tree = tree
        self.n_metagenome_reads = n_metagenome_reads
        self.avg_read_length = avg_read_length
        
    def __init__(self, args, database_controller, index):
        self.output = args.output_file
        self.sample_id_key = args.sample_id_key
        self.sample_id = args.sample_id
        self.index = index
        self.database_controller = database_controller
        self.n_metagenome_reads = None
        self.fraction_mapped_reads = None
        self.mapped_reads = None
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
                    taxpathsh = '|'.join([re.sub(r'^[a-z]__', '', name) if '_unclassified' not in name else '' for name in clade.split('|')])
                    outf.write( '\t'.join( [ leaf_taxid, rank, taxid, taxpathsh, str(relab*self.fraction_mapped_reads) ] ) + '\n' )
                    
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
    
    def report_results(self, tree, n_metagenome_reads, avg_read_length):
        """Reports the MetaPhlAn results"""
        super().report_results(tree, n_metagenome_reads, avg_read_length)
        self.get_mapped_fraction()
        if self.cami_output:
            self.report_cami_output()
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
                    outf.write( "\t".join(["UNCLASSIFIED", "-1", str(round((1-self.fraction_mapped_reads)*100,5)),""]) + "\n" )                   
                clade2abundance = self.get_clade2abundance()
                for clade, values in clade2abundance.items():
                    taxid, relab = values
                    if (clade, taxid) in self.database_controller.database_pkl['merged_taxon'] and not self.use_group_representative:
                        add_repr = '{}'.format(','.join( [ n[0] for n in self.database_controller.database_pkl['merged_taxon'][(clade, taxid)]] ))
                        outf.write( "\t".join( [clade,  taxid, str(relab*self.fraction_mapped_reads), add_repr ] ) + "\n" )
                    else:                        
                        outf.write( "\t".join( [clade,  taxid, str(relab*self.fraction_mapped_reads)] ) + "\n" )                            
            if add_repr is not None:
                warning("The metagenome profile contains clades that represent multiple species merged into a single representant. "
                        "An additional column listing the merged species is added to the MetaPhlAn output.")
            
    def __init__(self, args, database_controller, index):
        super().__init__(args, database_controller, index)
        self.tax_lev = args.tax_lev
        self.cami_output = args.CAMI_format_output
        self.use_group_representative = args.use_group_representative
        self.unclassified_estimation = args.unclassified_estimation
        
        
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
    
    def report_results(self, tree, n_metagenome_reads, avg_read_length):
        """Reports the MetaPhlAn results"""
        super().report_results(tree, n_metagenome_reads, avg_read_length)
        self.get_mapped_fraction()
        out_stream = open(self.output,"w") if self.output else sys.stdout
        with out_stream as outf:
            self.write_common_headers(outf)
            outf.write( "#Estimated reads mapped to known clades: {}\n".format(int(self.mapped_reads)))
            outf.write( "\t".join( ["#clade_name","clade_taxid","relative_abundance","coverage","estimated_number_of_reads_from_the_clade" ]) +"\n" )
            if self.unclassified_estimation:
                outf.write( "\t".join(["UNCLASSIFIED", "-1", str(round((1-self.fraction_mapped_reads)*100,5)),"-", str(self.n_metagenome_reads - self.mapped_reads)]) + "\n" )                   
            clade2abundance = self.get_clade2abundance()
            for clade, values in clade2abundance.items():
                taxid, relab, coverage, nreads = values
                outf.write( "\t".join( [clade,  taxid, str(relab*self.fraction_mapped_reads), str(coverage), str(nreads)] ) + "\n" )         
            
    def __init__(self, args, database_controller, index):
        super().__init__(args, database_controller, index)
        self.tax_lev = args.tax_lev
        self.unclassified_estimation = args.unclassified_estimation
    
class CladeProfilesAnalysis(MetaphlanAnalysis):
    def report_results(self, tree, n_metagenome_reads, avg_read_length):
        """Reports the MetaPhlAn results"""
        super().report_results(tree, n_metagenome_reads, avg_read_length)
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
    def report_results(self, tree, n_metagenome_reads, avg_read_length):
        """Reports the MetaPhlAn results"""
        super().report_results(tree, n_metagenome_reads, avg_read_length)
        clade2profiles = self.tree.clade_profiles()
        out_stream = open(self.output,"w") if self.output else sys.stdout
        with out_stream as outf:
            self.write_common_headers(outf)     
            for v in clade2profiles.values():
                outf.write( "\n".join(["\t".join([str(a),str(b/float(self.n_metagenome_reads)) if self.n_metagenome_reads else str(b)])
                                for a,b in v if b > 0.0]) + "\n" )
            
    def __init__(self, args, database_controller, index):
        super().__init__(args, database_controller, index)
    
class MarkerPresenceTableAnalysis(MetaphlanAnalysis):
    def report_results(self, tree, n_metagenome_reads, avg_read_length):
        """Reports the MetaPhlAn results"""
        super().report_results(tree, n_metagenome_reads, avg_read_length)
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
                       self.stat, self.stat_q, self.perc_nonzero, self.avoid_disqm, self.avg_read_length)
            
    def separate_reads2markers(self, reads2markers):
        """Separates the viral hits from the reads to markers dictionary

        Args:
            reads2markers (dict): full reads to markers dictionary

        Returns:
            (dict, dict): tuple of dictionaries for the SGBs/EUKs and viral contigs
        """
        return {r: m for r, m in reads2markers.items() if ('SGB' in m or 'EUK' in m) and not 'VDB' in m}, {r: m for r, m in reads2markers.items() if 'VDB' in m and not ('SGB' in m or 'EUK' in m)}

    
    def subsample_reads(self, reads2markers):
        """Subsamples the reads using the mapping results 

        Args:
            reads2markers (dict): full reads to markers dictionary

        Returns:
            dict: the subsampled reads to markers dictionary
        """
        if self.subsampling >= self.n_metagenome_reads:
            warning("WARNING: The specified subsampling ({}) is equal or higher than the original number of reads ({}). Subsampling will be skipped.".format(self.subsampling, self.n_metagenome_reads), init_new_line=True)
        elif self.subsampling < 10000:
            warning("The specified subsampling ({}) is below the recommended minimum of 10,000 reads.".format(self.subsampling), init_new_line=True)
        else:
            reads2markers =  dict(sorted(reads2markers.items()))
            if self.subsampling_seed.lower() != 'random':
                random.seed(int(self.subsampling_seed))
            reads2filtmarkers = {}
            sgb_reads2markers, viral_reads2markers = self.separate_reads2markers(reads2markers)            
            n_sgb_mapped_reads = int((len(sgb_reads2markers) * self.subsampling) / self.n_metagenome_reads)
            reads2filtmarkers = { r:sgb_reads2markers[r] for r in random.sample(list(sgb_reads2markers.keys()), n_sgb_mapped_reads) }         
            n_viral_mapped_reads = int((len(viral_reads2markers) * self.subsampling) / self.n_metagenome_reads)
            reads2filtmarkers.update({ r:viral_reads2markers[r] for r in random.sample(list(viral_reads2markers.keys()), n_viral_mapped_reads) })            
            reads2markers = reads2filtmarkers
            sgb_reads2markers.clear()
            viral_reads2markers.clear()
            self.n_metagenome_reads = self.subsampling
        return reads2markers

    def parse_mapping(self):
        """Parses the mapping results into a marker to reads dictionary

        Returns:
            dict: marker to reads dictionary
        """
        self.n_metagenome_reads, self.avg_read_length, reads2markers = self.mapping_controller.get_reads2markers()   
        self.build_taxonomy_tree()     
        if self.subsampling is not None:
            reads2markers = self.subsample_reads(reads2markers)
        elif self.n_metagenome_reads < 10000:
            warning("The number of reads in the sample ({}) is below the recommended minimum of 10,000 reads.".format(self.n_metagenome_reads))            
        markers2reads = defdict(set)   
        for r, m in reads2markers.items():
            markers2reads[m].add(r)
        self.add_reads_to_tree(markers2reads)
        if self.no_map:
            os.remove(self.mapping_controller.inp)
    
    def add_reads_to_tree(self, markers2reads):
        """Adds reads mapping to the markers to the tree structure

        Args:
            markers2reads (dict): dictionary with the markers to reads information
        """
        for marker,reads in sorted(markers2reads.items(), key=lambda pars: pars[0]):
            if marker not in self.tree.markers2lens:
                continue
            self.tree.add_reads( marker, len(reads),
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
        if self.install:
            self.index=self.database_controller.check_and_install_database()
            info('The database has been installed ({})'.format(self.index), stderr=True, exit=True)
        
        #self.database_controller.check_database()    
        if self.input_type in ['fastq', 'fasta']:            
             self.mapping_controller.run_mapping()        
        self.parse_mapping()
        self.metaphlan_analysis.report_results(self.tree, self.n_metagenome_reads, self.avg_read_length)
        if self.profile_vsc:
            self.vsc_controller.run_analysis()

    def __init__(self, args):
        
        self.verbose = args.verbose
        self.database_controller = MetaphlanDatabaseController(args)
        self.index = self.database_controller.resolve_index()
        #here should be the code choosing the mapping controller in the future
        self.mapping_controller = Bowtie2Controller(args, self.index)
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
        self.subsampling_seed = args.subsampling_seed
        self.install = args.install
        self.offline = args.offline
        self.force_download = args.force_download
        self.tree = None
        self.n_metagenome_reads = None
        self.avg_read_length = None
        self.profile_vsc = args.profile_vsc
        if args.profile_vsc:
            self.vsc_controller = VSCController(args, self.mapping_controller, 
                                                  self.database_controller)


DEFAULT_DB_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "metaphlan_databases")
DEFAULT_DB_FOLDER = os.environ.get('METAPHLAN_DB_DIR', DEFAULT_DB_FOLDER)
DB_URL = 'http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases'

def read_params():
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
        "* an intermediary mapping file of the metagenome generated by a previous MetaPhlAn run \n"
        "If the input file is missing, the script assumes that the input is provided using the standard \n"
        "input, or named pipes.\n"
        "IMPORTANT: the type of input needs to be specified with --input_type")
    g = p.add_argument_group('Required arguments')
    arg = g.add_argument
    input_type_choices = ['fastq', 'fasta', 'bowtie2out', 'sam']
    arg('--input_type', choices=input_type_choices, help="set whether the input is the FASTA file of metagenomic reads or \n"
        "the SAM file of the mapping of the reads against the MetaPhlAn db.\n")
    g = p.add_argument_group('Mapping arguments')
    arg = g.add_argument
    arg('--force', action='store_true',
        help="Force profiling of the input file by removing the bowtie2out file")
    arg('--bowtie2db', metavar="METAPHLAN_BOWTIE2_DB", type=str, default=DEFAULT_DB_FOLDER,
        help=("Folder containing the MetaPhlAn database. You can specify the location by exporting the DEFAULT_DB_FOLDER variable in the shell."
              "[default "+DEFAULT_DB_FOLDER+"]\n"))
    arg('-x', '--index', type=str, default='latest',
        help=("Specify the id of the database version to use. "
              "If \"latest\", MetaPhlAn will get the latest version.\n"
              "If an index name is provided, MetaPhlAn will try to use it, if available, and skip the online check.\n"
              "If the database files are not found on the local MetaPhlAn installation they\n"
              "will be automatically downloaded [default latest]\n"))
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
    arg('--bowtie2out', metavar="FILE_NAME", type=str, default=None,
        help="The file for saving the output of BowTie2")
    arg('--min_mapq_val', type=int, default=5,
        help="Minimum mapping quality value (MAPQ) [default 5]")
    arg('--no_map', action='store_true',
        help="Avoid storing the --bowtie2out map file")
    arg('--tmp_dir', metavar="", default=None, type=str,
        help="The folder used to store temporary files [default is the OS "
             "dependent tmp dir]")
    g = p.add_argument_group('Post-mapping arguments')
    arg = g.add_argument
    stat_choices = ['avg_g', 'avg_l', 'tavg_g',
                    'tavg_l', 'wavg_g', 'wavg_l', 'med']
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
    arg('--unclassified_estimation', action='store_true',
        help="Scale relative abundances to the number of reads mapping to identified clades in order to estimate unclassified taxa\n")
    arg('--biom', '--biom_output_file',  metavar="biom_output", type=str, default=None,
        help="If requesting biom file output: The name of the output file in biom format \n")
    arg('--mdelim', '--metadata_delimiter_char',  metavar="mdelim", type=str, default="|",
        help="Delimiter for bug metadata: - defaults to pipe. e.g. the pipe in k__Bacteria|p__Proteobacteria \n")
    g = p.add_argument_group('Viral Sequence Clusters Analisys')
    arg = g.add_argument
    arg("--profile_vsc", action="store_true",help="Add this parameter to profile Viruses with VSCs approach.")
    arg("--vsc_out", help="Path to the VSCs breadth-of-coverage output file", default="mp3_viruses.csv")
    arg("--vsc_breadth", help="Minimum Breadth of Coverage for a Viral Group to be reported.\n"
    "Default is 0.75 (at least 75 percent breadth to report)", default=0.75,type=float)
    g = p.add_argument_group('Other arguments')
    arg = g.add_argument
    arg('--nproc', metavar="N", type=int, default=4,
        help="The number of CPUs to use for parallelizing the mapping [default 4]")
    arg('--subsampling', type=int, default=None,
        help="Specify the number of reads to be considered from the input metagenomes [default None]")
    arg('--mapping_subsampling', action='store_true',
        help="If used, the subsamping will be done on the mapping results instead of on the reads.")
    arg('--subsampling_seed', type=str, default='1992',
        help="Random seed to use in the selection of the subsampled reads. Choose \"random\r for a random behaviour")
    arg('--subsampling_output', type=str, default=None,
        help="The output file for the subsampled reads. If not specified the subsampled reads will not be saved.")
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
    if not (args.subsampling_seed.lower() == 'random' or args.subsampling_seed.isdigit()):
        error("The --subsampling_seed parameter is not accepted. It should contain an integer number or \"random\".", exit=True)
    if args.inp and ',' in args.inp and not args.bowtie2out:
        error("--bowtie2out needs to be specified when multiple FASTQ or FASTA files (comma separated) are provided", exit=True)
    if args.bowtie2out and os.path.exists(args.bowtie2out) and not args.force and not args.profile_vsc:
        error("BowTie2 output file detected: {}\n. Please use it as input or remove it if you want to re-perform the BowTie2 run".format(args.bowtie2out), exit=True)
    if not (args.subsampling_seed.lower() == 'random' or args.subsampling_seed.isdigit()):
        error('The --subsampling_seed parameter is not accepted. It should contain an integer number or \"random\"', exit=True) 
    if args.input_type == 'sam' and not args.nreads:
        error('The --nreads parameter must be specified when using input files in SAM format', exit=True)
    if args.input_type not in ['fasta', 'fastq'] and args.no_map:
        error('The --no_map parameter can only be used with FASTA or FASTQ input formats', exit=True)
    if args.CAMI_format_output and args.t != 'rel_ab':
        error('The --CAMI_format_output parameter can only be used with the default analysis type (rel_ab)', exit=True)               
    if args.force and os.path.exists(args.bowtie2out):
        os.remove(args.bowtie2out)
        warning("Previous Bowtie2 output file has been removed from: {}".format(args.bowtie2out)) 
    if args.profile_vsc and args.input_type == 'bowtie2out':
        error("The Viral Sequence Clusters mode requires fasta or sam input!", init_new_line = True, exit = True)
    if args.profile_vsc and (args.input_type == 'fasta' or  args.input_type == 'fastq') and not args.samout:
        error("The Viral Sequence Clusters mode with fasta files requires to specify a SAM output file with the -s parameter", init_new_line = True, exit = True)



def main():
    t0 = time.time()
    args = read_params()
    if args.verbose:
        info("Start MetaPhlAn execution", stderr=True)
    check_params(args)
    metaphlan_runner = Metaphlan(args)
    metaphlan_runner.run_metaphlan()
    exec_time = time.time() - t0    
    if args.verbose:
        info("Finish MetaPhlAn execution ({} seconds)".format(round(exec_time, 2)), stderr=True)


if __name__ == '__main__':
    main()
