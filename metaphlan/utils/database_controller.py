__author__ = 'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it)'
__version__ = '4.0.4'
__date__ = '17 Jan 2023'


import os
import pickle
import bz2
import bz2
import hashlib
import stat
import sys
import tarfile
import time
import zipfile
import urllib.request
from glob import glob, iglob
import subprocess as subp
from Bio import SeqIO
try:
    from .external_exec import generate_markers_fasta
    from .util_fun import info, error, warning, byte_to_megabyte
except ImportError:
    from external_exec import generate_markers_fasta
    from util_fun import info, error, warning, byte_to_megabyte

#init bad code

class ReportHook():
    """ReportHook class"""
    def report(self, blocknum, block_size, total_size):
        """Prints download progress message"""
        if blocknum == 0:
            self.start_time = time.time()
            if total_size > 0:
                sys.stderr.write("Downloading file of size: {:.2f} MB\n"
                                .format(byte_to_megabyte(total_size)))
        else:
            total_downloaded = blocknum * block_size
            status = "{:3.2f} MB ".format(byte_to_megabyte(total_downloaded))
            if total_size > 0:
                percent_downloaded = total_downloaded * 100.0 / total_size
                # use carriage return plus sys.stderr to overwrite stderr
                download_rate = total_downloaded / (time.time() - self.start_time)
                estimated_time = (total_size - total_downloaded) / download_rate
                estimated_minutes = int(estimated_time / 60.0)
                estimated_seconds = estimated_time - estimated_minutes * 60.0
                status += ("{:3.2f} %  {:5.2f} MB/sec {:2.0f} min {:2.0f} sec "
                        .format(percent_downloaded,
                                byte_to_megabyte(download_rate),
                                estimated_minutes, estimated_seconds))
            status += "        \r"
            sys.stderr.write(status)

    def __init__(self):
        self.start_time = time.time()
        
class MetaphlanDatabaseController:      
    """MetaphlanDatabaseController class"""
    def download(self, url, download_file, force=False):
        """Download a file from a url"""
        if not os.path.isfile(download_file) or force:
            try:
                info("\nDownloading " + url, stderr=True)
                urllib.request.urlretrieve(url, download_file, reporthook=ReportHook().report)
            except EnvironmentError:
                warning("\nWarning: Unable to download " + url, exit=True)
        else:
            warning("\nFile {} already present!".format(download_file), stderr=True)

    def download_unpack_tar(self, download_file_name, folder, bowtie2_build, nproc):
        """Download the url to the file and decompress into the folder"""

        # Create the folder if it does not already exist
        if not os.path.isdir(folder):
            try:
                os.makedirs(folder)
            except EnvironmentError:
                error("Unable to create folder for database install: " + folder, exit=True)

        # Check the directory permissions
        if not os.access(folder, os.W_OK):
            error("The directory is not writeable: {}. Please modify the permissions.".format(folder), exit=True)

        #local path of the tarfile and md5file
        tar_file = os.path.join(folder, download_file_name + ".tar")
        md5_file = os.path.join(folder, download_file_name + ".md5")

        #Download the list of all the files in the FPT    
        url_tar_file = "http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/{}.tar".format(download_file_name)
        url_md5_file = "http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/{}.md5".format(download_file_name)

        # download tar and MD5 checksum
        self.download(url_tar_file, tar_file)
        self.download(url_md5_file, md5_file)
        md5_md5 = None
        md5_tar = None
        if os.path.isfile(md5_file):
            with open(md5_file) as f:
                for row in f:
                    md5_md5 = row.strip().split(' ')[0]
        else:
            error('File "{}" not found!\n'.format(md5_file), exit=True)

        # compute MD5 of .tar.bz2
        if os.path.isfile(tar_file):
            hash_md5 = hashlib.md5()
            with open(tar_file, "rb") as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    hash_md5.update(chunk)
            md5_tar = hash_md5.hexdigest()[:32]
        else:
            error('File "{}" not found!\n'.format(tar_file), exit=True)
        if (md5_tar is None) or (md5_md5 is None):
            error("MD5 checksums not found, something went wrong!", exit=True)

        # compare checksums
        if md5_tar != md5_md5:
            error("MD5 checksums do not correspond! If this happens again, you should remove the database files and "
                    "rerun MetaPhlAn so they are re-downloaded", exit=True)

        # untar
        try:
            tarfile_handle = tarfile.open(tar_file)
            tarfile_handle.extractall(path=folder)
            tarfile_handle.close()
        except EnvironmentError:
            error("Warning: Unable to extract {}.\n".format(tar_file), exit=True)

        # uncompress sequences
        for bz2_file in iglob(os.path.join(folder, download_file_name + "_*.fna.bz2")):
            fna_file = bz2_file[:-4]
            if not os.path.isfile(fna_file):
                info('\n\nDecompressing {} into {}'.format(bz2_file, fna_file), stderr=True)
                with open(fna_file, 'wb') as fna_h, \
                    bz2.BZ2File(bz2_file, 'rb') as bz2_h:
                    for data in iter(lambda: bz2_h.read(100 * 1024), b''):
                        fna_h.write(data)
            os.remove(bz2_file)

        # join fasta files
        info('\n\nJoining FASTA databases', stderr = True)
        with open(os.path.join(folder, download_file_name + ".fna"), 'w') as fna_h:
            for fna_file in iglob(os.path.join(folder, download_file_name + "_*.fna")):
                with open(fna_file, 'r') as fna_r:
                    for line in fna_r:
                        fna_h.write(line)
        fna_file = os.path.join(folder, download_file_name + ".fna")

        # build bowtie2 indexes
        if not glob(os.path.join(folder, download_file_name + "*.bt2l")):
            bt2_base = os.path.join(folder, download_file_name)
            bt2_cmd = [bowtie2_build, '--quiet']
            if nproc > 1:
                bt2_build_output = subp.check_output([bowtie2_build, '--usage'], stderr=subp.STDOUT)
                if 'threads' in str(bt2_build_output):
                    bt2_cmd += ['--threads', str(nproc)]
            bt2_cmd += ['-f', fna_file, bt2_base]
            info('\nBuilding Bowtie2 indexes', stderr = True)

            try:
                subp.check_call(bt2_cmd)
            except Exception as e:
                error("Fatal error running '{}'\nError message: '{}'\n\n".format(' '.join(bt2_cmd), e), exit=True)
        try:
            for bt2 in glob(os.path.join(folder, download_file_name + "*.bt2l")):
                os.chmod(bt2, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)  # change permissions to 664
        except PermissionError as e:
            error('Cannot change permission for {}. Make sure the files are readable.'.format(os.path.join(folder, download_file_name + "*.bt2l")), exit=True)

        # remove all the individual FASTA files but ViralDB
        info('Removing uncompressed databases', stderr=True)
        os.remove(fna_file)
        for fna_file in iglob(os.path.join(folder, download_file_name + "_*.fna")):
            if not fna_file.endswith('_VSG.fna'):
                os.remove(fna_file)
        
    def download_unpack_zip(self, url,download_file_name,folder,software_name):
        """Download the url to the file and decompress into the folder"""
        if not os.access(folder, os.W_OK):
            error("The directory is not writeable: {} . Please modify the permissions.".format(folder), exit=True)        
        download_file=os.path.join(folder, download_file_name)        
        self.download(url, download_file, True)        
        error_during_extract=False        
        try:
            zipfile_handle=zipfile.ZipFile(download_file)
            zipfile_handle.extractall(path=folder)
            zipfile_handle.close()
        except EnvironmentError:
            warning("Unable to extract {}.".format(software_name))
            error_during_extract=True            
        if not error_during_extract:
            try:
                os.unlink(download_file)
            except EnvironmentError:
                warning("Unable to remove the temp download: " + download_file)

    def resolve_latest_database(self, bowtie2_db,mpa_latest_url, force=False, offline=False):
        """Resolves the latest MetaPhlAn database"""
        if not offline and os.path.exists(os.path.join(bowtie2_db,'mpa_latest')):
            ctime_latest_db = int(os.path.getctime(os.path.join(bowtie2_db,'mpa_latest')))
            if int(time.time()) - ctime_latest_db > 31536000:         #1 year in epoch
                os.rename(os.path.join(bowtie2_db,'mpa_latest'),os.path.join(bowtie2_db,'mpa_previous'))
                self.download(mpa_latest_url, os.path.join(bowtie2_db,'mpa_latest'), force=True)
        if not os.path.exists(os.path.join(bowtie2_db,'mpa_latest') or force):
            if offline:
                error("Database cannot be downloaded with the --offline option activated", exit=True)    
            self.download(mpa_latest_url, os.path.join(bowtie2_db,'mpa_latest'))
        with open(os.path.join(bowtie2_db,'mpa_latest')) as mpa_latest:
            latest_db_version = [line.strip() for line in mpa_latest if not line.startswith('#')]
        return ''.join(latest_db_version)

    def check_and_install_database(self, index, bowtie2_db, bowtie2_build, nproc, force_redownload_latest, offline):
        """Checks and install the MetaPhlAn database"""
        # Create the folder if it does not already exist
        if not os.path.isdir(bowtie2_db):
            try:
                os.makedirs(bowtie2_db)
            except EnvironmentError:
                error("Unable to create folder for database install: " + bowtie2_db, exit=True)
        if index != 'latest' and len(glob(os.path.join(bowtie2_db, "*{}*".format(index)))) >= 6:
            return index
        try:
            if not offline and urllib.request.urlopen("http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_latest").getcode() != 200:
                pass
        except:
            warning('It seems that you do not have Internet access.')
            if os.path.exists(os.path.join(bowtie2_db,'mpa_latest')):
                warning('Cannot connect to the database server. The latest available local database will be used.')
                with open(os.path.join(bowtie2_db,'mpa_latest')) as mpa_latest:
                    latest_db_version = [line.strip() for line in mpa_latest if not line.startswith('#')]
            else:
                error("""Cannot find a local database. Please run MetaPhlAn using option "-x <database_name>".
                You can download the MetaPhlAn database from \n {} 
                    """.format('http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases'), exit=True)

        #try downloading from the segatalab website. 
        if index == 'latest':
            mpa_latest = 'http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_latest'
            index = self.resolve_latest_database(bowtie2_db, mpa_latest, force_redownload_latest, offline)        
        if not offline and os.path.exists(os.path.join(bowtie2_db,'mpa_previous')):
            with open(os.path.join(bowtie2_db,'mpa_previous')) as mpa_previous:
                previous_db_version = ''.join([line.strip() for line in mpa_previous if not line.startswith('#')])
        
            if index != previous_db_version:
                choice = ''
                while choice.upper() not in ['Y','N']:
                    choice = input('A newer version of the database ({}) is available. Do you want to download it and replace the current one ({})?\t[Y/N]'.format(index, previous_db_version))

                if choice.upper() == 'N':
                    os.rename(os.path.join(bowtie2_db,'mpa_previous'),os.path.join(bowtie2_db,'mpa_latest'))
                    index = previous_db_version                    
        if len(glob(os.path.join(bowtie2_db, "*{}*".format(index)))) >= 7:
            return index
        if offline:
            error("Database cannot be downloaded with the --offline option activated", exit=True)
            
        # download the tar archive and decompress
        info("\nDownloading MetaPhlAn database\nPlease note due to "
                        "the size this might take a few minutes", stderr=True)
        self.download_unpack_tar(index, bowtie2_db, bowtie2_build, nproc)
        info("\nDownload complete", stderr=True)
        return index

# end of bad code

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
            set: the list of markers from the input clades
        """
        self.load_database()
        return set((marker for marker in self.database_pkl['markers']
                    if self.database_pkl['markers'][marker]['clade'] in clades))

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
