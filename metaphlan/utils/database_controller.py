
import os
import bz2
import time
import sys
import hashlib
import urllib.request
import tarfile
from glob import glob, iglob
import pickle as pkl
import subprocess as subp
import stat
import pandas as pd
from Bio import SeqIO

try:
    from .external_exec import generate_markers_fasta
    from .util_fun import info, error, warning, byte_to_megabyte
except ImportError:
    from external_exec import generate_markers_fasta
    from util_fun import info, error, warning, byte_to_megabyte

DB_URL = 'http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases'


class MetaphlanDatabaseController(): 
    """MetaphlanDatabaseController class"""

    def set_db_dir(self, value):
        """Sets the clade path of MetaPhlAn database

        Args:
            value (str): the path to the MetaPhlAn database
        """
        self.db_dir = value

    def set_index(self, value):
        """Sets the clade path of MetaPhlAn database

        Args:
            value (str): the database index
        """
        self.index = value

    def report(self, blocknum, block_size, total_size):
        """Prints the download progress message
        
        Args:   
            blocknum (int): the block number
            block_size (int): the block size
            total_size (int): the total
        """
        if blocknum == 0:
            self.start_time = time.time()
            if total_size > 0:
                info("Downloading file of size: {:.4f} MB".format(byte_to_megabyte(total_size)), init_new_line = True)
                status = "        \r"
                sys.stderr.write(status)
        else:
            total_downloaded = blocknum * block_size
            status = "{:3.2f} MB ".format(byte_to_megabyte(total_downloaded))

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
        """Calculates the md5 of .tar.bz2 or reads the md5 file
        
        Args:   
            file_path (str): the path to the .tar.bz2 file
            md5 (bool): whether to calculate the md5 or read the md5 file
        """ 
        info('Checking md5 of {}'.format(file_path), init_new_line = True)
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
            error('File "{}" not found!'.format(file_path), init_new_line = True, exit = True)


    def download(self, url, download_file, force=False, exit=True):
        """Download a file from a url
        
        Args: 
            url (str): the url of the file to download
            download_file (str): the path to the file to download
            force (bool, optional): whether to force download the file. Defaults to False.
        """
        if not os.path.isfile(download_file) or force:
            try:
                info("Downloading " + url, init_new_line = True)
                urllib.request.urlretrieve(url, download_file, reporthook=self.report)
            except EnvironmentError as e:
                error('EnvironmentError "{}"\n Unable to download {}'.format(e, url), init_new_line = True, exit = exit)
        else:
            warning("File {} already present!".format(download_file), init_new_line = True)


    def download_and_untar(self, download_file_name, folder, origin):
        """Download a file and untar it
        
        Args:   
            download_file_name (str): the name of the file to download
            folder (str): the path to the folder to untar the file
            origin (str): the url of the file to download
        """
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
            error("MD5 checksums do not correspond!"
                  "You should remove the database files and rerun MetaPhlAn "
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
        if not os.path.isdir(self.db_dir):
            try:
                os.makedirs(self.db_dir)
            except EnvironmentError as e:
                error('EnvironmentError "{}"\n Unable to create folder for database install: {}'.format(e, self.db_dir), exit = True)

        # Check the directory permissions
        if not os.access(self.db_dir, os.W_OK):
            error("The directory is not writable: {}\n Please modify the permissions.".format(self.db_dir), exit = True)
            
        info('Downloading and uncompressing bowtie2 indexes', init_new_line = True)
        self.download_and_untar("{}_bt2".format(self.index), self.db_dir, os.path.join(DB_URL,"bowtie2_indexes"))
        info('Downloading and uncompressing additional files', init_new_line = True)
        self.download_and_untar(self.index, self.db_dir, DB_URL)

        self.download(os.path.join(DB_URL, self.index+'.nwk'), os.path.join(self.db_dir,self.index+'.nwk'), exit = False)

        # uncompress sequences
        for bz2_file in iglob(os.path.join(self.db_dir, self.index + "_*.fna.bz2")):
            fna_file = bz2_file[:-4]

            if not os.path.isfile(fna_file):
                info('Decompressing {} into {}'.format(bz2_file, fna_file), init_new_line = True)

                with open(fna_file, 'wb') as fna_h, \
                    bz2.BZ2File(bz2_file, 'rb') as bz2_h:
                    for data in iter(lambda: bz2_h.read(100 * 1024), b''):
                        fna_h.write(data)
            os.remove(bz2_file)  
 
    def check_database(self):
        """Check if all the database files are present and the database installed
        
        Returns:    
            bool: True if the database is installed, False otherwise"""
        if os.path.isdir(self.db_dir) and len(glob(os.path.join(self.db_dir, "*{}*bt2l".format(self.index)))) >= 6:
            if self.verbose:
                info('Bowtie2 indexes found', init_new_line = True)
            if os.path.exists(os.path.join(self.db_dir, self.index + "_VSG.fna")) and os.path.exists(os.path.join(self.db_dir, self.index + "_VINFO.csv")):
                if self.verbose:
                    info('ViralDB files found', init_new_line = True)
                if os.path.isfile(os.path.join(self.db_dir, "{}.pkl".format(self.index))):
                    if self.verbose:
                        info('Pickle file found', init_new_line = True)
                    mpa_pkl = os.path.join(self.db_dir, "{}.pkl".format(self.index))

                    with bz2.BZ2File(mpa_pkl, 'r') as handle:
                        self.database_pkl = pkl.load(handle)
                    return True
                else:
                    if self.verbose:
                        info('Pickle file not found ({})'.format(os.path.join(self.db_dir, "{}.pkl".format(self.index))), init_new_line = True)
            else:
                if self.verbose:
                    info('ViralDB files not found ({}, {})'.format(os.path.join(self.db_dir, self.index + "_VSG.fna"),os.path.join(self.db_dir, self.index + "_VINFO.csv")), init_new_line = True)
        else:
            if self.verbose:
                info('Bowtie2 indexes not found ({})'.format(os.path.join(self.db_dir, "*{}*bt2l".format(self.index))), init_new_line = True)
        return False


    def check_folder_exists(self):
        if not os.path.isdir(self.db_dir):
            try:
                os.makedirs(self.db_dir)
            except EnvironmentError as e:
                error('EnvironmentError "{}"\n Unable to create folder for database install: {}'.format(e, self.db_dir), exit = True)

    def install_database(self):
        """Install the database"""       

        self.check_folder_exists()
        self.download_unpack_tar()
        self.prepare_indexes()

    def check_and_install_database(self):
        """Check if the database is installed and install it if not
        
        Returns:
            str: the index of the database"""
        # check if the database is already present locally 
        if not self.force_download:
            if self.check_database():
                return self.index
        
        # not enough database files present locally and offline option is on
        if self.offline:
            error("The database cannot be downloaded with the --offline option activated and some database files for {} are missing in {}".format(self.index, self.db_dir), init_new_line = True, exit = True)

        # database not present, download and install
        info("MetaPhlAn database not present or partially present in {}. \n Downloading database\n Please note due to the size this might take a few minutes.".format(self.db_dir), init_new_line = True)
        self.install_database()      
        info("Download complete.", init_new_line = True)

        # check if the database is installed
        if self.check_database():
            return self.index
        else:
            error('Database installation failed. Please try again the installation in --verbose mode', init_new_line = True, exit = True)
    
    def get_index(self):
        """Get the index of the database"""
        return self.index
    

    def resolve_index(self):
        """Find out what is the index of the latest mpa DB available online or locally
        
        Returns:    
            str: the index of the latest mpa DB available online or locally
        """
        if self.index == 'latest':
            if not self.offline:
                # check internet connection
                try:
                    self.check_folder_exists()
                    if urllib.request.urlopen(os.path.join(DB_URL,'mpa_latest')).getcode() == 200:
                        pass
                except EnvironmentError as e:
                    warning('It seems that you do not have Internet access.', init_new_line = True)
                    # if you do not have internet access
                    if os.path.exists(os.path.join(self.db_dir,'mpa_latest')):
                        with open(os.path.join(self.db_dir,'mpa_latest')) as mpa_latest:
                            latest_db_version = ''.join([line.strip() for line in mpa_latest if not line.startswith('#')])
                            self.index = latest_db_version
                            warning('Cannot connect to the database server. The latest available local database will be used.'.format(self.index), init_new_line = True)
                            return self.index
                    else:
                        error('Cannot find a local database. Please run MetaPhlAn using option "-x <database_name>". '
                              'You can download the MetaPhlAn database from: \n {}'.format(DB_URL), init_new_line = True, exit = True)

                # download latest if not available locally
                if not os.path.exists(os.path.join(self.db_dir, 'mpa_latest')):
                    self.download(os.path.join(DB_URL, 'mpa_latest'), os.path.join(self.db_dir, 'mpa_latest'), force=True)

                else:
                    # if available, check how old it is. If too old, download a new mpa_latest
                    ctime_latest_db = int(os.path.getctime(os.path.join(self.db_dir, 'mpa_latest')))
                    if int(time.time()) - ctime_latest_db > 31536000:         #1 year in epoch
                        os.rename(os.path.join(self.db_dir, 'mpa_latest'),os.path.join(self.db_dir, 'mpa_previous'))
                        self.download(os.path.join(DB_URL, 'mpa_latest'), os.path.join(self.db_dir, 'mpa_latest'), force=True)

                        # if mpa_previous present, make the user choose to proceed with it or newer version
                        if not self.force_download:         
                            with open(os.path.join(self.db_dir,'mpa_previous')) as mpa_previous:
                                previous_db_version = ''.join([line.strip() for line in mpa_previous if not line.startswith('#')])
                            with open(os.path.join(self.db_dir, 'mpa_latest')) as mpa_latest:
                                latest_db_version = ''.join([line.strip() for line in mpa_latest if not line.startswith('#')])
                                
                            choice = ''
                            while choice.upper() not in ['Y','N']:
                                choice = input('A newer version of the database ({}) is available. Do you want to download it and replace the current one ({})?\t[Y/N]'.format(self.index, previous_db_version))

                            # if not, rename mpa_previous to mpa_latest and use it
                            if choice.upper() == 'N':
                                os.rename(os.path.join(self.db_dir,'mpa_previous'),os.path.join(self.db_dir,'mpa_latest')) 
            else:
                if not os.path.exists(os.path.join(self.db_dir, 'mpa_latest')):
                    error("Database cannot be downloaded with the --offline option activated and no existing database was detected in {}".format(self.db_dir), init_new_line = True, exit = True) 
        
            with open(os.path.join(self.db_dir, 'mpa_latest')) as mpa_latest:
                latest_db_version = ''.join([line.strip() for line in mpa_latest if not line.startswith('#')])
                self.index = latest_db_version

        return self.index

    def prepare_indexes(self):
        """Prepare for building indexes"""
        if len(glob(os.path.join(self.db_dir, self.index + "*.fna"))) > 1 and not glob(os.path.join(self.db_dir, self.index + ".fna")):
            info('Joining FASTA databases', init_new_line = True )
            if not os.path.exists(os.path.join(self.db_dir, self.index + "_VSG.fna")) and self.profile_vsc:
                error('Viral markers are missing. Please try to re-download the database', init_new_line = True, exit=True)
            with open(os.path.join(self.db_dir, self.index + ".fna"), 'w') as fna_h:
                for fna_file in iglob(os.path.join(self.db_dir, self.index + "_*.fna")):
                    with open(fna_file, 'r') as fna_r:
                        for line in fna_r:
                            fna_h.write(line)

        # remove partial FASTA file except for ViralDB
        for fna_file in iglob(os.path.join(self.db_dir, self.index + "_*.fna")):
            if not fna_file.endswith('_VSG.fna') and not fna_file.endswith('{}.fna'.format(self.index)):
                info('Removing uncompressed databases', init_new_line = True)
                os.remove(fna_file)
        
        # check bowtie2
        try:
            subp.check_call([self.bowtie2_exe, "-h"], stdout=subp.DEVNULL)
        except Exception as e:
            if self.long_reads:
                warning('OSError: "{}"\nFatal error running BowTie2 at {}. You can ignore this if only mapping with minimap2'.format(e,self.bowtie2_exe), init_new_line = True)
                return
            else:
                error('OSError: "{}"\nFatal error running BowTie2 at {}. Please check BowTie2 installation and path\n'.format(e,self.bowtie2_exe), exit=True)

        # check bowtie2 indexes, if not present, build them.
        if not glob(os.path.join(self.db_dir, self.index + "*.bt2l")):
            self.build_bwt_indexes()
        else:
            try:
                subp.check_call([self.bowtie2_exe+'-inspect', '-n', os.path.join(self.db_dir, self.index)], stdout=subp.DEVNULL, stderr=subp.DEVNULL)
            except Exception as e:
                warning('Downloaded indexes are not compatible with the installed version of Bowtie2', init_new_line = True)
                info('Building indexes from the FASTA files', init_new_line = True)
                for btw_file in iglob(os.path.join(self.db_dir, self.index + "*.bt2l")):
                    os.remove(btw_file)
                self.build_bwt_indexes()
        try:
            for bt2 in glob(os.path.join(self.db_dir, self.index + "*.bt2l")):
                os.chmod(bt2, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)  # change permissions to 664
        except PermissionError as e:
            error('PermissionError: "{}"\nCannot change permission for {}. Make sure the files are readable.'.format(e, os.path.join(self.db_dir, self.index + "*.bt2l")))

    def build_bwt_indexes(self):
        """Build BowTie indexes"""
        fna_file = os.path.join(self.db_dir, self.index + ".fna")
        
        bt2_base = os.path.join(self.db_dir, self.index)
        bt2_cmd = [self.bowtie2_build, '--quiet']

        if self.nproc > 1:
            bt2_build_output = subp.check_output([self.bowtie2_build, '--usage'], stderr=subp.STDOUT)

            if 'threads' in str(bt2_build_output):
                bt2_cmd += ['--threads', str(self.nproc)]

        bt2_cmd += ['-f', fna_file, bt2_base]

        if self.verbose:
            info('Building Bowtie2 indexes', init_new_line = True)

        try:
            subp.check_call(bt2_cmd)
        except Exception as e:
            error("Fatal error running '{}'\nError message: '{}'\n\n".format(' '.join(bt2_cmd), e), exit = True)

        if self.verbose:
            info('Removing fasta file: {}'.format(fna_file), init_new_line = True)

        os.remove(fna_file)

    def __init__(self, args):
        self.verbose = args.verbose
        self.index = args.index
        self.db_dir = args.db_dir
        self.nproc = args.nproc
        self.force_download = args.force_download
        self.offline = args.offline
        self.database_pkl = None
        self.bowtie2_exe = args.bowtie2_exe if args.bowtie2_exe else 'bowtie2'
        self.bowtie2_build = args.bowtie2_build
        self.long_reads = args.long_reads
        self.profile_vsc = args.profile_vsc

class StrainphlanDatabaseController():
    def load_database(self, verbose=True):
        """Loads the MetaPhlAn PKL database"""
        if self.database_pkl is None:
            if verbose:
                info('Loading MetaPhlAn {} database...'.format(self.get_database_name()))
            with bz2.open(self.database, 'rb') as f:
                self.database_pkl = pkl.load(f)
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
        with bz2.open(os.path.join(self.db_controller_script_folder, '{}_size.txt.bz2'.format(self.get_database_name())), 'rt') as rf:
            sgb2size = {line.strip().split('\t')[0]: int(
                line.strip().split('\t')[1]) for line in rf}
        return sgb2size
    
    def __init__(self, database, db_dir=None):
        """

        Args:
            database (pathlib.Path|str):
            db_dir (pathlib.Path|str, optional): the path to the database folder. Defaults to None.
        """
        self.db_controller_script_folder = os.path.dirname(os.path.abspath(__file__))
        if db_dir is None:
            self.default_db_folder = os.path.join(
                self.db_controller_script_folder, "..", "metaphlan_databases")
            self.default_db_folder = os.environ.get(
                'METAPHLAN_DB_DIR', self.default_db_folder)
        else:
            self.default_db_folder = db_dir
        self.database = self.resolve_database(database)
        self.database_pkl = None