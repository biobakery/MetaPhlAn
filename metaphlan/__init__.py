import bz2
import hashlib
import os
import re
import stat
import subprocess as subp
import sys
import tarfile
import time
import zipfile
from glob import glob, iglob
import urllib.request

def remove_prefix(text):
        return re.sub(r'^[a-z]__', '', text)

def read_and_split(ofn):
    return (l.decode('utf-8').strip().split('\t') for l in ofn)

def read_and_split_line(line):
    return line.decode('utf-8').strip().split('\t')

def plain_read_and_split(ofn):
    return (l.strip().split('\t') for l in ofn)

def plain_read_and_split_line(l):
    return l.strip().split('\t')

def mybytes(val):
    return bytes(val, encoding='utf-8')

def byte_to_megabyte(byte):
    """
    Convert byte value to megabyte
    """

    return byte / (1024.0**2)

class ReportHook():
    update_interval = 0.5
    def __init__(self):
        self.start_time = time.time()
        self.last_update = time.time()

    def report(self, blocknum, block_size, total_size):
        """
        Print download progress message
        """
        now = time.time()
        if blocknum == 0:
            self.start_time = time.time()
            if total_size > 0:
                sys.stderr.write("Downloading file of size: {:.2f} MB\n"
                                 .format(byte_to_megabyte(total_size)))
        elif (now - self.last_update) > self.update_interval:
            total_downloaded = blocknum * block_size
            status = "{:3.2f} MB ".format(byte_to_megabyte(total_downloaded))

            if total_size > 0:
                percent_downloaded = total_downloaded * 100.0 / total_size
                # use carriage return plus sys.stderr to overwrite stderr
                download_rate = total_downloaded / (now - self.start_time)
                estimated_time = (total_size - total_downloaded) / download_rate
                estimated_minutes = int(estimated_time / 60.0)
                estimated_seconds = estimated_time - estimated_minutes * 60.0
                status += ("{:3.2f} %  {:5.2f} MB/sec {:2.0f} min {:2.0f} sec "
                           .format(percent_downloaded,
                                   byte_to_megabyte(download_rate),
                                   estimated_minutes, estimated_seconds))

            status += "        \r"
            sys.stderr.write(status)
            self.last_update = now

def download(url, download_file, force=False):
    """
    Download a file from a url
    """

    if not os.path.isfile(download_file) or force:
        try:
            sys.stderr.write("\nDownloading " + url + "\n")
            file, headers = urllib.request.urlretrieve(url, download_file,
                                        reporthook=ReportHook().report)
        except EnvironmentError:
            sys.stderr.write("\nWarning: Unable to download " + url + "\n")
    else:
        sys.stderr.write("\nFile {} already present!\n".format(download_file))
        
        
def download_and_untar(download_file_name, folder, origin):
        #local path of the tarfile and md5file
    tar_file = os.path.join(folder, download_file_name + ".tar")
    md5_file = os.path.join(folder, download_file_name + ".md5")
    #Download the list of all the files in the FPT    
    url_tar_file = "{}/{}.tar".format(origin, download_file_name)
    url_md5_file = "{}/{}.md5".format(origin, download_file_name)
    # download tar and MD5 checksum
    download(url_tar_file, tar_file)
    download(url_md5_file, md5_file)
    md5_md5 = None
    md5_tar = None
    if os.path.isfile(md5_file):
        with open(md5_file) as f:
            for row in f:
                md5_md5 = row.strip().split(' ')[0]
    else:
        sys.stderr.write('File "{}" not found!\n'.format(md5_file))
    # compute MD5 of .tar.bz2
    if os.path.isfile(tar_file):
        hash_md5 = hashlib.md5()

        with open(tar_file, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)

        md5_tar = hash_md5.hexdigest()[:32]
    else:
        sys.stderr.write('File "{}" not found!\n'.format(tar_file))
    if (md5_tar is None) or (md5_md5 is None):
        sys.exit("MD5 checksums not found, something went wrong!")
    # compare checksums
    if md5_tar != md5_md5:
        sys.exit("MD5 checksums do not correspond! If this happens again, you should remove the database files and "
                 "rerun MetaPhlAn so they are re-downloaded")
    # untar
    try:
        tarfile_handle = tarfile.open(tar_file)
        tarfile_handle.extractall(path=folder)
        tarfile_handle.close()
        os.remove(tar_file)
        os.remove(md5_file)
    except EnvironmentError:
        sys.stderr.write("Warning: Unable to extract {}.\n".format(tar_file))

def download_unpack_tar(download_file_name, folder, bowtie2_build, nproc, use_zenodo):
    """
    Download the url to the file and decompress into the folder
    """

    # Create the folder if it does not already exist
    if not os.path.isdir(folder):
        try:
            os.makedirs(folder)
        except EnvironmentError:
            sys.exit("ERROR: Unable to create folder for database install: " + folder)

    # Check the directory permissions
    if not os.access(folder, os.W_OK):
        sys.exit("ERROR: The directory is not writeable: " + folder + ". "
                 "Please modify the permissions.")
        
    sys.stderr.write('\n\Downloading and uncompressing indexes\n')
    download_and_untar("{}_bt2".format(download_file_name), folder, "http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes")
    sys.stderr.write('\nDownloading and uncompressing additional files\n')
    download_and_untar(download_file_name, folder, "http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases")

    # uncompress sequences
    for bz2_file in iglob(os.path.join(folder, download_file_name + "_*.fna.bz2")):
        fna_file = bz2_file[:-4]

        if not os.path.isfile(fna_file):
            sys.stderr.write('\n\nDecompressing {} into {}\n'.format(bz2_file, fna_file))

            with open(fna_file, 'wb') as fna_h, \
                bz2.BZ2File(bz2_file, 'rb') as bz2_h:
                for data in iter(lambda: bz2_h.read(100 * 1024), b''):
                    fna_h.write(data)
        os.remove(bz2_file)  

    # build bowtie2 indexes
    if not glob(os.path.join(folder, download_file_name + "*.bt2l")):
        build_bwt_indexes(folder, download_file_name, bowtie2_build, nproc)
    else:
        try:
            subp.check_call(['bowtie2-inspect', '-n', os.path.join(folder, download_file_name)], stdout=subp.DEVNULL, stderr=subp.DEVNULL)
        except Exception as e:
            sys.stderr.write('Downloaded indexes are not compatible with the installed verion of Bowtie2\n')
            sys.stderr.write('Building indexes from the FASTA files\n')
            for btw_file in iglob(os.path.join(folder, download_file_name + "*.bt2l")):
                os.remove(btw_file)
            build_bwt_indexes(folder, download_file_name, bowtie2_build, nproc)

    try:
        for bt2 in glob(os.path.join(folder, download_file_name + "*.bt2l")):
            os.chmod(bt2, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)  # change permissions to 664
    except PermissionError as e:
        sys.stderr.write('Cannot change permission for {}. Make sure the files are readable.'.format(os.path.join(folder, download_file_name + "*.bt2l")))

    sys.stderr.write('Removing uncompressed databases\n')
    #os.remove(fna_file)

    # remove all the individual FASTA files but ViralDB
    for fna_file in iglob(os.path.join(folder, download_file_name + "_*.fna")):
        if not fna_file.endswith('_VSG.fna'):
            os.remove(fna_file)
            
   
def build_bwt_indexes(folder, download_file_name, bowtie2_build, nproc):
    sys.stderr.write('\n\nJoining FASTA databases\n')
    with open(os.path.join(folder, download_file_name + ".fna"), 'w') as fna_h:
        for fna_file in iglob(os.path.join(folder, download_file_name + "_*.fna")):
            with open(fna_file, 'r') as fna_r:
                for line in fna_r:
                    fna_h.write(line)
    fna_file = os.path.join(folder, download_file_name + ".fna")
    
    bt2_base = os.path.join(folder, download_file_name)
    bt2_cmd = [bowtie2_build, '--quiet']

    if nproc > 1:
        bt2_build_output = subp.check_output([bowtie2_build, '--usage'], stderr=subp.STDOUT)

        if 'threads' in str(bt2_build_output):
            bt2_cmd += ['--threads', str(nproc)]

    bt2_cmd += ['-f', fna_file, bt2_base]

    sys.stderr.write('\nBuilding Bowtie2 indexes\n')

    try:
        subp.check_call(bt2_cmd)
    except Exception as e:
        sys.stderr.write("Fatal error running '{}'\nError message: '{}'\n\n".format(' '.join(bt2_cmd), e))
        sys.exit(1)
    
def download_unpack_zip(url,download_file_name,folder,software_name):
    """
    Download the url to the file and decompress into the folder
    """
    
    # Check for write permission to the target folder
    if not os.access(folder, os.W_OK):
        print("WARNING: The directory is not writeable: "+
            folder + " . Please modify the permissions.")
    
    download_file=os.path.join(folder, download_file_name)
    
    download(url, download_file, True)
    
    error_during_extract=False
    
    try:
        zipfile_handle=zipfile.ZipFile(download_file)
        zipfile_handle.extractall(path=folder)
        zipfile_handle.close()
    except EnvironmentError:
        print("WARNING: Unable to extract "+software_name+".")
        error_during_extract=True
        
    if not error_during_extract:
        try:
            os.unlink(download_file)
        except EnvironmentError:
            print("WARNING: Unable to remove the temp download: " + download_file)

def resolve_latest_database(bowtie2_db,mpa_latest_url, force=False, offline=False):
    if not offline and os.path.exists(os.path.join(bowtie2_db,'mpa_latest')):
        ctime_latest_db = int(os.path.getctime(os.path.join(bowtie2_db,'mpa_latest')))
        if int(time.time()) - ctime_latest_db > 31536000:         #1 year in epoch
            os.rename(os.path.join(bowtie2_db,'mpa_latest'),os.path.join(bowtie2_db,'mpa_previous'))
            download(mpa_latest_url, os.path.join(bowtie2_db,'mpa_latest'), force=True)

    if not os.path.exists(os.path.join(bowtie2_db,'mpa_latest') or force):
        if offline:
            print("Database cannot be downloaded with the --offline option activated")
            sys.exit()        
        download(mpa_latest_url, os.path.join(bowtie2_db,'mpa_latest'))

    with open(os.path.join(bowtie2_db,'mpa_latest')) as mpa_latest:
        latest_db_version = [line.strip() for line in mpa_latest if not line.startswith('#')]
    
    return ''.join(latest_db_version)

def check_and_install_database(index, bowtie2_db, bowtie2_build, nproc, force_redownload_latest, offline):
    # Create the folder if it does not already exist
    if not os.path.isdir(bowtie2_db):
        try:
            os.makedirs(bowtie2_db)
        except EnvironmentError:
            sys.exit("ERROR: Unable to create folder for database install: " + bowtie2_db)

    if index != 'latest' and len(glob(os.path.join(bowtie2_db, "*{}*".format(index)))) >= 6:
        return index
    
    use_zenodo = False
    try:
        if not offline and urllib.request.urlopen("http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_latest").getcode() != 200:
            # use_zenodo = True
            pass
    except:
        print('WARNING: It seems that you do not have Internet access.')
        if os.path.exists(os.path.join(bowtie2_db,'mpa_latest')):
            print('WARNING: Cannot connect to the database server. The latest available local database will be used.')
            with open(os.path.join(bowtie2_db,'mpa_latest')) as mpa_latest:
                latest_db_version = [line.strip() for line in mpa_latest if not line.startswith('#')]
        else:
            print("""ERROR: Cannot find a local database. Please run MetaPhlAn using option "-x <database_name>".
            You can download the MetaPhlAn database from \n {} 
                  """.format('http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases'))
            sys.exit()

    #try downloading from the segatalab website. If fails, use zenodo
    if index == 'latest':
        mpa_latest = 'http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_latest'
        index = resolve_latest_database(bowtie2_db, mpa_latest, force_redownload_latest, offline)
    
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
        print("Database cannot be downloaded with the --offline option activated")
        sys.exit()
    # download the tar archive and decompress
    sys.stderr.write("\nDownloading MetaPhlAn database\nPlease note due to "
                     "the size this might take a few minutes\n")
    download_unpack_tar(index, bowtie2_db, bowtie2_build, nproc, use_zenodo)
    sys.stderr.write("\nDownload complete\n")
    return index
