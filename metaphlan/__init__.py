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
from glob import glob
from urllib.request import urlretrieve

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
    def __init__(self):
        self.start_time = time.time()

    def report(self, blocknum, block_size, total_size):
        """
        Print download progress message
        """

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

# set the location of the database download url
DATABASE_DOWNLOAD = "https://www.dropbox.com/sh/7qze7m7g9fe2xjg/AADHWzATSQcI0CNFD0sk7MAga"
ZENODO_DATABASE_DOWNLOAD = "https://zenodo.org/record/3957592"
DBX_FILE_LIST = "https://www.dropbox.com/sh/7qze7m7g9fe2xjg/AAA4XDP85WHon_eHvztxkamTa/file_list.txt?dl=1"

def download(url, download_file, force=False):
    """
    Download a file from a url
    """

    if not os.path.isfile(download_file) or force:
        try:
            sys.stderr.write("\nDownloading " + url + "\n")
            file, headers = urlretrieve(url, download_file,
                                        reporthook=ReportHook().report)
        except EnvironmentError:
            sys.stderr.write("\nWarning: Unable to download " + url + "\n")
    else:
        sys.stderr.write("\nFile {} already present!\n".format(download_file))

def download_unpack_tar(FILE_LIST, download_file_name, folder, bowtie2_build, nproc, use_zenodo):
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

    #local path of the tarfile and md5file
    tar_file = os.path.join(folder, download_file_name + ".tar")
    md5_file = os.path.join(folder, download_file_name + ".md5")

    #Download the list of all the files in the Dropbox folder
    if not use_zenodo:
        list_file_path = os.path.join(folder, "file_list.txt")
        if not os.path.exists(list_file_path):
            download(FILE_LIST, list_file_path)

        if os.path.isfile(list_file_path):
            with open(list_file_path) as f:
                ls_f = dict( [row.strip().split() for row in f])
            url_tar_file = ls_f[download_file_name + ".tar"]
            url_md5_file = ls_f[download_file_name + ".md5"]
    else:
        url_tar_file = "https://zenodo.org/record/3957592/files/{}.tar?download=1".format(download_file_name)
        url_md5_file = "https://zenodo.org/record/3957592/files/{}.md5?download=1".format(download_file_name)

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
    except EnvironmentError:
        sys.stderr.write("Warning: Unable to extract {}.\n".format(tar_file))

    # uncompress sequences
    bz2_file = os.path.join(folder, download_file_name + ".fna.bz2")
    fna_file = os.path.join(folder, download_file_name + ".fna")

    if not os.path.isfile(fna_file):
        sys.stderr.write('\n\nDecompressing {} into {}\n'.format(bz2_file, fna_file))

        with open(fna_file, 'wb') as fna_h, \
            bz2.BZ2File(bz2_file, 'rb') as bz2_h:
            for data in iter(lambda: bz2_h.read(100 * 1024), b''):
                fna_h.write(data)

    # build bowtie2 indexes
    if not glob(os.path.join(folder, download_file_name + "*.bt2")):
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

    for bt2 in glob(os.path.join(folder, download_file_name + "*.bt2")):
        os.chmod(bt2, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)  # change permissions to 664

    sys.stderr.write('Removing uncompress database {}\n'.format(fna_file))
    os.remove(fna_file)

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

def resolve_latest_database(bowtie2_db,mpa_latest_url, force=False):
    if os.path.exists(os.path.join(bowtie2_db,'mpa_latest')):
        ctime_latest_db = int(os.path.getctime(os.path.join(bowtie2_db,'mpa_latest')))
        if int(time.time()) - ctime_latest_db > 31536000:         #1 year in epoch
            os.rename(os.path.join(bowtie2_db,'mpa_latest'),os.path.join(bowtie2_db,'mpa_previous'))
            download(mpa_latest_url, os.path.join(bowtie2_db,'mpa_latest'), force=True)

    if not os.path.exists(os.path.join(bowtie2_db,'mpa_latest') or force):
        download(mpa_latest_url, os.path.join(bowtie2_db,'mpa_latest'))

    with open(os.path.join(bowtie2_db,'mpa_latest')) as mpa_latest:
        latest_db_version = [line.strip() for line in mpa_latest if not line.startswith('#')]
    
    return ''.join(latest_db_version)

def check_and_install_database(index, bowtie2_db, bowtie2_build, nproc, force_redownload_latest):
    # Create the folder if it does not already exist
    if not os.path.isdir(bowtie2_db):
        try:
            os.makedirs(bowtie2_db)
        except EnvironmentError:
            sys.exit("ERROR: Unable to create folder for database install: " + bowtie2_db)

    if index != 'latest' and len(glob(os.path.join(bowtie2_db, "*{}*".format(index)))) >= 6:
        return index

    list_file_path = os.path.join(bowtie2_db, "file_list.txt")
    #try downloading from Dropbox
    try:
        if not os.path.exists(list_file_path):
            download(DBX_FILE_LIST, list_file_path)

        if os.path.isfile(list_file_path):
            with open(list_file_path) as f:
                ls_f = dict( [row.strip().split() for row in f])
        use_zenodo = False
    except: #If fails, use zenodo
        ls_f = {'mpa_lates' : 'https://zenodo.org/record/3957592/files/mpa_latest?download=1' }
        use_zenodo = True

    """ Check if the database is installed, if not download and install """
    if index == 'latest':
        index = resolve_latest_database(bowtie2_db, ls_f['mpa_latest'], force_redownload_latest)

    if os.path.exists(os.path.join(bowtie2_db,'mpa_previous')):
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

    # download the tar archive and decompress
    sys.stderr.write("\nDownloading MetaPhlAn database\nPlease note due to "
                     "the size this might take a few minutes\n")
    download_unpack_tar(DBX_FILE_LIST, index, bowtie2_db, bowtie2_build, nproc, use_zenodo)
    sys.stderr.write("\nDownload complete\n")
    return index