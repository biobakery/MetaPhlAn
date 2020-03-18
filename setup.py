import setuptools
from setuptools.command.install import install
from io import open
import sys, shutil, os, zipfile, tarfile, subprocess, tempfile, re, time
from urllib.request import urlretrieve

install_requires = ['numpy', 'h5py', 'biom-format', 'biopython', 'pandas', 'scipy', 'requests', 'dendropy', 'pysam']
setup_requires = ['numpy', 'cython']

if sys.version_info[0] < 3:
    sys.stdout.write('MetaPhlAn requires Python 3 or higher. Please update you Python installation')

def byte_to_megabyte(byte):
    """
    Convert byte value to megabyte
    """
    
    return byte / (1024.0**2)

class ReportHook():
    def __init__(self):
        self.start_time=time.time()
        
    def report(self, blocknum, block_size, total_size):
        """
        Print download progress message
        """
        
        if blocknum == 0:
            self.start_time=time.time()
            if total_size > 0:
                print("Downloading file of size: " + "{:.2f}".format(byte_to_megabyte(total_size)) + " MB\n")
        else:
            total_downloaded=blocknum*block_size
            status = "{:3.2f} MB ".format(byte_to_megabyte(total_downloaded))
                    
            if total_size > 0:
                percent_downloaded=total_downloaded * 100.0 / total_size
                # use carriage return plus sys.stdout to overwrite stdout
                try:
                    download_rate=total_downloaded/(time.time()-self.start_time)
                    estimated_time=(total_size-total_downloaded)/download_rate
                except ZeroDivisionError:
                    download_rate=0
                    estimated_time=0
                estimated_minutes=int(estimated_time/60.0)
                estimated_seconds=estimated_time-estimated_minutes*60.0
                status +="{:3.2f}".format(percent_downloaded) + " %  " + \
                    "{:5.2f}".format(byte_to_megabyte(download_rate)) + " MB/sec " + \
                    "{:2.0f}".format(estimated_minutes) + " min " + \
                    "{:2.0f}".format(estimated_seconds) + " sec "
            status+="        \r"
            sys.stdout.write(status)

def download(url, download_file):
    """
    Download a file from a url
    """

    try:
        print("Downloading "+url)
        file, headers = urlretrieve(url,download_file,reporthook=ReportHook().report)
        # print final return to start new line of stdout
        print("\n")
    except EnvironmentError:
        print("WARNING: Unable to download "+url)
    
def download_unpack_zip(url,download_file_name,folder,software_name):
    """
    Download the url to the file and decompress into the folder
    """
    
    # Check for write permission to the target folder
    if not os.access(folder, os.W_OK):
        print("WARNING: The directory is not writeable: "+
            folder + " . Please modify the permissions.")
    
    download_file=os.path.join(folder, download_file_name)
    
    download(url, download_file)
    
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
       
def download_unpack_tar(url,download_file_name,folder,software_name):
    """
    Download the url to the file and decompress into the folder
    """
    
    # Check for write permission to the target folder
    if not os.access(folder, os.W_OK):
        print("WARNING: The directory is not writeable: "+
            folder + " . Please modify the permissions.")
    
    download_file=os.path.join(folder, download_file_name)
    
    download(url, download_file)
    
    error_during_extract=False
    
    try:
        tarfile_handle=tarfile.open(download_file)
        tarfile_handle.extractall(path=folder)
        tarfile_handle.close()
    except (EnvironmentError,tarfile.ReadError):
        print("WARNING: Unable to extract "+software_name+".")
        error_during_extract=True
        
    if not error_during_extract:
        try:
            os.unlink(download_file)
        except EnvironmentError:
            print("WARNING: Unable to remove the temp download: " + download_file)

class Install(install):
    def run(self):  
        print('Installing cmseq')
        cmseq_repo = "https://github.com/SegataLab/cmseq/archive/master.zip"
        cmseq_path = os.path.dirname(os.path.abspath(__file__))
        # cmseq_path = os.path.join(self.install_base, "lib", "python{}.{}".format(sys.version_info[0], sys.version_info[1]), "site-packages", "metaphlan","utils", 'cmseq')
        if os.path.exists(cmseq_path):
            os.rmdir(cmseq_path)
        else:
            download_unpack_zip(cmseq_repo, "cmseq.zip", 'work/metaphlan/utils/cmseq', "cmseq")
        install.run(self)

setuptools.setup(
    name='MetaPhlAn',
    version='3.0',
    author='Francesco Beghini',
    author_email='francesco.beghini@unitn.it',
    url='http://github.com/biobakery/MetaPhlAn/',
    license='LICENSE.txt',
    packages=setuptools.find_packages(),
    package_data = { 'MetaPhlAn' : [
        'metaphlan/utils/*'
    ]},
    cmdclass={
        'biocondainstall': Install
    },
    entry_points={
        'console_scripts': [
            'metaphlan = metaphlan.metaphlan:main',
            'strainphlan = metaphlan.strainphlan:main',

            'add_metadata_tree.py = metaphlan.utils.add_metadata_tree:main',
            'extract_markers.py = metaphlan.utils.extract_markers:main',
            'merge_metaphlan_tables.py  = metaphlan.utils.merge_metaphlan_tables:main',
            'plot_tree_graphlan.py = metaphlan.utils.plot_tree_graphlan:main',
            'read_fastx.py = metaphlan.utils.read_fastx:main',
            'sample2markers.py = metaphlan.utils.sample2markers:main',
            'strain_transmission.py = metaphlan.utils.strain_transmission:main',
        ]
    },
    description='MetaPhlAn is a computational tool for profiling the composition of microbial communities (Bacteria, Archaea and Eukaryotes) from metagenomic shotgun sequencing data (i.e. not 16S) with species-level. With the newly added StrainPhlAn module, it is now possible to perform accurate strain-level microbial profiling.',
    long_description=open('README.md').read(),
    install_requires=install_requires,
    setup_requires = setup_requires
)