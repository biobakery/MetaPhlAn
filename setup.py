import setuptools
from setuptools.command.install import install
from io import open
import sys, shutil, os, zipfile, tarfile, subprocess, tempfile, re, time
from urllib.request import urlretrieve

install_requires = ['numpy', 'h5py', 'biom-format', 'biopython', 'pandas', 'scipy', 'requests', 'dendropy', 'pysam', 'phylophlan'],

if sys.version_info[0] < 3:
    sys.stdout.write('MetaPhlAn requires Python 3 or higher. Please update you Python installation')

setuptools.setup(
    name='MetaPhlAn',
    version='4.2.2',
    author='Claudia Mengoni',
    author_email='claudia.mengoni@unitn.it',
    url='http://github.com/biobakery/MetaPhlAn/',
    license='LICENSE.txt',
    packages=setuptools.find_packages(),
    package_data = { 'metaphlan' : [
        'metaphlan_databases/*.txt',
        'utils/*',
        'utils/treeshrink/*',
        'utils/treeshrink/scripts/*',
        'utils/treeshrink/R_scripts/*',
        'utils/treeshrink/Rlib/BMS/*',
        'utils/treeshrink/Rlib/BMS/*/*',
    ] },    
    entry_points={
        'console_scripts': [
            'metaphlan = metaphlan.metaphlan:main',
            'strainphlan = metaphlan.strainphlan:main',
            'add_metadata_tree.py = metaphlan.utils.add_metadata_tree:main',
            'extract_markers.py = metaphlan.utils.extract_markers:main',
            'merge_metaphlan_tables.py  = metaphlan.utils.merge_metaphlan_tables:main',
            'merge_vsc_tables.py  = metaphlan.utils.merge_vsc_tables:main',
            'plot_tree_graphlan.py = metaphlan.utils.plot_tree_graphlan:main',
            'read_fastx.py = metaphlan.utils.read_fastx:main',
            'sample2markers.py = metaphlan.utils.sample2markers:main',
            'strain_transmission.py = metaphlan.utils.strain_transmission:main',
            'sgb_to_gtdb_profile.py = metaphlan.utils.sgb_to_gtdb_profile:main',
            'metaphlan2krona.py = metaphlan.utils.metaphlan2krona:main',
            'run_treeshrink.py = metaphlan.utils.treeshrink.run_treeshrink:main',
            'treeshrink.py = metaphlan.utils.treeshrink.treeshrink:main',
            'create_toy_database.py = metaphlan.utils.create_toy_database:main',
            'fix_relab_mpa4.py = metaphlan.utils.fix_relab_mpa4:main',
        ]
    },
    description='MetaPhlAn is a computational tool for profiling the composition of microbial communities (Bacteria, Archaea and Eukaryotes) from metagenomic shotgun sequencing data (i.e. not 16S) with species-level. With the newly added StrainPhlAn module, it is now possible to perform accurate strain-level microbial profiling.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    install_requires=install_requires
)
