from setuptools import setup, find_packages


with open('README.md', 'r') as f:
    long_description = f.read()

setup(name='metaphlan2',
    version='2.6.1',
    description='MetaPhlAn2 installer for QIIME 2 plugin',
    long_description=long_description,
    author='Francesco Asnicar',
    author_email='f.asnicar@unitn.it',
    license='',
    url='http://segatalab.cibio.unitn.it/tools/metaphlan2/',
    keywords='',
    packages=find_packages(),
    package_data={'metaphlan2': ['db_v20/*']},
    install_requires=['numpy', 'biom-format'],
    entry_points={'qiime2.plugins': ['metaphlan2=metaphlan2.plugin_setup:plugin']},
    zip_safe=False
)
