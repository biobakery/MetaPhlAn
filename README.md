[TOC]

# MetaPhlAn2: Metagenomic Phylogenetic Analysis

## What's new

* New MetaPhlAn2 marker genes extracted with a newer version of ChocoPhlAn based on UniRef
* Calculation of metagenome size for improved estimation of reads mapped to a given clade
* Inclusion of NCBI taxonomy ID in the ouput file
* 



## Description
MetaPhlAn2 is a computational tool for profiling the composition of microbial communities (Bacteria, Archaea and Eukaryotes) from metagenomic shotgun sequencing data (i.e. not 16S) with species-level. With the newly added StrainPhlAn module, it is now possible to perform accurate strain-level microbial profiling.

MetaPhlAn2 relies on ~1M unique clade-specific marker genes (the latest marker information file `mpa_v25_CHOCOPhlAn_201901_marker_info.txt.bz2` can be found in the Download page [here](https://bitbucket.org/biobakery/metaphlan2/downloads/mpa_v25_CHOCOPhlAn_201901_marker_info.txt.bz2)) identified from ~100,000 reference genomes (~99,500 bacterial and archaeal and ~500 eukaryotic), allowing:

* unambiguous taxonomic assignments;
* accurate estimation of organismal relative abundance;
* species-level resolution for bacteria, archaea, eukaryotes and viruses;
* strain identification and tracking
* orders of magnitude speedups compared to existing methods.
* metagenomic strain-level population genomics

If you use MetaPhlAn version 1, please cite:

[**Metagenomic microbial community profiling using unique clade-specific marker genes.**](https://www.nature.com/articles/nmeth.2066) *Nicola Segata, Levi Waldron, Annalisa Ballarini, Vagheesh Narasimhan, Olivier Jousson, &  Curtis Huttenhower*. Nature Methods 9, 811-814 (2012)

If you use MetaPhlAn2, please cite:

[**MetaPhlAn2 for enhanced metagenomic taxonomic profiling.**](http://www.nature.com/nmeth/journal/v12/n10/pdf/nmeth.3589.pdf) *Duy Tin Truong, Eric A Franzosa, Timothy L Tickle, Matthias Scholz, George Weingart, Edoardo Pasolli, Adrian Tett, Curtis Huttenhower & Nicola Segata*. Nature Methods 12, 902-903 (2015)

If you use StrainPhlAn, please cite the MetaPhlAn2 paper and the following StrainPhlAn paper:

[**Microbial strain-level population structure and genetic diversity from metagenomes.**](http://genome.cshlp.org/content/27/4/626.full.pdf) *Duy Tin Truong, Adrian Tett, Edoardo Pasolli, Curtis Huttenhower, & Nicola Segata*. Genome Research 27:626-638 (2017)

-------------

## MetaPhlAn2 and StrainPhlAn tutorials and resources

In addition to the information on this page, you can refer to the following additional resources.

* The [MetaPhlAn2 tutorial on bioBakery](https://bitbucket.org/biobakery/biobakery/wiki/metaphlan2).

* The [StrainPhlAn tutorial on bioBakery](https://bitbucket.org/biobakery/biobakery/wiki/strainphlan).

* The MetaPhlAn2 and StrainPhlAn [Google Group](http://groups.google.com/forum/#!forum/metaphlan-users) ([metaphlan-users@googlegroups.com](mailto:metaphlan-users@googlegroups.com))

* Related tools including [PanPhlAn](https://bitbucket.org/CibioCM/panphlan/src) (and its [tutorial](https://bitbucket.org/CibioCM/panphlan/wiki/Home)), [GraPhlAn](https://bitbucket.org/nsegata/graphlan/wiki/Home) (and it [tutorial](https://bitbucket.org/biobakery/biobakery/wiki/graphlan)), [PhyloPhlAn2](https://bitbucket.org/nsegata/phylophlan/wiki/Home) (and its [tutorial](https://bitbucket.org/biobakery/biobakery/wiki/phylophlan)), [HUMAnN2](https://bitbucket.org/biobakery/humann2/wiki/Home) (and its [tutorial](https://bitbucket.org/biobakery/biobakery/wiki/humann2)).

* The related [bioBakery workflows](https://bitbucket.org/biobakery/biobakery/wiki/biobakery_workflows)

-------------

## Pre-requisites

MetaPhlAn2 requires *python 3* with argparse, tempfile, [numpy](http://www.numpy.org/), and [Biopython](https://biopython.org/) libraries installed 
(apart for numpy and Biopython, the others are usually installed together with the python distribution). 

MetaPhlAn2 requires the `read_fastx.py` script to be present in the system path, if not found MetaPhlAn2 will try to locate it in the folder containing the `metaphlan2.py` 
script under `utils/read_fastx.py`.
In case you moved the `metaphlan2.py` script, please export the `read_fastx.py` script in your PATH bash variable.

**If you provide the SAM output of [BowTie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) as input, there are no additional prerequisite.**

* If you would like to use the BowTie2 integrated in MetaPhlAn2, you need to have BowTie2 version 2.0.0 or higher and perl installed (bowtie2 needs to be in the system path with execute _and_ read permission)

* If you want to produce the output as "biom" file you also need [biom](http://biom-format.org/) installed

* MetaPhlAn2 is integrated with advanced heatmap plotting with [hclust2](https://bitbucket.org/nsegata/hclust2) and cladogram visualization with [GraPhlAn](https://bitbucket.org/nsegata/graphlan/wiki/Home). If you use such visualization tool please refer to their prerequisites. 

## Installation
MetaPhlAn2 2.9 can be obtained by 

**direct download** from [Bitbucket](https://bitbucket.org/biobakery/metaphlan2/get/default.zip)


or **cloning the repository** using the following command ``$ hg clone https://bitbucket.org/biobakery/metaphlan2``

<!-- Through **Bioconda**

``$ conda install metaphlan2``

Through **Docker**

``$ docker pull segatalab/metaphlan2``  -->

MetaPhlAn2 needs the clade markers and the database to be downloaded locally.  

```
#!bash
$ metaphlan2.py --install 
```

By default, the latest MetaPhlAn2 database is downloaded and built. You can download a specific version with the `--index` parameter

```
#!bash
$ metaphlan2.py --install --index v25_CHOCOPhlAn_201901
```

--------------------------

## Basic Usage

This section presents some basic usages of MetaPhlAn2, for more advanced usages, please see at [its wiki](https://bitbucket.org/biobakery/biobakery/wiki/metaphlan2).

We assume here that ``metaphlan2.py`` is in the system path and that ``mpa_dir`` bash variable contains the main MetaPhlAn folder. You can set this two variables moving to your MetaPhlAn2 local folder and type:

```
#!bash
$ export PATH=`pwd`:$PATH
$ export mpa_dir=`pwd`
```

Here is the basic example to profile a metagenome from raw reads (requires BowTie2 in the system path with execution and read permissions, Perl installed). 

```
#!bash
$ metaphlan2.py metagenome.fastq --input_type fastq -o profiled_metagenome.txt
```

It is highly recommended to save the intermediate BowTie2 output for re-running MetaPhlAn2 extremely quickly (`--bowtie2out`), and use multiple CPUs (`--nproc`) if available:

```
#!bash
$ metaphlan2.py metagenome.fastq --bowtie2out metagenome.bowtie2.bz2 --nproc 5 --input_type fastq -o profiled_metagenome.txt
```


If you already mapped your metagenome against the marker DB (using a previous MetaPhlAn2 run), you can obtain the results in few seconds by using the previously saved `--bowtie2out` file and specifying the input (`--input_type bowtie2out`):

```
#!bash
$ metaphlan2.py metagenome.bowtie2.bz2 --nproc 5 --input_type bowtie2out -o profiled_metagenome.txt
```

`bowtie2out` files generated with MetaPhlAn2 versions below 2.9 are not compatibile. Starting from MetaPhlAn2 2.9, the BowTie2 ouput now includes the size of the profiled metagenome.
If you want to re-run MetaPhlAn2 using these file you should provide the metagenome size via `--nreads`

```
#!bash
$ metaphlan2.py metagenome.bowtie2.bz2 --nproc 5 --input_type bowtie2out --nreads 520000 -o profiled_metagenome.txt
```

You can also provide an externally BowTie2-mapped SAM if you specify this format with `--input_type`. Two steps here: first map your metagenome with BowTie2 and then feed MetaPhlAn2 with the obtained sam:

```
#!bash
$ bowtie2 --sam-no-hd --sam-no-sq --no-unal --very-sensitive -S metagenome.sam -x metaphlan_databases/mpa_v25_CHOCOPhlAn_201901 -U metagenome.fastq
$ metaphlan2.py metagenome.sam --input_type sam -o profiled_metagenome.txt
```

MetaPhlAn 2 can also natively **handle paired-end metagenomes** (but does not use the paired-end information), and, more generally, metagenomes stored in multiple files (but you need to specify the --bowtie2out parameter):

```
#!bash
$ metaphlan2.py metagenome_1.fastq,metagenome_2.fastq --bowtie2out metagenome.bowtie2.bz2 --nproc 5 --input_type fastq -o profiled_metagenome.txt
```

You can provide the specific database version with `--index`. By default MetaPhlAn2 is run with `--index latest`: the latest version of the database is used; if it is not available, MetaPhlAn2 will try to download from BitBucket.


For advanced options and other analysis types (such as strain tracking) please refer to the full command-line options.

## Full command-line options
  INCLUDE COMMAND LINE


<!-- ## Docker image

MetaPhlAn2 Docker image is based on [Biobox](http://bioboxes.org/about/). For the execution, MetaPhlAn2 Docker image requires a YAML file with the inputs defined as follows:

```
#!yaml
version: 1.0.0
arguments:
    - fastq:
        - type: fastq
          value: /bbx/mnt/input/fastq1.fastq.gz
        - type: fastq
          value: /bbx/mnt/input/fastq2.fastq
```

Each FASTQ sample defined in the ```fastq``` section is a single metagenomic sample.

The YAML file **MUST** be named ```biobox.yaml``` and saved in the input folder (e.g. ```input```) 

To run the docker image, execute the following command:


```
#!bash
docker run --volume="${PWD}/input:/bbx/mnt/input:ro" --volume="${PWD}/output:/bbx/mnt/output:rw" segatalab/metaphlan2
```

where:

- ```input``` is the folder with the input files and *biobox.yaml*. ```input``` **MUST** be mounted to ```/bbx/mnt/input```

- ```output``` is the output folder in which MetaPhlAn2 output is saved. Like ```input```, ```output``` **MUST** be mounted to ```/bbx/mnt/output```

More folders can be mounted by specifing the mount point via the ```--volume``` parameter. Sym-linking the input files in the input folder will result in a error since Docker will not be able to reach the link source. -->


## Utility Scripts

MetaPhlAn's repository features a few utility scripts to aid in manipulation of sample output and its visualization. These scripts can be found under the ``utils`` folder in the metaphlan2 directory.

### Merging Tables

The script **merge\_metaphlan\_tables.py** allows to combine MetaPhlAn2 output from several samples to be merged into one table Bugs (rows) vs Samples (columns) with the table enlisting the relative normalized abundances per sample per bug.

To merge multiple output files, run the script as below

```
#!bash
$ python utils/merge_metaphlan_tables.py metaphlan_output1.txt metaphlan_output2.txt metaphlan_output3.txt output/merged_abundance_table.txt
```

Wildcards can be used as needed:

```
#!bash
$ python utils/merge_metaphlan_tables.py metaphlan_output*.txt  output/merged_abundance_table.txt
```

**There is no limit to how many files you can merge.**



## Heatmap Visualization

Hclust2
The script **metaphlan\_hclust\_heatmap.py** allows to visualize the MetaPhlAn results in the form of a hierarchically-clustered heatmap. To generate the heatmap for a merged MetaPhlAn output table (as described above), please run the script as below.

### GraPhlAn Visualization

The tutorial of using GraPhlAn can be found from [the MetaPhlAn2 wiki](https://bitbucket.org/biobakery/biobakery/wiki/metaphlan2).


## Customizing the database
In order to add a marker to the database, the user needs the following steps:

* Reconstruct the marker sequences (in fasta format) from the MetaPhlAn2 bowtie2 database by:

```
#!bash
bowtie2-inspect metaphlan2/databases/mpa_v20_m200 > metaphlan2/markers.fasta
```

* Add the marker sequence stored in a file new_marker.fasta to the marker set:

```
#!bash
cat new_marker.fasta >> metaphlan2/markers.fasta
```

* Rebuild the bowtie2 database:

```
#!bash
mkdir metaphlan2/db_v21/mpa_v21_m200
bowtie2-build metaphlan2/markers.fasta metaphlan2/db_v21/mpa_v21_m200
```

* Assume that the new marker was extracted from genome1, genome2. Update the taxonomy file from python console as follows:

```
#!python

import cPickle as pickle
import bz2

db = pickle.load(bz2.BZ2File('databases/mpa_v20_m200.pkl', 'r'))

# Add the taxonomy of the new genomes
db['taxonomy']['taxonomy of genome1'] = length of genome1
db['taxonomy']['taxonomy of genome2'] = length of genome2

# Add the information of the new marker as the other markers
db['markers'][new_marker_name] = {
                                   'clade': the clade that the marker belongs to,
                                   'ext': {the name of the first external genome where the marker appears, 
                                           the name of the second external genome where the marker appears, 
                                          },
                                   'len': length of the marker,
                                   'score': score of the marker,
                                   'taxon': the taxon of the marker}
# To see an example, try to print the first marker information:
# print db['markers'].items()[0]

# Save the new mpa_pkl file
ofile = bz2.BZ2File('metaphlan2/db_v21/mpa_v21_m200.pkl', 'w')
pickle.dump(db, ofile, pickle.HIGHEST_PROTOCOL)
ofile.close()
```

* To use the new database, switch to metaphlan2/db_v21 instead of metaphlan2/db\_v20 when running metaphlan2.py with option "--mpa\_pkl".


## Metagenomic strain-level population genomics

StrainPhlAn is a computational tool for tracking individual strains across large set of samples. **The input** of StrainPhlAn is a set of metagenomic samples and for each species, **the output** is a multiple sequence alignment (MSA) file of all species strains reconstructed directly from the samples. From this MSA, StrainPhlAn calls RAxML (or other phylogenetic tree builders) to build the phylogenetic tree showing the strain evolution of the sample strains. 
For each sample, StrainPhlAn extracts the strain of a specific species by merging and concatenating all reads mapped against that species markers in the MetaPhlAn2 database.

In detail, let us start from a toy example with 6 HMP gut metagenomic samples (SRS055982-subjectID\_638754422, SRS022137-subjectID\_638754422, SRS019161-subjectID\_763496533, SRS013951-subjectID\_763496533, SRS014613-subjectID\_763840445, SRS064276-subjectID\_763840445) from 3 three subjects (each was sampled at two time points) and one *Bacteroides caccae* genome G000273725. 
**We would like to**:

* extract the *Bacteroides caccae* strains from these samples and compare them with the reference genome in a phylogenetic tree.
* know how many snps between those strains and the reference genome.

Running StrainPhlAn on these samples, we will obtain the *Bacteroides caccae* phylogentic tree and its multiple sequence alignment in the following figure (produced with [ete2](http://etetoolkit.org/) and [Jalview](http://www.jalview.org/)):

![tree_alignment.png](https://bitbucket.org/repo/rM969K/images/476974413-tree_alignment.png)

We can see that the strains from the same subject are grouped together. The tree also highlights that the strains from subject "763840445" (red color) do not change between two sampling time points whereas the strains from the other subjects have slightly evolved. From the tree, we also know that the strains of subject "763496533" is closer to the reference genome than those of the others. 
In addition, the table below shows the number of snps between the sample strains and the reference genome based on the strain alignment returned by MetaPhlAn\_Strainer.

![snp_distance.png](https://bitbucket.org/repo/rM969K/images/1771497600-snp_distance.png)

In the next sections, we will illustrate step by step how to run MetaPhlAn\_Strainer on this toy example to reproduce the above figures.

### Pre-requisites
StrainPhlAn requires *python 2.7* and the libraries [pysam](http://pysam.readthedocs.org/en/latest/) (tested on **version 0.8.3**), [biopython](http://biopython.org/wiki/Main_Page), [msgpack](https://pypi.python.org/pypi/msgpack-python), [pandas](https://pandas.pydata.org) (tested on **version 0.22**), [numpy](http://www.numpy.org/) (tested on **version 1.14.2**) and [scipy](https://www.scipy.org) (tested on **version 1.0.0**), [dendropy](https://pythonhosted.org/DendroPy/) (tested on version **3.12.0**). Besides, StrainPhlAn also needs the following programs in the executable path:

* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for mapping reads against the marker database.
* [MUSCLE](http://www.drive5.com/muscle/) for the alignment step.
* [samtools, bcftools and vcfutils.pl](http://samtools.sourceforge.net/) which can be downloaded from [here](https://github.com/samtools) for building consensus markers. Note that vcfutils.pl is included in bcftools and **StrainPhlAn only works with samtools version 0.1.19** as samtools has changed the output format after this version.
* [blast+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) for adding reference genomes to the phylogenetic tree (blastn and makeblastdb commands)
* [raxmlHPC and raxmlHPC-PTHREADS-SSE3](http://sco.h-its.org/exelixis/web/software/raxml/index.html) for building the phylogenetic trees.

All dependence binaries on Linux 64 bit can be downloaded in the folder "bin" from [this link](https://www.dropbox.com/sh/m4na8wefp53j8ej/AABA3yVsG26TbB0t1cnBS9-Ra?dl=0).

The script files in folder "strainphlan_src" should be changed to executable mode by:

```
#!python
chmod +x strainphlan_src/*.py
```

and add to the executable path:

```
#!python
export PATH=$PATH:$(pwd -P)/strainphlan_src
```

### Usage

Let's reproduce the toy example result in the introduction section. Note that all the commands to run the below steps are in the "strainer_tutorial/step?*.sh" files (? corresponding to the step number). All the below steps are excuted under the "strainer_tutorial" folder.
The steps are as follows:

Step 1. Download 6 HMP gut metagenomic samples, the metadata.txt file and one reference genome from the folder "fastqs" and "reference_genomes" in [this link](https://www.dropbox.com/sh/m4na8wefp53j8ej/AABA3yVsG26TbB0t1cnBS9-Ra?dl=0) and put these folders under the "strainer_tutorial" folder.

Step 2. Obtain the sam files from these samples by mapping them against MetaPhlAn2 database:

This step will run MetaPhlAn2 to map all metagenomic samples against the MetaPhlAn2 marker database and produce the sam files (\*.sam.bz2).
Each sam file (in SAM format) corresponding to each sample contains the reads mapped against the marker database of MetaPhlAn2.
The commands to run are:

```
#!python
mkdir -p sams
for f in $(ls fastqs/*.bz2)
do
    echo "Running metaphlan2 on ${f}"
    bn=$(basename ${f} | cut -d . -f 1)
    tar xjfO ${f} | ../metaphlan2.py --bowtie2db ../databases/mpa_v20_m200 --mpa_pkl ../databases/mpa_v20_m200.pkl --input_type multifastq --nproc 10 -s sams/${bn}.sam.bz2 --bowtie2out sams/${bn}.bowtie2_out.bz2 -o sams/${bn}.profile
done
```

After this step, you will have a folder "sams" containing the sam files (\*.sam.bz2) and other MetaPhlAn2 output files. 
This step will take around 270 minutes. If you want to skip this step, you can download the sam files from the folder "sams" in [this link](https://www.dropbox.com/sh/m4na8wefp53j8ej/AABA3yVsG26TbB0t1cnBS9-Ra?dl=0).

Step 3. Produce the consensus-marker files which are the input for StrainPhlAn:

For each sample, this step will reconstruct all species strains found in it and store them in a marker file (\*.markers). Those strains are referred as *sample-reconstructed strains*. Additional details in generating consensus sequences can be found [here](http://samtools.sourceforge.net/mpileup.shtml).
The commands to run are:


```
#!python
mkdir -p consensus_markers
cwd=$(pwd -P)
export PATH=${cwd}/../strainphlan_src:${PATH}
python ../strainphlan_src/sample2markers.py --ifn_samples sams/*.sam.bz2 --input_type sam --output_dir consensus_markers --nprocs 10 &> consensus_markers/log.txt
```

The result is the same if you want run several sample2markers.py scripts in parallel with each run for a sample (this maybe useful for some cluster-system settings).
After this step, you will have a folder "consensus_markers" containing all sample-marker files (\*.markers).
This steps will take around 44 minutes.  If you want to skip this step, you can download the consensus marker files from the folder "consensus_markers" in [this link](https://www.dropbox.com/sh/m4na8wefp53j8ej/AABA3yVsG26TbB0t1cnBS9-Ra?dl=0).

Step 4. Extract the markers of *Bacteroides\_caccae* from MetaPhlAn2 database (to add its reference genome later):

This step will extract the markers of *Bacteroides_caccae* in the database and then StrainPhlAn will identify the sequences in the reference genomes that are closet to them (in the next step by using blast). Those will be concatenated and referred as *reference-genome-reconstructed strains*. 
The commands to run are:

```
#!python
mkdir -p db_markers
bowtie2-inspect ../databases/mpa_v20_m200 > db_markers/all_markers.fasta
python ../strainphlan_src/extract_markers.py --mpa_pkl ../databases/mpa_v20_m200.pkl --ifn_markers db_markers/all_markers.fasta --clade s__Bacteroides_caccae --ofn_markers db_markers/s__Bacteroides_caccae.markers.fasta
```

Note that the "all\_markers.fasta" file consists can be reused for extracting other reference genomes. 
After this step, you should have two files in folder "db\_markers": "all\_markers.fasta" containing all marker sequences, and "s\_\_Bacteroides\_caccae.markers.fasta" containing the markers of *Bacteroides caccae*.
This step will take around 1 minute and can skipped if you do not need to add the reference genomes to the phylogenetic tree. Those markers can be found in the folder "db\_markers" in [this link](https://www.dropbox.com/sh/m4na8wefp53j8ej/AABA3yVsG26TbB0t1cnBS9-Ra?dl=0).

Before building the trees, we should get the list of all clades detected from the samples and save them in the "output/clades.txt" file by the following command:

```
#!python
python ../strainphlan.py --mpa_pkl ../databases/mpa_v20_m200.pkl --ifn_samples consensus_markers/*.markers --output_dir output --nprocs_main 10 --print_clades_only > output/clades.txt
```

The clade names in the output file "clades.txt" will be used for the next step.

Step 5. Build the multiple sequence alignment and phylogenetic tree:

This step will align and clean the *sample-reconstructed strains* (stored in the marker files produced in step 3) and *reference-genome-reconstructed strains* (extracted based on the database markers in step 4) to produce a multiple sequence alignment (MSA) and store it in the file "clade_name.fasta". From this MSA file, StrainPhlAn will call RAxML to build the phylogenetic tree.
Note that: all marker files (\*.markers) **must be used together** as the input for the strainphlan.py script because StrainPhlAn needs to align all of the strains at once.

The commands to run are:

```
#!python
mkdir -p output
python ../strainphlan.py --mpa_pkl ../databases/mpa_v20_m200.pkl --ifn_samples consensus_markers/*.markers --ifn_markers db_markers/s__Bacteroides_caccae.markers.fasta --ifn_ref_genomes reference_genomes/G000273725.fna.bz2 --output_dir output --nprocs_main 10 --clades s__Bacteroides_caccae &> output/log_full.txt
```

This step will take around 2 minutes. After this step, you will find the tree "output/RAxML\_bestTree.s\_\_Bacteroides\_caccae.tree". All the output files can be found in the folder "output" in [this link](https://www.dropbox.com/sh/m4na8wefp53j8ej/AABA3yVsG26TbB0t1cnBS9-Ra?dl=0).
You can view it by [Archaeopteryx](https://sites.google.com/site/cmzmasek/home/software/archaeopteryx) or any other viewers.

By default, if you do not specify reference genomes (by --ifn\_ref\_genomes) and any specific clade (by --clades), strainphlan.py will build the phylogenetic trees for all species that it can detect.

In order to add the metadata, we also provide a script called "add\_metadata\_tree.py" which can be used as follows:

```
#!python
python ../strainphlan_src/add_metadata_tree.py --ifn_trees output/RAxML_bestTree.s__Bacteroides_caccae.tree --ifn_metadatas fastqs/metadata.txt --metadatas subjectID
```

The script "add\_metadata\_tree.py" can accept multiple metadata files (space separated, wild card can also be used) and multiple trees. A metadata file is a tab separated file where the first row is the meta-headers, and the following rows contain the metadata for each sample. Multiple metadata files are used in the case where your samples come from more than one dataset and you do not want to merge the metadata files.
For more details of using "add\_metadata\_tree.py", please see its help (with option "-h").
An example of a metadata file is the "fastqs/metadata.txt" file with the below content:

```
#!python
sampleID        subjectID
SRS055982       638754422
SRS022137       638754422
SRS019161       763496533
SRS013951       763496533
SRS014613       763840445
SRS064276       763840445
G000273725  ReferenceGenomes
```

Note that "sampleID" is a compulsory field. 

After adding the metadata, you will obtain the tree files "*.tree.metadata" with metadata and view them by [Archaeopteryx](https://sites.google.com/site/cmzmasek/home/software/archaeopteryx) as in the previous step.

If you have installed [graphlan](https://bitbucket.org/nsegata/graphlan/wiki/Home), you can plot the tree with the command:


```
#!python
python ../strainphlan_src/plot_tree_graphlan.py --ifn_tree output/RAxML_bestTree.s__Bacteroides_caccae.tree.metadata --colorized_metadata subjectID
```

and obtain the following figure (output/RAxML\_bestTree.s\_\_Bacteroides\_caccae.tree.metadata.png):

![RAxML_bestTree.s__Bacteroides_caccae.tree.metadata.png](https://bitbucket.org/repo/rM969K/images/1574126761-RAxML_bestTree.s__Bacteroides_caccae.tree.metadata.png)

Step 6. If you want to remove the samples with high-probability of containing multiple strains, you can rebuild the tree by removing the multiple strains:

```
#!python
python ../strainphlan_src/build_tree_single_strain.py --ifn_alignments output/s__Bacteroides_caccae.fasta --nprocs 10 --log_ofn output/build_tree_single_strain.log
python ../strainphlan_src/add_metadata_tree.py --ifn_trees output/RAxML_bestTree.s__Bacteroides_caccae.remove_multiple_strains.tree --ifn_metadatas fastqs/metadata.txt --metadatas subjectID
```

You will obtain the refined tree "output/RAxML\_bestTree.s\_\_Bacteroides\_caccae.remove\_multiple\_strains.tree.metadata". This tree can be found in the folder "output" in [this link](https://www.dropbox.com/sh/m4na8wefp53j8ej/AABA3yVsG26TbB0t1cnBS9-Ra?dl=0).

### Some useful options
All option details can be viewed by strainphlan.py help:

```
#!python
python ../strainphlan.py -h
```

The default setting can be stringent for some cases where you have very few samples left in the phylogenetic tree. You can relax some parameters to add more samples back:

1. *marker\_in\_clade*: In each sample, the clades with the percentage of present markers less than this threshold are removed. Default "0.8". You can set this parameter to "0.5" to add some more samples.
2. *sample\_in\_marker*: If the percentage of samples that a marker present in is less than this threhold, that marker is removed. Default "0.8". You can set this parameter to "0.5" to add some more samples.
3. *N\_in\_marker*: The consensus markers with the percentage of N nucleotides greater than this threshold are removed. Default "0.2". You can set this parameter to "0.5" to add some more samples.
4. *gap\_in\_sample*: The samples with full sequences concatenated from all markers and having the percentage of gaps greater than this threshold will be removed. Default 0.2. You can set this parameter to "0.5" to add some more samples.
5. *relaxed\_parameters*: use this option to automatically set the above parameters to add some more samples by accepting some more gaps, Ns, etc. This option is equivalent to set: marker\_in\_clade=0.5, sample\_in\_marker=0.5, N\_in\_marker=0.5, gap\_in\_sample=0.5. Default "False".
6. *relaxed\_parameters2*: use this option to add more samples by accepting some noise. This is equivalent to set marker\_in\_clade=0.2, sample\_in\_marker=0.2, N\_in\_marker=0.8, gap\_in\_sample=0.8. Default "False".

### Some other useful output files
In the output folder, you can find the following files:

1. clade_name.fasta: the alignment file of all metagenomic strains.
3. *.marker_pos: this file shows the starting position of each marker in the strains.
3. *.info: this file shows the general information like the total length of the concatenated markers (full sequence length), number of used markers, etc.
4. *.polymorphic: this file shows the statistics on the polymorphic site, where "sample" is the sample name, "percentage\_of\_polymorphic_sites" is the percentage of sites that are suspected to be polymorphic, "avg\_freq" is the average frequency of the dominant alleles on all polymorphic sites, "avg\_coverage" is the average coverage at all polymorphic sites.
