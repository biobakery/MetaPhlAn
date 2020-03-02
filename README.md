# MetaPhlAn: Metagenomic Phylogenetic Analysis

## What's new in version 3
* New MetaPhlAn marker genes extracted with a newer version of ChocoPhlAn based on UniRef
* Estimation of metagenome composed by unknown microbes with parameter `--unknown_estimation`
* Automatic retrieval and installation of the latest MetaPhlAn database  with parameter `--index latest`
* Virus profiling with `--add_viruses`
* Calculation of metagenome size for improved estimation of reads mapped to a given clade
* Inclusion of NCBI taxonomy ID in the ouput file
* CAMI (Taxonomic) Profiling Output Format included
* Removal of reads with low MAPQ values
-------------

## Description
MetaPhlAn is a computational tool for profiling the composition of microbial communities (Bacteria, Archaea and Eukaryotes) from metagenomic shotgun sequencing data (i.e. not 16S) with species-level. With the newly added StrainPhlAn module, it is now possible to perform accurate strain-level microbial profiling.

MetaPhlAn relies on ~1.1M unique clade-specific marker genes (the latest marker information file `mpa_v296_CHOCOPhlAn_201901_marker_info.txt.bz2` can be found  [here](https://www.dropbox.com/sh/7qze7m7g9fe2xjg/AAAv_ShZiz7pNTaT_YONJTF7a/mpa_v296_CHOCOPhlAn_201901_marker_info.txt.bz2?dl=1)) identified from ~100,000 reference genomes (~99,500 bacterial and archaeal and ~500 eukaryotic), allowing:

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

* The [MetaPhlAn tutorial on bioBakery](https://github.com/biobakery/biobakery/wiki/metaphlan2).

* The [StrainPhlAn tutorial on bioBakery](https://github.com/biobakery/biobakery/wiki/strainphlan).

* The [MetaPhlAn](https://forum.biobakery.org/c/Microbial-community-profiling/MetaPhlAn/) and [StrainPhlAn](https://forum.biobakery.org/c/Microbial-community-profiling/StrainPhlAn/) Discourse forum.

* Related tools including [PanPhlAn](https://github.com/segatalab/panphlan) (and its [tutorial](https://github.com/segatalab/panphlan/wiki/Home)), [GraPhlAn](https://github.com/segatalab/graphlan) (and it [tutorial](https://github.com/biobakery/biobakery/wiki/graphlan)), [PhyloPhlAn2](https://github.com/biobakery/phylophlan) (and its [tutorial](https://github.com/biobakery/biobakery/wiki/phylophlan)), [HUMAnN2](https://github.com/biobakery/humann/) (and its [tutorial](https://github.com/biobakery/biobakery/wiki/humann2)).

* The related [bioBakery workflows](https://github.com/biobakery/biobakery/wiki/biobakery_workflows)

-------------

## Pre-requisites

MetaPhlAn requires *python 3* or newer with [numpy](http://www.numpy.org/), and [Biopython](https://biopython.org/) libraries installed.

**If you provide the SAM output of [BowTie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) as input, there are no additional prerequisite.**

* If you would like to use the BowTie2 integrated in MetaPhlAn, you need to have BowTie2 version 2.0.0 or higher and perl installed (bowtie2 needs to be in the system path with execute _and_ read permission)

* If you want to produce the output as "biom" file you also need [biom](http://biom-format.org/) installed

* MetaPhlAn is integrated with advanced heatmap plotting with [hclust2](https://github.com/segatalab/hclust2) and cladogram visualization with [GraPhlAn](https://github.com/segatalab/graphlan). If you use such visualization tool please refer to their prerequisites. 

## Installation
The best way to install MetaPhlAn is through conda:

```
$ conda install -c bioconda metaphlan2
```

It is recommended to create an isolated conda environment and install MetaPhlAn2 into it. 

```
$ conda create --name mpa2 -c bioconda metaphlan2
```

This allow to have the correct version of all the dependencies isolated from the system's python installation.

Before using MetaPhlAn, you should activate the `mpa` environment:

```
$ conda activate mpa2
```

You can also install and run MetaPhlAn2 through Docker

```
$ docker pull quay.io/biocontainers/metaphlan2:2.9.21
```

Alternatively, you can **manually download** from [GitBub](https://github.com/biobakery/MetaPhlAn/archive/3.0_alpha.zip) or **clone the repository** using the following command ``$ git clone https://github.com/biobakery/MetaPhlAn.git``.

If you choose this way, **you'll need to install manually all the dependencies!**


MetaPhlAn needs the clade markers and the database to be downloaded locally. To obtain them:

```
$ metaphlan.py --install 
```

By default, the latest MetaPhlAn2 database is downloaded and built. You can download a specific version with the `--index` parameter

```
$ metaphlan.py --install --index mpa_v296_CHOCOPhlAn_201901
```

--------------------------

## Basic Usage

This section presents some basic usages of MetaPhlAn, for more advanced usages, please see at [its wiki](https://github.com/biobakery/biobakery/wiki/metaphlan2).

We assume here that ``metaphlan.py`` is in the system path. Here is the basic example to profile a metagenome from raw reads (requires BowTie2 in the system path with execution and read permissions, Perl installed). 

```
$ metaphlan.py metagenome.fastq --input_type fastq -o profiled_metagenome.txt
```

### Starting from version 3, MetaPhlAn estimates the fraction of the metagenome composed by microbes that are unknown. The relative abundance profile is scaled according the percentage of reads mapping to a known clade. 

### To obtain the old behavior in version 3 you can set the flag "--no-unknown-estimation"

It is highly recommended to save the intermediate BowTie2 output for re-running MetaPhlAn extremely quickly (`--bowtie2out`), and use multiple CPUs (`--nproc`) if available:

```
$ metaphlan.py metagenome.fastq --bowtie2out metagenome.bowtie2.bz2 --nproc 5 --input_type fastq -o profiled_metagenome.txt
```


If you already mapped your metagenome against the marker DB (using a previous MetaPhlAn run), you can obtain the results in few seconds by using the previously saved `--bowtie2out` file and specifying the input (`--input_type bowtie2out`):

```
$ metaphlan.py metagenome.bowtie2.bz2 --nproc 5 --input_type bowtie2out -o profiled_metagenome.txt
```

`bowtie2out` files generated with MetaPhlAn versions below 3.0 are not compatibile. Starting from MetaPhlAn 3.0, the BowTie2 ouput now includes the size of the profiled metagenome.
If you want to re-run MetaPhlAn using these file you should provide the metagenome size via `--nreads`

```
$ metaphlan.py metagenome.bowtie2.bz2 --nproc 5 --input_type bowtie2out --nreads 520000 -o profiled_metagenome.txt
```

You can also provide an externally BowTie2-mapped SAM if you specify this format with `--input_type`. Two steps here: first map your metagenome with BowTie2 and then feed MetaPhlAn with the obtained sam:

```
$ bowtie2 --sam-no-hd --sam-no-sq --no-unal --very-sensitive -S metagenome.sam -x metaphlan_databases/mpa_v295_CHOCOPhlAn_201901 -U metagenome.fastq
$ metaphlan.py metagenome.sam --input_type sam -o profiled_metagenome.txt
```

MetaPhlAn can also natively **handle paired-end metagenomes** (but does not use the paired-end information), and, more generally, metagenomes stored in multiple files (but you need to specify the --bowtie2out parameter):

```
$ metaphlan.py metagenome_1.fastq,metagenome_2.fastq --bowtie2out metagenome.bowtie2.bz2 --nproc 5 --input_type fastq -o profiled_metagenome.txt
```

You can provide the specific database version with `--index`. 

By default MetaPhlAn is run with `--index latest`: the latest version of the database is used; if it is not available, MetaPhlAn will try to download from BitBucket.


For advanced options and other analysis types (such as strain tracking) please refer to the full command-line options.

## Full command-line options

## Utility Scripts

MetaPhlAn's repository features a few utility scripts to aid in manipulation of sample output and its visualization. These scripts can be found under the ``utils`` folder in the MetaPhlAn directory.

### Merging Tables

The script **merge\_metaphlan\_tables.py** allows to combine MetaPhlAn output from several samples to be merged into one table Bugs (rows) vs Samples (columns) with the table enlisting the relative normalized abundances per sample per bug.

To merge multiple output files, run the script as below

```
$ python utils/merge_metaphlan_tables.py metaphlan_output1.txt metaphlan_output2.txt metaphlan_output3.txt output/merged_abundance_table.txt
```

Wildcards can be used as needed:

```
$ python utils/merge_metaphlan_tables.py metaphlan_output*.txt  output/merged_abundance_table.txt
```

**Output files can be merged only if the profiling was performed with the same version of the MetaPhlAn database.**

**There is no limit to how many files you can merge.**



### Heatmap Visualization

The hclust2 script generates a hierarchically-clustered heatmap from MetaPhlAn abundance profiles. To generate the heatmap for a merged MetaPhlAn output table (as described above), you need to run the script as below:

```
hclust2.py \
  -i HMP.species.txt \
  -o HMP.sqrt_scale.png \
  --skip_rows 1 \
  --ftop 50 \
  --f_dist_f correlation \
  --s_dist_f braycurtis \
  --cell_aspect_ratio 9 \
  -s --fperc 99 \
  --flabel_size 4 \
  --metadata_rows 2,3,4 \
  --legend_file HMP.sqrt_scale.legend.png \
  --max_flabel_len 100 \
  --metadata_height 0.075 \
  --minv 0.01 \
  --no_slabels \
  --dpi 300 \
  --slinkage complete
```

### GraPhlAn Visualization

The tutorial of using GraPhlAn can be found from [the MetaPhlAn2 wiki](https://github.com/biobakery/biobakery/wiki/metaphlan2).


## Customizing the database
In order to add a marker to the database, the user needs the following steps:

* Reconstruct the marker sequences (in fasta format) from the MetaPhlAn2 BowTie2 database by:

```
bowtie2-inspect metaphlan_databases/mpa_v296_CHOCOPhlAn_201901 > metaphlan_databases/mpa_v296_CHOCOPhlAn_201901_markers.fasta
```

* Add the marker sequence stored in a file new_marker.fasta to the marker set:

```
cat new_marker.fasta >> metaphlan_databases/mpa_v296_CHOCOPhlAn_201901_markers.fasta
```

* Rebuild the bowtie2 database:

```
bowtie2-build metaphlan_databases/mpa_v296_CHOCOPhlAn_201901_markers.fasta metaphlan_databases/mpa_v26_CHOCOPhlAn_NEW
```

* Assume that the new marker was extracted from genome1, genome2. Update the taxonomy file from the Python console as follows:

```
import pickle
import bz2

db = pickle.load(bz2.open('metaphlan_databases/mpa_v296_CHOCOPhlAn_201901.pkl', 'r'))

# Add the taxonomy of the new genomes
db['taxonomy']['taxonomy of genome1'] = ('NCBI taxonomy id of genome1', length of genome1)
db['taxonomy']['taxonomy of genome2'] = ('NCBI taxonomy id of genome1', length of genome2)

# Add the information of the new marker as the other markers
db['markers'][new_marker_name] = {
                                   'clade': the clade that the marker belongs to,
                                   'ext': {the GCA of the first external genome where the marker appears,
                                           the GCA of the second external genome where the marker appears,
                                          },
                                   'len': length of the marker,
                                   'taxon': the taxon of the marker
                                }
                                   
# To see an example, try to print the first marker information:
# print db['markers'].items()[0]

# Save the new mpa_pkl file
with bz2.BZ2File('metaphlan_databases/mpa_v296_CHOCOPhlAn_NEW.pkl', 'w') as ofile:
    pickle.dump(db, ofile, pickle.HIGHEST_PROTOCOL)

```

* To use the new database, run metaphlan.py with option "--index mpa_v296_CHOCOPhlAn_NEW".


# StrainPhlAn 3.0: metagenomic strain-level population genomics

[TOC]

## Description

StrainPhlAn 3.0 is a computational tool for tracking individual strains across large set of samples. **The input** of StrainPhlAn 3.0 is a set of metagenomic samples and for each species, **the output** is a multiple sequence alignment (MSA) file of all species strains reconstructed directly from the samples. From this MSA, StrainPhlAn 3.0 calls PhyloPhlAn2 (http://segatalab.cibio.unitn.it/tools/phylophlan/index.html) to build the phylogenetic tree showing the strain evolution of the sample strains. 
For each sample, StrainPhlAn 3.0 extracts the strain of a specific species by merging and concatenating all reads mapped against that species markers in the MetaPhlAn2 database.

In detail, let us start from a toy example with 6 HMP gut metagenomic samples (SRS055982-subjectID\_638754422, SRS022137-subjectID\_638754422, SRS019161-subjectID\_763496533, SRS013951-subjectID\_763496533, SRS014613-subjectID\_763840445, SRS064276-subjectID\_763840445) from 3 three subjects (each was sampled at two time points) and one *Bacteroides caccae* genome G000273725. 
**We would like to**:

* extract the *Bacteroides caccae* strains from these samples and compare them with the reference genome in a phylogenetic tree.
* know how many snps between those strains and the reference genome.

Running StrainPhlAn 3.0 on these samples, we will obtain the *Bacteroides caccae* phylogentic tree and its multiple sequence alignment in the following figure (produced with [ete2](http://etetoolkit.org/) and [Jalview](http://www.jalview.org/)):

![tree_alignment.png](https://bitbucket.org/repo/EgaXA7G/images/108375392-tree_alignment.png)

We can see that the strains from the same subject are grouped together. The tree also highlights that the strains from each subject did not evolv between the two sampling time points. From the tree, we also know that the strains of subject "763496533" is closer to the reference genome than those of the others. 
In addition, the table below shows the number of snps between the sample strains and the reference genome based on the strain alignment returned by StrainPhlAn 3.0.

![svn_distance.PNG](https://bitbucket.org/repo/EgaXA7G/images/1682440489-svn_distance.PNG)

In the next sections, we will illustrate step by step how to run StrainPhlAn 3.0 on this toy example to reproduce the above figures.

### Pre-requisites
StrainPhlAn 3.0 requires *python 3* and the libraries [biopython](http://biopython.org/wiki/Main_Page) (tested on **version 1.73**), [msgpack](https://pypi.python.org/pypi/msgpack-python) (tested on **version 0.15.2**), [numpy](http://www.numpy.org/) (tested on **version 1.16.2**), [dendropy](https://pythonhosted.org/DendroPy/) (tested on **version 4.4.0**), [pandas](https://pandas.pydata.org/) (tested on **version 0.24.2**) and [pysam](https://pypi.org/project/pysam/) (tested on **version 0.15.2**). Besides, StrainPhlAn 3.0 also needs the following programs in the executable path:

* [samtools](http://samtools.sourceforge.net/) (tested on **version 1.9**) which can be downloaded from [here](https://github.com/samtools) for processing the consensus markers. 
* [blast+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (tested on **version 2.9**) for adding reference genomes to the phylogenetic tree (blastn and makeblastdb commands)
* [MAFFT](http://mafft.cbrc.jp/alignment/software/) (tested on **version 7.271**) for performing the MSA.
* [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/) (tested on **version 8.2.4**) for building the phylogenetic trees.

If MetaPhlAn2 was installed through conda and you have the `mpa` activated, all the pre-requisites are satisfied.


The script files in folder "strainphlan" should be changed to executable mode by:

```
#!python
chmod +x strainphlan/*.py
```

and add to the executable path:

```
#!python
export PATH=$PATH:$(pwd -P)/strainphlan
```

### Usage

Let's reproduce the toy example result in the introduction section. The steps are as follows:

Step 1. Download 6 HMP gut metagenomic samples, the metadata.txt file and one reference genome from the folder "fastqs" and "reference_genomes" in [this link](https://www.dropbox.com/sh/rq3xqm12h2amn6q/AAD81ckPfT_metLDJcl5YJGNa?dl=0) and put these folders under the "strainer_tutorial" folder.

Step 2. Obtain the sam files from these samples by mapping them against MetaPhlAn2 database:

This step will run MetaPhlAn2 to map all metagenomic samples against the MetaPhlAn2 marker database and produce the sam files (\*.sam.bz2).
Each sam file (in SAM format) corresponding to each sample contains the reads mapped against the marker database of MetaPhlAn2.
The commands to run are:

```
#!python
mkdir -p sams
mkdir -p bowtie2
mkdir -p profiles
for f in fastq/SRS*
do
    echo "Running metaphlan2 on ${f}"
    bn=$(basename ${f})
    python metaphlan2.py ${f} --index mpa_v294_CHOCOPhlAn_201901 --input_type multifastq -s sams/${bn}.sam.bz2 --bowtie2out bowtie2/${bn}.bowtie2.bz2 -o profiles/profiled_${bn}.txt 
done
```

After this step, you will have a folder "sams" containing the sam files (\*.sam.bz2) and other MetaPhlAn2 output files in the "bowtie2" and "profiles" folders. 
If you want to skip this step, you can download the sam files from the folder "sams" in [this link](https://www.dropbox.com/sh/rq3xqm12h2amn6q/AAD81ckPfT_metLDJcl5YJGNa?dl=0).

Step 3. Produce the consensus-marker files which are the input for StrainPhlAn 3.0:

For each sample, this step will reconstruct all species strains found in it and store them in a pickle file (\*.pkl). Those strains are referred as *sample-reconstructed strains*. 
The commands to run are:


```
#!python
mkdir -p consensus_markers
python sample_to_markers.py -i sams/*.sam.bz2 -o consensus_markers -n 8
```

The result is the same if you want run several sample_to_markers.py scripts in parallel with each run for a sample (this maybe useful for some cluster-system settings).
After this step, you will have a folder "consensus_markers" containing all sample-marker files (\*.pkl).
If you want to skip this step, you can download the consensus marker files from the folder "consensus_markers" in [this link](https://www.dropbox.com/sh/rq3xqm12h2amn6q/AAD81ckPfT_metLDJcl5YJGNa?dl=0).

Step 4. Extract the markers of *Bacteroides\_caccae* from MetaPhlAn2 database (to add its reference genome later):

This step will extract the markers of *Bacteroides_caccae* in the database and then StrainPhlAn 3.0 will identify the sequences in the reference genomes that are closet to them (in the next step by using blast). Those will be concatenated and referred as *reference-genome-reconstructed strains*. 
The commands to run are:

```
#!python
mkdir -p db_markers
python extract_markers.py -d ../metaphlan_databases/mpa_v294_CHOCOPhlAn_201901.pkl -c s__Bacteroides_caccae -o db_markers/
```

After this step, you should have one file in folder "db\_markers": "s\_\_Bacteroides\_caccae.markers.fna" containing the markers of *Bacteroides caccae*.
Those markers can be found in the folder "db\_markers" in [this link](https://www.dropbox.com/sh/rq3xqm12h2amn6q/AAD81ckPfT_metLDJcl5YJGNa?dl=0).


Step 5. Build the multiple sequence alignment and phylogenetic tree:

This step will filtered the selected clade markers based on their presence in the *sample-reconstructed strains* (stored in the marker files produced in step 3) and *reference-genomes* (if specified). Also the *sample-reconstructed strains* and *reference-genomes* will be filtered based on the presence of the selected clade markers. From this filtered markers and samples, StrainPhlAn 3.0 will call PhyloPhlAn2 to produce a multiple sequence alignment (MSA) and to build the phylogenetic tree.

The commands to run are:

```
#!python
mkdir -p output
python strainphlan.py -s consensus_markers/*.pkl -m db_markers/s__Bacteroides_caccae.fna -r reference_genomes/G000273725.fna.bz2 -o output -n 8 -c s__Bacteroides_caccae --phylophlan_mode accurate --mutation_rates
```

After this step, you will find the tree "output/RAxML_bestTree.s\_\_Bacteroides\_caccae.StrainPhlAn3.tre". All the output files can be found in the folder "output" in [this link](https://www.dropbox.com/sh/rq3xqm12h2amn6q/AAD81ckPfT_metLDJcl5YJGNa?dl=0).
You can view it by [Archaeopteryx](https://sites.google.com/site/cmzmasek/home/software/archaeopteryx) or any other viewers.


In order to add the metadata, we also provide a script called "add\_metadata\_tree.py" which can be used as follows:

```
#!python
python add_metadata_tree.py -t output/RAxML_bestTree.s__Bacteroides_caccae.StrainPhlAn3.tre -f fastqs/metadata.txt -m subjectID --string_to_remove .fastq.bz2
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

After adding the metadata, you will obtain the tree files "*.tre.metadata" with metadata and view them by [Archaeopteryx](https://sites.google.com/site/cmzmasek/home/software/archaeopteryx) as in the previous step.

If you have installed [graphlan](https://bitbucket.org/nsegata/graphlan/wiki/Home), you can plot the tree with the command:


```
#!python
python plot_tree_graphlan.py -t output/RAxML_bestTree.s__Bacteroides_caccae.StrainPhlAn3.tre.metadata -m subjectID
```

and obtain the following figure (output/RAxML_bestTree.s\_\_Bacteroides\_caccae.StrainPhlAn3.tre.metadata.png):
![RAxML_bestTree.s__Bacteroides_caccae.StrainPhlAn3.tre.metadata.png](https://bitbucket.org/repo/EgaXA7G/images/2601372061-RAxML_bestTree.s__Bacteroides_caccae.StrainPhlAn3.tre.metadata.png)

Note that this Script must be executed using Python2.

## Full command-line options


```
usage: strainphlan.py [-h] [-d DATABASE] [-m CLADE_MARKERS]
                       [-s SAMPLES [SAMPLES ...]]
                       [-r REFERENCES [REFERENCES ...]] [-c CLADE]
                       [-o OUTPUT_DIR] [-n NPROCS]
                       [--secondary_samples SECONDARY_SAMPLES [SECONDARY_SAMPLES ...]]
                       [--secondary_references SECONDARY_REFERENCES [SECONDARY_REFERENCES ...]]
                       [--trim_sequences TRIM_SEQUENCES]
                       [--marker_in_n_samples MARKER_IN_N_SAMPLES]
                       [--sample_with_n_markers SAMPLE_WITH_N_MARKERS]
                       [--secondary_sample_with_n_markers SECONDARY_SAMPLE_WITH_N_MARKERS]
                       [--phylophlan_configuration PHYLOPHLAN_CONFIGURATION]
                       [--mutation_rates]

optional arguments:
  -h, --help            show this help message and exit
  -d DATABASE, --database DATABASE
                        The input MetaPhlAn2.9 database
  -m CLADE_MARKERS, --clade_markers CLADE_MARKERS
                        The clade markers as FASTA file
  -s SAMPLES [SAMPLES ...], --samples SAMPLES [SAMPLES ...]
                        The reconstructed markers for each sample
  -r REFERENCES [REFERENCES ...], --references REFERENCES [REFERENCES ...]
                        The reference genomes
  -c CLADE, --clade CLADE
                        The clade to investigate
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        The output directory
  -n NPROCS, --nprocs NPROCS
                        The number of threads to use
  --secondary_samples SECONDARY_SAMPLES [SECONDARY_SAMPLES ...]
                        The reconstructed markers for each secondary sample
  --secondary_references SECONDARY_REFERENCES [SECONDARY_REFERENCES ...]
                        The secondary reference genomes
  --trim_sequences TRIM_SEQUENCES
                        The number of bases to remove from both ends when trimming 
                        markers. Default 50
  --marker_in_n_samples MARKER_IN_N_SAMPLES
                        Theshold defining the minimum percentage of samples to
                        keep a marker. Default 80 (%)
  --sample_with_n_markers SAMPLE_WITH_N_MARKERS
                        Threshold defining the minimun number of markers to
                        keep a sample. Default 20
  --secondary_sample_with_n_markers SECONDARY_SAMPLE_WITH_N_MARKERS
                        Threshold defining the minimun number of markers to
                        keep a secondary sample. Default 20
  --phylophlan_configuration PHYLOPHLAN_CONFIGURATION
                        The PhyloPhlAn configuration file
  --phylophlan_mode PHYLOPHLAN_MODE
                        The precision of the phylogenetic analysis {fast,
                        normal [default], accurate}
  --mutation_rates      If specified will produced a mutation rates table for
                        each of the aligned markers and a summary table for
                        the concatenated MSA. This operation can take long
                        time to finish  
  --print_clades_only   If specified only print the potential clades and stop
                        without building any tree
```

### Some other useful output files
In the output folder, you can find the following files:

1. clade_name.info: this file shows the general information like the total length of the concatenated markers (full sequence length), number of used markers, etc.
2. clade_name.polymorphic: this file shows the statistics on the polymorphic site, where "sample" is the sample name, "percentage\_of\_polymorphic_sites" is the percentage of sites that are suspected to be polymorphic, "avg\_freq" is the average frequency of the dominant alleles on all polymorphic sites, "avg\_coverage" is the average coverage at all polymorphic sites.
3. clade_name.StrainPhlAn3_concatenated.aln: the alignment file of all metagenomic strains.
5. clade_name.mutation: this file contains a mutation rates summary table for the concatenated MSA. This file will be generated if *--mutation_rates* param is specified.
4. clade_name_mutation_rates: this folder contains the mutation rates table for each of the aligned markers. This folder will be generated if *--mutation_rates* param is specified.
