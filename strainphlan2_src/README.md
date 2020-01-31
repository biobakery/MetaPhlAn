# StrainPhlAn2: metagenomic strain-level population genomics

[TOC]

## Description

StrainPhlAn2 is a computational tool for tracking individual strains across large set of samples. **The input** of StrainPhlAn2 is a set of metagenomic samples and for each species, **the output** is a multiple sequence alignment (MSA) file of all species strains reconstructed directly from the samples. From this MSA, StrainPhlAn2 calls PhyloPhlAn2 (http://segatalab.cibio.unitn.it/tools/phylophlan/index.html) to build the phylogenetic tree showing the strain evolution of the sample strains. 
For each sample, StrainPhlAn2 extracts the strain of a specific species by merging and concatenating all reads mapped against that species markers in the MetaPhlAn2 database.

In detail, let us start from a toy example with 6 HMP gut metagenomic samples (SRS055982-subjectID\_638754422, SRS022137-subjectID\_638754422, SRS019161-subjectID\_763496533, SRS013951-subjectID\_763496533, SRS014613-subjectID\_763840445, SRS064276-subjectID\_763840445) from 3 three subjects (each was sampled at two time points) and one *Bacteroides caccae* genome G000273725. 
**We would like to**:

* extract the *Bacteroides caccae* strains from these samples and compare them with the reference genome in a phylogenetic tree.
* know how many snps between those strains and the reference genome.

Running StrainPhlAn2 on these samples, we will obtain the *Bacteroides caccae* phylogentic tree and its multiple sequence alignment in the following figure (produced with [ete2](http://etetoolkit.org/) and [Jalview](http://www.jalview.org/)):

![tree_alignment.png](https://bitbucket.org/repo/EgaXA7G/images/108375392-tree_alignment.png)

We can see that the strains from the same subject are grouped together. The tree also highlights that the strains from each subject did not evolv between the two sampling time points. From the tree, we also know that the strains of subject "763496533" is closer to the reference genome than those of the others. 
In addition, the table below shows the number of snps between the sample strains and the reference genome based on the strain alignment returned by StrainPhlAn2.

![svn_distance.PNG](https://bitbucket.org/repo/EgaXA7G/images/1682440489-svn_distance.PNG)

In the next sections, we will illustrate step by step how to run StrainPhlAn2 on this toy example to reproduce the above figures.

### Pre-requisites
StrainPhlAn2 requires *python 3* and the libraries [biopython](http://biopython.org/wiki/Main_Page) (tested on **version 1.73**), [msgpack](https://pypi.python.org/pypi/msgpack-python) (tested on **version 0.15.2**), [numpy](http://www.numpy.org/) (tested on **version 1.16.2**), [dendropy](https://pythonhosted.org/DendroPy/) (tested on **version 4.4.0**), [pandas](https://pandas.pydata.org/) (tested on **version 0.24.2**) and [pysam](https://pypi.org/project/pysam/) (tested on **version 0.15.2**). Besides, StrainPhlAn2 also needs the following programs in the executable path:

* [samtools](http://samtools.sourceforge.net/) (tested on **version 1.9**) which can be downloaded from [here](https://github.com/samtools) for processing the consensus markers. 
* [blast+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (tested on **version 2.9**) for adding reference genomes to the phylogenetic tree (blastn and makeblastdb commands)
* [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/) (tested on **version 8.2.4**) for building the phylogenetic trees.

If MetaPhlAn2 was installed through conda and you have the `mpa` activated, all the pre-requisites are satisfied.


The script files in folder "strainphlan_src2" should be changed to executable mode by:

```
#!python
chmod +x strainphlan2_src/*.py
```

and add to the executable path:

```
#!python
export PATH=$PATH:$(pwd -P)/strainphlan_src2
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

Step 3. Produce the consensus-marker files which are the input for StrainPhlAn2:

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

This step will extract the markers of *Bacteroides_caccae* in the database and then StrainPhlAn2 will identify the sequences in the reference genomes that are closet to them (in the next step by using blast). Those will be concatenated and referred as *reference-genome-reconstructed strains*. 
The commands to run are:

```
#!python
mkdir -p db_markers
python extract_markers.py -d ../metaphlan_databases/mpa_v294_CHOCOPhlAn_201901.pkl -c s__Bacteroides_caccae -o db_markers/
```

After this step, you should have one file in folder "db\_markers": "s\_\_Bacteroides\_caccae.markers.fna" containing the markers of *Bacteroides caccae*.
Those markers can be found in the folder "db\_markers" in [this link](https://www.dropbox.com/sh/rq3xqm12h2amn6q/AAD81ckPfT_metLDJcl5YJGNa?dl=0).


Step 5. Build the multiple sequence alignment and phylogenetic tree:

This step will filtered the selected clade markers based on their presence in the *sample-reconstructed strains* (stored in the marker files produced in step 3) and *reference-genomes* (if specified). Also the *sample-reconstructed strains* and *reference-genomes* will be filtered based on the presence of the selected clade markers. From this filtered markers and samples, StrainPhlAn2 will call PhyloPhlAn2 to produce a multiple sequence alignment (MSA) and to build the phylogenetic tree.

The commands to run are:

```
#!python
mkdir -p output
python strainphlan2.py -s consensus_markers/*.pkl -m db_markers/s__Bacteroides_caccae.fna -r reference_genomes/G000273725.fna.bz2 -o output -n 8 -c s__Bacteroides_caccae --phylophlan_mode accurate --mutation_rates
```

After this step, you will find the tree "output/RAxML_bestTree.s\_\_Bacteroides\_caccae.StrainPhlAn2.tre". All the output files can be found in the folder "output" in [this link](https://www.dropbox.com/sh/rq3xqm12h2amn6q/AAD81ckPfT_metLDJcl5YJGNa?dl=0).
You can view it by [Archaeopteryx](https://sites.google.com/site/cmzmasek/home/software/archaeopteryx) or any other viewers.


In order to add the metadata, we also provide a script called "add\_metadata\_tree.py" which can be used as follows:

```
#!python
python add_metadata_tree.py -t output/RAxML_bestTree.s__Bacteroides_caccae.StrainPhlAn2.tre -f fastqs/metadata.txt -m subjectID
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
python plot_tree_graphlan.py -t output/RAxML_bestTree.s__Bacteroides_caccae.StrainPhlAn2.tre.metadata -m subjectID
```

and obtain the following figure (output/RAxML_bestTree.s\_\_Bacteroides\_caccae.StrainPhlAn2.tre.metadata.png):
![RAxML_bestTree.s__Bacteroides_caccae.StrainPhlAn2.tre.metadata.png](https://bitbucket.org/repo/EgaXA7G/images/2601372061-RAxML_bestTree.s__Bacteroides_caccae.StrainPhlAn2.tre.metadata.png)

Note that this Script must be executed using Python2.

## Full command-line options


```
usage: strainphlan2.py [-h] [-d DATABASE] [-m CLADE_MARKERS]
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
                        keep a marker. Default 80
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
                        without building any tree.
```

### Some other useful output files
In the output folder, you can find the following files:

1. clade_name.info: this file shows the general information like the total length of the concatenated markers (full sequence length), number of used markers, etc.
2. clade_name.polymorphic: this file shows the statistics on the polymorphic site, where "sample" is the sample name, "percentage\_of\_polymorphic_sites" is the percentage of sites that are suspected to be polymorphic, "avg\_freq" is the average frequency of the dominant alleles on all polymorphic sites, "avg\_coverage" is the average coverage at all polymorphic sites.
3. clade_name.StrainPhlAn2_concatenated.aln: the alignment file of all metagenomic strains.
5. clade_name.mutation: this file contains a mutation rates summary table for the concatenated MSA. This file will be generated if *--mutation_rates* param is specified.
4. clade_name_mutation_rates: this folder contains the mutation rates table for each of the aligned markers. This folder will be generated if *--mutation_rates* param is specified.