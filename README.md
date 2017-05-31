[TOC]

#**MetaPhlAn 2: Metagenomic Phylogenetic Analysis**#

AUTHORS: Duy Tin Truong (duytin.truong@unitn.it), Nicola Segata (nicola.segata@unitn.it)

##**Description**##
MetaPhlAn is a computational tool for profiling the composition of microbial communities (Bacteria, Archaea, Eukaryotes and Viruses) from metagenomic shotgun sequencing data with species level resolution. From version 2.0, MetaPhlAn is also able to identify specific strains (in the not-so-frequent cases in which the sample contains a previously sequenced strains) and to track strains across samples for all species.

MetaPhlAn 2 relies on ~1M unique clade-specific marker genes ([the marker information file can be found at src/utils/markers_info.txt.bz2 or here](https://bitbucket.org/biobakery/metaphlan2/src/473a41eba501df5f750da032d4f04b38db98dde1/utils/markers_info.txt.bz2?at=default)) identified from ~17,000 reference genomes (~13,500 bacterial and archaeal, ~3,500 viral, and ~110 eukaryotic), allowing:

* unambiguous taxonomic assignments;
* accurate estimation of organismal relative abundance;
* species-level resolution for bacteria, archaea, eukaryotes and viruses;
* strain identification and tracking
* orders of magnitude speedups compared to existing methods.
* metagenomic strain-level population genomics

If you use this software, please cite :

[**MetaPhlAn2 for enhanced metagenomic taxonomic profiling.**](http://www.nature.com/nmeth/journal/v12/n10/pdf/nmeth.3589.pdf)
 *Duy Tin Truong, Eric A Franzosa, Timothy L Tickle, Matthias Scholz, George Weingart, Edoardo Pasolli, Adrian Tett, Curtis Huttenhower & Nicola Segata*. 
Nature Methods 12, 902–903 (2015)

-------------

##**Pre-requisites**##

MetaPhlAn requires *python 2.7* or higher with argparse, tempfile and [*numpy*](http://www.numpy.org/) libraries installed 
  (apart for numpy they are usually installed together with the python distribution). 
  Python3 is also now supported.

**If you provide the SAM output of [BowTie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) as input, there are no additional prerequisite.**

* If you would like to use the BowTie2 integrated in MetaPhlAn, you need to have BowTie2 version 2.0.0 or higher and perl installed (bowtie2 needs to be in the system path with execute _and_ read permission)

* If you use the "utils/metaphlan_hclust_heatmap.py" script to plot and hierarchical cluster multiple MetaPhlAn-profiled samples you will also need the following python libraries: [matplotlib](http://matplotlib.org/index.html), [scipy](http://www.scipy.org/), [pylab](http://wiki.scipy.org/PyLab) (if not installed together with MatPlotLib).

* If you want to produce the output as "biom" file you also need [biom](http://biom-format.org/) installed

* MetaPhlAn is not tightly integrated with advanced heatmap plotting with [hclust2](https://bitbucket.org/nsegata/hclust2) and cladogram visualization with [GraPhlAn](https://bitbucket.org/nsegata/graphlan/wiki/Home). If you use such visualization tool please refer to their prerequisites. 

----------------------

##**Installation**##

MetaPhlAn 2.0 can be obtained by either

[Downloading MetaPhlAn v2.0](https://bitbucket.org/biobakery/metaphlan2/get/default.zip)  

**OR**

Cloning the repository via the following commands
``$ hg clone https://bitbucket.org/biobakery/metaphlan2``

--------------------------


##**Basic Usage**##

This section presents some basic usages of MetaPhlAn2, for more advanced usages, please see at [its wiki](https://bitbucket.org/biobakery/biobakery/wiki/metaphlan2).

We assume here that ``metaphlan2.py`` is in the system path and that ``mpa_dir`` bash variable contains the main MetaPhlAn folder. You can set this two variables moving to your MetaPhlAn2 local folder and type:
```
#!cmd
$ export PATH=`pwd`:$PATH
$ export mpa_dir=`pwd`
```

Here is the basic example to profile a metagenome from raw reads (requires BowTie2 in the system path with execution and read permissions, Perl installed). 

```
#!cmd
$ metaphlan2.py metagenome.fastq --input_type fastq > profiled_metagenome.txt
```

It is highly recommended to save the intermediate BowTie2 output for re-running MetaPhlAn extremely quickly (--bowtie2out), and use multiple CPUs (--nproc) if available:

```
#!cmd
$ metaphlan2.py metagenome.fastq --bowtie2out metagenome.bowtie2.bz2 --nproc 5 --input_type fastq > profiled_metagenome.txt
```

If you already mapped your metagenome against the marker DB (using a previous  MetaPhlAn run), you can obtain the results in few seconds by using the previously saved --bowtie2out file and specifying the input (--input_type bowtie2out):

```
#!cmd
$ metaphlan2.py metagenome.bowtie2.bz2 --nproc 5 --input_type bowtie2out > profiled_metagenome.txt
```

You can also provide an externally BowTie2-mapped SAM if you specify this format with --input_type. Two steps here: first map your metagenome with BowTie2 and then feed MetaPhlAn2 with the obtained sam:

```
#!cmd
$ bowtie2 --sam-no-hd --sam-no-sq --no-unal --very-sensitive -S metagenome.sam -x ${mpa_dir}/db_v20/mpa_v20_m200 -U metagenome.fastq
$ metaphlan2.py metagenome.sam --input_type sam > profiled_metagenome.txt
```

In order to make MetaPhlAn 2 easily compatible with complex metagenomic pipeline, there are now multiple alternative ways to pass the input:

```
#!cmd
$ cat metagenome.fastq | metaphlan2.py --input_type fastq > profiled_metagenome.txt
```

```
#!cmd
$ tar xjf metagenome.tar.bz2 --to-stdout | metaphlan2.py --input_type fastq --bowtie2db ${mpa_dir}/db_v20/mpa_v20_m200 > profiled_metagenome.txt
```

```
#!cmd
$ metaphlan2.py --input_type fastq < metagenome.fastq > profiled_metagenome.txt
```

```
#!cmd
$ metaphlan2.py --input_type fastq <(bzcat metagenome.fastq.bz2) > profiled_metagenome.txt
```

```
#!cmd
$ metaphlan2.py --input_type fastq <(zcat metagenome_1.fastq.gz metagenome_2.fastq.gz) > profiled_metagenome.txt
```

MetaPhlAn 2 can also natively **handle paired-end metagenomes** (but does not use the paired-end information), and, more generally, metagenomes stored in multiple files (but you need to specify the --bowtie2out parameter):

```
#!cmd
$ metaphlan2.py metagenome_1.fastq,metagenome_2.fastq --bowtie2out metagenome.bowtie2.bz2 --nproc 5 --input_type fastq > profiled_metagenome.txt
```

For advanced options and other analysis types (such as strain tracking) please refer to the full command-line options.

