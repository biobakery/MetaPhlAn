# MetaPhlAn: Metagenomic Phylogenetic Analysis
![Build MetaPhlAn package](https://github.com/biobakery/MetaPhlAn/workflows/Build%20MetaPhlAn%20package/badge.svg?branch=3.0)
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

MetaPhlAn relies on ~1.1M unique clade-specific marker genes (the latest marker information file `mpa_v30_CHOCOPhlAn_201901_marker_info.txt.bz2` can be found  [here](https://www.dropbox.com/sh/7qze7m7g9fe2xjg/AAAlyQITZuUCtBUJxpxhIroIa/mpa_v30_CHOCOPhlAn_201901_marker_info.txt.bz2?dl=1)) identified from ~100,000 reference genomes (~99,500 bacterial and archaeal and ~500 eukaryotic), allowing:

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

## Installation
The best way to install MetaPhlAn is through conda via the Bioconda channel. If you have not configured you Anaconda installation in order to fetch packages from Bioconda, please follow [these steps](https://bioconda.github.io/user/install.html#set-up-channels) in order to setup the channels.

You can install MetaPhlAn by running

```
$ conda install -c bioconda metaphlan
```

For installing it from the source code and for further installation instructions, please see the Wiki at the [Installation](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-3.0#installation) paragraph.

-------------

## MetaPhlAn2 and StrainPhlAn tutorials and resources

In addition to the information on this page, you can refer to the following additional resources.

* The [MetaPhlAn tutorial](https://github.com/biobakery/MetaPhlAn/wiki).

* The [StrainPhlAn tutorial](https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-3.0).

* The [MetaPhlAn](https://github.com/biobakery/biobakery/wiki/metaphlan2) and [StrainPhlAn](https://github.com/biobakery/biobakery/wiki/strainphlan2) wikis on bioBakery.

* The [MetaPhlAn](https://forum.biobakery.org/c/Microbial-community-profiling/MetaPhlAn/) and [StrainPhlAn](https://forum.biobakery.org/c/Microbial-community-profiling/StrainPhlAn/) Discourse forum.

* Related tools including [PanPhlAn](https://github.com/segatalab/panphlan) (and its [tutorial](https://github.com/segatalab/panphlan/wiki/Home)), [GraPhlAn](https://github.com/segatalab/graphlan) (and it [tutorial](https://github.com/biobakery/biobakery/wiki/graphlan)), [PhyloPhlAn2](https://github.com/biobakery/phylophlan) (and its [tutorial](https://github.com/biobakery/biobakery/wiki/phylophlan)), [HUMAnN2](https://github.com/biobakery/humann/) (and its [tutorial](https://github.com/biobakery/biobakery/wiki/humann2)).

* The related [bioBakery workflows](https://github.com/biobakery/biobakery/wiki/biobakery_workflows)
