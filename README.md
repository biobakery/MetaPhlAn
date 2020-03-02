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

