[TOC]

# MetaPhlAn2: Metagenomic Phylogenetic Analysis

## What's new

* Automatic retrieval and installation of the latest MetaPhlAn2 database
* New MetaPhlAn2 marker genes extracted with a newer version of ChocoPhlAn based on UniRef
* Calculation of metagenome size for improved estimation of reads mapped to a given clade
* Inclusion of NCBI taxonomy ID in the ouput file

-------------

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

MetaPhlAn2 requires *python 2.7* or newer with argparse, tempfile, [numpy](http://www.numpy.org/), and [Biopython](https://biopython.org/) libraries installed 
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


MetaPhlAn2 needs the clade markers and the database to be downloaded locally. To obtain them:

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

You can provide the specific database version with `--index`. 

By default MetaPhlAn2 is run with `--index latest`: the latest version of the database is used; if it is not available, MetaPhlAn2 will try to download from BitBucket.


For advanced options and other analysis types (such as strain tracking) please refer to the full command-line options.

## Full command-line options
```
usage: metaphlan2.py --input_type
                     {fastq,fasta,multifasta,multifastq,bowtie2out,sam}
                     [--mpa_pkl MPA_PKL] [--bowtie2db METAPHLAN_BOWTIE2_DB]
                     [-x INDEX] [--bt2_ps BowTie2 presets]
                     [--bowtie2_exe BOWTIE2_EXE]
                     [--bowtie2_build BOWTIE2_BUILD] [--bowtie2out FILE_NAME]
                     [--no_map] [--tmp_dir] [--tax_lev TAXONOMIC_LEVEL]
                     [--min_cu_len] [--min_alignment_len]
                     [--ignore_eukaryotes] [--ignore_bacteria]
                     [--ignore_archaea] [--stat_q] [--perc_nonzero]
                     [--ignore_markers IGNORE_MARKERS] [--avoid_disqm]
                     [--stat] [-t ANALYSIS TYPE] [--nreads NUMBER_OF_READS]
                     [--pres_th PRESENCE_THRESHOLD] [--clade] [--min_ab]
                     [-o output file] [--sample_id_key name]
                     [--sample_id value] [-s sam_output_file]
                     [--legacy-output] [--no-unknown-estimation]
                     [--biom biom_output] [--mdelim mdelim] [--nproc N]
                     [--install] [--force_download]
                     [--read_min_len READ_MIN_LEN] [-v] [-h]
                     [INPUT_FILE] [OUTPUT_FILE]

DESCRIPTION
 MetaPhlAn version 2.9 (15 May 2019): 
 METAgenomic PHyLogenetic ANalysis for metagenomic taxonomic profiling.

AUTHORS: Nicola Segata (nicola.segata@unitn.it), Duy Tin Truong, Francesco Asnicar (f.asnicar@unitn.it), Francesco Beghini (francesco.beghini@unitn.it)

COMMON COMMANDS

 We assume here that metaphlan2.py is in the system path and that mpa_dir bash variable contains the
 main MetaPhlAn folder. Also BowTie2 should be in the system path with execution and read
 permissions, and Perl should be installed)

========== MetaPhlAn 2 clade-abundance estimation ================= 

The basic usage of MetaPhlAn 2 consists in the identification of the clades (from phyla to species and 
strains in particular cases) present in the metagenome obtained from a microbiome sample and their 
relative abundance. This correspond to the default analysis type (-t rel_ab).

*  Profiling a metagenome from raw reads:
$ metaphlan2.py metagenome.fastq --input_type fastq -o profiled_metagenome.txt

*  You can take advantage of multiple CPUs and save the intermediate BowTie2 output for re-running
   MetaPhlAn extremely quickly:
$ metaphlan2.py metagenome.fastq --bowtie2out metagenome.bowtie2.bz2 --nproc 5 --input_type fastq -o profiled_metagenome.txt

*  If you already mapped your metagenome against the marker DB (using a previous MetaPhlAn run), you
   can obtain the results in few seconds by using the previously saved --bowtie2out file and 
   specifying the input (--input_type bowtie2out):
$ metaphlan2.py metagenome.bowtie2.bz2 --nproc 5 --input_type bowtie2out -o profiled_metagenome.txt

*  bowtie2out files generated with MetaPhlAn2 versions below 2.9 are not compatibile.
   Starting from MetaPhlAn2 2.9, the BowTie2 ouput now includes the size of the profiled metagenome.
   If you want to re-run MetaPhlAn2 using these file you should provide the metagenome size via --nreads:
$ metaphlan2.py metagenome.bowtie2.bz2 --nproc 5 --input_type bowtie2out --nreads 520000 -o profiled_metagenome.txt

*  You can also provide an externally BowTie2-mapped SAM if you specify this format with 
   --input_type. Two steps: first apply BowTie2 and then feed MetaPhlAn2 with the obtained sam:
$ bowtie2 --sam-no-hd --sam-no-sq --no-unal --very-sensitive -S metagenome.sam -x ${mpa_dir}/metaphlan_databases/mpa_v25_CHOCOPhlAn_201901 -U metagenome.fastq
$ metaphlan2.py metagenome.sam --input_type sam -o profiled_metagenome.txt

*  We can also natively handle paired-end metagenomes, and, more generally, metagenomes stored in 
  multiple files (but you need to specify the --bowtie2out parameter):
$ metaphlan2.py metagenome_1.fastq,metagenome_2.fastq --bowtie2out metagenome.bowtie2.bz2 --nproc 5 --input_type fastq

------------------------------------------------------------------- 
 

========== Marker level analysis ============================ 

MetaPhlAn 2 introduces the capability of charachterizing organisms at the strain level using non
aggregated marker information. Such capability comes with several slightly different flavours and 
are a way to perform strain tracking and comparison across multiple samples.
Usually, MetaPhlAn 2 is first ran with the default -t to profile the species present in
the community, and then a strain-level profiling can be performed to zoom-in into specific species
of interest. This operation can be performed quickly as it exploits the --bowtie2out intermediate 
file saved during the execution of the default analysis type.

*  The following command will output the abundance of each marker with a RPK (reads per kilo-base) 
   higher 0.0. (we are assuming that metagenome_outfmt.bz2 has been generated before as 
   shown above).
$ metaphlan2.py -t marker_ab_table metagenome_outfmt.bz2 --input_type bowtie2out -o marker_abundance_table.txt
   The obtained RPK can be optionally normalized by the total number of reads in the metagenome 
   to guarantee fair comparisons of abundances across samples. The number of reads in the metagenome
   needs to be passed with the '--nreads' argument

*  The list of markers present in the sample can be obtained with '-t marker_pres_table'
$ metaphlan2.py -t marker_pres_table metagenome_outfmt.bz2 --input_type bowtie2out -o marker_abundance_table.txt
   The --pres_th argument (default 1.0) set the minimum RPK value to consider a marker present

*  The list '-t clade_profiles' analysis type reports the same information of '-t marker_ab_table'
   but the markers are reported on a clade-by-clade basis.
$ metaphlan2.py -t clade_profiles metagenome_outfmt.bz2 --input_type bowtie2out -o marker_abundance_table.txt

*  Finally, to obtain all markers present for a specific clade and all its subclades, the 
   '-t clade_specific_strain_tracker' should be used. For example, the following command
   is reporting the presence/absence of the markers for the B. fragulis species and its strains
   the optional argument --min_ab specifies the minimum clade abundance for reporting the markers

$ metaphlan2.py -t clade_specific_strain_tracker --clade s__Bacteroides_fragilis metagenome_outfmt.bz2 --input_type bowtie2out -o marker_abundance_table.txt

------------------------------------------------------------------- 

positional arguments:
  INPUT_FILE            the input file can be:
                        * a fastq file containing metagenomic reads
                        OR
                        * a BowTie2 produced SAM file. 
                        OR
                        * an intermediary mapping file of the metagenome generated by a previous MetaPhlAn run 
                        If the input file is missing, the script assumes that the input is provided using the standard 
                        input, or named pipes.
                        IMPORTANT: the type of input needs to be specified with --input_type
  OUTPUT_FILE           the tab-separated output file of the predicted taxon relative abundances 
                        [stdout if not present]

Required arguments:
  --input_type {fastq,fasta,multifasta,multifastq,bowtie2out,sam}
                        set whether the input is the multifasta file of metagenomic reads or 
                        the SAM file of the mapping of the reads against the MetaPhlAn db.
                        [default 'automatic', i.e. the script will try to guess the input format]

Mapping arguments:
  --mpa_pkl MPA_PKL     The metadata pickled MetaPhlAn file [deprecated]
  --bowtie2db METAPHLAN_BOWTIE2_DB
                        The BowTie2 database file of the MetaPhlAn database. Used if --input_type is fastq, fasta, multifasta, or multifastq [default ${mpa_dir}/metaphlan_databases]
  -x INDEX, --index INDEX
                        Specify the id of the database version to use. If the database
                        files are not found on the local MetaPhlAn2 installation they
                        will be automatically downloaded [default latest]
  --bt2_ps BowTie2 presets
                        Presets options for BowTie2 (applied only when a multifasta file is provided)
                        The choices enabled in MetaPhlAn are:
                         * sensitive
                         * very-sensitive
                         * sensitive-local
                         * very-sensitive-local
                        [default very-sensitive]
  --bowtie2_exe BOWTIE2_EXE
                        Full path and name of the BowTie2 executable. This option allowsMetaPhlAn to reach the executable even when it is not in the system PATH or the system PATH is unreachable
  --bowtie2_build BOWTIE2_BUILD
                        Full path to the bowtie2-build command to use, deafult assumes that 'bowtie2-build is present in the system path
  --bowtie2out FILE_NAME
                        The file for saving the output of BowTie2
  --no_map              Avoid storing the --bowtie2out map file
  --tmp_dir             The folder used to store temporary files [default is the OS dependent tmp dir]

Post-mapping arguments:
  --tax_lev TAXONOMIC_LEVEL
                        The taxonomic level for the relative abundance output:
                        'a' : all taxonomic levels
                        'k' : kingdoms
                        'p' : phyla only
                        'c' : classes only
                        'o' : orders only
                        'f' : families only
                        'g' : genera only
                        's' : species only
                        [default 'a']
  --min_cu_len          minimum total nucleotide length for the markers in a clade for
                        estimating the abundance without considering sub-clade abundances
                        [default 2000]
  --min_alignment_len   The sam records for aligned reads with the longest subalignment
                        length smaller than this threshold will be discarded.
                        [default None]
  --ignore_eukaryotes   Do not profile eukaryotic organisms
  --ignore_bacteria     Do not profile bacterial organisms
  --ignore_archaea      Do not profile archeal organisms
  --stat_q              Quantile value for the robust average
                        [default 0.1]
  --perc_nonzero        Percentage of markers with a non zero relative abundance for misidentify a species
                        [default 0.33]
  --ignore_markers IGNORE_MARKERS
                        File containing a list of markers to ignore. 
  --avoid_disqm         Deactivate the procedure of disambiguating the quasi-markers based on the 
                        marker abundance pattern found in the sample. It is generally recommended 
                        to keep the disambiguation procedure in order to minimize false positives
  --stat                EXPERIMENTAL! Statistical approach for converting marker abundances into clade abundances
                        'avg_g'  : clade global (i.e. normalizing all markers together) average
                        'avg_l'  : average of length-normalized marker counts
                        'tavg_g' : truncated clade global average at --stat_q quantile
                        'tavg_l' : trunated average of length-normalized marker counts (at --stat_q)
                        'wavg_g' : winsorized clade global average (at --stat_q)
                        'wavg_l' : winsorized average of length-normalized marker counts (at --stat_q)
                        'med'    : median of length-normalized marker counts
                        [default tavg_g]

Additional analysis types and arguments:
  -t ANALYSIS TYPE      Type of analysis to perform: 
                         * rel_ab: profiling a metagenomes in terms of relative abundances
                         * rel_ab_w_read_stats: profiling a metagenomes in terms of relative abundances and estimate the number of reads comming from each clade.
                         * reads_map: mapping from reads to clades (only reads hitting a marker)
                         * clade_profiles: normalized marker counts for clades with at least a non-null marker
                         * marker_ab_table: normalized marker counts (only when > 0.0 and normalized by metagenome size if --nreads is specified)
                         * marker_counts: non-normalized marker counts [use with extreme caution]
                         * marker_pres_table: list of markers present in the sample (threshold at 1.0 if not differently specified with --pres_th
                        [default 'rel_ab']
  --nreads NUMBER_OF_READS
                        The total number of reads in the original metagenome. It is used only when 
                        -t marker_table is specified for normalizing the length-normalized counts 
                        with the metagenome size as well. No normalization applied if --nreads is not 
                        specified
  --pres_th PRESENCE_THRESHOLD
                        Threshold for calling a marker present by the -t marker_pres_table option
  --clade               The clade for clade_specific_strain_tracker analysis
  --min_ab              The minimum percentage abundace for the clade in the clade_specific_strain_tracker analysis

Output arguments:
  -o output file, --output_file output file
                        The output file (if not specified as positional argument)
  --sample_id_key name  Specify the sample ID key for this analysis. Defaults to '#SampleID'.
  --sample_id value     Specify the sample ID for this analysis. Defaults to 'Metaphlan2_Analysis'.
  -s sam_output_file, --samout sam_output_file
                        The sam output file
  --legacy-output       Old two columns output
  --no-unknown-estimation
                        Ignore estimation of reads mapping to unkwnown clades
  --biom biom_output, --biom_output_file biom_output
                        If requesting biom file output: The name of the output file in biom format 
  --mdelim mdelim, --metadata_delimiter_char mdelim
                        Delimiter for bug metadata: - defaults to pipe. e.g. the pipe in k__Bacteria|p__Proteobacteria 

Other arguments:
  --nproc N             The number of CPUs to use for parallelizing the mapping [default 4]
  --install             Only checks if the MetaPhlAn2 DB is installed and installs it if not. All other parameters are ignored.
  --force_download      Force the re-download of the latest MetaPhlAn2 database.
  --read_min_len READ_MIN_LEN
                        Specify the minimum length of the reads to be considered when parsing the input file with 'read_fastx.py' script, default value is 70
  -v, --version         Prints the current MetaPhlAn version and exit
  -h, --help            show this help message and exit
```

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

**Output files can be merged only if the profiling was performed with the same version of the MetaPhlAn2 database.**

**There is no limit to how many files you can merge.**



### Heatmap Visualization

The hclust2 script generates a hierarchically-clustered heatmap from MetaPhlAn2 abundance profiles. To generate the heatmap for a merged MetaPhlAn2 output table (as described above), you need to run the script as below:

```
#!bash
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

The tutorial of using GraPhlAn can be found from [the MetaPhlAn2 wiki](https://bitbucket.org/biobakery/biobakery/wiki/metaphlan2).


## Customizing the database
In order to add a marker to the database, the user needs the following steps:

* Reconstruct the marker sequences (in fasta format) from the MetaPhlAn2 BowTie2 database by:

```
#!bash
bowtie2-inspect metaphlan_databases/mpa_v25_CHOCOPhlAn_201901 > metaphlan_databases/mpa_v25_CHOCOPhlAn_201901_markers.fasta
```

* Add the marker sequence stored in a file new_marker.fasta to the marker set:

```
#!bash
cat new_marker.fasta >> metaphlan_databases/mpa_v25_CHOCOPhlAn_201901_markers.fasta
```

* Rebuild the bowtie2 database:

```
#!bash
bowtie2-build metaphlan_databases/mpa_v25_CHOCOPhlAn_201901_markers.fasta metaphlan_databases/mpa_v25_CHOCOPhlAn_NEW
```

* Assume that the new marker was extracted from genome1, genome2. Update the taxonomy file from the Python console as follows:

```
#!python

import pickle
import bz2

db = pickle.load(bz2.open('metaphlan_databases/mpa_v25_CHOCOPhlAn_201901.pkl', 'r'))

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
with bz2.BZ2File('metaphlan_databases/mpa_v25_CHOCOPhlAn_NEW.pkl', 'w') as ofile:
    pickle.dump(db, ofile, pickle.HIGHEST_PROTOCOL)

```

* To use the new database, run metaphlan2.py with option "--index v25_CHOCOPhlAn_NEW".