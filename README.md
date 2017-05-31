[TOC]

# MetaPhlAn 2: Metagenomic Phylogenetic Analysis

AUTHORS: Duy Tin Truong (duytin.truong@unitn.it), Nicola Segata (nicola.segata@unitn.it)

## Description
MetaPhlAn is a computational tool for profiling the composition of microbial communities (Bacteria, Archaea, Eukaryotes and Viruses) from metagenomic shotgun sequencing data with species level resolution. From version 2.0, MetaPhlAn is also able to identify specific strains (in the not-so-frequent cases in which the sample contains a previously sequenced strains) and to track strains across samples for all species.

MetaPhlAn 2 relies on ~1M unique clade-specific marker genes ([the marker information file can be found at src/utils/markers_info.txt.bz2 or here](https://bitbucket.org/biobakery/metaphlan2/src/473a41eba501df5f750da032d4f04b38db98dde1/utils/markers_info.txt.bz2?at=default)) identified from ~17,000 reference genomes (~13,500 bacterial and archaeal, ~3,500 viral, and ~110 eukaryotic), allowing:

* unambiguous taxonomic assignments;
* accurate estimation of organismal relative abundance;
* species-level resolution for bacteria, archaea, eukaryotes and viruses;
* strain identification and tracking
* orders of magnitude speedups compared to existing methods.
* metagenomic strain-level population genomics

If you use this software, please cite:

[**MetaPhlAn2 for enhanced metagenomic taxonomic profiling.**](http://www.nature.com/nmeth/journal/v12/n10/pdf/nmeth.3589.pdf) *Duy Tin Truong, Eric A Franzosa, Timothy L Tickle, Matthias Scholz, George Weingart, Edoardo Pasolli, Adrian Tett, Curtis Huttenhower & Nicola Segata*. Nature Methods 12, 902-903 (2015)

-------------

## Pre-requisites

MetaPhlAn requires *python 2.7* or higher with argparse, tempfile and [*numpy*](http://www.numpy.org/) libraries installed 
  (apart for numpy they are usually installed together with the python distribution). 
  Python3 is also now supported.

**If you provide the SAM output of [BowTie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) as input, there are no additional prerequisite.**

* If you would like to use the BowTie2 integrated in MetaPhlAn, you need to have BowTie2 version 2.0.0 or higher and perl installed (bowtie2 needs to be in the system path with execute _and_ read permission)

* If you use the "utils/metaphlan_hclust_heatmap.py" script to plot and hierarchical cluster multiple MetaPhlAn-profiled samples you will also need the following python libraries: [matplotlib](http://matplotlib.org/index.html), [scipy](http://www.scipy.org/), [pylab](http://wiki.scipy.org/PyLab) (if not installed together with MatPlotLib).

* If you want to produce the output as "biom" file you also need [biom](http://biom-format.org/) installed

* MetaPhlAn is not tightly integrated with advanced heatmap plotting with [hclust2](https://bitbucket.org/nsegata/hclust2) and cladogram visualization with [GraPhlAn](https://bitbucket.org/nsegata/graphlan/wiki/Home). If you use such visualization tool please refer to their prerequisites. 

----------------------

## Installation

MetaPhlAn 2.0 can be obtained by either

[Downloading MetaPhlAn v2.0](https://bitbucket.org/biobakery/metaphlan2/get/default.zip)  

**OR**

Cloning the repository via the following commands
``$ hg clone https://bitbucket.org/biobakery/metaphlan2``

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
$ metaphlan2.py metagenome.fastq --input_type fastq > profiled_metagenome.txt
```

It is highly recommended to save the intermediate BowTie2 output for re-running MetaPhlAn extremely quickly (--bowtie2out), and use multiple CPUs (--nproc) if available:

```
#!bash
$ metaphlan2.py metagenome.fastq --bowtie2out metagenome.bowtie2.bz2 --nproc 5 --input_type fastq > profiled_metagenome.txt
```

If you already mapped your metagenome against the marker DB (using a previous  MetaPhlAn run), you can obtain the results in few seconds by using the previously saved --bowtie2out file and specifying the input (--input_type bowtie2out):

```
#!bash
$ metaphlan2.py metagenome.bowtie2.bz2 --nproc 5 --input_type bowtie2out > profiled_metagenome.txt
```

You can also provide an externally BowTie2-mapped SAM if you specify this format with --input_type. Two steps here: first map your metagenome with BowTie2 and then feed MetaPhlAn2 with the obtained sam:

```
#!bash
$ bowtie2 --sam-no-hd --sam-no-sq --no-unal --very-sensitive -S metagenome.sam -x ${mpa_dir}/db_v20/mpa_v20_m200 -U metagenome.fastq
$ metaphlan2.py metagenome.sam --input_type sam > profiled_metagenome.txt
```

In order to make MetaPhlAn 2 easily compatible with complex metagenomic pipeline, there are now multiple alternative ways to pass the input:

```
#!bash
$ cat metagenome.fastq | metaphlan2.py --input_type fastq > profiled_metagenome.txt
```

```
#!bash
$ tar xjf metagenome.tar.bz2 --to-stdout | metaphlan2.py --input_type fastq --bowtie2db ${mpa_dir}/db_v20/mpa_v20_m200 > profiled_metagenome.txt
```

```
#!bash
$ metaphlan2.py --input_type fastq < metagenome.fastq > profiled_metagenome.txt
```

```
#!bash
$ metaphlan2.py --input_type fastq <(bzcat metagenome.fastq.bz2) > profiled_metagenome.txt
```

```
#!bash
$ metaphlan2.py --input_type fastq <(zcat metagenome_1.fastq.gz metagenome_2.fastq.gz) > profiled_metagenome.txt
```

MetaPhlAn 2 can also natively **handle paired-end metagenomes** (but does not use the paired-end information), and, more generally, metagenomes stored in multiple files (but you need to specify the --bowtie2out parameter):

```
#!bash
$ metaphlan2.py metagenome_1.fastq,metagenome_2.fastq --bowtie2out metagenome.bowtie2.bz2 --nproc 5 --input_type fastq > profiled_metagenome.txt
```

For advanced options and other analysis types (such as strain tracking) please refer to the full command-line options.

## Full command-line options


```
usage: metaphlan2.py --input_type
                     {fastq,fasta,multifasta,multifastq,bowtie2out,sam}
                     [--mpa_pkl MPA_PKL] [--bowtie2db METAPHLAN_BOWTIE2_DB]
                     [--bt2_ps BowTie2 presets] [--bowtie2_exe BOWTIE2_EXE]
                     [--bowtie2out FILE_NAME] [--no_map] [--tmp_dir]
                     [--tax_lev TAXONOMIC_LEVEL] [--min_cu_len]
                     [--min_alignment_len] [--ignore_viruses]
                     [--ignore_eukaryotes] [--ignore_bacteria]
                     [--ignore_archaea] [--stat_q]
                     [--ignore_markers IGNORE_MARKERS] [--avoid_disqm]
                     [--stat] [-t ANALYSIS TYPE] [--nreads NUMBER_OF_READS]
                     [--pres_th PRESENCE_THRESHOLD] [--clade] [--min_ab] [-h]
                     [-o output file] [--sample_id_key name]
                     [--sample_id value] [-s sam_output_file]
                     [--biom biom_output] [--mdelim mdelim] [--nproc N] [-v]
                     [INPUT_FILE] [OUTPUT_FILE]

DESCRIPTION
 MetaPhlAn version 2.1.0 (28 April 2015): 
 METAgenomic PHyLogenetic ANalysis for metagenomic taxonomic profiling.

AUTHORS: Nicola Segata (nicola.segata@unitn.it), Duy Tin Truong (duytin.truong@unitn.it)

COMMON COMMANDS

 We assume here that metaphlan2.py is in the system path and that mpa_dir bash variable contains the
 main MetaPhlAn folder. Also BowTie2 should be in the system path with execution and read
 permissions, and Perl should be installed.

========== MetaPhlAn 2 clade-abundance estimation ================= 

The basic usage of MetaPhlAn 2 consists in the identification of the clades (from phyla to species and 
strains in particular cases) present in the metagenome obtained from a microbiome sample and their 
relative abundance. This correspond to the default analysis type (--analysis_type rel_ab).

*  Profiling a metagenome from raw reads:
$ metaphlan2.py metagenome.fastq --input_type fastq

*  You can take advantage of multiple CPUs and save the intermediate BowTie2 output for re-running
   MetaPhlAn extremely quickly:
$ metaphlan2.py metagenome.fastq --bowtie2out metagenome.bowtie2.bz2 --nproc 5 --input_type fastq

*  If you already mapped your metagenome against the marker DB (using a previous MetaPhlAn run), you
   can obtain the results in few seconds by using the previously saved --bowtie2out file and 
   specifying the input (--input_type bowtie2out):
$ metaphlan2.py metagenome.bowtie2.bz2 --nproc 5 --input_type bowtie2out

*  You can also provide an externally BowTie2-mapped SAM if you specify this format with 
   --input_type. Two steps: first apply BowTie2 and then feed MetaPhlAn2 with the obtained sam:
$ bowtie2 --sam-no-hd --sam-no-sq --no-unal --very-sensitive -S metagenome.sam -x ${mpa_dir}/db_v20/mpa_v20_m200 -U metagenome.fastq
$ metaphlan2.py metagenome.sam --input_type sam > profiled_metagenome.txt

*  Multiple alternative ways to pass the input are also available:
$ cat metagenome.fastq | metaphlan2.py --input_type fastq 
$ tar xjf metagenome.tar.bz2 --to-stdout | metaphlan2.py --input_type fastq 
$ metaphlan2.py --input_type fastq < metagenome.fastq
$ metaphlan2.py --input_type fastq <(bzcat metagenome.fastq.bz2)
$ metaphlan2.py --input_type fastq <(zcat metagenome_1.fastq.gz metagenome_2.fastq.gz)

*  We can also natively handle paired-end metagenomes, and, more generally, metagenomes stored in 
  multiple files (but you need to specify the --bowtie2out parameter):
$ metaphlan2.py metagenome_1.fastq,metagenome_2.fastq --bowtie2out metagenome.bowtie2.bz2 --nproc 5 --input_type fastq

------------------------------------------------------------------- 
 

========== MetaPhlAn 2 strain tracking ============================ 

MetaPhlAn 2 introduces the capability of charaterizing organisms at the strain level using non
aggregated marker information. Such capability comes with several slightly different flavours and 
are a way to perform strain tracking and comparison across multiple samples.
Usually, MetaPhlAn 2 is first ran with the default --analysis_type to profile the species present in
the community, and then a strain-level profiling can be performed to zoom-in into specific species
of interest. This operation can be performed quickly as it exploits the --bowtie2out intermediate 
file saved during the execution of the default analysis type.

*  The following command will output the abundance of each marker with a RPK (reads per kil-base) 
   higher 0.0. (we are assuming that metagenome_outfmt.bz2 has been generated before as 
   shown above).
$ metaphlan2.py -t marker_ab_table metagenome_outfmt.bz2 --input_type bowtie2out > marker_abundance_table.txt
   The obtained RPK can be optionally normalized by the total number of reads in the metagenome 
   to guarantee fair comparisons of abundances across samples. The number of reads in the metagenome
   needs to be passed with the '--nreads' argument

*  The list of markers present in the sample can be obtained with '-t marker_pres_table'
$ metaphlan2.py -t marker_pres_table metagenome_outfmt.bz2 --input_type bowtie2out > marker_abundance_table.txt
   The --pres_th argument (default 1.0) set the minimum RPK value to consider a marker present

*  The list '-t clade_profiles' analysis type reports the same information of '-t marker_ab_table'
   but the markers are reported on a clade-by-clade basis.
$ metaphlan2.py -t clade_profiles metagenome_outfmt.bz2 --input_type bowtie2out > marker_abundance_table.txt

*  Finally, to obtain all markers present for a specific clade and all its subclades, the 
   '-t clade_specific_strain_tracker' should be used. For example, the following command
   is reporting the presence/absence of the markers for the B. fragulis species and its strains
$ metaphlan2.py -t clade_specific_strain_tracker --clade s__Bacteroides_fragilis metagenome_outfmt.bz2 db_v20/mpa_v20_m200.pkl --input_type bowtie2out > marker_abundance_table.txt
   the optional argument --min_ab specifies the minimum clade abundance for reporting the markers

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
  --mpa_pkl MPA_PKL     the metadata pickled MetaPhlAn file
  --input_type {fastq,fasta,multifasta,multifastq,bowtie2out,sam}
                        set whether the input is the multifasta file of metagenomic reads or 
                        the SAM file of the mapping of the reads against the MetaPhlAn db.
                        [default 'automatic', i.e. the script will try to guess the input format]

Mapping arguments:
  --bowtie2db METAPHLAN_BOWTIE2_DB
                        The BowTie2 database file of the MetaPhlAn database. 
                        Used if --input_type is fastq, fasta, multifasta, or multifastq
  --bt2_ps BowTie2 presets
                        presets options for BowTie2 (applied only when a multifasta file is provided)
                        The choices enabled in MetaPhlAn are:
                         * sensitive
                         * very-sensitive
                         * sensitive-local
                         * very-sensitive-local
                        [default very-sensitive]
  --bowtie2_exe BOWTIE2_EXE
                        Full path and name of the BowTie2 executable. This option allows 
                        MetaPhlAn to reach the executable even when it is not in the system 
                        PATH or the system PATH is unreachable
  --bowtie2out FILE_NAME
                        The file for saving the output of BowTie2
  --no_map              Avoid storing the --bowtie2out map file
  --tmp_dir             the folder used to store temporary files 
                        [default is the OS dependent tmp dir]

Post-mapping arguments:
  --tax_lev TAXONOMIC_LEVEL
                        The taxonomic level for the relative abundance output:
                        'a' : all taxonomic levels
                        'k' : kingdoms (Bacteria and Archaea) only
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
  --ignore_viruses      Do not profile viral organisms
  --ignore_eukaryotes   Do not profile eukaryotic organisms
  --ignore_bacteria     Do not profile bacterial organisms
  --ignore_archaea      Do not profile archeal organisms
  --stat_q              Quantile value for the robust average
                        [default 0.1]
  --ignore_markers IGNORE_MARKERS
                        File containing a list of markers to ignore. 
  --avoid_disqm         Descrivate the procedure of disambiguating the quasi-markers based on the 
                        marker abundance pattern found in the sample. It is generally recommended 
                        too keep the disambiguation procedure in order to minimize false positives
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
  -h, --help            show this help message and exit

Output arguments:
  -o output file, --output_file output file
                        The output file (if not specified as positional argument)
  --sample_id_key name  Specify the sample ID key for this analysis. Defaults to '#SampleID'.
  --sample_id value     Specify the sample ID for this analysis. Defaults to 'Metaphlan2_Analysis'.
  -s sam_output_file, --samout sam_output_file
                        The sam output file
  --biom biom_output, --biom_output_file biom_output
                        If requesting biom file output: The name of the output file in biom format 
  --mdelim mdelim, --metadata_delimiter_char mdelim
                        Delimiter for bug metadata: - defaults to pipe. e.g. the pipe in k__Bacteria|p__Proteobacteria 

Other arguments:
  --nproc N             The number of CPUs to use for parallelizing the mapping
                        [default 1, i.e. no parallelism]
  -v, --version         Prints the current MetaPhlAn version and exit


```

## Utility Scripts

MetaPhlAn's repository features a few utility scripts to aid in manipulation of sample output and its visualization. These scripts can be found under the ``utils`` folder in the metaphlan2 directory.

### Merging Tables

The script **merge\_metaphlan\_tables.py** allows to combine MetaPhlAn output from several samples to be merged into one table Bugs (rows) vs Samples (columns) with the table enlisting the relative normalized abundances per sample per bug.

To merge multiple output files, run the script as below

```
#!bash
$ python utils/merge_metaphlan_tables.py metaphlan_output1.txt metaphlan_output2.txt metaphlan_output3.txt > output/merged_abundance_table.txt
```

Wildcards can be used as needed:

```
#!bash
$ python utils/merge_metaphlan_tables.py metaphlan_output*.txt > output/merged_abundance_table.txt
```

**There is no limit to how many files you can merge.**

## Heatmap Visualization

The script **metaphlan\_hclust\_heatmap.py** allows to visualize the MetaPhlAn results in the form of a hierarchically-clustered heatmap. To generate the heatmap for a merged MetaPhlAn output table (as described above), please run the script as below.

```
#!bash
$ python utils/metaphlan_hclust_heatmap.py -c bbcry --top 25 --minv 0.1 -s log --in output/merged_abundance_table.txt --out output_images/abundance_heatmap.png
```

For detailed command-line instructions, please refer to below:


```
$ utils/metaphlan_hclust_heatmap.py -h
usage: metaphlan_hclust_heatmap.py [-h] --in INPUT_FILE --out OUTPUT_FILE
                                   [-m {single,complete,average,weighted,centroid,median,ward}]
                                   [-d {euclidean,minkowski,cityblock,seuclidean,sqeuclidean,cosine,correlation,hamming,jaccard,chebyshev,canberra,braycurtis,mahalanobis,yule,matching,dice,kulsinski,rogerstanimoto,russellrao,sokalmichener,sokalsneath,wminkowski,ward}]
                                   [-f {euclidean,minkowski,cityblock,seuclidean,sqeuclidean,cosine,correlation,hamming,jaccard,chebyshev,canberra,braycurtis,mahalanobis,yule,matching,dice,kulsinski,rogerstanimoto,russellrao,sokalmichener,sokalsneath,wminkowski,ward}]
                                   [-s scale norm] [-x X] [-y Y] [--minv MINV]
                                   [--maxv max value]
                                   [--tax_lev TAXONOMIC_LEVEL] [--perc PERC]
                                   [--top TOP] [--sdend_h SDEND_H]
                                   [--fdend_w FDEND_W] [--cm_h CM_H]
                                   [--cm_ticks label for ticks of the colormap]
                                   [--font_size FONT_SIZE]
                                   [--clust_line_w CLUST_LINE_W]
                                   [-c {Accent,Blues,BrBG,BuGn,BuPu,Dark2,GnBu,Greens,Greys,OrRd,Oranges,PRGn,Paired,Pastel1,Pastel2,PiYG,PuBu,PuBuGn,PuOr,PuRd,Purples,RdBu,RdGy,RdPu,RdYlBu,RdYlGn,Reds,Set1,Set2,Set3,Spectral,YlGn,YlGnBu,YlOrBr,YlOrRd,afmhot,autumn,binary,bone,brg,bwr,cool,copper,flag,gist_earth,gist_gray,gist_heat,gist_ncar,gist_rainbow,gist_stern,gist_yarg,gnuplot,gnuplot2,gray,hot,hsv,jet,ocean,pink,prism,rainbow,seismic,spectral,spring,summer,terrain,winter,bbcyr,bbcry}]

This scripts generates heatmaps with hierarchical clustering of both samples
and microbial clades. The script can also subsample the number of clades to                                                                                                                                                                                                    
display based on the their nth percentile abundance value in each sample                                                                                                                                                                                                       
                                                                                                                                                                                                                                                                               
optional arguments:                                                                                                                                                                                                                                                            
  -h, --help            show this help message and exit                                                                                                                                                                                                                        
  --in INPUT_FILE       The input file of microbial relative abundances. This                                                                                                                                                                                                  
                        file is typically obtained with the                                                                                                                                                                                                                    
                        "utils/merge_metaphlan_tables.py"                                                                                                                                                                                                                      
  --out OUTPUT_FILE     The output image. The extension of the file determines                                                                                                                                                                                                 
                        the image format. png, pdf, and svg are the preferred                                                                                                                                                                                                  
                        format                                                                                                                                                                                                                                                 
  -m {single,complete,average,weighted,centroid,median,ward}                                                                                                                                                                                                                   
                        The hierarchical clustering method, default is                                                                                                                                                                                                         
                        "average"                                                                                                                                                                                                                                              
  -d {euclidean,minkowski,cityblock,seuclidean,sqeuclidean,cosine,correlation,hamming,jaccard,chebyshev,canberra,braycurtis,mahalanobis,yule,matching,dice,kulsinski,rogerstanimoto,russellrao,sokalmichener,sokalsneath,wminkowski,ward}                                      
                        The distance function for samples. Default is                                                                                                                                                                                                          
                        "braycurtis"                                                                                                                                                                                                                                           
  -f {euclidean,minkowski,cityblock,seuclidean,sqeuclidean,cosine,correlation,hamming,jaccard,chebyshev,canberra,braycurtis,mahalanobis,yule,matching,dice,kulsinski,rogerstanimoto,russellrao,sokalmichener,sokalsneath,wminkowski,ward}                                      
                        The distance function for microbes. Default is                                                                                                                                                                                                         
                        "correlation"                                                                                                                                                                                                                                          
  -s scale norm                                                                                                                                                                                                                                                                
  -x X                  Width of heatmap cells. Automatically set, this option                                                                                                                                                                                                 
                        should not be necessary unless for very large heatmaps                                                                                                                                                                                                 
  -y Y                  Height of heatmap cells. Automatically set, this                                                                                                                                                                                                       
                        option should not be necessary unless for very large                                                                                                                                                                                                   
                        heatmaps                                                                                                                                                                                                                                               
  --minv MINV           Minimum value to display. Default is 0.0, values                                                                                                                                                                                                       
                        around 0.001 are also reasonable                                                                                                                                                                                                                       
  --maxv max value      Maximum value to display. Default is maximum value                                                                                                                                                                                                     
                        present, can be set e.g. to 100 to display the full
                        scale
  --tax_lev TAXONOMIC_LEVEL
                        The taxonomic level to display: 'a' : all taxonomic
                        levels 'k' : kingdoms (Bacteria and Archaea) only 'p'
                        : phyla only 'c' : classes only 'o' : orders only 'f'
                        : families only 'g' : genera only 's' : species only
                        [default 's']
  --perc PERC           Percentile to be used for ordering the microbes in
                        order to select with --top the most abundant microbes
                        only. Default is 90
  --top TOP             Display the --top most abundant microbes only
                        (ordering based on --perc)
  --sdend_h SDEND_H     Set the height of the sample dendrogram. Default is
                        0.1
  --fdend_w FDEND_W     Set the width of the microbes dendrogram. Default is
                        0.1
  --cm_h CM_H           Set the height of the colormap. Default = 0.03
  --cm_ticks label for ticks of the colormap
  --font_size FONT_SIZE
                        Set label font sizes. Default is 7
  --clust_line_w CLUST_LINE_W
                        Set the line width for the dendrograms
  -c {Accent,Blues,BrBG,BuGn,BuPu,Dark2,GnBu,Greens,Greys,OrRd,Oranges,PRGn,Paired,Pastel1,Pastel2,PiYG,PuBu,PuBuGn,PuOr,PuRd,Purples,RdBu,RdGy,RdPu,RdYlBu,RdYlGn,Reds,Set1,Set2,Set3,Spectral,YlGn,YlGnBu,YlOrBr,YlOrRd,afmhot,autumn,binary,bone,brg,bwr,cool,copper,flag,gist_earth,gist_gray,gist_heat,gist_ncar,gist_rainbow,gist_stern,gist_yarg,gnuplot,gnuplot2,gray,hot,hsv,jet,ocean,pink,prism,rainbow,seismic,spectral,spring,summer,terrain,winter,bbcyr,bbcry}
                        Set the colormap. Default is "jet".
```

### GraPhlAn Visualization

The tutorial of using GraPhlAn can be found from [the MetaPhlAn2 wiki](https://bitbucket.org/biobakery/biobakery/wiki/metaphlan2).


## Customizing the database
In order to add a marker to the database, the user needs the following steps:

* Reconstruct the marker sequences (in fasta format) from the MetaPhlAn2 bowtie2 database by:

```
#!bash
bowtie2-inspect metaphlan2/db_v20/mpa_v20_m200 > metaphlan2/markers.fasta
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

db = pickle.load(bz2.BZ2File('db_v20/mpa_v20_m200.pkl', 'r'))

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
StrainPhlAn requires *python 2.7* and the libraries [pysam](http://pysam.readthedocs.org/en/latest/) (tested on **version 0.8.3**), [biopython](http://biopython.org/wiki/Main_Page), [msgpack](https://pypi.python.org/pypi/msgpack-python) and [numpy](http://www.numpy.org/), [dendropy](https://pythonhosted.org/DendroPy/) (tested on version **3.12.0**). Besides, StrainPhlAn also needs the following programs in the executable path:

* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for mapping reads against the marker database.
* [MUSCLE](http://www.drive5.com/muscle/) for the alignment step.
* [samtools, bcftools and vcfutils.pl](http://samtools.sourceforge.net/) which can be downloaded from [here](https://github.com/samtools) for building consensus markers. Note that vcfutils.pl is included in bcftools and **StrainPhlAn only works with samtools version 0.1.19** as samtools has changed the output format after this version.
* [blastn](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) for adding reference genomes to the phylogenetic tree.
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
    tar xjfO ${f} | ../metaphlan2.py --bowtie2db ../db_v20/mpa_v20_m200 --mpa_pkl ../db_v20/mpa_v20_m200.pkl --input_type multifastq --nproc 10 -s sams/${bn}.sam.bz2 --bowtie2out sams/${bn}.bowtie2_out.bz2 -o sams/${bn}.profile
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
bowtie2-inspect ../db_v20/mpa_v20_m200 > db_markers/all_markers.fasta
python ../strainphlan_src/extract_markers.py --mpa_pkl ../db_v20/mpa_v20_m200.pkl --ifn_markers db_markers/all_markers.fasta --clade s__Bacteroides_caccae --ofn_markers db_markers/s__Bacteroides_caccae.markers.fasta
```

Note that the "all\_markers.fasta" file consists can be reused for extracting other reference genomes. 
After this step, you should have two files in folder "db\_markers": "all\_markers.fasta" containing all marker sequences, and "s\_\_Bacteroides\_caccae.markers.fasta" containing the markers of *Bacteroides caccae*.
This step will take around 1 minute and can skipped if you do not need to add the reference genomes to the phylogenetic tree. Those markers can be found in the folder "db\_markers" in [this link](https://www.dropbox.com/sh/m4na8wefp53j8ej/AABA3yVsG26TbB0t1cnBS9-Ra?dl=0).

Before building the trees, we should get the list of all clades detected from the samples and save them in the "output/clades.txt" file by the following command:

```
#!python
python ../strainphlan.py --mpa_pkl ../db_v20/mpa_v20_m200.pkl --ifn_samples consensus_markers/*.markers --output_dir output --nprocs_main 10 --print_clades_only > output/clades.txt
```

The clade names in the output file "clades.txt" will be used for the next step.

Step 5. Build the multiple sequence alignment and phylogenetic tree:

This step will align and clean the *sample-reconstructed strains* (stored in the marker files produced in step 3) and *reference-genome-reconstructed strains* (extracted based on the database markers in step 4) to produce a multiple sequence alignment (MSA) and store it in the file "clade_name.fasta". From this MSA file, StrainPhlAn will call RAxML to build the phylogenetic tree.
Note that: all marker files (\*.markers) **must be used together** as the input for the strainphlan.py script because StrainPhlAn needs to align all of the strains at once.

The commands to run are:

```
#!python
mkdir -p output
python ../strainphlan.py --mpa_pkl ../db_v20/mpa_v20_m200.pkl --ifn_samples consensus_markers/*.markers --ifn_markers db_markers/s__Bacteroides_caccae.markers.fasta --ifn_ref_genomes reference_genomes/G000273725.fna.bz2 --output_dir output --nprocs_main 10 --clades s__Bacteroides_caccae &> output/log_full.txt
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