[TOC]

#**MetaPhlAn 2.0: Metagenomic Phylogenetic Analysis**#

AUTHORS: Nicola Segata (nicola.segata@unitn.it)

##**Description**##
MetaPhlAn is a computational tool for profiling the composition of microbial communities (Bacteria, Archaea, Eukaryotes and Viruses) from metagenomic shotgun sequencing data with species level resolution. From version 2.0 MetaPhlAn is also able to identify specific strains (in the not-so-frequent cases in which the sample contains a previously sequenced strains) and to track strains across samples for all species.

MetaPhlAn 2.0 relies on ~1M unique clade-specific marker genes identified from ~17,000 reference genomes (~13,500 bacterial and archaeal, ~3,500 viral, and ~110 eukaryotic), allowing:

* unambiguous taxonomic assignments;
* accurate estimation of organismal relative abundance;
* species-level resolution for bacteria, archaea, eukaryotes and viruses;
* strain identification and tracking
* orders of magnitude speedups compared to existing methods.

If you use this software, please cite :

[**Metagenomic microbial community profiling using unique clade-specific marker genes.**](http://www.nature.com/nmeth/journal/v9/n8/full/nmeth.2066.html)
 *Nicola Segata, Levi Waldron, Annalisa Ballarini, Vagheesh Narasimhan, Olivier Jousson, Curtis Huttenhower*. 
Nature Methods, 8, 811–814, 2012

-------------

##**Pre-requisites**##

* MetaPhlAn requires python 2.7 or higher with argparse, tempfile and numpy libraries installed 
  (apart for numpy they are usually installed together with the python distribution). 
  Python3 is not supported yet.

**If you provide the output of BLASTN or BowTie2 as input, there are no additional prerequisite.**

If you would like to use the BowTie2 integrated in MetaPhlAn, you need to have:

* BowTie2 version 2.0.0 or higher and perl (bowtie2 needs to be in the system path with execute _and_ read permission)

If you use the "utils/metaphlan_hclust_heatmap.py" script to plot and hierarchial cluster multiple metaphlan-profiled samples you will also need the following python libraries: matplotlib, scipy, pylab 

----------------------

##**Installation**##

MetaPhlAn v2.0 can be obtained by either

* [Downloading MetaPhlAn v2.0](metaphlan2/downloads)

**OR**

* Cloning the repository via the following commands
* * ``$ hg clone https://bitbucket.org/biobakery/metaphlan2``

--------------------------


##**Basic Usage**##

* Profiling a metagenome from raw reads (requires BowTie2 in the system path with execution and read permissions, Perl installed, and the BowTie2 marker DB provided with MetaPhlAn). Here, mpa.pkl is the marker metadata file provided with the MetaPhlAn package:

* * ``$ metaphlan2.py metagenome.fastq --mpa_pkl mpa.pkl --bowtie2db bowtie2db/mpa``

* You can take advantage of multiple CPUs and you can save the intermediate BowTie2 output for re-running MetaPhlAn extremely quickly:

* * ``metaphlan2.py metagenome.fastq --mpa_pkl mpa.pkl --bowtie2db bowtie2db/mpa --nproc 5 --bowtie2out metagenome.bt2out.bz2``

* If you already mapped your metagenome against the marker DB (using a previous  MetaPhlAn run, you can obtain the results in few seconds:

* * ``metaphlan2.py --input_type bowtie2out --mpa_pkl mpa.pkl metagenome.bowtie2out.bz2``

* Standard SAM file can be converted into bowtie2out using the command below, and then run through MetaPhlAn as described above.

* * ``cat file.sam | cut -f 1,3 | grep -v "*" > file.bowtie2out.txt``

* The metagenome can also be passed from the standard input but it is necessary to specify the input format explicitly:

* * ``tar xjf metagenome.tar.bz2 --to-stdout | metaphlan2.py --input_type multifastq --mpa_pkl mpa.pkl --bowtie2db bowtie2db/mpa``

* Also the pre-computed BowTie2 output can be provided with a pipe (again specifying the input type): 

* * ``metaphlan2.py --input_type bowtie2out --mpa_pkl mpa.pkl < metagenome.bowtie2out.txt > profiling_output.txt``

* You can also set advanced options for the BowTie2 step selecting the preset option 
  among 'sensitive','very-sensitive','sensitive-local','very-sensitive-local' 
  (valid for metagenome as input only):

* * ``metaphlan2.py --bt2_ps very-sensitive-local --mpa_pkl mpa.pkl metagenome.fasta``



##**Full command-line options**##


```
#!cmd

usage: metaphlan2.py [-h] [-v] [--mpa_pkl] [--stat] [-t ANALYSIS TYPE]
                     [--tax_lev TAXONOMIC_LEVEL] [--nreads NUMBER_OF_READS]
                     [--pres_th PRESENCE_THRESHOLD]
                     [--bowtie2db METAPHLAN_BOWTIE2_DB]
                     [--bt2_ps BowTie2 presets] [--tmp_dir] [--min_cu_len]
                     [--input_type {automatic,multifasta,multifastq,bowtie2out,sam}]
                     [--ignore_viruses] [--ignore_eukaryotes]
                     [--ignore_bacteria] [--ignore_archaea] [--stat_q]
                     [--avoid_disqm] [--bowtie2_exe BOWTIE2_EXE]
                     [--bowtie2out FILE_NAME] [--no_map] [-o output file]
                     [--nproc N] [--biom biom_output] [--mdelim mdelim]
                     [INPUT_FILE] [OUTPUT_FILE]

positional arguments:
  INPUT_FILE            the input file can be:
                        * a multi-fasta file containing metagenomic reads
                        OR
                        * a NCBI BLAST output file (-outfmt 6 format) of the metagenome against the MetaPhlAn database. 
                        OR
                        * a BowTie2 output file of the metagenome generated by a previous MetaPhlAn run 
                        The software will recognize the format automatically.
                        If the input file is missing, the script assumes that the input is provided using the standard 
                        input, and the input format has to be specified with --input_type
  OUTPUT_FILE           the tab-separated output file of the predicted taxon relative abundances 
                        [stdout if not present]

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Prints the current MetaPhlAn version and exit
  --mpa_pkl             the metadata pickled MetaPhlAn file
  --stat                EXPERIMENTAL! Statistical approach for converting marker abundances into clade abundances
                        'avg_g'  : clade global (i.e. normalizing all markers together) average
                        'avg_l'  : average of length-normalized marker counts
                        'tavg_g' : truncated clade global average at --stat_q quantile
                        'tavg_l' : trunated average of length-normalized marker counts (at --stat_q)
                        'wavg_g' : winsorized clade global average (at --stat_q)
                        'wavg_l' : winsorized average of length-normalized marker counts (at --stat_q)
                        'med'    : median of length-normalized marker counts
                        [default tavg_g]
  -t ANALYSIS TYPE      Type of analysis to perform: 
                         * rel_ab: profiling a metagenomes in terms of relative abundances
                         * reads_map: mapping from reads to clades (only reads hitting a marker)
                         * clade_profiles: normalized marker counts for clades with at least a non-null marker
                         * marker_ab_table: normalized marker counts (only when > 0.0 and normalized by metagenome size if --nreads is specified)
                         * marker_pres_table: list of markers present in the sample (threshold at 1.0 if not differently specified with --pres_th
                        [default 'rel_ab']
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
  --nreads NUMBER_OF_READS
                        The total number of reads in the original metagenome. It is used only when 
                        -t marker_table is specified for normalizing the length-normalized counts 
                        with the metagenome size as well. No normalization applied if --nreads is not 
                        specified
  --pres_th PRESENCE_THRESHOLD
                        Threshold for calling a marker present by the -t marker_pres_table option
  --bowtie2db METAPHLAN_BOWTIE2_DB
                        The BowTie2 database file of the MetaPhlAn database 
  --bt2_ps BowTie2 presets
                        presets options for BowTie2 (applied only when a multifasta file is provided)
                        The choices enabled in MetaPhlAn are:
                         * sensitive
                         * very-sensitive
                         * sensitive-local
                         * very-sensitive-local
                        [default very-sensitive]
  --tmp_dir             the folder used to store temporary files 
                        [default is the OS dependent tmp dir]
  --min_cu_len          minimum total nucleotide length for the markers in a clade for
                        estimating the abundance without considering sub-clade abundances
                        [default 2000]
  --input_type {automatic,multifasta,multifastq,bowtie2out,sam}
                        set whet`er the input is the multifasta file of metagenomic reads or 
                        the SAM file of the mapping of the reads against the MetaPhlAn db.
                        [default 'automatic', i.e. the script will try to guess the input format]
  --ignore_viruses      Do not profile viral organisms
  --ignore_eukaryotes   Do not profile eukaryotic organisms
  --ignore_bacteria     Do not profile bacterial organisms
  --ignore_archaea      Do not profile archeal organisms
  --stat_q              Quantile value for the robust average
                        [default 0.1]
  --avoid_disqm         Descrivate the procedure of disambiguating the quasi-markers based on the 
                        marker abundance pattern found in the sample. It is generally recommended 
                        too keep the disambiguation procedure in order to minimize false positives
  --bowtie2_exe BOWTIE2_EXE
                        Full path and name of the BowTie2 executable. This option allows 
                        MetaPhlAn to reach the executable even when it is not in the system 
                        PATH or the system PATH is unreachable
  --bowtie2out FILE_NAME
                        The file for saving the output of BowTie2
  --no_map              Avoid storing the --bowtie2out map file
  -o output file, --output_file output file
                        The output file (if not specified as positional argument)
  --nproc N             The number of CPUs to use for parallelizing the mapping
                        [default 1, i.e. no parallelism]
  --biom biom_output, --biom_output_file biom_output
                        If requesting biom file output: The name of the output file in biom format 
  --mdelim mdelim, --metadata_delimiter_char mdelim
                        Delimiter for bug metadata: - defaults to pipe. e.g. the pipe in k__Bacteria|p__Proteobacteria
```

##**Utility Scripts**##

MetaPhlAn's repository features a few utility scripts to aid in manipulation of sample output and its visualization. These scripts can be found under the ``utils`` folder in the metaphlan2 directory.

###**Merging Tables**###

The script **merge_metaphlan_tables.py** allows to combine MetaPhlAn output from several samples to be merged into one table Bugs (rows) vs Samples (columns) with the table enlisting the relative normalized abundances per sample per bug.

To merge multiple output files, run the script as below

``$ python utils/merge_metaphlan_tables.py metaphlan_output1.txt metaphlan_output2.txt metaphlan_output3.txt > output/merged_abundance_table.txt``

**There is no limit to how many files you can merge.**

##**Heatmap Visualization**##

The script **metaphlan_hclust_heatmap.py** allows to visualize the MetaPhlAn results in the form of a hierarchically-clustered heatmap. To generate the heatmap for a merged metaphlan output table (as described above), please run the script as below.

``$ python utils/metaphlan_hclust_heatmap.py -c bbcry --top 25 --minv 0.1 -s log --in output/merged_abundance_table.txt --out output_images/abundance_heatmap.png``

For detailed command-line instructions, please refer to below:


```
#!

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

###**GraPhlAn Visualization**###

Utilities also features **metaphlan2graphlan.py** script that provides a way to automatically create the two input files to create a [GraPhlAn](http://huttenhower.sph.harvard.edu/graphlan) cladogram. To convert the MetaPhlAn output into input for GraPhlAn, please run the following script.

``$ python utils/metaphlan2graphlan2.py merged_abundance_table.txt --tree_file merged_table.tree --annot_file merged_table.annot``

The above script will create the files **merged_table.tree** and **merged_table.annot** which you can then provide to GraPhlAn to create the cladogram.

For details, please refer to GraPhlAn's documentation.