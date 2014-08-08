============================================
MetaPhlAn 2.0: Metagenomic Phylogenetic Analysis
============================================

AUTHORS: Nicola Segata (nicola.segata@unitn.it)

DESCRIPTION
 MetaPhlAn version 2.0.0 beta2 (12 July 2014): METAgenomic PHyLogenetic ANalysis for
 taxonomic classification of metagenomic reads.

MetaPhlAn is a computational tool for profiling the composition of microbial communities (Bacteria, Archaea, Eukaryotes and Viruses) from metagenomic shotgun sequencing data. 
MetaPhlAn relies on unique clade-specific marker genes identified from ~17,000 reference genomes (~13,500 bacterial and archaeal, ~3,500 viral, and ~110 eukaryotic), allowing:
- orders of magnitude speedups compared to existing methods;
- unambiguous taxonomic assignments;
- accurate estimation of organismal relative abundance;
- species-level resolution for bacteria, archaea, eukaryotes and viruses.

If you use this software, please cite our paper: ****
"Metagenomic microbial community profiling using unique clade-specific marker genes" 
Nicola Segata, Levi Waldron, Annalisa Ballarini, Vagheesh Narasimhan, Olivier Jousson, Curtis Huttenhower. 
Nature Methods, in press

-------------
PREREQUISITES
-------------
* MetaPhlAn requires python 2.7 or higher with argparse, tempfile and numpy libraries installed 
  (apart for numpy they are usually installed together with the python distribution). 
  Python3 is not supported yet.

If you provide the output of BLASTN or BowTie2 as input, there are no additional prerequisite.

If you would like to use the BowTie2 integrated in MetaPhlAn, you need to have:
* BowTie2 version 2.0.0 or higher and perl (bowtie2 needs to be in the system path with execute _and_ read permission)

If you use the "utils/metaphlan_hclust_heatmap.py" script to plot and hierarchial cluster multiple metaphlan-profiled samples you will also need the following python libraries: matplotlib, scipy, pylab 

-----------
BASIC USAGE
-----------

COMMON COMMANDS

* Profiling a metagenome from raw reads (requires BowTie2 in the system path 
  with execution and read permissions, Perl installed, and the BowTie2 marker DB 
  provided with MetaPhlAn):
metaphlan2.py metagenome.fastq --mpa_pkl mpa.pkl --bowtie2db bowtie2db/mpa
  mpa.pkl is the marker metadata file provided with the MetaPhlAn package
  Although not optimal, also reads in fasta format can be used. 

* You can take advantage of multiple CPUs and you can save the intermediate BowTie2
  output
 for re-running MetaPhlAn extremely quickly:
metaphlan.py metagenome.fastq --mpa_pkl mpa.pkl --bowtie2db bowtie2db/mpa --nproc 5 --bowtie2out metagenome.bt2out.bz2

* If you already mapped your metagenome against the marker DB (using a previous 
  MetaPhlAn run, you can obtain the results in few seconds:
metaphlan2.py --input_type bowtie2out --mpa_pkl mpa.pkl metagenome.bowtie2out.bz2
  (notice that 'bowtie2out' file is automatically bzip2 compressed/uncompressed) 
  a standard SAM file as follows: 
cat file.sam | cut -f 1,3 | grep -v "*" > file.bowtie2out.txt

* The metagenome can also be passed from the standard input but 
  it is necessary to specify the input format explicitly:
tar xjf metagenome.tar.bz2 --to-stdout | metaphlan.py --input_type multifastq --mpa_pkl mpa.pkl --bowtie2db bowtie2db/mpa

* Also the pre-computed BowTie2 output can be provided with a pipe (again 
  specifying the input type): 
metaphlan2.py --input_type bowtie2out --mpa_pkl mpa.pkl < metagenome.bowtie2out.txt > profiling_output.txt

* You can also set advanced options for the BowTie2 step selecting the preset option 
  among 'sensitive','very-sensitive','sensitive-local','very-sensitive-local' 
  (valid for metagenome as input only):
metaphlan2.py --bt2_ps very-sensitive-local --mpa_pkl mpa.pkl metagenome.fasta


-------------------------
FULL COMMAND LINE OPTIONS
-------------------------

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
