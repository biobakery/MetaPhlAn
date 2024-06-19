## Version 4.1.1 (Mar 11th, 2024)
### Database updates
* We just released the new vJun23_202403 and vOct22_202403 databases
  * Same SGBs as for the vJun23_202307 and vOct22_202212 versions, respectively, but the NCBI taxonomy assignment has been fixed to keep the taxa consistent across the MetaPhlAn taxonomic tree, allowing accurate relative abundance estimation also at higher taxonomic levels
### New features
* [MetaPhlAn] The new `fix_relab_mpa4.py` script enables to fix errors in the relative abundances in profiles generated with previous databases
* [MetaPhlAn] Implementation of the option `--subsampling_paired [N_PAIRED_READS]` to subsample paired-end input reads. It needs to be used in conjunction with `-1 [FORWARD_READS_FILE]` and `-2 [REVERSE_READS_FILE]`
### Fixes
* [MetaPhlAn] Fixed a bug that would halt MetaPhlAn execution when the option `--profile_vsc` was used but had no viral hits
* [MetaPhlAn] Fixed a bug that would halt MetaPhlAn execution when the number of reads to map was zero
* [StrainPhlAn] Fixed a bug in the new implementation (since v4.1) of `–-print_clades_only`
  <br/>

## Version 4.1.0 (Feb 20th, 2024)
### Database updates
* We just released the new vJun23_202307 database
  * Addition of ~45k reference genomes from NCBI
  * Addition of ~50k MAGs from ocean, ~40k MAGs from soil, ~30k MAGs from domestic animals and non-human primates, ~4k MAGs from giant turtles, ~7.5k MAGs from skin microbiome, ~20k MAGs from dental plaque, ~15k MAGs from Asian populations, ~2.7k MAGs from ancient and modern Bolivians and other small datasets from diverse sources
  * Expansion of the markers database with 36,822 SGBs (6,272 more SGBs than in vOct22)
* Inclusion of the new Viral Sequence Clusters (VSCs) database
  * Containing 3,944 VSCs clustered into 1,345 Viral Sequence Groups (VSGs).
  * Including a total of 45,872 representative VSGs sequences.
  * Each cluster/group is labeled as known (kVSG) or unknown (uVSG) depending on the presence of at least a viral RefSeq reference genome within the cluster/group.
### New features
* [MetaPhlAn] The new `--profile_vsc` parameter (together with  `--vsc_out` and  `--vsc_breadth`) enables the profiling of viral sequence clusters.
* [MetaPhlAn] The `--subsampling` now subsamples the FASTQ files and not the mapping results
* [MetaPhlAn] The new `--mapping_subsampling` parameter enables the previous mapping subsampling behaviour
* [MetaPhlAn] The new `--subsampling_output` parameter enables to save the subsampled FASTQ file
* [MetaPhlAn] The new `create_toy_database.py` script enables the custom filtering of the MetaPhlAn databases
### Changed features
* [MetaPhlAn] The average read length is included in the output header with the -t rel_ab_w_read_stats parameter
* [StrainPhlAn] Quasi-markers behaviour in line with that of MetaPhlAn
* [StrainPhlAn] sample2markers.py output is now in JSON format
* [StrainPhlAn] Simplified sample and marker filtering parameters, integrated with primary/secondary samples
* [StrainPhlAn] Faster inference of small and medium phylogenies
* [StrainPhlAn] Faster execution of the parameter `–-print_clades_only`
<br/>

## Version 4.0.6 (Mar 1st, 2023)
### Changed features
* [MetaPhlAn] The GTDB taxonomic assignment for the vOct22 database is now available.

<br/>

## Version 4.0.5 (Feb 23rd, 2023)
### Database updates
* We just released the new vOct22 database
* Addition of ~200k new genomes
* 3,580 more SGBs than the vJan21
* 2,548 genomes considered reference genomes in vJan21 were relabelled as MAGs in NCBI -> 1,550 kSGBs in vJan21 are now uSGBs in vOct22
* Removed redundant reference genomes from the vJan21 genomic database using a MASH distance threshold at 0.1%
* Local reclustering to improve SGB definitions of oversized or too-close SGBs
* Improved GGB and FGB definitions by reclustering SGB centroids from scratch
* Improved phylum assignment of SGBs with no reference genomes at FGB level using MASH distances on amino acids to find the closest kSGB
### Changed features
* [StrainPhlAn] Improved StrainPhlAn's speed when running with the --print_clades_only option
### Missing features
* [MetaPhlAn] The GTDB taxonomic assignment for the vOct22 database is not available yet (expected release: end of Feb 2023)
* [MetaPhlAn] The phylogenetic tree of life for the vOct22 database is not available yet (expected release: TBD)

<br/>

## Version 4.0.4 (Jan 17th, 2023)
### Changed features
* [MetaPhlAn] Download of the pre-computed Bowtie2 database is now the default option during installation
* [StrainPhlAn] Improved StrainPhlAn's sample2makers.py script performance and speed
### Fixes
* [StrainPhlAn] Fixes error when using --abs_n_samples_threshold in the PhyloPhlAn call

<br/>

## Version 4.0.3 (Oct 24th, 2022)
### Changed features
* [MetaPhlAn] Removal of the NCBI taxID from the merged profiles produced by the `merge_metaphlan_profiles.py` script
* [StrainPhlAn] Improved StrainPhlAn's performance in the markers/samples filtering step
### Fixes
* [MetaPhlAn] `-t rel_ab_w_read_stats` now produces the reads stats also at the SGB level
* [MetaPhlAn] Fixes overstimation of reads aligned to known clades
* [MetaPhlAn] Fixes error when not providing the number of reads using SAM files as input
* [StrainPhlAn] Fixes `No markers were found for the clade` error while executing StrainPhlAn without providing the clade markers FASTA file

<br/>

## Version 4.0.2 (Sep 22nd, 2022)
### New features
* [MetaPhlAn] The new `--subsampling` parameter allows reads' subsampling on the flight
* [MetaPhlAn] The new `--subsampling_seed` parameter enables a deterministic or randomized subsampling of the reads
* [MetaPhlAn] The new `--gtdb_profiles` of the `merge_metaphlan_profiles.py` allows the merge of GTDB-based MetaPhlAn profiles
* [StrainPhlAn] The new `--breadth_thres` parameter allows StrainPhlAn to filter the consensus markers sequences after the execution of `sample2markers.py`
* [StrainPhlAn] Interactive selection of the available SGBs when the clade is specified at the species level
* [StrainPhlAn] The new `--non_interactive` parameter disables user interaction when running StrainPhlAn 
* [StrainPhlAn] The new `--abs_n_markers_thres` and `--abs_n_samples_thres` parameters enables the specification of the samples/markers filtering thresholds in absolute numbers 
* [StrainPhlAn] The new `--treeshrink` parameter enables StrainPhlAn to run TreeShrink for outlier removal in the tree 
* [StrainPhlAn] Addition of the `VallesColomerM_2022_Jan21_thresholds.tsv` for compatibility with the mpa_vJan21 database
* [StrainPhlAn] The new `--clades` parameter enables `sample2markers.py` to restrict the reconstruction of markers to the specified clades

### Changed features
* [StrainPhlAn] The `-c` parameter of the `extract_markers.py` script now allows the specification of multiple clades
* [StrainPhlAn] The `--print_clades_only` parameter now produces an output `print_clades_only.tsv` report
* [StrainPhlAn] Compatibility with clade markers compressed in bz2 format
* [StrainPhlAn] The `strain_transmission.py` script now uses by the default the `VallesColomerM_2022_Jan21_thresholds.tsv` thresholds
### Fixes
* [MetaPhlAn] `metaphlan2krona.py` and `hclust2` have been added to the bioconda recipe

<br/>

## Version 4.0.1 (Aug 25th, 2022)
### New features
* [MetaPhlAn] The new `--offline` parameter stops MetaPhlAn from automatically checking for updates
### Changed features
* [StrainPhlAn] Improved StrainPhlAn's gaps management with the newest version of PhyloPhlAn (version 3.0.3)
* [StrainPhlAn] Improved set of colors for the `plot_tree_graphlan.py script`
### Fixes
* [MetaPhlAn] Fixes `KeyError: 't'` error when running MetaPhlAn with the `--CAMI_format_output` parameter

<br/>

## Version 4.0.0 (Aug 22nd, 2022)
### Database updates
* Adoption of the species-level genome bins system (SGBs)
* New MetaPhlAn marker genes extracted identified from ~1M microbial genomes
* Ability to profile 21,978 known (kSGBs) and 4,992 unknown (uSGBs) microbial species
* Better representation of, not only the human gut microbiome but also many other animal and ecological environments
### New features
* [MetaPhlAn] Compatibility with MetaPhlAn 3 databases with parameter `--mpa3`
### Changed features
* [MetaPhlAn] Estimation of metagenome composed by microbes not included in the database with parameter `--unclassified_estimation`

<br/>


## Version 3.1.0 (Jul 25th, 2022)
### Database updates
* 433 low-quality species were removed from the MetaPhlAn 3.1 marker database and 2,680 species were added (for a new total of 15,766; a 17% increase)
* Marker genes for a subset of existing bioBakery 3 species were also revised
* Most existing bioBakery 3 species pangenomes were updated with revised or expanded gene content
### Changed features
* [MetaPhlAn] MetaPhlAn 3.1 software has been updated to work with revised marker database

<br/>


## Version 3.0
* New MetaPhlAn marker genes extracted with a newer version of ChocoPhlAn based on UniRef 
* Estimation of metagenome composed by unknown microbes with parameter `--unknown_estimation`
* Automatic retrieval and installation of the latest MetaPhlAn database  with parameter `--index latest`
* Virus profiling with `--add_viruses`
* Calculation of metagenome size for improved estimation of reads mapped to a given clade
* Inclusion of NCBI taxonomy ID in the ouput file
* CAMI (Taxonomic) Profiling Output Format included
* Removal of reads with low MAPQ values
## Version 2.2.0
- added option "marker_counts" (by Nicola)
## Version 2.1.0
- added min_alignment_len option to filter out short alignments in local mode. For long reads (>150) it is now recommended to use local mapping together with "--min_alignment_len 100" to filter out very short alignments. (by Tin)
- added "--samout" option to store the mapping file in SAM format (the SAM will be compressed if the extension of the specified output file ends with ".bz2") (by Tin)
- fix: MetaPhlAn2 now ingores about ~300 markers that were a-specific (thanks to Eric)
## Version 2.0.0
- fix: Biom >= 2.0.0 has the clade IDs second and the sample ids third'
- added extract_markers.py
- fix: #5; revamp biom generation; set clade IDs as enumeration
- added utils/metaphlan2krona.py
