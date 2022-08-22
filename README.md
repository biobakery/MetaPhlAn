# MetaPhlAn: Metagenomic Phylogenetic Analysis
## What's new in version 4
* Adoption of the species-level genome bins system (SGBs, http://segatalab.cibio.unitn.it/data/Pasolli_et_al.html)
* New MetaPhlAn marker genes extracted identified from ~1M microbial genomes
* Ability to profile 21,978 known (kSGBs) and 4,992 unknown (uSGBs) microbial species
* Better representation of, not only the human gut microbiome but also many other animal and ecological environments
* Estimation of metagenome composed by microbes not included in the database with parameter `--unclassified_estimation`
* Compatibility with MetaPhlAn 3 databases with parameter `--mpa3`
-------------

## Description
MetaPhlAn is a computational tool for profiling the composition of microbial communities (Bacteria, Archaea and Eukaryotes) from metagenomic shotgun sequencing data (i.e. not 16S) with species-level. With StrainPhlAn, it is possible to perform accurate strain-level microbial profiling.

MetaPhlAn 4 relies on ~5.1M unique clade-specific marker genes identified from ~1M microbial genomes (~236,600 references and 771,500 metagenomic assembled genomes) spanning 26,970 species-level genome bins (SGBs, http://segatalab.cibio.unitn.it/data/Pasolli_et_al.html), 4,992 of them taxonomically unidentified at the species level, allowing:

* unambiguous taxonomic assignments;
* an accurate estimation of organismal relative abundance;
* SGB-level resolution for bacteria, archaea and eukaryotes;
* strain identification and tracking
* orders of magnitude speedups compared to existing methods.
* metagenomic strain-level population genomics

If you use MetaPhlAn, please cite:

[**Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3**](https://elifesciences.org/articles/65088) *Francesco Beghini, Lauren J McIver, Aitor Blanco-Míguez, Leonard Dubois, Francesco Asnicar, Sagun Maharjan, Ana Mailyan, Paolo Manghi, Matthias Scholz, Andrew Maltez Thomas, Mireia Valles-Colomer, George Weingart, Yancong Zhang, Moreno Zolfo, Curtis Huttenhower, Eric A Franzosa, Nicola Segata*. eLife (2021)

If you use StrainPhlAn, please cite the MetaPhlAn paper and the following StrainPhlAn paper:

[**Microbial strain-level population structure and genetic diversity from metagenomes.**](http://genome.cshlp.org/content/27/4/626.full.pdf) *Duy Tin Truong, Adrian Tett, Edoardo Pasolli, Curtis Huttenhower, & Nicola Segata*. Genome Research 27:626-638 (2017)


-------------

## MetaPhlAn and StrainPhlAn tutorials and resources

In addition to the information on this page, you can refer to the following additional resources.

* The [MetaPhlAn tutorial](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4).

* The [StrainPhlAn tutorial](https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4).

* Related tools including [PanPhlAn](https://github.com/segatalab/panphlan) (and its [tutorial](https://github.com/segatalab/panphlan/wiki/Home)), [GraPhlAn](https://github.com/segatalab/graphlan) (and it [tutorial](https://github.com/biobakery/biobakery/wiki/graphlan)), [PhyloPhlAn 3](https://github.com/biobakery/phylophlan) (and its [tutorial](https://github.com/biobakery/biobakery/wiki/phylophlan)), [HUMAnN](https://github.com/biobakery/humann/) (and its [tutorial](https://github.com/biobakery/biobakery/wiki/humann2)).

* The related [bioBakery workflows](https://github.com/biobakery/biobakery/wiki/biobakery_workflows)
