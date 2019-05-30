#!/bin/bash
mkdir -p db_markers
bowtie2-inspect ../metaphlan_databases/mpa_v29_CHOCOPhlAn_201901 > db_markers/all_markers.fasta
python2 ../strainphlan_src/extract_markers.py \
        --mpa_pkl ../metaphlan_databases/mpa_v29_CHOCOPhlAn_201901.pkl \
        --ifn_markers db_markers/all_markers.fasta \
        --clade s__Bacteroides_caccae \
        --ofn_markers db_markers/s__Bacteroides_caccae.markers.fasta
