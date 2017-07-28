#!/bin/bash
mkdir -p db_markers
bowtie2-inspect ../db_v20/mpa_v20_m200 > db_markers/all_markers.fasta
python ../strainphlan_src/extract_markers.py --mpa_pkl ../db_v20/mpa_v20_m200.pkl --ifn_markers db_markers/all_markers.fasta --clade s__Bacteroides_caccae --ofn_markers db_markers/s__Bacteroides_caccae.markers.fasta
