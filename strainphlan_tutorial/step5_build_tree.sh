#!/bin/bash
mkdir -p output
python ../strainphlan.py --mpa_pkl ../db_v20/mpa_v20_m200.pkl --ifn_samples consensus_markers/*.markers --ifn_markers db_markers/s__Bacteroides_caccae.markers.fasta --ifn_ref_genomes reference_genomes/G000273725.fna.bz2 --output_dir output --nprocs_main 10 --clades s__Bacteroides_caccae &> output/log.txt
python ../strainphlan_src/add_metadata_tree.py --ifn_trees output/RAxML_bestTree.s__Bacteroides_caccae.tree --ifn_metadatas fastqs/metadata.txt --metadatas subjectID
