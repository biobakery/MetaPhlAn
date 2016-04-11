#!/bin/bash
python ../mpa3src/build_tree_single_strain.py --ifn_alignments output/s__Bacteroides_caccae.fasta --nprocs 10 --log_ofn output/build_tree_single_strain.log
python ../mpa3src/add_metadata.py --ifn_trees output/RAxML_bestTree.s__Bacteroides_caccae.remove_multiple_strains.tree --ifn_metadatas fastqs/metadata.txt --metadatas subjectID
