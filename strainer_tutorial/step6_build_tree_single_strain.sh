#!/bin/bash
python ../strainer_src/build_tree_single_strain.py --ifn_alignments output/s__Bacteroides_caccae.fasta --nprocs 10 --log_ofn output/build_tree_single_strain.log
python ../strainer_src/add_metadata_tree.py --ifn_trees output/RAxML_bestTree.s__Bacteroides_caccae.remove_multiple_strains.tree --ifn_metadatas fastqs/metadata.txt --metadatas subjectID
