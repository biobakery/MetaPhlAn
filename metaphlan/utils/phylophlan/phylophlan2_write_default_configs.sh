#!/bin/bash


# This script assumes that PhyloPhlAn2 is installed and the commands are available
# in the command line, if not, the first commented row is the execution of the
# phylophlan2_write_config_file.py script from the local PhyloPhlAn2 folder.

# supermatrix_nt.cfg
# ./phylophlan2_write_config_file.py -o supermatrix_nt.cfg \
phylophlan2_write_config_file.py -o supermatrix_nt.cfg \
    -d n \
    --db_dna makeblastdb \
    --map_dna blastn \
    --msa mafft \
    --trim trimal \
    --tree1 fasttree \
    --tree2 raxml \
    --overwrite \
    --verbose

# supertree_nt.cfg
# ./phylophlan2_write_config_file.py -o supertree_nt.cfg \
phylophlan2_write_config_file.py -o supertree_nt.cfg \
    -d n \
    --db_dna makeblastdb \
    --map_dna blastn \
    --msa mafft \
    --trim trimal \
    --gene_tree1 fasttree \
    --gene_tree2 raxml \
    --tree1 astral \
    --overwrite \
    --verbose

# supermatrix_aa.cfg
# ./phylophlan2_write_config_file.py -o supermatrix_aa.cfg \
phylophlan2_write_config_file.py -o supermatrix_aa.cfg \
    -d a \
    --db_aa diamond \
    --map_dna diamond \
    --map_aa diamond \
    --msa mafft \
    --trim trimal \
    --tree1 fasttree \
    --tree2 raxml \
    --overwrite \
    --verbose

# supertree_aa.cfg
# ./phylophlan2_write_config_file.py -o supertree_aa.cfg \
phylophlan2_write_config_file.py -o supertree_aa.cfg \
    -d a \
    --db_aa diamond \
    --map_dna diamond \
    --map_aa diamond \
    --msa mafft \
    --trim trimal \
    --gene_tree1 fasttree \
    --gene_tree2 raxml \
    --tree1 astral \
    --overwrite \
    --verbose
