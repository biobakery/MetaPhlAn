#!/bin/bash
mkdir -p sams
for f in $(ls fastqs/*.bz2)
do
    echo "Running metaphlan2 on ${f}"
    bn=$(basename ${f} | cut -d . -f 1)
    tar xjfO ${f} | ../metaphlan2.py --bowtie2db ../db_v20/mpa_v20_m200 --mpa_pkl ../db_v20/mpa_v20_m200.pkl --input_type multifastq --nproc 10 -s sams/${bn}.sam.bz2 --bowtie2out sams/${bn}.bowtie2_out.bz2 -o sams/${bn}.profile
done
