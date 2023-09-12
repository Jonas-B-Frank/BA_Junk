#! /usr/bin/bash

snvs_with_segments=../j/$1
onek1k=../onek1k/onek1k_subset_columns.tsv
dest_extra_columns=../snvs_with_segments_extra_columns
dest_joined=../joined_mmml_onek1k

for file in $snvs_with_segments/*; do
    pid=$(echo $file | tr -dc '0-9')

    echo "Add columns on PID ${pid}"
    python3 add_columns.py -f $file -o ${dest_extra_columns}/extra_columns_${pid}.vcf.gz
    
    # echo "Intersect on PID ${pid}"
    # intersectBed -wa -wb -a ${dest_extra_columns}/extra_columns_${pid}.vcf.gz -b $onek1k | gzip > ${dest_joined}/mmmml_onek1k_${pid}.vcf.gz

    echo "Done for PID ${pid}"
done