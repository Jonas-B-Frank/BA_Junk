#! /usr/bin/bash

snvs=../snv_files
segments=../segment_files
dest=../snvs_with_segments

for snv_file in $snvs/*; do
    pid=$(echo $snv_file | tr -dc '0-9')

    cnv_file="${segments}/segments_${pid}.txt"
    if test -f "$cnv_file"; then
        echo "Intersect on PID ${pid}"
        intersectBed -wa -wb -a $snv_file -b $cnv_file | gzip > ${dest}/snvs_with_segments_${pid}.vcf.gz
    else 
        echo "No segments file for PID ${pid}"
    fi
done