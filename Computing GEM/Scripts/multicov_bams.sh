#!/bin/bash
for bam in `ls *.bam | sed 's/.bam//' `
do
bedtools multicov -D -q 255 -bams ${bam}.bam -bed mm10.canonical.bed \
| cut -f4,7 > ${bam}.tmpcnt & 
done
