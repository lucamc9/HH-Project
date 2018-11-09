#!/bin/bash
for bam in `ls *.bam | sed 's/.bam//' `
do 
sort ${bam}.tmpcnt \
| perl agg.pl > ${bam}_count.txt &
done
wait
