#!/bin/bash

# Get all bam names in directory
bam_list=""
for bam in *.bam
do 
 bam_list="$bam_list $bam"
done

# Execute python script
python -W ignore get_termination_sites.py --bed mm10_final_wn.bed --utr mm10_3utr_wn.bed --name test --bam $bam_list

