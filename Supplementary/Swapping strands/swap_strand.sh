#!/bin/bash

# A simple shell script to swap the erroneous strands assigned to the genes for the correct ones
# Luca G.Mc s1442231

# Located in the directory must be a BED file with erroneous strands "mm10.canonical.bed" and 
# the one with correct strands "mm10_3utr.bed" but including exons


while IFS=$'\t' read -r -a gene
do
 found=($(grep "${gene[3]}" -m1 mm10_3utr.bed))
 len=${#found[@]}
 if [ "$len" = "6" ]; then 
    echo -e "${gene[0]} \t ${gene[1]} \t ${gene[2]} \t ${gene[3]} \t ${found[5]}" >> mm10_final.bed
 else
    echo -e "${gene[0]} \t ${gene[1]} \t ${gene[2]} \t ${gene[3]} \t ${found[4]}" >> mm10_final.bed
 fi
 echo "found gene: ${gene[3]}"
done < mm10.canonical.bed
