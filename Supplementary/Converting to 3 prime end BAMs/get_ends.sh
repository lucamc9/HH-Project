#!/bin/bash
# This script takes a reference genome and converts all .bam files in the directory
# into alternate .bam files with only those reads in the region bounded by the
# end of transcript and 1000 base pairs 

# Luca G.Mc s1442231
 
mkdir end_transcripts

for bam in *.bam
do
 i=0
 j=0
 mkdir end_transcripts/tss
 mkdir end_transcripts/tss/sub_tss

 while IFS=$'\t' read -r -a gene
 do
 
 # Create subfolders of max 1000 bam files and merge, to avoid merge's max limit of 1024 files
 if [ "$i" = "1000" ]; then
  samtools merge end_transcripts/tss/merged_$j.bam end_transcripts/tss/sub_tss/*.bam
  rm -rf end_transcripts/tss/sub_tss
  mkdir end_transcripts/tss/sub_tss
  i=0
  j=$((j + 1))
  echo "merged sub_tss"
 fi 

 # Trim elements from white spaces

 chr="$(echo -e "${gene[0]}" | sed -e 's/[[:space:]]*$//')"
 loc1="$(echo -e "${gene[1]}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
 loc2="$(echo -e "${gene[2]}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
 name=${gene[3]}
 strand="$(echo -e "${gene[4]}" | sed -e 's/^[[:space:]]*//')"
 colon=":"
 slash="-"
 
 if [ "$strand" = "+" ]; then
  echo "$bam: ${gene[3]} $strand"
  i=$((i + 1))
  # Construct region <END-1000, END>
  end=$loc2
  begin=$((end - 1000))
  region=$chr$colon$begin$slash$end
  echo "$region"
  samtools view -b -1 -o end_transcripts/tss/sub_tss/ts_$i.bam $bam $region

 else 
  echo "$bam: ${gene[3]} $strand"
  i=$((i + 1))
  # Construct region <BEGIN, BEGIN+1000>
  begin=$loc1
  end=$((begin + 1000))
  region=$chr$colon$begin$slash$end
  echo "$region"
  samtools view -b -1 -o end_transcripts/tss/sub_tss/ts_$i.bam $bam $region
 
 fi 
 done < mm10_final.bed

 # Merge all files in tss
 samtools merge end_transcripts/end_$bam end_transcripts/tss/*.bam
 rm -rf end_transcripts/tss

done
 

