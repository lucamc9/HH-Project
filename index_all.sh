#!/bin/bash

for x in *.bam
do 
echo "Indexing: $x"
samtools index $x
done
