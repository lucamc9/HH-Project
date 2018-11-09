#!/bin/bash
awk '$0 !~ /#/{arr[$1]=arr[$1] " " $2}END{for(i in arr)print i,arr[i]}' *count.txt | tr ' '  '\t' | cut -f1,3- > matrix.tbl
