#!/bin/bash
cd /Users/apetenkaya/anaconda2/pkgs/salmon-0.8.2-1/bin/
for fn in /Users/apetenkaya/Desktop/input_fastq/logs/*;
do
filename=$(basename -- "$fn")
salmon quant -i /Users/apetenkaya/Desktop/human_trans_index -l A \
-r "${fn}" --fldMean 250 --fldSD 25 --writeUnmappedNames -p 8 \
-o "/Users/apetenkaya/Desktop/${filename}_quant"
done