#!/bin/bash
cd /salmon-0.8.2-1/bin/
for fn in /input_fastq/logs/*;
do
filename=$(basename -- "$fn")
salmon quant -i /human_trans_index -l A \
-r "${fn}" --fldMean 250 --fldSD 25 --writeUnmappedNames -p 8 \
-o "/${filename}_quant"
done
