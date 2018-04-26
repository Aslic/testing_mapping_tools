#!/bin/bash

for fn in /Users/apetenkaya/Desktop/simulated_reads1/logs/*;
do
	samp=`basename ${fn}`
	kallisto quant -i /Users/apetenkaya/Desktop/kallisto/human_transcripts.idx -o "/Users/apetenkaya/Desktop/${samp}_kallisto" --single -b 100 -l 250 -s 25 "${fn}"
done

