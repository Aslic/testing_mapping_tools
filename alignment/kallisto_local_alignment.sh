#!/bin/bash

for fn in /logs/*;
do
	samp=`basename ${fn}`
	kallisto quant -i /kallisto/human_transcripts.idx -o "/Desktop/${samp}_kallisto" --single -b 100 -l 250 -s 25 "${fn}"
done

