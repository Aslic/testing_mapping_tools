#!/bin/bash

for fn in simulated_reads1/fasta/*.fasta;
do
samp=`basename ${fn}`
echo "${samp}"
echo "${fn}"
	perl /fasta_to_fastq.pl "${fn}" > "simulated_reads1/logs/converted_${samp}.fq"
done
