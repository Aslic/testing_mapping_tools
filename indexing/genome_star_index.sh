#!/bin/bash
#SBATCH --partition=interactive --time=05:00:00 --output=star_genome_index.txt
#SBATCH -N 1
#SBATCH --mem=80gb

cd /genomeref/
module load rnastar

STAR \
--runThreadN 32 --runMode genomeGenerate --genomeDir /genomeref/ \
--genomeFastaFiles /genomeref/Homo_sapiens.GRCh38.dna.primary_assembly.fa


