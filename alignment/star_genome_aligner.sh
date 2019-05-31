#!/bin/bash
#SBATCH --partition=interactive --time=08:00:00 --output=genome_output_star.txt
#SBATCH -N 1
#SBATCH --mem=80gb

cd /logs/
module load rnastar
for fn in /logs/*.fq;
do
filename=$(basename -- "$fn")
STAR \
--runThreadN 20 --twopassMode Basic --runMode alignReads --genomeDir /genomeref/ \
--sjdbGTFfile /genomeref/Homo_sapiens.GRCh38.92.chr.gtf   --sjdbOverhang 100 \
--readFilesIn "${fn}" \
--outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix "/genomeAligned/${filename}"
done
