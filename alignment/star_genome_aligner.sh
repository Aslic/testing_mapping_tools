#!/bin/bash
#SBATCH --partition=interactive --time=08:00:00 --output=genome_output_star.txt
#SBATCH -N 1
#SBATCH --mem=80gb

cd /scratch/apetenka/logs/
module load rnastar
for fn in /scratch/apetenka/logs/*.fq;
do
filename=$(basename -- "$fn")
STAR \
--runThreadN 20 --twopassMode Basic --runMode alignReads --genomeDir /scratch/apetenka/genomeref/ \
--sjdbGTFfile /scratch/apetenka/genomeref/Homo_sapiens.GRCh38.92.chr.gtf   --sjdbOverhang 100 \
--readFilesIn "${fn}" \
--outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix "/scratch/apetenka/genomeAligned/${filename}"
done
