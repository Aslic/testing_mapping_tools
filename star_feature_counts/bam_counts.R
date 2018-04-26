source("https://bioconductor.org/biocLite.R")
#biocLite("Rsamtools")
#biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#biocLite("BiocParallel")
#biocLite("GenomicAlignments")
#biocLite("DESeq2")
set.seed(1234)
#library(BiocParallel)
#library(GenomicAlignments)
library(Rsamtools)
#library(stringr)
#library(magrittr)
#library(DESeq2)
#library(plyr)
library(Rsubread)

################feature count of genome aligned samples##############
setwd("~/Desktop/genomeAligned")
filenames<-list.files(pattern = ".bam")

#use Rsubread for exon counting and then grouping them by transcript id
gtffile <- "/Users/apetenkaya/Desktop/reference_files/Homo_sapiens.GRCh38.92.chr.gtf"
fc <- featureCounts(files=filenames, 
                    annot.ext=gtffile, 
                    isGTFAnnotationFile=TRUE,
                    isPairedEnd=FALSE,
                    GTF.featureType="exon",
                    GTF.attrType="transcript_id",
                    allowMultiOverlap = FALSE,
                    countMultiMappingReads = TRUE)

genome_counts<-fc$counts
genome_counts<-genome_counts[!rowSums(genome_counts)==0,]
setwd("~/Desktop")
write.csv(genome_counts,"genome_al_star_counts.csv")

#use Rsubread for exon counting and then grouping them by gene id

fc1 <- featureCounts(files=filenames, 
                    annot.ext=gtffile, 
                    isGTFAnnotationFile=TRUE,
                    isPairedEnd=FALSE,
                    GTF.featureType="exon",
                    GTF.attrType="gene_id",
                    allowMultiOverlap = FALSE,
                    countMultiMappingReads = TRUE)

genome_counts<-fc1$counts
genome_counts<-genome_counts[!rowSums(genome_counts)==0,]
setwd("~/Desktop")
write.csv(genome_counts,"genome_al_star_gene_counts.csv")

