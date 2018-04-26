
library(Rsamtools)
library(Rsubread)

################feature count of genome aligned samples##############
setwd("~/Desktop/genomeAligned")
filenames<-list.files(pattern = ".bam")


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

