library(Biostrings)
library(reshape2)
library(plyr)
library(stringr)

fasta_counter<-function(directory){
  setwd(directory)
  files<-list.files(pattern = "fasta")
  all_files<-list()
  for (i in 1:length(files)){
    filepath<-paste0(directory,"/",files[i])
    fasta_file<- filepath
    fasta = readDNAStringSet(fasta_file)
    trans<-fasta@ranges@NAMES
    trans_expanded<-colsplit(trans,"\\|",c("read","id_1","ref","ncbi_id","name"))
    trans_expanded$ncbi_root<-gsub(".*[_]([^.]+)[.].*", "\\1", trans_expanded$ncbi_id)
    trans_expanded$ncbi_root<-paste0(substr(trans_expanded$ncbi_id,1,3),trans_expanded$ncbi_root)
    summary<-ddply(trans_expanded,"ncbi_root",summarise, n=length(ncbi_root))
    print(summary(summary))
    print(i)
    print(filepath)
    all_files[[i]]<-summary
  }
  return(all_files)
}

k<-fasta_counter("~/Desktop/simulated_reads1/fasta")
all<-do.call(cbind,k)
all_files<-all[ ,-grep("ncbi_root",colnames(all))]
final<-cbind(all$ncbi_root,all_files)
names(final)<-c("id",list.files(pattern = "fasta"))
setwd("~/Desktop")
write.csv(final,"ground_truth_counts.csv",row.names = FALSE)

