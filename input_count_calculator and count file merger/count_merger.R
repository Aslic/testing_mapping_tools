setwd("/kallisto/kallisto_ch22_bs100/")#give the directory where the read folders are stored
dirs<-list.files()
names<-vector()
files<-list()
for(i in 1:length(dirs)){
  setwd(dirs[i])
  print (getwd())
  names[i]<-dirs[i]
  files[[i]]<-read.table("abundance.tsv",header = TRUE,stringsAsFactors = FALSE,sep="\t")#give the file name
  print(summary(files[[i]]))
  setwd("../")
}

all_files<-do.call(cbind,files)
read_mat<-all_files[,grepl("est_counts",names(all_files))] #give the header where the counts are
read_mat<-cbind(all_files[,1],read_mat)
names(read_mat)<-c("id",names)
informative<-read_mat[rowSums(read_mat[2:ncol(read_mat)])!=0,]
setwd("/Users/apetenkaya/Desktop/")
write.csv(informative,"kallistobs_counts.csv")
