args <- commandArgs(trailingOnly = TRUE)
#resultpath <- as.character(args[3])
#setwd(resultpath)

#set size of promoter region relative to the TSS.
up<-as.numeric(args[2])
down<-as.numeric(args[3])

gc<-read.delim("/mnt/work1/users/lupiengroup/Paul/resources/gencodev19.txt",stringsAsFactors=F)
genes<-read.delim(as.character(args[1]),header=F,stringsAsFactors=F)


a1<-gc[which(gc$name2 %in% genes[,1]),]
a2<-data.frame(Loc=character(nrow(a1)),Strand=a1$strand,Name=character(nrow(a1)), stringsAsFactors=F)

for (i in 1:nrow(a1)){
  if (a1$strand[i]=="+"){
    a2$Loc[i]<-paste(a1$chrom[i],a1$txStart[i]-up,a1$txStart[i]+down,sep="\t")
    a2$Name[i]<-paste(a1$name2[i],"_",a1$X.name[i],sep="")
  }
  else if (a1$strand[i]=="-"){
    a2$Loc[i]<-paste(a1$chrom[i],a1$txEnd[i]-down,a1$txEnd[i]+up,sep="\t")
    a2$Name[i]<-paste(a1$name2[i],"_",a1$X.name[i],sep="")
  }
}

write.table(a2,as.character(args[4]),sep="\t",col.names=F,row.names=F,quote=F)
