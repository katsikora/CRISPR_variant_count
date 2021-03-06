#run in R-3.3.1
#set working directory
wdir<-commandArgs(trailingOnly=TRUE)[1]
#system(paste0('mkdir -p ',wdir)) #for debugging
setwd(wdir)
message(sprintf("working directory is %s",getwd()))

bedIN<-commandArgs(trailingOnly=TRUE)[2]
message(sprintf("Processing %s",bedIN))

require(data.table)

bi<-fread(bedIN,header=FALSE,sep="\t")
colnames(bi)<-c("CHR","START","END","ReadName","MAPQ","Strand","CIGAR")
bi$Length<-bi$END-bi$START

tgtInt<-commandArgs(trailingOnly=TRUE)[3]
tgt.chrom<-gsub(":.+","",tgtInt)
tgt.start<-gsub(".+:","",gsub("-.+","",tgtInt))
tgt.end<-gsub(".+-","",tgtInt)
tg.len<-as.numeric(tgt.end)-as.numeric(tgt.start)

bi.filt<-bi[(bi$Length>=(tg.len-1) | bi$Length<=(tg.len+1)) & bi$CHR==tgt.chrom,]
ts<-bi.filt$ReadName
write.table(ts,file=gsub(".sorted.bed",".whitelist.txt",bedIN),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)



