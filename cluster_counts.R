#run in R-3.3.1
#set working directory
wdir<-commandArgs(trailingOnly=TRUE)[1]
#system(paste0('mkdir -p ',wdir)) #for debugging
setwd(wdir)
message(sprintf("working directory is %s",getwd()))

fx<-commandArgs(trailingOnly=TRUE)[2]
message(sprintf("Processing %s",fx))

fxdir<-dir(fx,pattern="*.whitelist.dedup.fasta",full.names=TRUE)
fxshort<-gsub(".whitelist.dedup.fasta","",basename(fxdir))

require(seqinr)

fxL<-vector("list",length(fxdir))
names(fxL)<-fxshort
ctL<-vector("list",length(fxdir))

for(i in seq_along(fxdir)){
    fi<-read.fasta(fxdir[i],seqtype="DNA",as.string=TRUE,forceDNAtolower=FALSE)
    fxL[[i]]<-fi
    fn<-names(fi)
    fd<-data.frame(fn,fxshort[i])
    colnames(fd)<-c("CluID","SampleID")
    ctL[[i]]<-fd

}

fxdat<-as.data.frame(do.call(rbind,ctL),stringsAsFactors=FALSE)
fxdat$Cluster<-as.numeric(gsub("\\-.+","",fxdat$CluID))
fxdat$Count<-as.numeric(gsub(".+\\-","",fxdat$CluID))

write.table(fxdat,file="Cluster.counts.txt",sep="\t",quote=FALSE,row.names=FALSE)

##produce upset plot

require(UpSetR)
require(ggplot2)

png("Upset.plot.png",height=800,width=1200)
upset(fromList(fxL),nsets=length(fxL),nintersects=NA, order.by = "freq",matrix.color="darkslategrey",main.bar.color="darkslategrey",mainbar.y.label = "Amplicon Seq Intersections", sets.x.label = "Variants per sample")
dev.off()
