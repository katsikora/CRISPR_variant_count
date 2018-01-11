#run in R-3.3.1
#set working directory
wdir<-commandArgs(trailingOnly=TRUE)[1]
#system(paste0('mkdir -p ',wdir)) #for debugging
setwd(wdir)
message(sprintf("working directory is %s",getwd()))

fx<-commandArgs(trailingOnly=TRUE)[2]
message(sprintf("Processing %s",fx))

require(data.table)
require(seqinr)
require(dplyr)

fxdir<-dir(fx,pattern="*.whitelist.fasta",full.names=TRUE)
fxshort<-gsub(".whitelist.fasta","",basename(fxdir))

#####start from bed file and cluster CIGAR strings
bedIN<-dir(sub("filtered_sequences","bams",basename(fxdir)),pattern=".sorted.bed",full.names=TRUE)
bedshort<-gsub(".sorted.bed","",basename(bedIN))
cluLF<-vector("list",length(bedIN))



for(i in seq_along(bedL)){

    bi<-fread(bedIN[i],header=FALSE,sep="\t")
    colnames(bi)<-c("CHR","START","END","ReadName","MAPQ","Strand","CIGAR")
    bi$Length<-bi$END-bi$START
    bi.filt<-bi[bi$Length==names(which.max(table(bi$Length[grep("S",bi$CIGAR,invert=TRUE)]))),]
    bi.filt2<-bi.filt[bi.filt$CHR %in% "X",]
    bi.filt3<-bi.filt2[grep("S",bi.filt2$CIGAR,invert=TRUE),]
    bi.filt4<-bi.filt3[grep("H",bi.filt2$CIGAR,invert=TRUE),]
    tabcig<-table(bi.filt4$CIGAR)
    bi.filt5<-bi.filt4[bi.filt4$CIGAR %in% names(tabcig[as.numeric(tabcig)>=20]),]
    tabi<-as.data.frame(table(bi.filt5$CIGAR),stringsAsFactors=FALSE)
    tabi$SampleID<-bedshort[i]
    tabi$Seq<-"NA"
    tabi$SeqAb<-"NA"
    
    fi<-read.fasta(fxdir[i],seqtype="DNA",as.string=TRUE,forceDNAtolower=FALSE)
    cu<-unique(tabi$Var1)
    for(k in seq_along(cu)){
        rnk<-bi.filt5$ReadName[bi.filt5$CIGAR %in% cu[k]]
        fal<-fi[names(fi) %in% rnk]
        tabseq<-table(unlist(fal))
        tabfilt<-tabseq[tabseq>10]
        if(length(tabfilt)>0){
            tabi$Seq[k]<-paste(names(tabfilt),collapse=";")
            tabi$SeqAb[k]<-paste(as.character(tabfilt),collapse=";")
        }
    }
    cluLF[[i]]<-tabi
    print(paste0(i,"_processed"))
}

save(cluLF,file="cluLF.RData")

cludat<-rbindlist(cluLF)
save(cludat,file="cludat.RData")
write.table(cludat,file="Cluster.counts.txt",row.names=FALSE,sep="\t",quote=FALSE)

subclu<-unlist(lapply(strsplit(cludat$SeqAb,split=";"),length))

cigrv<-rep(cludat$Var1,subclu)
freqrv<-rep(cludat$Freq,subclu)
sidrv<-rep(cludat$SampleID,subclu)
seqrv<-unlist(lapply(cludat$Seq,function(X)strsplit(X,split=";")))
abrv<-unlist(lapply(cludat$SeqAb,function(X)strsplit(X,split=";")))

cludat2<-data.table(cigrv,freqrv,sidrv,seqrv,abrv)
colnames(cludat2)<-c("CIGAR","Freq","SampleID","Seq","SeqAb")
cludat2<-cludat2[!(cludat2$Seq=="NA"),]
cludat2$CluID<-paste(cludat2$CIGAR,cludat2$Seq,sep="_")
save(cludat2,file="cludat2.RData")

clu2sum<-as.data.frame(table(cludat2$SampleID),stringsAsFactors=FALSE)
write.table(clu2sum,file="Variant.counts.txt",row.names=FALSE,sep="\t",quote=FALSE)

##produce upset plot

require(UpSetR)
require(ggplot2)

cluL2<-vector("list",length(unique(cludat2$SampleID)))
names(cluL2)<-unique(cludat2$SampleID)

for(i in seq_along(cluL2)){
    cluL2[[i]]<-cludat2$CluID[cludat2$SampleID %in% unique(cludat2$SampleID)[i]]
    print(paste0(i,"_processed"))

}


png("Upset.plot.clufilt.CIGAR.png",height=800,width=1200)
upset(fromList(cluL2),nsets=length(cluL2),nintersects=NA, order.by = "freq",matrix.color="darkslategrey",main.bar.color="darkslategrey",mainbar.y.label = "Amplicon Seq Intersections", sets.x.label = "Variants per sample")
dev.off()
