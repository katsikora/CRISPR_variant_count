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
fxshort<-gsub("_prin.whitelist.fasta","",basename(fxdir))

#####start from bed file and cluster CIGAR strings
bedIN<-dir(sub("filtered_sequences","bams",fx),pattern=".sorted.bed",full.names=TRUE)
bedshort<-gsub("_prin.sorted.bed","",basename(bedIN))
cluLF<-vector("list",length(bedIN))



for(i in seq_along(cluLF)){

    bi<-fread(bedIN[i],header=FALSE,sep="\t")
    colnames(bi)<-c("CHR","START","END","ReadName","MAPQ","Strand","CIGAR")
    bi$Length<-bi$END-bi$START
    bi.filt<-bi[bi$Length==names(which.max(table(bi$Length[grep("S",bi$CIGAR,invert=TRUE)]))),]
    bi.filt2<-bi.filt[bi.filt$CHR %in% "X",]
    bi.filt3<-bi.filt2[grep("S",bi.filt2$CIGAR,invert=TRUE),]
    bi.filt4<-bi.filt3[grep("H",bi.filt2$CIGAR,invert=TRUE),]
    tabi<-bi.filt4
    tabi$SampleID<-bedshort[i]
    tabi$Seq<-"NA"
        
    fi<-read.fasta(fxdir[i],seqtype="DNA",as.string=TRUE,forceDNAtolower=FALSE)
    fal<-unlist(fi[names(fi) %in% tabi$ReadName])
    tabi$Seq<-fal[match(tabi$ReadName,names(fal))]
    cluLF[[i]]<-tabi
    print(dim(tabi))
    print(paste0(i,"_processed"))
}

save(cluLF,file="cluLF.RData")

cludat<-rbindlist(cluLF)
save(cludat,file="cludat.RData")


##pattern match
require(stringr)
Seqloc<-str_locate(cludat$Seq, "GAAGGAGAT.+TGTGTGC")
Seqsub<-substr(cludat$Seq,start=Seqloc[,"start"],stop=Seqloc[,"end"]) ###NAs introduced where no substring found
cludat$Seqsub<-Seqsub
save(cludat,file="cludat.RData")

#summarize counts/cluster and filter for min 10 seqs
sumdat<-data.table(summarize(group_by(cludat,SampleID,Seqsub),Count=length(Seqsub)),stringsAsFactors=FALSE)
sumdat$Length<-str_length(sumdat$Seqsub)
##add a filtering step for length equal to WT or WT-1 -> 41 or 40nt
sumdat.filt<-sumdat[sumdat$Count>=10&!(is.na(sumdat$Seqsub))&sumdat$Length!=40&sumdat$Length!=41,]

save(sumdat.filt,file="sumdat.filt.RData")
write.table(sumdat.filt,file="Cluster.counts.txt",row.names=FALSE,sep="\t",quote=FALSE)


vardat<-as.data.frame(summarize(group_by(sumdat.filt,SampleID),NVar=length(Seqsub)),stringsAsFactors=FALSE)
write.table(vardat,file="Variant.counts.txt",row.names=FALSE,sep="\t",quote=FALSE)

##produce upset plot

require(UpSetR)
require(ggplot2)

cluL2<-vector("list",length(unique(sumdat.filt$SampleID)))
names(cluL2)<-unique(sumdat.filt$SampleID)

for(i in seq_along(cluL2)){
    cluL2[[i]]<-sumdat.filt$CluID[sumdat.filt$SampleID %in% unique(sumdat.filt$SampleID)[i]]
    print(paste0(i,"_processed"))

}


png("Upset.plot.clufilt.CIGAR.png",height=800,width=1200)
upset(fromList(cluL2),nsets=length(cluL2),nintersects=NA, order.by = "freq",matrix.color="darkslategrey",main.bar.color="darkslategrey",mainbar.y.label = "Amplicon Seq Intersections", sets.x.label = "Variants per sample")
dev.off()
