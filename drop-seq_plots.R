args = commandArgs(trailingOnly=TRUE)
FQZ<-as.character(strsplit(args, " ")[1])
FQZ
if (FQZ=="NULL") {
  message("Error, no fastq files specified, quitting now")
  q()
}

library("reshape2")

#FQZ="3t3sg4-2k_S1_L001_R1_001.fastq.gz"

PLOTS_FILE=paste(FQZ,".plots.pdf",sep="")
pdf(file=PLOTS_FILE)

#bc count dist
BC_FILENAME=paste(FQZ,".bc.tsv",sep="")
x<-read.table(BC_FILENAME)
colnames(x)=c("index","count","barcode")
plot(x$index,x$count,xlab="Index",ylab="Read count",main="Barcode representation")

#duplicates
DUPL_FILENAME=paste(FQZ,".dupl.txt",sep="")
x<-read.table(DUPL_FILENAME)
NROW=nrow(x)
x<-x[-NROW,]
tail(x)
x<-x[order(-x$V2),] 
x$dup=x$V2-x$V3
plot(x$V2,x$dup,xlab="Total reads",ylab="Duplicate reads",main="Duplicate read content")
plot(x$V2,x$dup/x$V2*100,xlab="Total reads",ylab="Percent duplicates (%)",main="Duplicate read proportion")

#STAR mapping
STARMAP_FILENAME=paste(FQZ,".starmapping.txt",sep="")
y<-read.table(STARMAP_FILENAME)
plot(log2(y$V2),y$V3,xlab="Total reads (log2)",ylab="Percent mapped (%)",main="Proportion mapped with STAR")

#rRNA contaminaton plot
RRNA_FILENAME=paste(FQZ,".rrna.counts.txt",sep="")
x<-read.table(RRNA_FILENAME)
x$rrna=x$V3/x$V2
NROW=nrow(x)
x<-x[-NROW,]
x<-x[order(-x$rrna),]
plot(x$rrna*100,xlab="Index",ylab="Percent mapped (%)",main="Proportion mapped to rRNA with BWA")

#kallisto mapping
KAL_COUNTS=paste(FQZ,".kalcounts.txt",sep="")
x<-read.table(KAL_COUNTS)
NROW=nrow(x)
x<-x[-NROW,]
x$TOT=x$V2+x$V3
x<-x[order(-x$TOT),]
plot(log2(x$V3),x$V4/x$V3*100,xlab="Total reads (log2)",ylab="Percent mapped (%)",main="Proportion mapped with Kallisto")

#kal mx gen
KAL_3COL=paste(FQZ,".kal.3col",sep="")
tmp<-read.table(KAL_3COL)
y<-as.matrix(acast(tmp, V2~V1, value.var="V3"))
MXFILE=paste(FQZ,".kal.mx",sep="")
write.table(y,file=MXFILE,sep="\t",quote=F)

#MDS plot original
yy<-scale(y)
plot(cmdscale(dist(t(yy))), main="MDS plot", xlab="Coordinate 1", ylab="Coordinate 2", type = "n")
text(cmdscale(dist(t(yy))), labels=colnames(yy), cex=0.5, ) 

#MDS plot filter low read barcodes out
yy<-y[,colSums(y)>=1000]
yy<-scale(yy)
plot(cmdscale(dist(t(yy))), main="MDS plot low reads removed", xlab="Coordinate 1", ylab="Coordinate 2", type = "n")
text(cmdscale(dist(t(yy))), labels=colnames(yy), cex=0.5, )

dev.off()

