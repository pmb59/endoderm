## Script to detect sequence-specific biases in ATAC-seq ATAC17.bam

#for all chrs

gc(reset=TRUE)

G <- read.table("GRCh38_15.chrom.sizes", head=F)

bam="ATAC17.bam"
chrN = as.character(G$V1) #"chr22"  #c("chr21","chr22")
chrL = as.numeric(G$V2)   #50818468 #c(46709983, 50818468)
#chrL=900000
library(Rsamtools)

j = 1 # A single bam file

for (i in 1:length(chrN) ) {

  which <- GRanges(seqnames = chrN[i], ranges = IRanges(c(1), c(chrL[i])))

  what <- c("rname", "strand", "pos", "qwidth")


  param <- ScanBamParam(which = which, what = what)
  SampleBam <- scanBam(file=bam[j], param=param)
  sapply(SampleBam[[1]], class);


  SamplelstPlus <- lapply(names(SampleBam[[1]]), function(elt) { do.call(c, unname(lapply(SampleBam, "[[", elt)))})
  names(SamplelstPlus) <- names(SampleBam[[1]])
  SampledfPlus <- do.call("DataFrame", SamplelstPlus)


  SampledfPlus$pos[which(SampledfPlus$strand==1)] <- SampledfPlus$pos[which(SampledfPlus$strand==1)] -1 + 4 #SEQ starts in -1, +4 bp ATAC

  SampledfPlus$qwidth <- 1

  SampledfPlus <- SampledfPlus[which(SampledfPlus$strand==1),]


  plus <- data.frame(chr=chrN[i], start=SampledfPlus$pos-3, end=SampledfPlus$pos+3, forward="forward", one=1, strand="+")  


  options(scipen=100)
  write.table(x=plus[which(plus$start>0 & plus$end<=chrL[i]),], file =paste(chrN[i], "pos_Tn5_forward.bed",sep="_"), append = FALSE, quote = F, sep = "\t", row.names = F,col.names = F)
  options(scipen=0)

  rm(SampledfPlus)

          
  SamplelstMinus <- lapply(names(SampleBam[[1]]), function(elt) { do.call(c, unname(lapply(SampleBam, "[[", elt)))})
  names(SamplelstMinus) <- names(SampleBam[[1]])
  SampledfMinus  <- do.call("DataFrame", SamplelstMinus )



  SampledfMinus$pos[which(SampledfMinus$strand==2)] <- SampledfMinus$pos[which(SampledfMinus$strand==2)] -1 + SampledfMinus$qwidth[which(SampledfMinus$strand==2)] -4 #- 5 #+  1, -5 bp ATAC
  SampledfMinus$qwidth <- 1

  SampledfMinus <- SampledfMinus[which(SampledfMinus$strand==2),]

  minus <- data.frame(chr=chrN[i], start=SampledfMinus$pos-3, end=SampledfMinus$pos+3, forward="reverse", one=1, strand="-")


  options(scipen=100)
  write.table(x=minus[which(minus$start>0 & minus$end<=chrL[i]),], file = paste(chrN[i], "pos_Tn5_reverse.bed", sep="_"), append = FALSE, quote = F, sep = "\t", row.names = F,col.names = F)
  options(scipen=0)

  rm(SampledfMinus)


}




