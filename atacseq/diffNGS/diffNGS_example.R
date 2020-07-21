# Get regions of open chromatin (ATAC-seq) changing significantly between condition 1 (C1) and condition 2 (C2) 

args<-commandArgs(TRUE)

#set working dir
setwd("../diffNGS")

### Parameters of the Analysis
  adjP <- 1e-4   # ( Bonferroni adj P-value cutoff )  
  AbsFc <- 2
  #Fl <- 750
  Nbasis <- 10
  Bins <- 31
  CNDS  <- c(as.character(args[2]),as.character(args[3]),as.character(args[4]),as.character(args[5]) )  
  bws   <- c(as.character(args[6]),as.character(args[7]),as.character(args[8]),as.character(args[9]) )  
  peaks <- as.character(args[1]) 
  
 
### Load the package or install if not present
  library(fda)
  library(ICSNP)
  library(genomation)
  library(GenomicRanges)
  library(RColorBrewer)

### Read bed files with location of aggregated peaks
  bed <- readGeneric(peaks, keep.all.metadata = FALSE) 
  head(bed)
  length(bed)

#Create data.frame structure to store results
results <- data.frame(region.chr=seqnames(bed), region.start=start(bed) ,region.end=end(bed) , condition_C2vsC1=  paste(rev(unique(CNDS)), collapse="_vs_"), avg.C2= rep(NA,length(bed)), avg.C1= rep(NA,length(bed)), pval=rep(NA,length(bed)),  Bonferroni.pval=rep(NA,length(bed)),  fc= rep(NA,length(bed)) ,  diff= rep(FALSE,length(bed)), diff.score= rep(NA,length(bed))     )
head(results)


source("diffNGS.R")  
x <- diffNGS(bedFile= peaks , headerBed=FALSE, bigwigs=bws, conditions=CNDS, pcs = 1, variation = 0.01, nbasis=Nbasis, NB=Bins)


source("processDiffNGS.R")
processDiffNGS(rawPvals = x$p.values, results=results, CNDS = CNDS, AbsFc = AbsFc, adjP = adjP, plotpdf=TRUE, Xmin=-10, Xmax=10, GTF="genes_ensembl_76_transcriptome-GRCh38_15.gtf", W=10e3 )


