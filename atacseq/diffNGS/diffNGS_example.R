rm(list=objects())
graphics.off()


# Get regions of open chromatin (atac-seq) changing significantly between condition 1 (C1) and condition 2 (C2) 
args<-commandArgs(TRUE)


#set working dir
#setwd("/Users/pm12/Desktop/EPIGENODE-manuscript/diffNGS")
setwd("/lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC2/diffNGS")

### Parameters of the Analysis
  adjP <- 1e-4  #1e-6  ( Bonferroni adj P-value cutoff )  
  AbsFc <- 2
  #Fl <- 750
  Nbasis <- 10
  Bins <- 31
  CNDS  <- c(as.character(args[2]),as.character(args[3]),as.character(args[4]),as.character(args[5]) )  #c(rep("EG1",2), rep("h12",2) )
  bws   <- c(as.character(args[6]),as.character(args[7]),as.character(args[8]),as.character(args[9]) )  #c( "ATAC1.norm.bw", "ATAC2.norm.bw", "ATAC9.norm.bw", "ATAC10.norm.bw" )   ### Normalized bigwig files
  peaks <- as.character(args[1])  #"temp.bed"
  
 
### Load the package or install if not present
  library(fda)
  library(ICSNP)
  library(genomation)
  library(GenomicRanges)
  library(RColorBrewer)

#  if (!require("fda"))      { install.packages("fda") } 
#  if (!require("ICSNP"))    { install.packages("ICSNP") } 
#  #BioC
#  source("http://bioconductor.org/biocLite.R")
#  if (!require("genomation"))     { biocLite("genomation")  } 
#  if (!require("GenomicRanges"))  { biocLite("GenomicRanges")  } 


### Read bed files with location of aggregated peaks
  bed <- readGeneric(peaks, keep.all.metadata = FALSE) #TRUE  #[1:REGIONS,]
  head(bed)
  length(bed)

#Create data.frame structure to store results
results <- data.frame(region.chr=seqnames(bed), region.start=start(bed) ,region.end=end(bed) , condition_C2vsC1=  paste(rev(unique(CNDS)), collapse="_vs_"), avg.C2= rep(NA,length(bed)), avg.C1= rep(NA,length(bed)), pval=rep(NA,length(bed)),  Bonferroni.pval=rep(NA,length(bed)),  fc= rep(NA,length(bed)) ,  diff= rep(FALSE,length(bed)), diff.score= rep(NA,length(bed))     )
head(results)


source("diffNGS.R")  # Modified function
x <- diffNGS(bedFile= peaks , headerBed=FALSE, bigwigs=bws, conditions=CNDS, pcs = 1, variation = 0.01, nbasis=Nbasis, NB=Bins)


source("processDiffNGS.R")
processDiffNGS(rawPvals = x$p.values, results=results, CNDS = CNDS, AbsFc = AbsFc, adjP = adjP, plotpdf=TRUE, Xmin=-10, Xmax=10, GTF="genes_ensembl_76_transcriptome-GRCh38_15.gtf", W=10e3 )




# end # --------------------------------------------------------------------------------




