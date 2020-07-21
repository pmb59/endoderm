#Modified for ATAC21/ATAC22 (36h EPIGENODE Phase II)

## EPIGENODE Project - ATAC-seq, 2016-06-18 ##
## FootprintMixture     ##
## Detection of TF footprints while controlling for Tn5-specific sequence bias ##

# FIMO should be run first (used 50k top-scored motifs) (stringent analysis) 
# For now will use TF_Information_all_motifs.txt (Directly determined or best inferred motif)
# Ideally, select all FIMO ocurrences after a certain cut-off (relaxed analysis)
# Move coordinates +4bp / -5bp
# Generate SeqBias.txt from ATAC17.bam>R script>SeqTools  (so far in chr22)

#graphics.off()
rm(list=objects())

args<-commandArgs(TRUE)

#Supress global warnings
options(warn=-1)

#Remove scientific notation (to avoid problems in Start/End genomic coordinates when using BEDtools)
#options(scipen=999)

library(grid)
library(Rsamtools)
library(seqLogo)
library(Biostrings)
library(png)
library(genomation)
library(GenomicRanges)
library(scales)
library(gtools);
#library(mixtools); #install.packages("mixtools_0.1.0.tar.gz") !
library(RColorBrewer)
library(tictoc)
source("/lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC2/FootprintMixture/MixtureModel.r");
#library(miscTools)
#library(sm)
#library(vioplot)

##THIS WAS TO SELECT ONLY A SUBSET OF ALL MOTIFS IN MY LAPTOP#
#setwd("/lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/fimo/fimo_occurrences_500k");  #Dir with Top 60k FIMO matches
setwd("/lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC/fimo/fimo_occurrences_500k")
##setwd("/lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/fimo/fimo_occurrences"); 
Stored <- unlist(strsplit(list.files()  ,split=".txt"))   #c("M1957_1.02", "M4532_1.02" , "M5943_1.02")
OfInterest <- as.character(read.csv2("/lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC/fimo/CISBP/TF_Information.txt", sep="\t")$Motif_ID)
#which are common in both list?
TF_motifs <- intersect(Stored , OfInterest)


InMyList <- which (TF_motifs == as.character(args[1])  )
if (length(InMyList ) >0 ) {
TF_motifs = as.character( args[1] );


#Bam alignments
Bam_ATAC  <- as.character( args[2] )    #c( "ATAC1" )
#Bam_ATAC  <-  c( "ATAC1", "ATAC2", "ATAC9", "ATAC10", "ATAC11" , "ATAC12",  "ATAC13", "ATAC14",  "ATAC15", "ATAC16", "ATAC17" )

##Bam_ATAC <- c("ATAC1ds_sorted","ATAC2ds_sorted","ATAC3ds_sorted","ATAC4ds_sorted","ATAC5ds_sorted","ATAC6ds_sorted","ATAC7ds_sorted","ATAC8ds_sorted","ATAC9ds_sorted","ATAC10ds_sorted","ATAC11ds_sorted","ATAC12ds_sorted","ATAC13ds_sorted","ATAC14ds_sorted","ATAC15ds_sorted","ATAC16ds_sorted") 
###Bam_ATAC <- c("ATAC1ds_sorted","ATAC2ds_sorted","ATAC9ds_sorted","ATAC10ds_sorted","ATAC11ds_sorted","ATAC12ds_sorted","ATAC13ds_sorted","ATAC14ds_sorted","ATAC15ds_sorted","ATAC16ds_sorted") 

for (q in 1:length(TF_motifs) ){

  
  for (ba in 1:length(Bam_ATAC) ){
    print(paste( paste("#Motif = ",q, sep="") ,paste("#bamFile = ", Bam_ATAC[ba] , sep="") ,sep=" ; ") )
    
    tic();
    #--------------------------------------------------------------------------------------------------
    #PARAMS:
    setwd("/lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC2/FootprintMixture/");  #working directory
    bamFile = Bam_ATAC[ba]   #"ATAC1-2"
    bamPath = "/lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC2/" #mergedBams/"  #"/Users/pm12/Desktop/EPIGENODE-manuscript/FootprintMixture3/"
    
    motif = TF_motifs[q]   #M5943_1.02"   #USF1 "M4509_1.02" #"M1957_1.02"  #"M5808_1.02" 
    Nmotifs = 10e3 # 50e3 as a stringet set (Yardimci et al)  / OR MINIMUM NUMBER OF MOTIFS TO ANALYZE
    NROWS = 1.0*Nmotifs #1.1* 60e3 #Keep only autosomal and sex-chrs matches
    EXT = 25  #35 #25
    chrom_Sizes="GRCh38_15.chrom.sizes" #chromosomes that we actually use
    controlSeqBias="Seq"  # "Flat" or "Seq"
    genome <- "/lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC2/Homo_sapiens.GRCh38_15.fa"
    database <- "/lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC/fimo/CISBP/TF_Information.txt" #TF_Information_all_motifs_plus.txt" 
    ##fimo_out <- "/Users/pm12/Desktop/EPIGENODE-manuscript/FootprintMixture3/fimo_out"
    fimo_txt <- "/lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC/fimo/fimo_occurrences_500k"
#fimo_txt <- "/lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/fimo/fimo_occurrences"
    #resultsPath="/lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/FootprintMixture/Footprints_plots"
    resultsPath="/lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC2/FootprintMixture/Footprints_plots"
    db <- read.csv(database, head=T, sep="\t")
    #head(db)
    #--------------------------------------------------------------------------------------------------
    
    
    # Get the NAME of the TF of interest (paste if more than 1)
    # To avoid FILENAME TOO LONG problem (in few cases), select the first 20 characters
    TF_Name <- substr(  x=paste( as.character(db$TF_Name[which(db$Motif_ID == motif)]) , collapse="_") , start=1, stop= 20 )
    print(paste(TF_Name, motif, sep=" - "))
    
    ## CIS-BP to MEME format for this motif (Farm, Done)
    ##( Rscript cisbp2meme.r)
    
    ## Run FIMO to get matches in the human genome (Farm, Done)
    #fimo.sh
    
    
    fimo_Matches = paste(fimo_txt,"/",motif,".txt",sep="")  #paste(motif,".fimo.txt", sep="")  #Motifs are ranked
    
    ## Peak Calling for DNase-seq using MACS (proved to work well)
    ##New MACS version
    ##1. To find enriched cutting sites such as some DNase-Seq datasets. In this case, all 5' ends of sequenced reads should be extended in both direction to smooth the pileup signals. If the wanted smoothing window is 200bps, then use '--nomodel --shift -100 --extsize 200'.=
    ##Name="NA18486"
    #setwd("/Users/pm12/Desktop/soft/python_local/bin")
    #system(paste( paste("./macs2 callpeak -f BAM -g hs -n ATAC1 --nomodel -q 0.01 --outdir /Users/pm12/Desktop/EPIGENODE-analysis/FootprintMixture2/results -t", bamFile ,sep=" "), " > /Users/pm12/Desktop/EPIGENODE-analysis/FootprintMixture2/results/ATAC1-info.txt", sep=" " ) , intern=FALSE)
    ##show.output.on.console = TRUE
    ##https://github.com/taoliu/MACS --> see paramameters recommende for DNAseq
    
    #Check if we have NROWS motifs in the file
    if(  length( (readLines(fimo_Matches)) )   >=  (NROWS -1) )   {
      print('More than 10k motif occurrences found...Starting Analysis...')
      
      #Filter Out Motif ocurrence out of chr1-22, X,Y
      fimo <- read.table(fimo_Matches, header = FALSE  )   #!!###   , nrows =NROWS)

#	#Assure integre values (no scientific notation)
#	fimo$V3 <- as.integer( fimo$V3 )
#	fimo$V4 <- as.integer( fimo$V4 )

	##########################################
	# check End > Start in FIMO input
	suspic <- which( fimo$V4 < fimo$V3 )
	if ( length(suspic)>0 ) {
		suspicTemp <- fimo$V3[suspic]
		fimo$V3[suspic] <- fimo$V4[suspic]
		fimo$V4[suspic] <- suspicTemp
		rm(suspicTemp)
	}
	rm(suspic) 

	# check End != Start in FIMO input
	suspic <- which( fimo$V4 != fimo$V3 )
	if ( length(suspic)>0 ) {
		fimo <- fimo[suspic , ]
	}
	rm(suspic)
	#########################################

      #read Valid chromosome space, and also check that END motif & START motif (Extended) are within Valid genomic ranges (hg38) 
      chrs <-  read.table(chrom_Sizes, head=FALSE) 
      chrs$V1 <- as.character(chrs$V1)
      
      valid <- c()
      fastaL=3
      for (j in 1:nrow(fimo)) {
        if ( length( which( chrs$V1 ==  as.character( fimo$V2[j]) )) >0     ){ 
		if(  (fimo$V3[j] - EXT - fastaL)  > 0  ){
			AuxChr <- which( chrs$V1 ==  as.character( fimo$V2[j]) );
			if(  (fimo$V4[j] + EXT + fastaL +1)  < chrs$V2[AuxChr]  ){
				valid <- c(valid, j) }  #print(j)   # & to && !!
			}
			rm(AuxChr)
		}
      }
      nrow(fimo)
      length(valid)
      
      
      #Store bed file corresponding to the matches (candidate TFBSs)- used the Valid
      fimo <- fimo[valid, ]; dim(fimo)
#MODIFY THIS TO GET ALL MOTIFS
#      #get only Nmotifs (check if the ones reported by FIMO are of enough number)
#      if (nrow(fimo) >= Nmotifs){
#        #ocurrences <- fimo[1:Nmotifs, ]
#	ocurrences <- fimo #[1:Nmotifs, ]
#        head(ocurrences)
#        dim(ocurrences)
#      }
#      if (nrow(fimo) < Nmotifs){
#        ocurrences <- fimo
#      }
ocurrences <- fimo      
      dim(ocurrences)
      
#	#Assure integre values (no scientific notation)
#	ocurrences$V3 <- as.integer( ocurrences$V3 )
#	ocurrences$V4 <- as.integer( ocurrences$V4 )
      
      
      #50k (Nmotifs) strongest matches
#options(scipen=999)
      write.table(x=ocurrences[,2:4], file = paste(motif,bamFile,".bed.data.txt",sep=""), append = FALSE, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
      write.table(x=ocurrences[,2:7], file = paste(motif,bamFile,".bed.genomation.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
#options(scipen=0)
     print('.bed.data.txt and .bed.genomation.txt... SAVED IN DISK')

      #Motif Length
      Lm <- nchar(x=as.character(ocurrences$V8[1])) #e.g., 19
      
      
      fastaOut <-paste(motif,bamFile,".fa", sep="") 
      fastaL=3
      if (controlSeqBias == "Seq") {
        #Store fasta file corresponding to the matches from the FIMO output
        #DNA sequences surrounding the candidate motif need to be +3bp wider than PadLen (EXT = 25 bp)
        ocurrences$V3 <- ocurrences$V3 - EXT - fastaL
        ocurrences$V4 <- ocurrences$V4 + EXT + fastaL +1
 #   #Assure integre values (no scientific notation)
    ocurrences$V3 <- as.integer( ocurrences$V3 )
    ocurrences$V4 <- as.integer( ocurrences$V4 )	

#options(scipen=999)
        write.table(x=ocurrences[,2:4], file = paste(motif,bamFile,".bed4fasta.bed",sep=""), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
#options(scipen=0)        
	#BEDtools
	Sys.sleep(5)
        fastaBed <-paste(motif,bamFile,".bed4fasta.bed", sep="") 
        system(paste("./bedtools getfasta -fi Homo_sapiens.GRCh38_15.fa -bed",fastaBed, "-fo",fastaOut, sep=" "))
      }
      
      print('FASTA file for the motifs...DONE')
      
      ## Capture and Save all the DNase I cuts in a text file
      #setwd("/Users/pm12/Desktop/EPIGENODE-analysis/FootprintMixture2")
      
      #Get ATAC insertions
      #Select Bam
      bam.file=paste( bamPath, bamFile,".bam", sep="")
      
      #footprints
      #footprints <- read.table(gzfile("/lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/msCentipede/M5744_1.02_motifs_GRCh38.txt.gz"), header=T)
      footprints <- readGeneric(paste(motif,bamFile,".bed.genomation.txt",sep=""), header=F, keep.all.metadata = TRUE, strand = 4)
      

      #Genomation
      #footprints = footprints[seqnames(footprints) == "chr21"]
      #footprints = footprints[order(-footprints$PwmScore )]
      #footprints = resize(footprints, width = 500, fix = "center")
      head(footprints)
      start(footprints) <- start(footprints) - EXT
      end(footprints)   <- end(footprints) + EXT
      nBins <- Lm+EXT+EXT
      sm <-  ScoreMatrixBin(target = bam.file, bin.num = nBins, windows = footprints, rpm=F, type="bam", strand.aware = TRUE, extend=1)  #footprints[1:1000,]
      #plot(colMeans(sm),type='l', col='blue')
      
      
      footprints2 <- footprints
      start(footprints2) <- start(footprints2) + 4 #Tn5
      end(footprints2)   <- end(footprints2)   +  4#    #Tn5
      head(footprints2)
      nBins2 <- Lm+EXT+EXT #+8 
      #Read Tn5 cuts
      smF <- ScoreMatrixBin(target = bam.file, bin.num = nBins2, windows = footprints2, rpm=F, type="bam", strand.aware = TRUE, extend=1, param = ScanBamParam(which=reduce(footprints2, ignore.strand=T), flag=scanBamFlag(isMinusStrand=FALSE)))  #footprints[1:1000,]
      ##points(colMeans( smF[,1:(nBins2-4)] ),type='l', col='green')
      #points(colMeans( smF),type='l', col='green')
      
      footprints3 <- footprints
      start(footprints3) <- start(footprints3) -4
      end(footprints3)   <- end(footprints3)   -4   #Tn5
      head(footprints3)
      nBins3 <- Lm+EXT+EXT #+8
      #Read Tn5 cuts
      smR <- ScoreMatrixBin(target = bam.file, bin.num = nBins3, windows = footprints3, rpm=F, type="bam", strand.aware = TRUE, extend=1, param = ScanBamParam(which=reduce(footprints3, ignore.strand=T), flag=scanBamFlag(isMinusStrand=TRUE)))  #footprints[1:1000,]
      ##points(colMeans(smR[,5:(nBins3)]),type='l', col='purple')
      #points(colMeans(smR),type='l', col='purple')
      
      nc=Lm+EXT+EXT 
      sm[,1:nc] <- smF   + smR  #smR[,5:(nBins2)] #smR[,1: nc]
      #sm[,1:nc] <- smF[,1:(nBins2-4)]   + smR[,5:(nBins3)] #smR[,5:(nBins2)] #smR[,1: nc]
      #points(colMeans(sm),type='l', col='black', lwd=2)
      
      #Save cut.data.txt
      dim(sm)
      write.table(x=sm, file = paste(motif,bamFile,".cut.data.txt",sep=""), append = FALSE, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
      
      
      L <- (end(footprints)[1])  - (start(footprints)[1] )

#      
#      pdf(file=paste(paste(TF_Name,motif,bamFile,sep="_"),"pdf",sep="."), height=7, width=8)
#      #par(mfrow=c(3,1))
#      #layout(matrix(c(1,2,4,5,3,6,7,3,6), 3, 3, byrow = TRUE))
##      layout(matrix( c(1,7,6,2,4,6,3,5,6), 3, 3, byrow = TRUE))
##layout(matrix( c(9,4,7,8,2,3,6,1,5), 3, 3, byrow = TRUE))
#layout(matrix( c(9,5,7,8,3,4,2,1,6), 3, 3, byrow = TRUE))
#
#      par(mar=c(5,4,4,3)+.1)
#      plot(colMeans(smF)+colMeans(smR), xaxt = "n", col=alpha("#8B008B",0.5),lty=1, type='l' , ylab='' , xlab="Distance to motif centre (bp)", cex.lab=1.4, main=TF_Name, lwd=2, frame=FALSE)
#axis(1, at=c(1,length(colMeans(smF))/2,length(colMeans(smF))), labels=c(-(L/2),0,L/2 ))
#   mtext("Mean Tn5 integrations",side=2,line=3,col="black") #  alpha("#8B008B",0.5)  )
#      #points(colMeans(smR), col=alpha("blue",0.5),lty=1 , type='l' , lwd=1.5)
#      #legend("topright", legend=c("+","-"),lty=c(1,1), col=c("blue","red"), bty='n' )
##      par(new = T)
#      #see MixtureModel.r
#      
      
      print('Starting FootprintMixture Analysis For MOTIF and SAMPLE....')
      ## Files OK. READY to run FootprintMixture ##
      
      #The number of DNase-seq reads that map to each coordinate surrounding the candidate binding sites. Typically the model uses a 25bp windows that surround the candidate binding site (upstream and downstream). This data can be stored as a matrix (a tab separated file), similar to example file (cut.data.txt).
      c <- read.table(paste(motif,bamFile,".cut.data.txt",sep=""));
      #The coordinates of candidate binding sites in bed format (bed.data.txt).
      b <- read.table(paste(motif,bamFile,".bed.data.txt",sep=""));
      #The coordinates of ChIP-seq peaks for identifying known binding sites. A bed file in narrowpeak format may be used, similar to example file(peak.data.txt). In the absence of known ChIP-seq peaks, coordinates of DNase hypersensitive sites may be used since these regions tend to be enriched for binding sites.
      p <- read.table("EPIGENODE_merged_ATAC_peaks_final.bed");
	
	
      #dim(c)
      #dim(b)
      #dim(p)
      
      #head(c)
      #head(b)
      #head(p)
      #source("MixtureModel.r")
        BuildSeqBiasBackground <- function(FastaName){
        BiasFile <- "SeqBias.txt";
        inFA <- FastaName;
        outSignal <-  paste("signal", TF_Name,motif,bamFile,".txt",sep="_") ;  #"signal.txt";
        system( paste("perl RebuildSignal.pl ", BiasFile, " ", inFA, " > ", outSignal, sep="" ) );
        signal <- read.table(outSignal);
        M <- signal / sum(signal);
	#pedro
	system ( paste("rm", paste("signal", TF_Name,motif,bamFile,".txt",sep="_"), sep= " ") );
        return(M);
        }


      #Seq
      m <- MultMMixture_Full(TF_Bed=b,Cuts=c,peakbed=p,Plot=F,PadLen=EXT,Collapse=T,k=2,ReturnPar=T,Fixed=T,Background=controlSeqBias,FastaName=fastaOut);
      
      
#      plotMeta(sm, xcoords = c(-L/2, L/2), lwd=1, cex.lab=1.4, line.col ="darkblue", ylab='Average Tn5 integrations', xlab="distance to motif centre", main=TF_Name) ##104E8B
#      abline(v=-L/2+EXT, lty=2, col="darkgray", lwd=1.5)
#      abline(v=-L/2+EXT+Lm, lty=2, col="darkgray", lwd=1.5)
#      #if ( mean(sm) <  0.1) {WI = 100}  #winsorize in heatMatrix might fail -> use (0,100)
#      #if ( mean(sm) >= 0.1) {WI = 99 }
#    WI=99  #100
#      heatMatrix(sm, order=TRUE, cex.lab=0.8,legend.name=c("Tn5 transposase integration (scaled)"),main=paste(nrow(sm),"regions",sep=" "), cex.main=0.8,xcoords =c(-L/2, L/2) ,winsorize=c(0,WI), col=colorRampPalette(c("white","darkblue"))(100) );  #    c("white","darkblue")) #main=TF_Name    "white","darkblue"


      #Footprint likelihood values (FLR)
      CUTOFF = 5  #10  #8  #10
      #A likelihood ratio of greater than 1 indicates the test result is associated with the footprint
      #A likelihood ratio less than 1 indicates that the result is associated with absence of the footprint
      trueF  <- which(m$llr >= CUTOFF)
      falseF <- which(m$llr < CUTOFF)
      nanF   <- which( is.nan(m$llr) == TRUE)

#      
#      #plot mean
#      xl <- -(L/2):(L/2) 
#      ##plot F and R strands
#      #plot(xl,colMeans(smF[trueF,]), type='l', lwd=1.5, col="blue",frame=FALSE, ylim=c(0,max(colMeans(smF[trueF,]), colMeans(smR[trueF,])) ), ylab="Tn5 transposase read-start sites", xlab="position", cex.lab=1.4, main=paste(motif,"Bound", sep=" ")  )
#      #points(xl,colMeans(smR[trueF,]), type='l', lwd=1.5, col="red")
#      #legend("top", col=c("blue","red"),legend=c("+","-"), lwd=1.5, bty="n") #x="top"
#      #plot(xl,colMeans(smF[falseF,]), type='l', lwd=1.5, col="blue", lty=2,frame=FALSE, ylim=c(0,max(colMeans(smF[trueF,]), colMeans(smR[trueF,])) ), ylab="Tn5 transposase read-start sites", xlab="position", cex.lab=1.4, main=paste(motif,"Unbound", sep=" ")  )
#      #points(xl,colMeans(smR[falseF,]), type='l', lwd=1.5, col="red", lty=2)
#      #legend("top", col=c("blue","red"),legend=c("+","-"), lwd=1.5, bty="n") #x="top"
#      
#      
#      bcol <- brewer.pal(n = 3, name = "Dark2")[c(3,1,2)]
#      #par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
#      #par(mar=c(6.4,4.1,4.1, 8.1), xpd=FALSE) # this is usually the default
## bcol <- brewer.pal(n = 11, name = "Spectral")[c(11,1,8)]
#
#      
#if( length(trueF) > 0 ) { 
#
#      par(mar=c(8.1,1.1,0,6.1), xpd=TRUE)
#   #YLIMPLOTF=c(min(colMeans(c[falseF,]) -  apply(t(c[falseF,]),1,sd)),  max(colMeans(c[trueF,]) +  apply(t(c[trueF,]),1,sd))  )
#      #plot( xl, colMeans(c), type='l', lwd=1.5, frame=FALSE, ylim=c(0,max(colMeans(c[trueF,]))), ylab="Tn5 transposase read-start sites", xlab="position", cex.lab=1.4, main=motif,lty=2) #ylim=c(0,3)
#      plot(xl,colMeans(c[trueF,]), type='l', lwd=1.5, col=bcol[1],frame=FALSE, ylim=c(0,max(colMeans(c[trueF,]))), ylab="Mean Tn5 integrations", xlab=paste("Distance to motif centre (bp)",TF_Name, sep="\n"), cex.lab=1.4, main="", mgp= c(3, 0.5, 0))
#      points(xl,colMeans(c[falseF,]), type='l', lwd=1.5, col=bcol[2])
#   #points(xl, colMeans(c[trueF,]) +  apply(t(c[trueF,]),1,sd), col=alpha(bcol[1],0.5), lwd=1  )
#   #points(xl, colMeans(c[trueF,]) -  apply(t(c[trueF,]),1,sd), col=alpha(bcol[1],0.5), lwd=1  )
#   #points(xl, colMeans(c[falseF,]) +  apply(t(c[falseF,]),1,sd), col=alpha(bcol[2],0.5), lwd=1  )
#   #points(xl, colMeans(c[falseF,]) -  apply(t(c[falseF,]),1,sd), col=alpha(bcol[2],0.5), lwd=1  )
#   mtext("Mean Tn5 integrations",side=2,line=3,col="black")
#      #abline(v=round(-Lm/2), col="gray",lty=2, lwd=1.5)
#      #abline(v=round(Lm/2), col="gray",lty=2, lwd=1.5)
#      segments(y0=0, y1=max(colMeans(c[trueF,])),  x0= -Lm/2, x1= -Lm/2, col="gray",lty=2, lwd=1.5)
#      segments(y0=0, y1=max(colMeans(c[trueF,])) , x0= Lm/2, x1= Lm/2, col="gray",lty=2, lwd=1.5)
#      
#      #par(mar=c(6.4,4.1,4.1, 8.1),xpd=TRUE)
#      #legend(EXT+Lm,1, col=c("black",bcol[1],bcol[2]),lty=c(2,1,1),legend=c("Mean", "Bound","Unbound"), lwd=1.5, bty="n", inset=c(90,0)) #x="top"
#      legend(EXT,max(colMeans(c[trueF,])), col=c(bcol[1],bcol[2]),legend=c("Bound","Unbound"), lwd=1.5, bty="n", inset=c(100,0)) #x="top"
#      #(EXT+Lm,1)
#      
#      
#      #boxplot
#      par(mar=c(5.1,6.1,4.1,7.1), xpd=TRUE) #outline =T
#      boxplot(m$llr[trueF],m$llr[falseF],outline =F, main=paste(length(trueF)+length(falseF), "matches\n of the PWM analysed",sep=" "), ylab="Footprint Likelihood Ratio", names=c(length(trueF), length(falseF)), boxwex = 0.6, horizontal = F, border =c(bcol[1],bcol[2]) , cex.lab=1.4, frame=F )
#      legend(x="topright", col=c(bcol[1],bcol[2]),legend=c("Bound","Unbound"), lwd=1.5, bty="n", inset=c(-0.9,0)) #x="top"
#      
#      #heatplot of the groups
#      #if ( mean(sm) <  0.25) {WI = 100}  #winsorize in heatMatrix might fail -> use (0,100) DONE above
#      #if ( mean(sm) >= 0.25) {WI = 99 }
##      heatMatrix(sm, group=list(Bound=trueF, Unbound=falseF, Unknown=nanF),group.col=c(bcol[1],bcol[2],bcol[3]), order=TRUE, legend.name="Tn5 transposase integration (scaled)",xlab="distance to motif centre", main=TF_Name, xcoords =c(-L/2, L/2) ,winsorize=c(0,WI), col=colorRampPalette(c("white","darkblue"))(100), cex.lab=0.8)
#      
#      
#      #library(png)
#      par(mar=c(0,4,8.5,9.5), xpd=TRUE)
#      ima <- readPNG(paste("/lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/fimo/CISBP/logos_all_motifs/",motif,"_fwd.png",sep=""))
#      
#      #Set up the plot area
#      plot(1:2, type='n', main="Plotting Over an Image", xlab="x", ylab="y",frame.plot=F, yaxt='n' , xaxt='n', ann=FALSE)
#      mtext=motif
#      #Get the plot information so the image will fill the plot box, and draw it
#      lim <- par()
#      rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4],main=motif)
#      #abline(v=1, col="gray",lty=2, lwd=1.5)
#      #abline(v=2, col="gray",lty=2, lwd=1.5)
#      segments(x0=1.1, y0=1.1, x1 = 1.48, y1 = 0.9,  col="gray",lty=1, lwd=1.5)
#      segments(x0=2.0, y0=1.1, x1 = 1.53, y1 = 0.9,  col="gray",lty=1, lwd=1.5)
#      
#} # END if( length(trueF) > 0 ) {       
#

      #dev.off()
    	
      #UPDATE April 2016
      #Save coordinates of BOUND  binding sites in bed format (bed.data.txt).
#setwd("/lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/FootprintMixture/Footprints_plots/");
setwd("/lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC2/FootprintMixture/Footprints_plots/")

      if (length(trueF)> 0 ) {
      write.table(x=cbind(b[trueF,],  m$llr[trueF], rep(motif, length( m$llr[trueF])) , rep(bamFile, length( m$llr[trueF])) , rep(TF_Name, length( m$llr[trueF]))    ),  file = paste(TF_Name,motif,bamFile,"Bound.LLR5.bed",sep="_"),   append = F, quote = F, sep = "\t", row.names = F, col.names = F);
      }
      if (length(falseF)> 0 ) {
      write.table(x=cbind(b[falseF,], m$llr[falseF] ), file = paste(TF_Name,motif,bamFile,"Unbound.LLR5.bed",sep="_"), append = F, quote = F, sep = "\t", row.names = F, col.names = F);
      }
      if (length(nanF)> 0 ) {
      write.table(x=cbind(b[nanF,], m$llr[nanF] ), file = paste(TF_Name,motif,bamFile,"Unknown.LLR5.bed",sep="_"), append = F, quote = F, sep = "\t", row.names = F, col.names = F);
      }
setwd("/lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC2/FootprintMixture/"); 

#
#if( length(trueF) > 0 ) { 
#      #Properties TF footprints
#      openess = ( sum(colMeans(c[trueF,])[1:(1+EXT)]) + 1 +  sum(colMeans(c[trueF,])[(1+EXT+Lm):(ncol(c)) ]) ) /2;  #sm-> c[trueF,]
#      print(openess)
#      directionality = (sum(colMeans(c[trueF,])[1:(1+EXT)]) +1) / (sum(colMeans(c[trueF,])[(1+EXT+Lm):(ncol(c)) ]) +1) ;
#      print(directionality)
#      depth = openess  /  ( (sum(colMeans(c[trueF,])[(1+EXT):(1+EXT+Lm)])  ) +1);
#      print(depth)
#}
#

#if( length(trueF) <= 0 ) { 
#	openess = NA
#	directionality = NA
#	depth = NA
#}
#
      #save results/plots in Results File
      write.table(x=t(c(TF_Name, motif,Lm, EXT, bamFile, nrow(sm), length(m$llr), length(trueF), length(falseF), length(nanF) )) , file = paste(resultsPath, paste(TF_Name,motif,bamFile,".fm",sep="_") ,sep="/"), append = F, quote = F, sep = "\t", row.names = F,  col.names = F)     
#system(paste("mv", paste(paste(TF_Name,motif,bamFile, sep="_"),"pdf",sep=".")  , resultsPath, sep=" ")  )
 
     
#
#if( length(trueF) > 0 ) {
#      
#	#Properties TF footprints(2)
#	openess_True <- c()
#	directionality_True <- c()
#	depth_True <- c()
#	for (di in 1:nrow(c[trueF,]) ){
#		
#		openess_True[di] <- ( sum(c[trueF[di],][1:(1+EXT)]) + 1 + sum(c[trueF[di],][(1+EXT+Lm):(ncol(c)) ]) ) /2; 	
#		directionality_True[di] <- (sum(c[trueF[di],][1:(1+EXT)])+1) / (sum(c[trueF[di],][(1+EXT+Lm):(ncol(c)) ])+1) ;
#		depth_True[di] <- openess_True[di]  / ( (sum(c[trueF[di],][(1+EXT):(1+EXT+Lm)])  )+1) ;
#	}
#
#	rm(di)
#	openess_False <- c()
#       directionality_False  <- c()
#        depth_False  <- c()
#        for (di in 1:nrow(c[falseF,]) ){
#
#                openess_False[di] <- ( sum(c[falseF[di],][1:(1+EXT)]) + 1 + sum(c[falseF[di],][(1+EXT+Lm):(ncol(c)) ]) ) /2;
#                directionality_False[di] <- (sum(c[falseF[di],][1:(1+EXT)]) +1) / (sum(c[falseF[di],][(1+EXT+Lm):(ncol(c)) ]) +1) ;
#                depth_False[di] <- openess_False[di]  / ( (sum(c[falseF[di],][(1+EXT):(1+EXT+Lm)])  ) +1) ;
#        }
#
#
#	par(mar=c(6.1,5.1,4.1,5.1), xpd=TRUE) 
#	openess_True <- log2(openess_True )
#	directionality_True <- log2(directionality_True )
#	depth_True <- log2(depth_True )
#	
#	openess_False <- log2(openess_False )
#	directionality_False <- log2( directionality_False  ) 
#	depth_False <- log2(depth_False  )
#
#
#	boxplot(openess_True,openess_False, directionality_True, directionality_False, depth_True, depth_False, las=2, ylab="log2(Value)", names=c('Bound','Unbound','Bound','Unbound','Bound','Unbound') , at=c(1,2,4,5,7,8), col=c("red","red","sienna","sienna","royalblue2","royalblue2") , cex.lab=1.4, frame=F , boxwex = 0.6, horizontal = F, outline=F, xlab="" )
#	
#	legend(x="top", col=c("red","sienna","royalblue2"),legend=c("Chromatin accessibility","Directional binding","Footprint depth"), lwd=1.5, bty="n", inset=c(0,-0.5)) 
#
##	par(mfrow=c(1,3), mar=c(7,5,7,5) )  # , xpd=TRUE )
##	boxplot(openess_True,openess_False, names=c('Bound','Unbound') , las=2, outline=T,xlab="",  ylab="Chromatin openess", boxwex = 0.8, horizontal = F, border =c(bcol[1],bcol[2]) , cex.lab=1.4, frame=F )
##	boxplot(directionality_True,directionality_False, names=c('Bound','Unbound') , las=2, outline=T,xlab="", ylab="Directionality", boxwex = 0.8, horizontal = F, border =c(bcol[1],bcol[2]) , cex.lab=1.4, frame=F )
##	boxplot(depth_True,depth_False, names=c('Bound','Unbound'), las=2, outline=T, xlab="", ylab="Footprint depth", boxwex = 0.8, horizontal = F, border =c(bcol[1],bcol[2]) , cex.lab=1.4, frame=F )
#
#} #END if( length(trueF) > 0 ) {
#
#	dev.off()
#	system(paste("mv", paste(paste(TF_Name,motif,bamFile, sep="_"),"pdf",sep=".")  , resultsPath, sep=" ")  )
#
#

      # remove intermediate files
#      system ( paste("rm", paste("signal", TF_Name,motif,bamFile,".txt",sep="_"), sep= " ") )
      system (  paste("rm", paste(motif,bamFile,".bed.data.txt",sep=""), sep= " "   ) )
      system (  paste("rm", paste(motif,bamFile,".bed.genomation.txt",sep=""), sep= " "   ) )
      system (  paste("rm", paste(motif,bamFile,".fa", sep="") , sep= " "   ) )
      system (  paste("rm", paste(motif,bamFile,".bed4fasta.bed",sep="") , sep= " "   ) )
      system (  paste("rm", paste(motif,bamFile,".cut.data.txt",sep="") , sep= " "   ) )
      #Remove PDF
#      system (  paste("rm",paste(TF_Name,motif,bamFile,".pdf",sep="_")  , sep= " "   )   )
    

#      #remove objects
#      if( length(trueF) > 0 ) {
#      rm(lim,ima, openess_True, directionality_True, depth_True, openess_False, directionality_False, depth_False  )
#      }
      rm(smF,smR,sm,trueF,valid,falseF, footprints,footprints2,footprints3,nanF,nBins,nBins2,nBins3 )
      rm(TF_Name,nc,m,L,b,c,p,fimo,ocurrences  ) #WI, xl


    } #end IF
    
    
    toc()
    
  } #end  TF_motifs FOR loop
  
} #end Bam_ATACC FOR loop

} #end of If TF_motif

#end###########################################################################
###############################################################################
###############################################################################
###############################################################################
# Automate the pipeline for >1 bam and >1 motif
# Plots for overall results
# - Scatterplot of LLRs for 1 TF in 1 bam vs another bam
# - Heatmaps....


