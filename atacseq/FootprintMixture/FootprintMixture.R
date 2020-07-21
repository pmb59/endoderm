# For ATAC21/ATAC22 (36h)

## FootprintMixture     ##
## Detection of TF footprints while controlling for Tn5-specific sequence bias ##

# FIMO should be run first (used top-scored motifs) 
# For now will use TF_Information_all_motifs.txt (Directly determined or best inferred motif)
# Move coordinates +4bp / -5bp
# Generate SeqBias.txt from ATAC17.bam >R script>SeqTools  


args<-commandArgs(TRUE)

#Supress global warnings
options(warn=-1)

library(grid)
library(Rsamtools)
library(seqLogo)
library(Biostrings)
library(png)
library(genomation)
library(GenomicRanges)
library(scales)
library(gtools);
library(RColorBrewer)
library(tictoc)
source("MixtureModel.r");

setwd("../fimo/fimo_occurrences_500k")
Stored <- unlist(strsplit(list.files()  ,split=".txt"))  
OfInterest <- as.character(read.csv2("../fimo/CISBP/TF_Information.txt", sep="\t")$Motif_ID)
#which are common in both list?
TF_motifs <- intersect(Stored , OfInterest)

InMyList <- which (TF_motifs == as.character(args[1])  )
if (length(InMyList ) >0 ) {
TF_motifs = as.character( args[1] );


#Bam alignments
Bam_ATAC  <- as.character( args[2] )   


for (q in 1:length(TF_motifs) ){

  for (ba in 1:length(Bam_ATAC) ){
    print(paste( paste("#Motif = ",q, sep="") ,paste("#bamFile = ", Bam_ATAC[ba] , sep="") ,sep=" ; ") )
    
    tic();
    #--------------------------------------------------------------------------------------------------
    #PARAMS:
    setwd("../FootprintMixture/");  #working directory
    bamFile = Bam_ATAC[ba]   
    bamPath = "/../" 
    
    motif = TF_motifs[q]   #M5943_1.02"   #USF1 "M4509_1.02" #"M1957_1.02"  #"M5808_1.02" 
    Nmotifs = 10e3   # Min number of motifs
    NROWS = 1.0*Nmotifs 
    EXT = 25  # bp extension
    chrom_Sizes="GRCh38_15.chrom.sizes" #chromosomes that we actually use
    controlSeqBias="Seq"  # "Flat" or "Seq"
    genome <- "Homo_sapiens.GRCh38_15.fa"
    database <- "./fimo/CISBP/TF_Information.txt" 
    fimo_txt <- "./fimo/fimo_occurrences_500k"
    resultsPath="./FootprintMixture/Footprints_plots"
    db <- read.csv(database, head=T, sep="\t")
    #--------------------------------------------------------------------------------------------------
    
    
    # Get the NAME of the TF of interest (paste if more than 1)
    # To avoid FILENAME TOO LONG problem (in few cases), select the first 20 characters
    TF_Name <- substr(  x=paste( as.character(db$TF_Name[which(db$Motif_ID == motif)]) , collapse="_") , start=1, stop= 20 )
    print(paste(TF_Name, motif, sep=" - "))
    
    ## CIS-BP to MEME format for this motif (Done)
    ##( Rscript cisbp2meme.r)
    
    ## Run FIMO to get matches in the human genome (Done)
    #fimo.sh
    
    
    fimo_Matches = paste(fimo_txt,"/",motif,".txt",sep="")   #Motifs are ranked
    
    
    # Check if we have NROWS motifs in the file
    if(  length( (readLines(fimo_Matches)) )   >=  (NROWS -1) )   {
      print('More than 10k motif occurrences found...Starting Analysis...')
      
      #Filter Out Motif ocurrence out of chr1-22, X,Y
      fimo <- read.table(fimo_Matches, header = FALSE  )  


	##########################################
	# check that End > Start in FIMO input
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
				valid <- c(valid, j) }  
			}
			rm(AuxChr)
		}
      }
      nrow(fimo)
      length(valid)
      
      
      #Store bed file corresponding to the matches (candidate TFBSs)- used the Valid
      fimo <- fimo[valid, ]; dim(fimo)

      ocurrences <- fimo      
      dim(ocurrences)
      


      write.table(x=ocurrences[,2:4], file = paste(motif,bamFile,".bed.data.txt",sep=""), append = FALSE, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
      write.table(x=ocurrences[,2:7], file = paste(motif,bamFile,".bed.genomation.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

      print('.bed.data.txt and .bed.genomation.txt... SAVED IN DISK')

      #Motif Length
      Lm <- nchar(x=as.character(ocurrences$V8[1]))  #e.g., 19
      
      
      fastaOut <-paste(motif,bamFile,".fa", sep="") 
      fastaL=3
      if (controlSeqBias == "Seq") {
        	# Store fasta file corresponding to the matches from the FIMO output
        	# DNA sequences surrounding the candidate motif need to be +3bp wider than PadLen (EXT = 25 bp)
        	ocurrences$V3 <- ocurrences$V3 - EXT - fastaL
        	ocurrences$V4 <- ocurrences$V4 + EXT + fastaL +1
    		# Assure integre values (no scientific notation)
    		ocurrences$V3 <- as.integer( ocurrences$V3 )
    		ocurrences$V4 <- as.integer( ocurrences$V4 )	
        	write.table(x=ocurrences[,2:4], file = paste(motif,bamFile,".bed4fasta.bed",sep=""), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)      
		# BEDtools
		Sys.sleep(5)
        	fastaBed <-paste(motif,bamFile,".bed4fasta.bed", sep="") 
        	system(paste("./bedtools getfasta -fi Homo_sapiens.GRCh38_15.fa -bed",fastaBed, "-fo",fastaOut, sep=" "))
      }
      
      print('FASTA file for the motifs...DONE')
      
      
      #Get TN5 insertions
      #Select Bam
      bam.file=paste( bamPath, bamFile,".bam", sep="")
      
      #footprints
      footprints <- readGeneric(paste(motif,bamFile,".bed.genomation.txt",sep=""), header=F, keep.all.metadata = TRUE, strand = 4)
      

      #Genomation
      head(footprints)
      start(footprints) <- start(footprints) - EXT
      end(footprints)   <- end(footprints) + EXT
      nBins <- Lm+EXT+EXT
      sm <-  ScoreMatrixBin(target = bam.file, bin.num = nBins, windows = footprints, rpm=F, type="bam", strand.aware = TRUE, extend=1) 
    
      
      footprints2 <- footprints
      start(footprints2) <- start(footprints2) + 4
      end(footprints2)   <- end(footprints2)   + 4
      head(footprints2)
      nBins2 <- Lm+EXT+EXT 
      #Read Tn5 cuts
      smF <- ScoreMatrixBin(target = bam.file, bin.num = nBins2, windows = footprints2, rpm=F, type="bam", strand.aware = TRUE, extend=1, param = ScanBamParam(which=reduce(footprints2, ignore.strand=T), flag=scanBamFlag(isMinusStrand=FALSE)))  
      
      footprints3 <- footprints
      start(footprints3) <- start(footprints3) -4
      end(footprints3)   <- end(footprints3)   -4   
      head(footprints3)
      nBins3 <- Lm+EXT+EXT 
      #Read Tn5 cuts
      smR <- ScoreMatrixBin(target = bam.file, bin.num = nBins3, windows = footprints3, rpm=F, type="bam", strand.aware = TRUE, extend=1, param = ScanBamParam(which=reduce(footprints3, ignore.strand=T), flag=scanBamFlag(isMinusStrand=TRUE)))  
      
      nc=Lm+EXT+EXT 
      sm[,1:nc] <- smF   + smR 
	    
      
      #Save cut.data.txt
      dim(sm)
      write.table(x=sm, file = paste(motif,bamFile,".cut.data.txt",sep=""), append = FALSE, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
      
      
      L <- (end(footprints)[1])  - (start(footprints)[1] )


      print('Starting FootprintMixture Analysis For MOTIF and SAMPLE....')
      ## Files OK. READY to run FootprintMixture ##
      
      #The number of ATAC-seq reads that map to each coordinate surrounding the candidate binding sites. Typically the model uses a 25bp windows that surround the candidate binding site (upstream and downstream). This data can be stored as a matrix (a tab separated file), similar to example file (cut.data.txt).
      c <- read.table(paste(motif,bamFile,".cut.data.txt",sep=""));
      #The coordinates of candidate binding sites in bed format (bed.data.txt).
      b <- read.table(paste(motif,bamFile,".bed.data.txt",sep=""));
      #The coordinates of ChIP-seq peaks for identifying known binding sites. A bed file in narrowpeak format may be used, similar to example file(peak.data.txt). In the absence of known ChIP-seq peaks, coordinates of DNase hypersensitive sites may be used since these regions tend to be enriched for binding sites.
      p <- read.table("EPIGENODE_merged_ATAC_peaks_final.bed");
	
      BuildSeqBiasBackground <- function(FastaName){
      BiasFile <- "SeqBias.txt";
      inFA <- FastaName;
      outSignal <-  paste("signal", TF_Name,motif,bamFile,".txt",sep="_") ;  #"signal.txt";
      system( paste("perl RebuildSignal.pl ", BiasFile, " ", inFA, " > ", outSignal, sep="" ) );
      signal <- read.table(outSignal);
      M <- signal / sum(signal);
      system ( paste("rm", paste("signal", TF_Name,motif,bamFile,".txt",sep="_"), sep= " ") );
      return(M);
      }

      #Seq
      m <- MultMMixture_Full(TF_Bed=b,Cuts=c,peakbed=p,Plot=F,PadLen=EXT,Collapse=T,k=2,ReturnPar=T,Fixed=T,Background=controlSeqBias,FastaName=fastaOut);
      

      # Footprint likelihood values (FLR)
      CUTOFF = 5  
      # A likelihood ratio of greater than 1 indicates the test result is associated with the footprint
      # A likelihood ratio less than 1 indicates that the result is associated with absence of the footprint
      trueF  <- which(m$llr >= CUTOFF)
      falseF <- which(m$llr < CUTOFF)
      nanF   <- which( is.nan(m$llr) == TRUE)


      # Save coordinates of BOUND binding sites in bed format (bed.data.txt)
      setwd("/../FootprintMixture/Footprints_plots/")

      if (length(trueF)> 0 ) {
      write.table(x=cbind(b[trueF,],  m$llr[trueF], rep(motif, length( m$llr[trueF])) , rep(bamFile, length( m$llr[trueF])) , rep(TF_Name, length( m$llr[trueF]))    ),  file = paste(TF_Name,motif,bamFile,"Bound.LLR5.bed",sep="_"),   append = F, quote = F, sep = "\t", row.names = F, col.names = F);
      }
      if (length(falseF)> 0 ) {
      write.table(x=cbind(b[falseF,], m$llr[falseF] ), file = paste(TF_Name,motif,bamFile,"Unbound.LLR5.bed",sep="_"), append = F, quote = F, sep = "\t", row.names = F, col.names = F);
      }
      if (length(nanF)> 0 ) {
      write.table(x=cbind(b[nanF,], m$llr[nanF] ), file = paste(TF_Name,motif,bamFile,"Unknown.LLR5.bed",sep="_"), append = F, quote = F, sep = "\t", row.names = F, col.names = F);
      }
      setwd("/../FootprintMixture/"); 


      #save results/plots in Results File
      write.table(x=t(c(TF_Name, motif,Lm, EXT, bamFile, nrow(sm), length(m$llr), length(trueF), length(falseF), length(nanF) )) , file = paste(resultsPath, paste(TF_Name,motif,bamFile,".fm",sep="_") ,sep="/"), append = F, quote = F, sep = "\t", row.names = F,  col.names = F)     


      # remove intermediate files
      system (  paste("rm", paste(motif,bamFile,".bed.data.txt",sep=""), sep= " "   ) )
      system (  paste("rm", paste(motif,bamFile,".bed.genomation.txt",sep=""), sep= " "   ) )
      system (  paste("rm", paste(motif,bamFile,".fa", sep="") , sep= " "   ) )
      system (  paste("rm", paste(motif,bamFile,".bed4fasta.bed",sep="") , sep= " "   ) )
      system (  paste("rm", paste(motif,bamFile,".cut.data.txt",sep="") , sep= " "   ) )
      rm(smF,smR,sm,trueF,valid,falseF, footprints,footprints2,footprints3,nanF,nBins,nBins2,nBins3 )
      rm(TF_Name,nc,m,L,b,c,p,fimo,ocurrences  ) 


    } 

    toc()
    
  } 
  
} 

} 

