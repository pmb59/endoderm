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
source("./FootprintMixture/MixtureModel.r");


setwd("../fimo/fimo_occurrences_500k")
Stored <- unlist(strsplit(list.files()  ,split=".txt"))  
OfInterest <- as.character(read.csv2("../fimo/CISBP/TF_Information.txt", sep="\t")$Motif_ID)
#which are common in both list?
TF_motifs <- intersect(Stored , OfInterest)


InMyList <- which (TF_motifs == as.character(args[1])  )
if (length(InMyList ) >0 ) {
TF_motifs = as.character( args[1] );


#Bam alignments
Bam_ATAC  <- as.character( args[2] )    #c( "ATAC1" )
#Bam_ATAC  <-  c( "ATAC1", "ATAC2", "ATAC9", "ATAC10", "ATAC11" , "ATAC12",  "ATAC13", "ATAC14",  "ATAC15", "ATAC16", "ATAC17" )


for (q in 1:length(TF_motifs) ){

  
  for (ba in 1:length(Bam_ATAC) ){
    print(paste( paste("#Motif = ",q, sep="") ,paste("#bamFile = ", Bam_ATAC[ba] , sep="") ,sep=" ; ") )
    
    tic();
    #--------------------------------------------------------------------------------------------------
    #PARAMS:
    setwd("../FootprintMixture/");  #working directory
    bamFile = Bam_ATAC[ba]   #"ATAC1-2"
    bamPath = "/../" 
    
    motif = TF_motifs[q]   #M5943_1.02"   #USF1 "M4509_1.02" #"M1957_1.02"  #"M5808_1.02" 
    Nmotifs = 10e3 # 50e3 as a stringet set (Yardimci et al)  / OR MINIMUM NUMBER OF MOTIFS TO ANALYZE
    NROWS = 1.0*Nmotifs #1.1* 60e3 #Keep only autosomal and sex-chrs matches
    EXT = 25  #35 #25
    chrom_Sizes="GRCh38_15.chrom.sizes" #chromosomes that we actually use
    controlSeqBias="Seq"  # "Flat" or "Seq"
    genome <- "./Homo_sapiens.GRCh38_15.fa"
    database <- "./fimo/CISBP/TF_Information.txt" #TF_Information_all_motifs_plus.txt" 
    fimo_txt <- "./fimo/fimo_occurrences_500k"
    resultsPath="./FootprintMixture/Footprints_plots"
    db <- read.csv(database, head=T, sep="\t")
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
      
      


      #Footprint likelihood values (FLR)
      CUTOFF = 5  #10  #8  #10
      #A likelihood ratio of greater than 1 indicates the test result is associated with the footprint
      #A likelihood ratio less than 1 indicates that the result is associated with absence of the footprint
      trueF  <- which(m$llr >= CUTOFF)
      falseF <- which(m$llr < CUTOFF)
      nanF   <- which( is.nan(m$llr) == TRUE)



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


      #save results/plots in Results File
      write.table(x=t(c(TF_Name, motif,Lm, EXT, bamFile, nrow(sm), length(m$llr), length(trueF), length(falseF), length(nanF) )) , file = paste(resultsPath, paste(TF_Name,motif,bamFile,".fm",sep="_") ,sep="/"), append = F, quote = F, sep = "\t", row.names = F,  col.names = F)     
#system(paste("mv", paste(paste(TF_Name,motif,bamFile, sep="_"),"pdf",sep=".")  , resultsPath, sep=" ")  )
 
     


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



