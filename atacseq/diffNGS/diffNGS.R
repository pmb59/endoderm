
################################################################################################################
#  diffNGS
#  http://github.com/pmb59/diffNGS
################################################################################################################
# Pedro Madrigal { pmb59 [at] cam.ac.uk }
# 
# This is a modified version of function narrowpeaksDiff.R of the Bioconductor package 
# NarrowPeaks: Shape-based Analysis of Variation in ChIP-seq using Functional PCA
# http://bioconductor.org/packages/devel/bioc/html/NarrowPeaks.html
#
# It quantifies differences between the shapes, and uses Hotelling's T2 tests on the functional principal 
# component scores to identify significant differences across conditions. An application of the package 
# for Arabidopsis datasets is described in Mateos, Madrigal, et al. (2015) Genome Biology: 16:31.
#
# High RAM mem is required if fdamatrix object is large #
################################################################################################################

diffNGS <- function(bedFile,  headerBed= TRUE, bigwigs , conditions , nbasis=50, pcs = 10, variation = 0.6, NB=20) {   
 
	fdamatrix <- list()  #To store signal profiles of all datasets

	#Extend peak centre Upstream and Dowstream 'flank' bp
	peaks <- readGeneric(bedFile, keep.all.metadata = FALSE)
	start(peaks) <- start(peaks)  # - flank
	end(peaks)   <- end(peaks)    # + flank
	# Parameters in genomation for signal extraction form the bigwig
	nBins = NB #(2*NB)+1  #(2*flank)+1
	scaleData = FALSE  #Signal is already normalized

	#Read bigwigs
	for (k in 1:length(bigwigs) ){
    
    		print(paste( paste("Reading bigWig File...", bigwigs[k], sep="") , Sys.time() , sep=" @ ")    )
	  	fdamatrix[[k]] <- matrix(0.0, ncol=NB, nrow= length(peaks) )    #ncol=1+2*NB
	  	fdamatrix[[k]]  <- ScoreMatrixBin(target = bigwigs[k], bin.num = nBins, windows = peaks, type="bigWig",rpm=F, bin.op="max" )  #bin.op="max"
    	  	# here we have the matrix with all the profiles
    	  	# correct non-numeric values in case of NAs
        	correctval <- 1e-3 #1e-2  # To avoid numerical problems
    	  	fdamatrix[[k]][which(is.na(fdamatrix[[k]])==TRUE)] <- correctval
    	  	fdamatrix[[k]][which(is.nan(fdamatrix[[k]])==TRUE)] <- correctval
    	  	fdamatrix[[k]][which(is.numeric(fdamatrix[[k]])==FALSE)] <- correctval
	}
  
  	#plot
  	Ylim =0
  	for (k in 1:length(bigwigs) ){ Ylim <- max(Ylim, max( colMeans(fdamatrix[[k]],na.rm=TRUE) )  )  } 
	
  	pdf(file=paste( paste(rev(unique(CNDS)),collapse="_vs_"), "mean", "pdf"  , sep='.')  )
  	par(mfrow=c(2,2))
  	cls <- sort(rep(brewer.pal(4,"Accent")[1:length(unique(CNDS))]  ,2 ))  
  	for (k in 1:length(bigwigs) ){ 
  		plot(colMeans(fdamatrix[[k]],na.rm=TRUE),type="l", ylim=c(0,Ylim), lwd=3, main=conditions[k], xlab=paste(paste("#bins =",nBins, sep=" "), "; scaled region", sep=" " )  , ylab="normalized signal", cex.lab=1, col=cls[k] )
  		mtext(bigwigs[k])
  	}
  	dev.off()
  	rm (Ylim)
  
	#Create list of data.frames for the p-values
	PVALS <- list()
	uniqueCond <- unique(conditions)
	for(j in 1:length(peaks) ){
	  	Mtemp <- matrix(NA, ncol=length(uniqueCond), nrow=length(uniqueCond))
		colnames(Mtemp) <- uniqueCond
		rownames(Mtemp) <- uniqueCond
		PVALS[[j]] <- as.data.frame(Mtemp)
		rm(Mtemp)
	}
	
	#Create B-spline basis (order 4) of the signals before FPCA:
	bspl <- create.bspline.basis(rangeval=c(-NB,NB),nbasis=nbasis, norder=4)   # Cubic B-splines
	argvalsBS <- -floor(NB/2):floor(NB/2)
	
	for(j in 1:length(peaks) ){
    		print(paste( paste("Analyzing peak region #",j,sep=""), Sys.time() , sep=" @ ")  )  
		tempMatrix <- matrix(0.0, ncol=NB, nrow= length(bigwigs) )                     
		for (m in 1:length(bigwigs)) {
      			if (sum (fdamatrix[[m]][j,]) == 0) {   tempMatrix[m,] <- runif(length( fdamatrix[[m]][j,] ), 0.0, 1e-3)   }  
      			if (sum (fdamatrix[[m]][j,]) != 0) {   tempMatrix[m,] <- fdamatrix[[m]][j,]  }
		}
    
		#to avoid numeric problems in pca
		if ( length(range(tempMatrix)) > 1   ){
    
		fdaData <- Data2fd(y=t(tempMatrix), argvals= argvalsBS, basisobj=bspl)		
		pc <- pca.fd(fdobj=fdaData, nharm = pcs, harmfdPar=fdPar(fdaData),centerfns = FALSE)
		
		# Select the PC scores for the components which amount >= 'variation'
		# print(pc$varprop)
		CS <- cumsum(pc$varprop)
		requiredPCs <- which(CS >= variation)[1]
		
		# Hotelling's T2 test between the conditions provided
		for (o in 1:length(uniqueCond)){
			for (p in 1:length(uniqueCond)){
				Xid <- which(conditions == uniqueCond[o] )
				Yid <- which(conditions == uniqueCond[p] )
				PVALS[[j]][o,p]  <- HotellingsT2(X=as.matrix(pc$scores[Xid ,1:requiredPCs],ncol=requiredPCs),Y=as.matrix(pc$scores[Yid,1:requiredPCs],ncol=requiredPCs), test="chi")$p.value
			}
		}
	}
		
	if ( length(range(tempMatrix)) < 2   ){
  		for (o in 1:length(uniqueCond)){ for (p in 1:length(uniqueCond)){  PVALS[[j]][o,p]  <- 1.0  } }
	}  
  

	# report the P-values in a list of data.frames (PVALS)		
	}		
	
	return(list(fdaprofiles=fdamatrix, p.values = PVALS ) ) # report signal and raw p-values 


}     


