

#Build "SeqBias.txt"(FootprintMixture) file for ATAC-seq
x<-read.table("GRCh38_Kmers.txt", head=F)  #relative freq of each 6-mer in the human genome - sum(x$V2)


#jellyfish output processing (y)
y<- x
jelly <- readLines("mer_counts_dumps.fa")


L<-strsplit(jelly, split=">")
kmer_jelly <- rep(0,length(L)/2)
freq_jelly <- rep(0,length(L)/2)
for (i in 1:length(L)) {
	if( length( L[[i]]) == 1  ) {kmer_jelly <-c(kmer_jelly,  L[[i]][1] )  }
	if( length( L[[i]]) == 2  ) {freq_jelly <-c(freq_jelly,  L[[i]][2] )  }
}
#length(freq_jelly)
#length(kmer_jelly)

y$V2 <- 0.0


for (i in 1:length(y$V1)) {
	ix <- which( kmer_jelly ==   as.character(y$V1[i]))
	if (length(ix) > 0){
		y$V2[i] <- freq_jelly[ix]
	}
}


N <- sum(as.numeric(y$V2)) 

y$V2 <-  ( as.numeric(y$V2) / N ) / x$V2    # 100 *
#See Yardimci et al, NAR


write.table(x=y, file = "SeqBiasTn5.txt", append = FALSE, quote = F, sep = " ", row.names = F, col.names = F) ;



