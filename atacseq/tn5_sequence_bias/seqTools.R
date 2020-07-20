
library(seqTools)

res <- countFastaKmers(filenames = "pos_Tn5.fa", k=6)


#Build "SeqBias.txt" file for ATAC-seq
x<-read.table("GRCh38_Kmers.txt", head=F)  #relative freq of each 6-mer in the human genome

N <- sum(res[,1]) 
res[,1] <-  ( res[,1] / N ) / x$V2    # 100 *
#See Yardimci et al, NAR


write.table(x=res, file = "SeqBiasTn5.txt", append = FALSE, quote = F, sep = " ", row.names = T, col.names = F) ;



