# to create 'gene_expression_FPKM.csv'

# Read counts
fc <- read.table(file='featureCountsEndo.txt' , head=TRUE, sep='\t')
head(fc)

fcounts <- fc[, 7:ncol(fc)]
dim(fcounts)


rownames(fcounts) <- fc$Geneid
head(fcounts)
dim(fcounts)

cinfo <- data.frame(SAMPLE=rep('s', ncol(fcounts) ) )
rownames(cinfo) <- colnames(fcounts )
cinfo


library("NOISeq")   
mylength <- as.vector(fc$Length)
names(mylength) <- fc$Geneid
dim(mylength)
head(mylength)  #dim(null)


mydata <- readData(data = fcounts, length = mylength, factors = cinfo)
mydata

myRPKM = rpkm(assayData(mydata)$exprs, long = mylength, k = 0, lc = 1)

#colnames(myRPKM ) <- paste( cinfo$celltype, colnames(myRPKM )  ,sep='.')
head(myRPKM)  # RPKMs
dim(myRPKM)

FPKM <- data.frame( gene_id = rownames(myRPKM), 
                    h0 = rowMeans(myRPKM[,1:3]) ,       
                    h12 = rowMeans(myRPKM[,4:6]) ,          
                    h24 = rowMeans(myRPKM[,7:9]) ,          
                    h36 = rowMeans(myRPKM[,10:12]) ,          
                    h48 = rowMeans(myRPKM[,13:15]) ,  
                    h60 = rowMeans(myRPKM[,16:18]) ,  
                    h72 = rowMeans(myRPKM[,19:21])    
                    )
head(FPKM)

# remove genes with 0 FPKM
dim(FPKM)
FPKM <- FPKM[ which(rowMeans(FPKM[,2:8]) != 0 ) , ]
dim(FPKM)




#save gene expression matrix
library("biomaRt")   
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl = useMart(host='http://oct2014.archive.ensembl.org',
                  biomart='ENSEMBL_MART_ENSEMBL',
                  dataset='hsapiens_gene_ensembl')

listAttributes(ensembl)

bm <- getBM(attributes=c('ensembl_gene_id','external_gene_name','gene_biotype' ),
            mart = ensembl)
colnames(bm) <- c('gene_id' ,'gene_name' ,'gene_biotype')
head(bm)
dim(bm)

x <- merge(x=FPKM, y=bm, by='gene_id')
head(x)


#save table with FPKMs (only genes with at least 1 count)
#write.csv(x, file='FPKM_endodiff.csv', row.names=FALSE)


#get prot_coding
x <- x[ which(x$gene_biotype == 'protein_coding'), ]
dim(x)
head(x)

x$gene_name <- as.character(x$gene_name)


# Select Genes of Interest
sel <- which(x$gene_name=="SOX2" | x$gene_name=="NANOG" | x$gene_name=="PRDM14" | x$gene_name=="POU5F1" |  x$gene_name=="OTX2" | 
			#early induced
			x$gene_name=="ID1" | x$gene_name=="ID2" | x$gene_name=="ID3" | x$gene_name=="ID4" | x$gene_name=="DACT1" | x$gene_name=="BAMBI" | x$gene_name=="NOTUM" | x$gene_name=="WLS" | x$gene_name=="LEFTY1" | 
			#mesendo
			x$gene_name=="EOMES" | x$gene_name=="T"  |x$gene_name=="MIXL1"|x$gene_name=="FOXH1" |x$gene_name=="CER1" |
			
			#endoderm
			 x$gene_name=="LZTS1" | x$gene_name=="SOX17" | x$gene_name=="GATA6" | x$gene_name=="FOXA2"| x$gene_name=="GSC"|x$gene_name=="LHX1"|x$gene_name=="HHEX"|x$gene_name=="TBX3"|x$gene_name=="CXCR4"|x$gene_name=="FOXC1"|x$gene_name=="GATA4"|
			#mesoderm
			x$gene_name=="HAND1" | x$gene_name=="RUNX1" |x$gene_name=="GATA1" |x$gene_name=="TBX5" |x$gene_name=="NKX2-5" |
			#ectoderm
			x$gene_name=="PAX6" | x$gene_name=="SOX1"|  x$gene_name=="OLIG3" | x$gene_name=="GBX2" |
			#Other 
			 x$gene_name=="CTCF" | x$gene_name=="TP53" | x$gene_name=="NRF1" | x$gene_name=="TFAP2A" | x$gene_name=="TFAP2C"|
       #HK
       x$gene_name=="PGK1" | x$gene_name=="GAPDH" | x$gene_name=="LDHA" |
      #AP-1 
			 x$gene_name=="JUN" | x$gene_name=="JUNB" | x$gene_name=="JUND" | x$gene_name=="FOS"
			 

 )
 
 x$gene_name[sel]


split <- rep('Endoderm differentiation',length(sel)) 

split[c(12,28, 36,44,49  )] <- "Pluri"
split[c( 29,32:34, 47 )] <- "Mesendo"
split[c(2,3,17,18,21,23,24,26, 30,35,  51)] <- "Endo"
split[c( 7,11,22  )] <- "HouseKeep"
split[c( 1,38, 42,45 )] <- "Ecto"
split[c( 5,8, 13,31, 46 )] <- "Meso"
split[c( 4,9,10,25, 27 )] <- "Other"
split[c( 20, 39,40,43 )] <- "AP1"
split[c( 6,14:16,19,37,41, 48, 50 )] <- "EarlyInd"


x2 <- x[sel,2:8]

y <- t(x2) 
colnames(y) <- x$gene_name[sel]
head(y)


#  New modification (July 2016)############
library('ColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(viridis)

pal= circlize::colorRamp2(c(0, 1, 2, 5, 7), viridis::viridis(5))

pdf('Fig1_heatmap.pdf', width=4, height=15)
Heatmap(t(log2(y+1) ), name = "log2(FPKM+1)" ,cluster_columns = F,cluster_rows= T , split=split, 
         rect_gp = gpar(col = "white", lwd = 1),
        col = pal
    )
dev.off()
        
        






