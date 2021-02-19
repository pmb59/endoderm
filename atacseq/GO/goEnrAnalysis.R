library(org.Hs.eg.db)
library(GO.db)
library(GOstats)
library(KEGG.db)
library(limma)
library(KEGG.db)
library(biomaRt)
library(ggplot2)
library(topGO)

FILE <- 'h12_vs_h0_FC_2_adjP_1e-04_diffNGS'    
regions <- 'DOWN'  # UP or DOWN
#======================================================================================

x <- read.csv(paste(FILE,"csv", sep='.'), header=T)

if (regions =='UP'){
  x <- x[ which(x$diff == 'TRUE.Up'), ]
  print(dim(x))
}
if (regions =='DOWN'){
  x <- x[ which(x$diff == 'TRUE.Down'), ]
  print(dim(x))
}

newX <- c()
for (i in 1:nrow(x)){
  print(i)
  temp <- strsplit( as.character(x$gene_names[i]), split=';')[[1]]
  if (length(temp) > 0 ){
    newX <- c(newX , temp)
  }
  rm(temp)
}


# Get gene names #
hgnc_symbol <- newX 

results <- data.frame(hgnc_symbol =hgnc_symbol ,One=1 )

ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="oct2014.archive.ensembl.org")
my_chr <- c(1:22, 'X', 'Y')

SELECTED <- as.character(results$hgnc_symbol )

my_entrez_gene <- getBM(attributes='entrezgene',
                        filters = 'hgnc_symbol',
                        values = SELECTED,
                        mart = ensembl)

universe_entrez_gene <- getBM(attributes='entrezgene',
                              filters = 'chromosome_name',
                              values = my_chr,
                              mart = ensembl)

# GOHyperGParams
params <- new('GOHyperGParams',
              geneIds=unique(my_entrez_gene),
              universeGeneIds=unique(universe_entrez_gene),
              ontology='BP',   
              pvalueCutoff=0.001,  
              conditional=FALSE,  
              testDirection='over',     
              annotation="org.Hs.eg.db"
)


hgOver <- hyperGTest(params)

result <- summary(hgOver)
head(result,10)

write.table(x=result,   file = paste(FILE,'GOstats_pvalueCutoff=0.001',regions,'txt',sep='.'),   append = FALSE, quote = F, sep = "\t", row.names = F, col.names = T)


