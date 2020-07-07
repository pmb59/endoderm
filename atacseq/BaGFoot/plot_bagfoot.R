setwd("/Users/pmb59/Desktop/unt")

rm( list=objects() )

library(hash)
library(data.table)
library(digest)
library(Cairo)
library(aplpack)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Mmusculus.UCSC.mm9)
library(bagfoot)

source('aux.R')

#read bagfoot output
dat <- read.table('T24_T48_On_consolidated_hotspot_footprint_depth_table.csv');
gen_bagplot_chisq(dat, dataname1='T24', dataname2='T48')


df <- read.csv("bagplot_cutcount_diff_total_footprinting_depth_T48-T24_qvalue_bagplot_output.csv")
head(df)
dim(df)

# first we filter the table to leave only genes differentially expressed during differentiation
f1 <- as.character(read.csv("/Users/pmb59/Desktop/unt/DESeq2.Diff_Gene.Expression/signif.-h0_vs_h12.csv")$gene_name)
f2 <- as.character(read.csv("/Users/pmb59/Desktop/unt/DESeq2.Diff_Gene.Expression/signif.-h12_vs_h24.csv")$gene_name)
f3 <- as.character(read.csv("/Users/pmb59/Desktop/unt/DESeq2.Diff_Gene.Expression/signif.-h24_vs_h36.csv")$gene_name)
f4 <- as.character(read.csv("/Users/pmb59/Desktop/unt/DESeq2.Diff_Gene.Expression/signif.-h36_vs_h48.csv")$gene_name)
f5 <- as.character(read.csv("/Users/pmb59/Desktop/unt/DESeq2.Diff_Gene.Expression/signif.-h48_vs_h60.csv")$gene_name)
f6 <- as.character(read.csv("/Users/pmb59/Desktop/unt/DESeq2.Diff_Gene.Expression/signif.-h60_vs_h72.csv")$gene_name)
GOI <- unique(c(f1,f2,f3,f4,f5,f6))
GOI


keep <- c()
for (i in 1:nrow(df)){
 temp<-  which( GOI  == as.character( df$name[i] )  ) 
  if(length(temp)>0){ keep <- c(keep,i) }
 rm(temp)
}

df <- df[unique(keep), ]
head(df)
dim(df)

#gen_bagplot(dat, dataname1='T24', dataname2='T48', factor=2.5, noFence=F)
#gen_bagplot_chisq(dat, dataname1='T24', dataname2='T48')


df <- df [ which(df$name!='ARID5A'), ]

df$name[ which(df$name!='EOMES'    & df$name!='MIXL1'   & df$name!='GATA6'  & df$name!='JUNB' 
               & df$name!='POU5F1' & df$name!='JUND'    & df$name!='FOS'    & df$name!='BACH2' 
               & df$name!='FOSB'   & df$name!='GATA2'   & df$name!='GATA3'  & df$name!='GATA4'
               & df$name!='GATA5'  & df$name!='FOSL2'   & df$name!='TFAP2C' & df$name!='TFAP2A' 
               & df$name!='CTCF'   & df$name!='FOXA2'   & df$name!='TP53' & df$name!='T' 
               & df$name!='FOSL1'  & df$name!='FOXO1'   & df$name!='FOXO3'
               & df$name!='FOXO4'  & df$name!='SMARCC1' & df$name!='NRF1') ] <- ''
#df$name[ which(df$pvalue>0.01) ] <- ''

library(ggplot2)
library(ggrepel)
# load the functions from this Gist:
devtools::source_gist("00772ccea2dd0b0f1745", filename = "000_geom_bag.r")
devtools::source_gist("00772ccea2dd0b0f1745", filename = "001_bag_functions.r")


ggplot(df, aes(deltahyp, deltafd , label=name)) +   #colour = category, fill = category
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_bag() +
  theme_minimal()  +
  geom_point(aes( color=pvalue   )) +  #, shape=category
  geom_text_repel( color='black',
                   segment.size = 0.2,force=5,
                    segment.color = "red"
  )+
  xlab( expression(paste('Difference in flanking chromatin accessibility (',Delta,'FA)' ) ) ) + 
  ylab(expression(paste('Difference in footprint depth (',-Delta,'FPD)' ) ) )+
  scale_color_gradient2(midpoint=0.4, low="tomato", mid="royalblue2",
                         high="darkgray", space ="Lab" , breaks=c(0, 0.25, 0.5, 0.75, 1), limits=c(0,1) )+
  ggtitle('Endoderm differentiation [Time 48h] - [Time 24h]')+
  labs(color='P-value') 
ggsave('bagfoot_T48_T24.pdf', height=6, width=7.6)
ggsave('bagfoot_T48_T24.tiff', height=6, width=7.6)


