#!/bin/bash

# First run this script:
# Rscript bias_1.r (bias_1.r)

##merge all .bed into a single file
#cat *pos_Tn5_*.bed > pos_Tn5.bed
##cat *pos_Tn5_forward.bed > pos_Tn5_forward.bed
##cat *pos_Tn5_reverse.bed > pos_Tn5_reverse.bed

chmod +x bedtools
genome="/lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/Homo_sapiens.GRCh38_15.fa"

#./bedtools getfasta -s -fi $genome -bed pos_Tn5.bed -fo pos_Tn5.fa
##Not necessary with jellyfish
##sed -i 's/^.*chr.*$/>chr/' pos_Tn5.fa

#./bedtools getfasta -s -fi $genome -bed pos_Tn5_forward.bed -fo pos_Tn5_forward.fa
#./bedtools getfasta -s -fi $genome -bed pos_Tn5_reverse.bed -fo pos_Tn5_reverse.fa
#sed -i 's/^.*chr.*$/>chr/' pos_Tn5_forward.fa
#sed -i 's/^.*chr.*$/>chr/' pos_Tn5_reverse.fa

DIRDATA="/lustre/scratch109/sanger/pm12/biasEPIGENODE"

#cat FASTAs

#Create SeqBias File
##cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin
##./Rscript /lustre/scratch109/sanger/pm12/biasEPIGENODE/seqTools.R
./jellyfish count -m 6 -s 3G -t 1 pos_Tn5.fa  # -C
./jellyfish dump mer_counts.jf > mer_counts_dumps.fa

#process jellyfish output
cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin
./Rscript /lustre/scratch109/sanger/pm12/biasEPIGENODE/jellyfish.R


