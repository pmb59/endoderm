#!/bin/bash

# Find differential footprints using Wellington_bootstrap (pyDNase-0.2.4)


#create directories
for i in h0_h12 h12_h24 h24_h36 h36_h48 h48_h72 ; do   
mkdir bootstrap_${i};
done

DATADIR="../mergedBams/"
OUTDIR="../wellington";
BED="../wellington/extended.bed"  # EPIGENODE_merged_ATAC_peaks.bed, extended 50bp each colum (start and end)


#Get differential footprints

#1
wellington_bootstrap.py -p 16 -A -fdr 0.05 --FDR_limit -10 -fp 4,30,1 -fdriter 100 ${DATADIR}ATAC1-2.bam ${DATADIR}ATAC9-10.bam ${BED} ${OUTDIR}/bootstrap_h0_h12/h0_only ${OUTDIR}/bootstrap_h0_h12/h12_only "
#2
wellington_bootstrap.py -p 16 -A -fdr 0.05 --FDR_limit -10 -fp 4,30,1 -fdriter 100 ${DATADIR}ATAC9-10.bam ${DATADIR}ATAC11-12.bam ${BED} ${OUTDIR}/bootstrap_h12_h24/h12_only ${OUTDIR}/bootstrap_h12_h24/h24_only "
#3
wellington_bootstrap.py -p 16 -A -fdr 0.05 --FDR_limit -10 -fp 4,30,1 -fdriter 100 ${DATADIR}ATAC11-12.bam ${DATADIR}ATAC21-22.bam ${BED} ${OUTDIR}/bootstrap_h24_h36/h24_only ${OUTDIR}/bootstrap_h24_h36/h36_only "
#4
wellington_bootstrap.py -p 16 -A -fdr 0.05 --FDR_limit -10 -fp 4,30,1 -fdriter 100 ${DATADIR}ATAC21-22.bam ${DATADIR}ATAC13-14.bam ${BED} ${OUTDIR}/bootstrap_h36_h48/h36_only ${OUTDIR}/bootstrap_h36_h48/h48_only "
#5
wellington_bootstrap.py -p 16 -A -fdr 0.05 --FDR_limit -10 -fp 4,30,1 -fdriter 100 ${DATADIR}ATAC13-14.bam ${DATADIR}ATAC15-16.bam ${BED} ${OUTDIR}/bootstrap_h48_h72/h48_only ${OUTDIR}/bootstrap_h48_h72/h72_only "


