#!/bin/bash

# pyDNase-0.2.4


OUTDIR="../wellington";
BED="../wellington/extended.bed"  # EPIGENODE_merged_ATAC_peaks.bed"  
#BED is extended 50bp each colum (start and end)

#Get footprints
for i in 21-22 ; do #1-2 9-10 11-12 13-14 15-16 ; do
bsub -o ${OUTDIR}/o.wellington${i} -e ${OUTDIR}/e.wellington${i} -sp 100 -G team170 -q long -n32 -R"select[mem>32000] rusage[mem=32000] span[hosts=1]" -M32000 " ./wellington_footprints.py -p 32 -A -fdr 0.1 --FDR_limit -4 --pv_cutoffs -4 -fp 4,30,1 -fdriter 500 --one_dimension -o atac_${i} ${BED} /lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC2/mergedBam/ATAC${i}.bam ${OUTDIR}/atac_${i} " ;
done

for i in 1-2 9-10 11-12 13-14 15-16 3-4 5-6 7-8 ; do
bsub -o ${OUTDIR}/o.wellington${i} -e ${OUTDIR}/e.wellington${i} -sp 100 -G team170 -q long -n32 -R"select[mem>32000] rusage[mem=32000] span[hosts=1]" -M32000 " ./wellington_footprints.py -p 32 -A -fdr 0.1 --FDR_limit -4 --pv_cutoffs -4 -fp 4,30,1 -fdriter 500 --one_dimension -o atac_${i} ${BED} /lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC/mergedBams/ATAC${i}.bam ${OUTDIR}/atac_${i} " ;
done

for i in 17 ; do
bsub -o ${OUTDIR}/o.wellington${i} -e ${OUTDIR}/e.wellington${i} -sp 100 -G team170 -q long -n32 -R"select[mem>32000] rusage[mem=32000] span[hosts=1]" -M32000 " ./wellington_footprints.py -p 32 -A -fdr 0.1 --FDR_limit -4 --pv_cutoffs -4 -fp 4,30,1 -fdriter 500 --one_dimension -o atac_${i} ${BED} /lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC/mergedBams/ATAC${i}.bam ${OUTDIR}/atac_${i} " ;
done





