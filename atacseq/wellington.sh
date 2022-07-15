#!/bin/bash

# pyDNase-0.2.4

OUTDIR="$(pwd)/wellington";
BED="$(pwd)/wellington/extended.bed"  # EPIGENODE_merged_ATAC_peaks.bed, extended 50bp each colum (start and end)

#Get Wellington footprints

for i in 1-2 9-10 11-12 13-14 15-16 21-22 17 ; do
    wellington_footprints.py -p 32 -A -fdr 0.1 --FDR_limit -4 --pv_cutoffs -4 -fp 4,30,1 -fdriter 500 --one_dimension -o atac_${i} ${BED} ../mergedBams/ATAC${i}.bam ${OUTDIR}/atac_${i}  ;
done

# Output files can be found at
# http://ngs.sanger.ac.uk/production/endoderm/

