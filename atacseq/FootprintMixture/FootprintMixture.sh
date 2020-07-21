#!/bin/bash

chmod +x bedtools

outDir="../Footprints_plots"

#motifs  in /fimo_occurrences_60k/
cd /../fimo/fimo_occurrences_500k

n=$(ls *txt | sed -e 's/\_..*$//' )

cd ./FootprintMixture


for SAMPLE in ATAC21 ATAC22 ; do   # samples at 36h

  for moFile in $(echo $n) ; do    # M1957  M0085 M0820 M0415  ; do

    Rscript FootprintMixture.R ${moFile}_1.02 ${SAMPLE} 

  done

done

