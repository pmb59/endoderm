#!/bin/bash

chmod +x bedtools

##outDir="/lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/FootprintMixture/Footprints_plots"
#outDir="/lustre/scratch110/sanger/pm12/Footprints_plots"
outDir="/lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC2/FootprintMixture/Footprints_plots"

#motifs  in /fimo_occurrences_60k/
#cd /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/fimo/fimo_occurrences_500k  #_60k ;
cd /lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC/fimo/fimo_occurrences_500k

n=$(ls *txt | sed -e 's/\_..*$//' )

#cd /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/FootprintMixture
cd /lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC2/FootprintMixture

#M3334=GATA6
#M4432=CTCF
#M1018=HLX

for SAMPLE in ATAC21 ATAC22 ; do   # ATAC21 ATAC22

for moFile in $(echo $n) ; do    # M1957  M0085 M0820 M0415  ; do

bsub -o ${outDir}/o.${moFile}_1.02_${SAMPLE} -e ${outDir}/e.${moFile}_1.02_${SAMPLE} -sp 100 -G team170 -q long -n1 -R"select[mem>64000] rusage[mem=64000] span[hosts=1]" -M64000 " cd /lustre/scratch117/cellgen/team204/pm12/scratch109/R/R-3.2.1/bin;  ./Rscript /lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC2/FootprintMixture/FootprintMixture.R ${moFile}_1.02 ${SAMPLE} "

done

done

#End here for ATAC 36h (EPIGENODE-ATAC2)
##############################
##############################

#64000 RAM

##############################
#### RESEND FAILED JOBS
#############################
outDir="/lustre/scratch110/sanger/pm12/Footprints_plots"

moFile='M6528'
for SAMPLE in ATAC1 ATAC2 ; do
bsub -o ${outDir}/o.${moFile}_1.02_${SAMPLE} -e ${outDir}/e.${moFile}_1.02_${SAMPLE} -sp 100 -G team170 -q yesterday -n1 -R"select[mem>64000] rusage[mem=64000] span[hosts=1]" -M64000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/FootprintMixture/FootprintMixture.R ${moFile}_1.02 ${SAMPLE}  "
done

moFile='M2947'
SAMPLE='ATAC2'
bsub -o ${outDir}/o.${moFile}_1.02_${SAMPLE} -e ${outDir}/e.${moFile}_1.02_${SAMPLE} -sp 100 -G team170 -q yesterday -n1 -R"select[mem>64000] rusage[mem=64000] span[hosts=1]" -M64000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/FootprintMixture/FootprintMixture.R ${moFile}_1.02 ${SAMPLE}  "

moFile='M5887'
SAMPLE='ATAC2'
bsub -o ${outDir}/o.${moFile}_1.02_${SAMPLE} -e ${outDir}/e.${moFile}_1.02_${SAMPLE} -sp 100 -G team170 -q yesterday -n1 -R"select[mem>64000] rusage[mem=64000] span[hosts=1]" -M64000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/FootprintMixture/FootprintMixture.R ${moFile}_1.02 ${SAMPLE}  "

moFile='M5763'
SAMPLE='ATAC2'
bsub -o ${outDir}/o.${moFile}_1.02_${SAMPLE} -e ${outDir}/e.${moFile}_1.02_${SAMPLE} -sp 100 -G team170 -q yesterday -n1 -R"select[mem>64000] rusage[mem=64000] span[hosts=1]" -M64000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/FootprintMixture/FootprintMixture.R ${moFile}_1.02 ${SAMPLE}  "


moFile='M6404'
SAMPLE='ATAC10'
bsub -o ${outDir}/o.${moFile}_1.02_${SAMPLE} -e ${outDir}/e.${moFile}_1.02_${SAMPLE} -sp 100 -G team170 -q yesterday -n1 -R"select[mem>64000] rusage[mem=64000] span[hosts=1]" -M64000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/FootprintMixture/FootprintMixture.R ${moFile}_1.02 ${SAMPLE}  "

moFile='M6198'
SAMPLE='ATAC11'
bsub -o ${outDir}/o.${moFile}_1.02_${SAMPLE} -e ${outDir}/e.${moFile}_1.02_${SAMPLE} -sp 100 -G team170 -q yesterday -n1 -R"select[mem>64000] rusage[mem=64000] span[hosts=1]" -M64000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/FootprintMixture/FootprintMixture.R ${moFile}_1.02 ${SAMPLE}  "

moFile='M2296'
SAMPLE='ATAC13'
bsub -o ${outDir}/o.${moFile}_1.02_${SAMPLE} -e ${outDir}/e.${moFile}_1.02_${SAMPLE} -sp 100 -G team170 -q yesterday -n1 -R"select[mem>64000] rusage[mem=64000] span[hosts=1]" -M64000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/FootprintMixture/FootprintMixture.R ${moFile}_1.02 ${SAMPLE}  "

moFile='M5519'
SAMPLE='ATAC14'
bsub -o ${outDir}/o.${moFile}_1.02_${SAMPLE} -e ${outDir}/e.${moFile}_1.02_${SAMPLE} -sp 100 -G team170 -q yesterday -n1 -R"select[mem>64000] rusage[mem=64000] span[hosts=1]" -M64000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/FootprintMixture/FootprintMixture.R ${moFile}_1.02 ${SAMPLE}  "

moFile='M3004'
SAMPLE='ATAC16'
bsub -o ${outDir}/o.${moFile}_1.02_${SAMPLE} -e ${outDir}/e.${moFile}_1.02_${SAMPLE} -sp 100 -G team170 -q yesterday -n1 -R"select[mem>64000] rusage[mem=64000] span[hosts=1]" -M64000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/FootprintMixture/FootprintMixture.R ${moFile}_1.02 ${SAMPLE}  "

moFile='M3009'
SAMPLE='ATAC16'
bsub -o ${outDir}/o.${moFile}_1.02_${SAMPLE} -e ${outDir}/e.${moFile}_1.02_${SAMPLE} -sp 100 -G team170 -q yesterday -n1 -R"select[mem>64000] rusage[mem=64000] span[hosts=1]" -M64000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/FootprintMixture/FootprintMixture.R ${moFile}_1.02 ${SAMPLE}  "

#done

