#!/bin/bash

# Differential ChIP-seq Analysis using DiffReps

################################################################################################
## merge input samples
################################################################################################
samtools merge INPUT1.bam SP51.bam SP52.bam
samtools merge INPUT2.bam SP106.bam SP112.bam

##INDEX BAMS:
samtools index INPUT1.bam
samtools index INPUT2.bam

################################################################################################
# Run DiffReps
################################################################################################
#H3K4me3
# 0_12;
sh analysis.sh SP1 SP2 SP3 SP4 INPUTSall ./0_12_H3K4me3/ 316 600 60 20 30 A
# 12_24;
sh analysis.sh SP3 SP4 SP5 SP6 INPUTSall ./12_24_H3K4me3/ 316 600 60 20 30 A
# 24_36;
sh analysis.sh SP5 SP6 SP113 X INPUTSall ./24_36_H3K4me3/ 316 600 60 20 30 B
# 36_48;
sh analysis.sh SP113 X SP7 SP8 INPUTSall ./36_48_H3K4me3/ 316 600 60 20 30 C
# 48_72;
sh analysis.sh SP7 SP8 SP9 SP10 INPUTSall ./48_72_H3K4me3/ 316 600 60 20 30 A

#H3K27ac
# 0_12;
sh analysis.sh SP21 SP22 SP23 SP24 INPUTSall ./0_12_H3K27ac/ 316 600 60 20 100 A
# 12_24;
sh analysis.sh SP23 SP24 SP25 SP26 INPUTSall ./12_24_H3K27ac/ 316 600 60 20 100 A
# 24_36;
sh analysis.sh SP25 SP26 SP115 X INPUTSall ./24_36_H3K27ac/ 316 600 60 20 100 B
# 36_48;
sh analysis.sh SP115 X SP27 SP28 INPUTSall ./36_48_H3K27ac/ 316 600 60 20 100 C
# 48_72;
sh analysis.sh SP27 SP28 SP29 SP30 INPUTSall ./48_72_H3K27ac/ 316 600 60 20 100 A


#H3K27me3 (36h did not work!)
# 0_12;
sh analysis.sh SP101 SP107 SP102 SP108 INPUTSall ./0_12_H3K27me3/ 316 2000 200 10 50 A
# 12_24;
sh analysis.sh SP102 SP108 SP103 SP109 INPUTSall ./12_24_H3K27me3/ 316 2000 200 10 50 A
# 24_48;
sh analysis.sh SP103 SP109 SP104 SP110 INPUTSall ./24_48_H3K27me3/ 316 2000 200 10 50 A
# 48_72;
sh analysis.sh SP104 SP110 SP105 SP111 INPUTSall ./48_72_H3K27me3/ 316 2000 200 10 50 A


#H3K4me1
# 0_12;
sh analysis.sh SP31 SP32 SP33 SP34 INPUTSall ./0_12_H3K4me1/ 316 6000 300 5 200 A
# 12_24;
sh analysis.sh SP33 SP34 SP35 SP36 INPUTSall ./12_24_H3K4me1/ 316 6000 300 5 200 A
# 24_36;
sh analysis.sh SP35 SP36 SP116 X INPUTSall ./24_36_H3K4me1/ 316 6000 300 5 200 B
# 36_48;
sh analysis.sh SP116 X SP37 SP38 INPUTSall ./36_48_H3K4me1/ 316 6000 300 5 200 C
# 48_72;
sh analysis.sh SP37 SP38 SP39 SP40 INPUTSall ./48_72_H3K4me1/ 316 6000 300 5 200 A


#H3K36me3
# 0_12;
sh analysis.sh SP41 SP42 SP43 SP44 INPUTSall ./0_12_H3K36me3/ 316 10000 200 2 200 A
# 12_24;
sh analysis.sh SP43 SP44 SP45 SP46 INPUTSall ./12_24_H3K36me3/ 316 10000 200 2 200 A
# 24_36;
sh analysis.sh SP45 SP46 SP117 X INPUTSall ./24_36_H3K36me3/ 316 10000 200 2 200 B
# 36_48;
sh analysis.sh SP117 X SP47 SP48 INPUTSall ./36_48_H3K36me3/ 316 10000 200 2 200 C
# 48_72;
sh analysis.sh SP47 SP48 SP49 SP50 INPUTSall ./48_72_H3K36me3/ 316 10000 200 2 200 A



