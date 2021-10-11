#!/bin/bash

#Peak calling for ChIP-seq data

input1="$(pwd)/peakranger/INPUT1.bam"
input2="$(pwd)/peakranger/INPUT2.bam"

# H3K4me3
for i in 1 2 3 4 5 6 7 8 9 10 ; do
peakranger ranger -d ${dir1}SP${i}.bam  -c ${input1} --format bam -l 316 -t 2 -b 200 -q 0.05 -o H3K4me3_SP${i} ;
done

# H3K27ac
for i in 21 22 23 24 25 26 27 28 29 30 ; do
peakranger ranger -d ${dir1}SP${i}.bam  -c ${input1} --format bam -l 316 -t 2 -b 200 -q 0.05 -o H3K27ac_SP${i} ;
done

# H3K27me3
for i in 101 102 103 104 105 ; do
peakranger ccat -d ${dir2}SP${i}.bam   -c ${input2} --format bam -l 316 -t 2 --win_size 1000 --win_step 100 --min_count 70 --min_score 7 -q 0.05 -o H3K27me3_SP${i} ;
done

# H3K27me3
for i in  107 108 109 110 111 ; do
peakranger ccat -d ${dir3}SP${i}.bam   -c ${input2} --format bam -l 316 -t 2 --win_size 1000 --win_step 100 --min_count 70 --min_score 7 -q 0.05 -o H3K27me3_SP${i} ;
done

# H3K4me1
for i in 31 32 33 34 35 36 37 38 39 40 ; do
peakranger bcp -d ${dir1}SP${i}.bam  -c ${input1} --format bam -l 316 -o H3K4me1_SP${i} ;
done

# H3K36me3
for i in 41 42 43 44 45 46 47 48 49 50   ; do
peakranger bcp -d ${dir1}SP${i}.bam  -c ${input1} --format bam -l 316 -o H3K36me3_SP${i} ;
done

