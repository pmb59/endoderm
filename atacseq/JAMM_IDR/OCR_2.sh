#!/bin/bash

#working dir
WDIR=""


#IDR analaysis for replicates
j=$1    #"1"
k=$2    #"2"



#------------------------------------------
# IDR analysis for all/filtered peaks (NEW SCRIPT)
#------------------------------------------
# Then truncate the number of peaks to the top 100k-125k. Using more than this simply increases the running time of the IDR analysis with no advantage.
# Select as much as peaks as possible
#
#Toy example:
#A=$(wc -l ${WDIR}jamm_${j}_${k}/rep1/peaks/filtered.peaks.narrowPeak);
#A=$(wc -l ${WDIR}jamm_${j}_${k}/rep1/peaks/filtered.peaks.narrowPeak);
#e.cram2bam_16424_1
# $(Rscript -e 'args<-commandArgs(TRUE); cat(min(length(readLines(args[1])) ,length(readLines(args[2])) ,4 ))' )


#A=$(wc -l ${WDIR}jamm_${j}_${k}/rep1/peaks/filtered.peaks.narrowPeak);
#B=$(wc -l ${WDIR}jamm_${j}_${k}/rep2/peaks/filtered.peaks.narrowPeak);
##A=${WDIR}"jamm_"${j}"_"${k}"/rep1/peaks/filtered.peaks.narrowPeak"
##B=${WDIR}"jamm_"${j}"_"${k}"/rep2/peaks/filtered.peaks.narrowPeak"

#echo $A
#echo $B;
# $(Rscript -e 'args<-commandArgs(TRUE); cat(min(args[1],args[2]), sep="\n")' $A $B)
# add min 300k
sort -k 7nr,7nr ${WDIR}jamm_${j}_${k}/rep1/peaks/filtered.peaks.narrowPeak | head -n $(Rscript -e 'args<-commandArgs(TRUE); cat(min(length(readLines(args[1])) ,length(readLines(args[2])), 4000000) )' ${WDIR}jamm_${j}_${k}/rep1/peaks/filtered.peaks.narrowPeak ${WDIR}jamm_${j}_${k}/rep2/peaks/filtered.peaks.narrowPeak ) | gzip -c > ${WDIR}jamm_${j}_${k}/rep1/peaks/Ranked.peaks.gz
sort -k 7nr,7nr ${WDIR}jamm_${j}_${k}/rep2/peaks/filtered.peaks.narrowPeak | head -n $(Rscript -e 'args<-commandArgs(TRUE); cat(min(length(readLines(args[1])) ,length(readLines(args[2])), 4000000) )' ${WDIR}jamm_${j}_${k}/rep1/peaks/filtered.peaks.narrowPeak ${WDIR}jamm_${j}_${k}/rep2/peaks/filtered.peaks.narrowPeak ) | gzip -c > ${WDIR}jamm_${j}_${k}/rep2/peaks/Ranked.peaks.gz

gunzip ${WDIR}jamm_${j}_${k}/rep1/peaks/Ranked.peaks.gz
gunzip ${WDIR}jamm_${j}_${k}/rep2/peaks/Ranked.peaks.gz
#head ${WDIR}jamm_${j}_${k}/rep1/peaks/Ranked.peaks
#head ${WDIR}jamm_${j}_${k}/rep2/peaks/Ranked.peaks



# IDR analysis  (https://sites.google.com/site/anshulkundaje/projects/idr)
#Typical usage for MACSv2 peak caller peaks:
#Rscript batch-consistency-analysis.r [peakfile1] [peakfile2] -1 [outfile.prefix] 0 F p.value
#cd /Users/pm12/Desktop/temp/IDR/idrCode

##modify manually genome_table.txt##
## do it just once  ###
#rm /lustre/scratch109/sanger/pm12/IDR/idrCode/genome_table.txt
#cp GRCh38_15.chrom.sizes /lustre/scratch109/sanger/pm12/IDR/idrCode/genome_table.txt
cd /lustre/scratch117/cellgen/team204/pm12/scratch109/IDR/idrCode
Rscript batch-consistency-analysis.r ${WDIR}jamm_${j}_${k}/rep1/peaks/Ranked.peaks ${WDIR}jamm_${j}_${k}/rep2/peaks/Ranked.peaks -1 ${WDIR}jamm_${j}_${k}/IDR 0 F signal.value

#remove intermediate files
rm ${WDIR}jamm_${j}_${k}/rep1/peaks/Ranked.peaks
rm ${WDIR}jamm_${j}_${k}/rep2/peaks/Ranked.peaks

#You can plot the IDR plots using
#This is used to plot the IDR plots and diagnostic plots for a single or multiple pairs of replicates.
#Rscript batch-consistency-plot.r [npairs] [output.prefix] [input.file.prefix1] [input.file.prefix2] [input.file.prefix3] ....

Rscript batch-consistency-plot.r 1 ${WDIR}jamm_${j}_${k}/IDRplot ${WDIR}jamm_${j}_${k}/IDR


cd $WDIR
##Peak the top n peaks in Replicate integration at IDR<0.05
N=$(grep "0.05" ${WDIR}jamm_${j}_${k}/IDR-npeaks-aboveIDR.txt | awk '{print $4}' )
echo $N
# Final Set of Peaks
sort -k 7nr,7nr ${WDIR}jamm_${j}_${k}/peaks/filtered.peaks.narrowPeak | head -n $N | gzip -c > ${WDIR}jamm_${j}_${k}/peaks_IDR_0.05_jamm_${j}_${k}.narrowPeak.gz
gunzip ${WDIR}jamm_${j}_${k}/peaks_IDR_0.05_jamm_${j}_${k}.narrowPeak.gz


#ENCODE blacklist filtering (Final List)
#chmod +x bedtools
#sort first? Not necessary
# ./bedtools intersect -wa -u -a ${WDIR}jamm_${j}_${k}/peaks_IDR_0.05_jamm_${j}_${k}.narrowPeak -b /lustre/scratch109/sanger/pm12/MapabilityConsensusExcludable/hg38_wgEncodeDacMapabilityConsensusExcludable.bed > ${WDIR}jamm_${j}_${k}/Final_peaks_IDR_0.05_jamm_${j}_${k}.narrowPeak




