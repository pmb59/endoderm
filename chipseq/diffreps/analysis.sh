#!/bin/bash

# Submit one job per comparison

export PATH=/software/perl-5.8.9/bin:$PATH

file1=$1
file2=$2
file3=$3
file4=$4
ctrl=$5
wd=$6
frac=$7
W=$8
ST=$9
NSD=${10}
#MODE=$11
GAP=${11}
#sample configuration: (A: 2_vs_2; B: 2_vs_1; C: 1_vs_2)
SampleN=${12}


#convert to bed
chmod +x bedtools

#A
if [ $SampleN == "A" ]
then

	./bedtools bamtobed -i ${file1}.bam > ${wd}${file1}.bed;
	./bedtools bamtobed -i ${file2}.bam > ${wd}${file2}.bed;
	./bedtools bamtobed -i ${file3}.bam > ${wd}${file3}.bed;
	./bedtools bamtobed -i ${file4}.bam > ${wd}${file4}.bed;
	./bedtools bamtobed -i ${ctrl}.bam  > ${wd}${ctrl}.bed;

	cd ${wd}
	#diffReps - Detect Differential Sites from ChIP-seq with Biological Replicates.
	drDir="../diffreps-master/bin/";
	chromFile="GRCh38_15.chrom.sizes";

	perl5.8.9 ${drDir}diffReps.pl --treatment ${wd}${file3}.bed ${wd}${file4}.bed --control ${wd}${file1}.bed ${wd}${file2}.bed --report ${wd}output_results.txt --chrlen $chromFile --btr ${wd}${ctrl}.bed --bco ${wd}${ctrl}.bed --window $W --step $ST --gap $GAP --nsd $NSD --meth gt --pval 0.000001 --frag $frac --nproc 16 --noanno --nohs

	#--mode $MODE --window $W # --gap 30

	#remove temporary bed files form ${wd}
	rm ${file1}.bed
	rm ${file2}.bed
	rm ${file3}.bed
	rm ${file4}.bed
	rm ${ctrl}.bed


#B
elif [ $SampleN == "B" ]
then

	./bedtools bamtobed -i ${file1}.bam > ${wd}${file1}.bed;
	./bedtools bamtobed -i ${file2}.bam > ${wd}${file2}.bed;
	./bedtools bamtobed -i ${file3}.bam > ${wd}${file3}.bed;
	#./bedtools bamtobed -i ${file4}.bam > ${wd}${file4}.bed;
	./bedtools bamtobed -i ${ctrl}.bam  > ${wd}${ctrl}.bed;

	cd ${wd}
	#diffReps - Detect Differential Sites from ChIP-seq with Biological Replicates.
	drDir="../diffreps-master/bin/";
	chromFile="GRCh38_15.chrom.sizes";

	perl5.8.9 ${drDir}diffReps.pl --treatment ${wd}${file3}.bed --control ${wd}${file1}.bed ${wd}${file2}.bed --report ${wd}output_results.txt --chrlen $chromFile --btr ${wd}${ctrl}.bed --bco ${wd}${ctrl}.bed --window $W --step $ST --gap $GAP --nsd $NSD --meth gt --pval 0.000001 --frag $frac --nproc 16 --noanno --nohs

	#--mode $MODE --window $W # --gap 30

	#remove temporary bed files form ${wd}
	rm ${file1}.bed
	rm ${file2}.bed
	rm ${file3}.bed
	#rm ${file4}.bed
	rm ${ctrl}.bed

else

#C
#if [ $SampleN == "C" ]
#then


	./bedtools bamtobed -i ${file1}.bam > ${wd}${file1}.bed;
	#./bedtools bamtobed -i ${file2}.bam > ${wd}${file2}.bed;
	./bedtools bamtobed -i ${file3}.bam > ${wd}${file3}.bed;
	./bedtools bamtobed -i ${file4}.bam > ${wd}${file4}.bed;
	./bedtools bamtobed -i ${ctrl}.bam  > ${wd}${ctrl}.bed;

	cd ${wd}
	#diffReps - Detect Differential Sites from ChIP-seq with Biological Replicates.
	drDir="../diffreps-master/bin/";
	chromFile="GRCh38_15.chrom.sizes";

	perl5.8.9 ${drDir}diffReps.pl --treatment ${wd}${file3}.bed ${wd}${file4}.bed --control ${wd}${file1}.bed --report ${wd}output_results.txt --chrlen $chromFile --btr ${wd}${ctrl}.bed --bco ${wd}${ctrl}.bed --window $W --step $ST --gap $GAP --nsd $NSD --meth gt --pval 0.000001 --frag $frac --nproc 16 --noanno --nohs

	#--mode $MODE --window $W # --gap 30

	#remove temporary bed files form ${wd}
	rm ${file1}.bed
	#rm ${file2}.bed
	rm ${file3}.bed
	rm ${file4}.bed
	rm ${ctrl}.bed

fi


