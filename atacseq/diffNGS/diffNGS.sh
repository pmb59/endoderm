#!/bin/bash


##Submit job
##bsub -o o.diffNGS_1  -e e.diffNGS_1  -sp 100 -G team170 -q normal -n1 -R"select[mem>2000] rusage[mem=2000]" -M2000 sh diffNGS.sh


chmod +x /lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC2/diffNGS/diffNGS_example.R

####EPIGENODE
bsub -o o.diffNGS_1  -e e.diffNGS_1  -sp 100 -G team170 -q basement -n2 -R"select[mem>20000] rusage[mem=20000] span[hosts=1]" -M20000 " cd /lustre/scratch117/cellgen/team204/pm12/scratch109/R/R-3.2.1/bin;  ./Rscript /lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC2/diffNGS/diffNGS_example.R EPIGENODE_merged_ATAC_peaks_final.bed EG1 EG1 h12 h12 ATAC1.norm.bw ATAC2.norm.bw ATAC9.norm.bw ATAC10.norm.bw "
bsub -o o.diffNGS_2  -e e.diffNGS_2  -sp 100 -G team170 -q basement -n2 -R"select[mem>20000] rusage[mem=20000] span[hosts=1]" -M20000 " cd /lustre/scratch117/cellgen/team204/pm12/scratch109/R/R-3.2.1/bin;  ./Rscript /lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC2/diffNGS/diffNGS_example.R EPIGENODE_merged_ATAC_peaks_final.bed h12 h12 h24 h24 ATAC9.norm.bw ATAC10.norm.bw ATAC11.norm.bw ATAC12.norm.bw "
bsub -o o.diffNGS_3  -e e.diffNGS_3  -sp 100 -G team170 -q basement -n2 -R"select[mem>20000] rusage[mem=20000] span[hosts=1]" -M20000 " cd /lustre/scratch117/cellgen/team204/pm12/scratch109/R/R-3.2.1/bin;  ./Rscript /lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC2/diffNGS/diffNGS_example.R EPIGENODE_merged_ATAC_peaks_final.bed h24 h24 h36 h36 ATAC11.norm.bw ATAC12.norm.bw ATAC21.norm.bw ATAC22.norm.bw "
bsub -o o.diffNGS_4  -e e.diffNGS_4  -sp 100 -G team170 -q basement -n2 -R"select[mem>20000] rusage[mem=20000] span[hosts=1]" -M20000 " cd /lustre/scratch117/cellgen/team204/pm12/scratch109/R/R-3.2.1/bin;  ./Rscript /lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC2/diffNGS/diffNGS_example.R EPIGENODE_merged_ATAC_peaks_final.bed h36 h36 h48 h48 ATAC21.norm.bw ATAC22.norm.bw ATAC13.norm.bw ATAC14.norm.bw "
bsub -o o.diffNGS_5  -e e.diffNGS_5  -sp 100 -G team170 -q basement -n2 -R"select[mem>20000] rusage[mem=20000] span[hosts=1]" -M20000 " cd /lustre/scratch117/cellgen/team204/pm12/scratch109/R/R-3.2.1/bin;  ./Rscript /lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC2/diffNGS/diffNGS_example.R EPIGENODE_merged_ATAC_peaks_final.bed h48 h48 h72 h72 ATAC13.norm.bw ATAC14.norm.bw ATAC15.norm.bw ATAC16.norm.bw "

#Plots
chmod +x /lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC2/diffNGS/Figure2a.R
bsub -o o.Figure2a -e e.Figure2a -sp 100 -G VallierGroup -q long -n1 -R"select[mem>48000] rusage[mem=48000] span[hosts=1]" -M48000 " cd /lustre/scratch117/cellgen/team204/pm12/scratch109/R/R-3.2.1/bin;  ./Rscript /lustre/scratch117/cellgen/team204/pm12/scratch109/EPIGENODE-ATAC2/diffNGS/Figure2a.R "


#### End

#bsub -o o.diffNGS_2  -e e.diffNGS_2  -sp 100 -G team170 -q basement -n2 -R"select[mem>20000] rusage[mem=20000] span[hosts=1]" -M20000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/diffNGS/diffNGS_example.R EPIGENODE_merged_ATAC_peaks_final.bed h12 h12 h24 h24 ATAC9.norm.bw ATAC10.norm.bw ATAC11.norm.bw ATAC12.norm.bw "
#bsub -o o.diffNGS_3  -e e.diffNGS_3  -sp 100 -G team170 -q basement -n2 -R"select[mem>20000] rusage[mem=20000] span[hosts=1]" -M20000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/diffNGS/diffNGS_example.R EPIGENODE_merged_ATAC_peaks_final.bed h48 h48 h72 h72 ATAC13.norm.bw ATAC14.norm.bw ATAC15.norm.bw ATAC16.norm.bw "
#bsub -o o.diffNGS_4  -e e.diffNGS_4  -sp 100 -G team170 -q basement -n2 -R"select[mem>20000] rusage[mem=20000] span[hosts=1]" -M20000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/diffNGS/diffNGS_example.R EPIGENODE_merged_ATAC_peaks_final.bed h24 h24 h48 h48 ATAC11.norm.bw ATAC12.norm.bw ATAC13.norm.bw ATAC14.norm.bw "
#bsub -o o.diffNGS_5  -e e.diffNGS_5  -sp 100 -G team170 -q basement -n2 -R"select[mem>20000] rusage[mem=20000] span[hosts=1]" -M20000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/diffNGS/diffNGS_example.R EPIGENODE_merged_ATAC_peaks_final.bed h48 h48 h72 h72 ATAC13.norm.bw ATAC14.norm.bw ATAC15.norm.bw ATAC16.norm.bw "


####Cell cycle
#bsub -o o.diffNGS_5  -e e.diffNGS_5  -sp 100 -G team170 -q basement -n2 -R"select[mem>20000] rusage[mem=20000]" -M20000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/diffNGS/diffNGS_example.R UD_merged_ATAC_peaks_final.bed EG1 EG1 LG1 LG1 ATAC1.norm.bw ATAC2.norm.bw ATAC3.norm.bw ATAC4.norm.bw "
#bsub -o o.diffNGS_6  -e e.diffNGS_6  -sp 100 -G team170 -q basement -n2 -R"select[mem>20000] rusage[mem=20000]" -M20000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/diffNGS/diffNGS_example.R UD_merged_ATAC_peaks_final.bed LG1 LG1 G1_S G1_S ATAC3.norm.bw ATAC4.norm.bw ATAC5.norm.bw ATAC6.norm.bw "
#bsub -o o.diffNGS_7  -e e.diffNGS_7  -sp 100 -G team170 -q basement -n2 -R"select[mem>20000] rusage[mem=20000]" -M20000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/diffNGS/diffNGS_example.R UD_merged_ATAC_peaks_final.bed G1_S G1_S S_G2_M S_G2_M ATAC5.norm.bw ATAC6.norm.bw ATAC7.norm.bw ATAC8.norm.bw "
#bsub -o o.diffNGS_8  -e e.diffNGS_8  -sp 100 -G team170 -q basement -n2 -R"select[mem>20000] rusage[mem=20000]" -M20000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/diffNGS/diffNGS_example.R UD_merged_ATAC_peaks_final.bed S_G2_M S_G2_M EG1 EG1 ATAC7.norm.bw ATAC8.norm.bw ATAC1.norm.bw ATAC2.norm.bw "

#####EPIGENODE-extra comparisons
#bsub -o o.diffNGS_9   -e e.diffNGS_9   -sp 100 -G team170 -q basement -n2 -R"select[mem>20000] rusage[mem=20000]" -M20000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/diffNGS/diffNGS_example.R EPIGENODE_merged_ATAC_peaks_final.bed EG1 EG1 h24 h24 ATAC1.norm.bw ATAC2.norm.bw ATAC11.norm.bw ATAC12.norm.bw "
#bsub -o o.diffNGS_10  -e e.diffNGS_10  -sp 100 -G team170 -q basement -n2 -R"select[mem>20000] rusage[mem=20000]" -M20000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/diffNGS/diffNGS_example.R EPIGENODE_merged_ATAC_peaks_final.bed EG1 EG1 h48 h48 ATAC1.norm.bw ATAC2.norm.bw ATAC13.norm.bw ATAC14.norm.bw "
#bsub -o o.diffNGS_11  -e e.diffNGS_11  -sp 100 -G team170 -q basement -n2 -R"select[mem>20000] rusage[mem=20000]" -M20000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/diffNGS/diffNGS_example.R EPIGENODE_merged_ATAC_peaks_final.bed EG1 EG1 h72 h72 ATAC1.norm.bw ATAC2.norm.bw ATAC15.norm.bw ATAC16.norm.bw "
#bsub -o o.diffNGS_12  -e e.diffNGS_12  -sp 100 -G team170 -q basement -n2 -R"select[mem>20000] rusage[mem=20000]" -M20000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/diffNGS/diffNGS_example.R EPIGENODE_merged_ATAC_peaks_final.bed h12 h12 h48 h48 ATAC9.norm.bw ATAC10.norm.bw ATAC13.norm.bw ATAC14.norm.bw "
#bsub -o o.diffNGS_13  -e e.diffNGS_13  -sp 100 -G team170 -q basement -n2 -R"select[mem>20000] rusage[mem=20000]" -M20000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/diffNGS/diffNGS_example.R EPIGENODE_merged_ATAC_peaks_final.bed h12 h12 h72 h72 ATAC9.norm.bw ATAC10.norm.bw ATAC15.norm.bw ATAC16.norm.bw "
#bsub -o o.diffNGS_14  -e e.diffNGS_14  -sp 100 -G team170 -q basement -n2 -R"select[mem>20000] rusage[mem=20000]" -M20000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/diffNGS/diffNGS_example.R EPIGENODE_merged_ATAC_peaks_final.bed h24 h24 h72 h72 ATAC11.norm.bw ATAC12.norm.bw ATAC15.norm.bw ATAC16.norm.bw "


#####LG1 EPIGENODE vs CELL CYCLE
#bsub -o o.diffNGS_15  -e e.diffNGS_15  -sp 100 -G team170 -q basement -n2 -R"select[mem>20000] rusage[mem=20000]" -M20000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/diffNGS/diffNGS_example.R EPIGENODE_merged_ATAC_peaks_final.bed LG1 LG1 h12 h12 ATAC3.norm.bw ATAC4.norm.bw ATAC9.norm.bw ATAC10.norm.bw "

##### ATAC1-2 (0h) vs Dnase-seq unsorted 0h
#bsub -o o.diffNGS_IP  -e e.diffNGS_IP  -sp 100 -G team170 -q basement -n2 -R"select[mem>20000] rusage[mem=20000]" -M20000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch110/sanger/pm12/EPIGENODE-ATAC/diffNGS/diffNGS_example.R  xx.bed EG1 EG1 DNase DNase ATAC1.norm.bw ATAC2.norm.bw H9_DNaseI_rep1.norm.bw  H9_DNaseI_rep2.norm.bw "



##Plots
#chmod +x /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/diffNGS/Figure2a.R
#bsub -o o.Figure2a -e e.Figure2a  -sp 100 -G team170 -q long -n1 -R"select[mem>48000] rusage[mem=48000]" -M48000 " cd /lustre/scratch109/sanger/pm12/R/R-3.2.1/bin;  ./Rscript /lustre/scratch109/sanger/pm12/EPIGENODE-ATAC/diffNGS/Figure2a.R "

