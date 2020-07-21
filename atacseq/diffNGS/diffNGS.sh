#!/bin/bash

Rscript diffNGS_example.R EPIGENODE_merged_ATAC_peaks_final.bed EG1 EG1 h12 h12 ATAC1.norm.bw ATAC2.norm.bw ATAC9.norm.bw ATAC10.norm.bw 
Rscript diffNGS_example.R EPIGENODE_merged_ATAC_peaks_final.bed h12 h12 h24 h24 ATAC9.norm.bw ATAC10.norm.bw ATAC11.norm.bw ATAC12.norm.bw 
Rscript diffNGS_example.R EPIGENODE_merged_ATAC_peaks_final.bed h24 h24 h36 h36 ATAC11.norm.bw ATAC12.norm.bw ATAC21.norm.bw ATAC22.norm.bw 
Rscript diffNGS_example.R EPIGENODE_merged_ATAC_peaks_final.bed h36 h36 h48 h48 ATAC21.norm.bw ATAC22.norm.bw ATAC13.norm.bw ATAC14.norm.bw 
Rscript diffNGS_example.R EPIGENODE_merged_ATAC_peaks_final.bed h48 h48 h72 h72 ATAC13.norm.bw ATAC14.norm.bw ATAC15.norm.bw ATAC16.norm.bw 

