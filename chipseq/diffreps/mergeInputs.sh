#!/bin/bash

samtools merge INPUTSall.bam INPUT1.bam INPUT2.bam SP118.bam
samtools index INPUTSall.bam
