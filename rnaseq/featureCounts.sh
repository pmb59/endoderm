#!/bin/bash

GTF='Homo_sapiens.GRCh38.76.gtf'

A=$(ls *.bam )

featureCounts -p -C -T 2 -t exon -g gene_id -a $GTF -o featureCountsEndo.txt $(echo $A)


