#!/bin/bash

# f=$(ls *.fa)
# for J in $(echo $f) ; do  # INVARIANT.fa
# sh homer.sh $J
# done


FILE=$1

WDIR="./HOMER/"
FASTA=${WDIR}${FILE}  
MOTIFS="/../MEME/cisbp_meme/all_CISBP.homer"


n=$(wc -l <  ${FASTA})

scrambleFasta.pl $FASTA > ${FASTA}.scramble  -n $(echo $n/2 | bc -l)

# score known motifs for enrichment)
homer2 known -p 1 -i $FASTA  -b ${FASTA}.scramble -m $MOTIFS -o ${FASTA}.known


