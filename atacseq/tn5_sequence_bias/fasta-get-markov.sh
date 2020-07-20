#!/bin/bash


cd ../meme/bin

chmod +x fasta-get-markov

./fasta-get-markov -m 5 Homo_sapiens.GRCh38_15.fa GRCh38_15_background_m5.txt

#Get freq of eah 6-mer in the human genome

tail -n 4097 GRCh38_15_background_m5.txt | head

tail -n 4096 GRCh38_15_background_m5.txt > GRCh38_Kmers.txt

