
mamba create -n intervene intervene -c bioconda
conda activate intervene
wget http://ngs.sanger.ac.uk/production/endoderm/hg38/diffNGS.Diff_Open_Chromatin_Regions.FC_2_adjP_1e-04.zip
unzip diffNGS.Diff_Open_Chromatin_Regions.FC_2_adjP_1e-04.zip

Rscript csv2bed.R 

intervene upset --input ./diffNGS.Diff_Open_Chromatin_Regions.FC_2_adjP_1e-04/*.bed \
    --dpi 700 \
    --figsize 7 5 \
    --mbcolor "#17202A" --sbcolor "#5D6D7E" \
    --names='0h - 12h','12h - 24h','24h - 36h','36h - 48h','48h - 72h'

