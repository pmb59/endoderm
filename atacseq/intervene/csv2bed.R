
setwd( "./diffNGS.Diff_Open_Chromatin_Regions.FC_2_adjP_1e-04" )

inputs <- list.files( path = ".", pattern = ".csv" )

for ( i in inputs ){
  
  readcsv <- read.csv( i, head=TRUE )

  write.table( readcsv[,1:3], file = paste0( i,'.bed' ), append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = FALSE)
}
