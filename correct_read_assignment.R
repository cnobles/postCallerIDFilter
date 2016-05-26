options(stringsAsFactors = FALSE)
#Load required dependancies
rPackages <- c("plyr", 
               "dplyr", 
               "GenomicRanges", 
               "Biostrings", 
               "igraph", 
               "argparse") 

stopifnot(all(sapply(rPackages, 
                     require, 
                     character.only=TRUE, 
                     quietly=TRUE, 
                     warn.conflicts=FALSE)))

setArguments <- function(){
  parser <- ArgumentParser(
    description = "Post intSiteCaller primerID filter for unique sites."
  )
  parser$add_argument("-d", default = getwd(), 
                      help = "postCallerIDData directory.")
  parser$add_argument("-c", "--codeDir", type="character", nargs=1,
                      help= "Code directory.")
  parser$add_argument("-s", "--specimen", type="character", nargs=1,
                      help= "Specimen to process.")
  
  arguments <- parser$parse_args()
  arguments
}

arguments <- setArguments()
print(arguments)

dataDir <- arguments$d
codeDir <- arguments$codeDir
specimen <- arguments$s

source(paste0(codeDir, "/utilities.R"))

load(paste0(dataDir, "/", specimen, "_prefilReads.RData"))

corrected_reads <- assign_sampleName_by_primerID(sites)

corrected_sites <- corrected_reads$sites_reassigned
reads_corrected <- corrected_reads$reassignment_frame

save(
  corrected_sites, 
  file = paste0(dataDir, "/", specimen, "_postfilReads.RData")
)

write.table(
  reads_corrected,
  file = paste0(dataDir, "/", specimen, "_reassignTable.tsv"),
  quote = FALSE,
)

