options(stringsAsFactors = FALSE)
suppressMessages(library("argparse"))
suppressMessages(library("pander"))

setArguments <- function(){
  parser <- ArgumentParser(
    description = "Post intSiteCaller primerID filter for unique sites."
  )
  parser$add_argument("-d", "--dataDir", default = getwd(), 
                      help = "postCallerIDData directory.")
  parser$add_argument("-c", "--codeDir", type="character", nargs=1,
                      help= "Code directory.")
  parser$add_argument("-s", "--specimen", type="character", nargs=1,
                      help= "Specimen to process.")
  
  arguments <- parser$parse_args()
  arguments
}

arguments <- setArguments()
pandoc.table(data.frame(
    "Variables" = paste0(names(arguments), ":"), 
    "Values" = unname(unlist(arguments))),
  justify = c("right", "left"))

#Load required dependancies
dependancies <- c("plyr", 
                  "dplyr", 
                  "GenomicRanges", 
                  "Biostrings", 
                  "igraph",
                  "devtools") 

null <- suppressMessages(sapply(dependancies, require, 
                                character.only=TRUE, 
                                quietly=TRUE, 
                                warn.conflicts=FALSE))
dependancies_present <- sapply(dependancies, function(package){
  package <- paste0("package:", package)
  logic <- package %in% search()
})

if(FALSE %in% dependancies_present){
  df <- data.frame(
    package = as.character(dependancies), 
    loaded = dependancies_present,
    row.names = NULL)
  pandoc.table(df, style = "grid", caption = "Loaded and Unloaded packages.")
  stop("\nLoad required packages. Check above for missing dependancies.")
}else{
  remove(dependancies, dependancies_present, null)
  message("Required packages loaded.")
}

#Assign variables and source
dataDir <- arguments$dataDir
codeDir <- arguments$codeDir
specimen <- arguments$specimen

source(paste0(codeDir, "/utilities.R"))

load(paste0(dataDir, "/prefilReads_", specimen, ".RData"))

corrected_reads <- assign_sampleName_by_primerID(sites)

corrected_sites <- corrected_reads$sites_reassigned
reads_corrected <- corrected_reads$reassignment_frame

save(corrected_sites, file = paste0(dataDir, "/filReads_", specimen, ".RData"))

write.table(
  reads_corrected,
  file = paste0(dataDir, "/reassignTable_", specimen, ".tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

