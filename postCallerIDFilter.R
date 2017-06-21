options(stringsAsFactors = FALSE)
suppressMessages(library("argparse"))
suppressMessages(library("pander"))

codeDir <- dirname(
  sub("--file=", 
      "", 
      grep("--file=", commandArgs(trailingOnly=FALSE), value=T)
      )
  )

#Set up and gather commandline arguments
setArguments <- function(){
  parser <- ArgumentParser(
    description = "Post intSiteCaller primerID filter for unique sites."
    )
  parser$add_argument("-d", "--analysisDir", default = getwd(), 
                      help = "Primary analysis directory.")
  parser$add_argument("-c", "--codeDir", type = "character", nargs=1,
                      default = codeDir,
                      help = "Code directory.")
  parser$add_argument("-r", "--refGenome", type = "character", nargs = 1,
                      default = "hg38", help = "Reference genome, i.e. hg38")
  
  arguments <- parser$parse_args()
  arguments
}

arguments <- setArguments()
pandoc.table(data.frame(
    "Variables" = paste0(names(arguments), ":"), 
    "Values" = unname(unlist(arguments))),
  justify = c("right", "left"))

primeDir <- arguments$analysisDir
codeDir <- arguments$codeDir

if(!"postCallerIDData" %in% list.files(primeDir)){
  system("mkdir postCallerIDData")
}

setwd(paste0(primeDir, "/postCallerIDData"))

#Load required dependancies
addDependancies <- c("dplyr", "GenomicRanges", "Biostrings", "igraph") 

addDependsLoaded <- suppressMessages(
  sapply(add_dependencies, require, character.only = TRUE))
if(!all(addDependsLoaded)){
  pandoc.table(addDependsLoaded, style = "grid")
  stop("Check dependancies.")
}else{
  remove(addDependencies, addDependsLoaded, null)
  pandoc.strong("Required packages loaded.")
}

#Load all data needed for analysis
source(paste0(codeDir, "/utilities.R"))
sampleInfo <- read.delim(paste0(primeDir, "/sampleInfo.tsv"))

if("refGenome" %in% colnames(sampleInfo)){
  sampleInfo <- sampleInfo[sampleInfo$refGenome == arguments$ref_genome,]
}

sampleInfo$specimen <- sapply(strsplit(sampleInfo$alias, "-"), "[[", 1)

message("Loading the following specimens:")
pandoc.list(unique(sampleInfo$specimen))

allSites <- lapply(
  sampleInfo$alias, load_intSiteCaller_data, dataType = "allSites", dataDir = primeDir)

primerIDs <- lapply(
  sampleInfo$alias, load_intSiteCaller_data, dataType = "primerIDData", dataDir = primeDir)

names(allSites) <- names(primerIDs) <- sampleInfo$alias
allSites <- allSites[sapply(allSites, length) > 0]
primerIDs <- primerIDs[sapply(primerIDs, length) > 0]

if(exists("allSites")){pandoc.strong("Unique sites loaded.")}
if(exists("primerIDs")){pandoc.strong("PrimerIDs loaded.")}

#Join primerIDs to read alignments
allSites <- do.call(c, lapply(1:length(allSites), function(i){
  allSites[[i]]
  }))

primerIDs <- do.call(c, lapply(1:length(primerIDs), function(i){
  primerIDs[[i]]
  }))

allSites$specimen <- sampleInfo[
  match(allSites$sampleName, sampleInfo$alias), 
  "specimen"
  ]

allSites$primerID <- primerIDs[names(allSites)]

mcols(allSites) <- mcols(allSites)[, c("specimen", "sampleName", "primerID")]

message("PrimerIDs merged to alignment information.")

#Split up data by specimen, identify crossover primerIDs
spSites <- split(allSites, allSites$specimen)

lapply(1:length(spSites), function(i){
  sites <- spSites[[i]]
  save(sites, file = paste0("prefilReads_", names(spSites[i]), ".RData"))
  })

lapply(names(spSites), function(specimen){
  bsub(jobName=sprintf("BushmanPostCallerProcessing_%s", specimen),
       maxmem=12000,
       logFile=paste0("processLog_", specimen, ".txt"),
       command=paste0("Rscript ", codeDir, "/correct_read_assignment.R ",
                      "-d ", primeDir, "/postCallerIDData ",
                      "-c ", codeDir, " ",
                      "-s ", specimen)
  )
})
