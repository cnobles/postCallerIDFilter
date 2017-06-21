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
  parser$add_argument("-o", "--outputDir", default = "postCallerIDData",
                      help = "Output directory, within the primary analysis directory. Default 'postCallerIDData'.")
  parser$add_argument("-c", "--codeDir", type = "character", nargs=1,
                      default = codeDir,
                      help = "Code directory.")
  parser$add_argument("-r", "--refGenome", type = "character", nargs = 1,
                      default = "hg38", help = "Reference genome, i.e. hg38")
  parser$add_argument("-p", "--process", type = "character", nargs=1,
                      default = "r-parallel",
                      help = "Parallel processing method, options include: serial, r-parallel, bsub.")
  parser$add_argument("--cores", type = "integer", nargs = 1,
                      default = 0, help = "Specify number of cores to use during parallel processing with r-parallel.")
  
  
  arguments <- parser$parse_args()
  arguments
}

args <- setArguments()
pandoc.table(data.frame(
    "Variables" = paste0(names(args), ":"), 
    "Values" = unname(unlist(args))),
  justify = c("right", "left"))

primeDir <- args$analysisDir
outputDir <- args$outputDir
codeDir <- args$codeDir

if(!outputDir %in% list.files(primeDir)){
  system(paste0("mkdir ", outputDir))
}

setwd(paste0(primeDir, "/", outputDir))

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

null <- lapply(1:length(spSites), function(i){
  sites <- spSites[[i]]
  save(sites, file = paste0("prefilReads_", names(spSites[i]), ".RData"))
  })

if(args$processing == "bsub"){
  null <- lapply(names(spSites), function(specimen){
    bsub(jobName = sprintf("BushmanPostCallerProcessing_%s", specimen),
         maxmem = 12000,
         logFile = paste0("processLog_", specimen, ".txt"),
         command = paste0("Rscript ", codeDir, "/correct_read_assignment.R ",
                          "-d ", primeDir, "/", args$outputDir, " ",
                          "-c ", codeDir, " ",
                          "-s ", specimen))
  })
}else if(args$processing == "r-parallel"){
  stopifnot(require("parallel"))
  
  buster <- makeCluster(args$cores)
  
  clusterExport(
    buster,
    varlist = list("args"))
  
  null <- parLapply(buster, names(spSites), function(specimen){
    library(stringr)
    
    cmd <- sprintf('Rscript %1$s/correct_read_assignment.R -d %2$s/%3$s -c %1$s -r %4$s -s %5$s',
                   args$codeDir, args$analysisDir, args$outputDir, args$refGenome, specimen)
    
    pander(sprintf("System call for processing: %1$s \n", specimen))
    pander(cmd)
    cmdOut <- system(cmd, intern = TRUE)
    pander(paste0(cmdOut, collapse = '\n'))
  })
  stopCluster(buster)
  
}else if(args$processing == "serial"){
  stop("Serial processing currently in development."
}
