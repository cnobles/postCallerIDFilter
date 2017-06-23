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
                      help = "Primary analysis directory path.")
  parser$add_argument("-o", "--outputDir", help = "Output directory path.")
  parser$add_argument("-c", "--codeDir", type = "character", nargs=1,
                      default = codeDir,
                      help = "Code directory.")
  parser$add_argument("-r", "--refGenome", type = "character", nargs = 1,
                      default = "hg38", help = "Reference genome, i.e. hg38")
  parser$add_argument("-p", "--process", type = "character", nargs=1,
                      default = "r-parallel",
                      help = "Parallel processing method, options include: serial, r-parallel, bsub.")
  parser$add_argument("--cores", type = "integer", nargs = 1,
                      default = 1, help = "Specify number of cores to use during parallel processing with r-parallel, default is 1.")
  
  
  arguments <- parser$parse_args()
  arguments
}

args <- setArguments()
pander("Post intSiteCaller primerID filtering\n")
pander(paste0("Date: ", Sys.Date(), "\n"))

pandoc.title("\nInput variables:")
pandoc.table(data.frame(
    "Variables" = paste0(names(args), ":"), 
    "Values" = unname(unlist(args))),
  justify = c("right", "left"),
  style = "simple")

primeDir <- args$analysisDir
outputDir <- args$outputDir
codeDir <- args$codeDir
if(!args$process %in% c("bsub", "r-parallel", "serial")){
  stop("Parallel processing method must be one of: 'bsub', 'r-parallel', or 'serial'.")
}


if(!dir.exists(outputDir)){
  system(paste0("mkdir ", outputDir))
}

#Load required dependancies
addDependencies <- c("dplyr", "GenomicRanges", "Biostrings", "igraph", "stringr") 

addDependsLoaded <- suppressMessages(
  sapply(addDependencies, require, character.only = TRUE))
if(!all(addDependsLoaded)){
  pandoc.table(addDependsLoaded, style = "simple")
  stop("Check dependancies.")
}else{
  remove(addDependencies, addDependsLoaded)
  pander("Required packages loaded.")
}

#Load all data needed for analysis
source(paste0(codeDir, "/utilities.R"))
sampleInfo <- read.delim(paste0(primeDir, "/sampleInfo.tsv"))

if("refGenome" %in% colnames(sampleInfo)){
  sampleInfo <- sampleInfo[sampleInfo$refGenome == args$refGenome,]
}

sampleInfo$specimen <- sapply(strsplit(sampleInfo$alias, "-"), "[[", 1)

message("\nLoading the following specimens:")
pandoc.table(
  data.frame("Specimens" = unique(sampleInfo$specimen)), 
  style = "simple",
  justify = "left")

allSites <- lapply(
  sampleInfo$alias, 
  load_intSiteCaller_data, 
  dataType = "allSites", 
  dataDir = primeDir)

primerIDs <- lapply(
  sampleInfo$alias, 
  load_intSiteCaller_data, 
  dataType = "primerIDData", 
  dataDir = primeDir)

names(allSites) <- names(primerIDs) <- sampleInfo$alias
allSites <- allSites[
  sapply(allSites, length) > 0 & sapply(allSites, class) == "GRanges"]
primerIDs <- primerIDs[sapply(primerIDs, length) > 0]

if(exists("allSites")){pander("Unique sites loaded.\n")}
if(sum(sapply(allSites, length))  == 0) stop("No unique sites found.")
if(exists("primerIDs")){pander("PrimerIDs loaded.\n")}

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

pander("PrimerIDs merged to alignment information.\n")

#Split up data by specimen, identify crossover primerIDs
spSites <- split(allSites, allSites$specimen)

null <- lapply(1:length(spSites), function(i){
  sites <- spSites[[i]]
  save(sites, file = paste0(outputDir,"/prefilReads_", names(spSites[i]), ".RData"))
  })

if(args$process == "bsub"){
  null <- lapply(names(spSites), function(specimen){
    bsub(jobName = sprintf("BushmanPostCallerProcessing_%s", specimen),
         maxmem = 12000,
         logFile = paste0("processLog_", specimen, ".txt"),
         command = paste0("Rscript ", codeDir, "/correct_read_assignment.R ",
                          "-d ", outputDir, " ",
                          "-c ", codeDir, " ",
                          "-s ", specimen))
  })
}else if(args$process == "r-parallel"){
  stopifnot(require("parallel"))
  
  buster <- makeCluster(args$cores)
  
  clusterExport(
    buster,
    varlist = list("args"))
  
  null <- parLapply(buster, names(spSites), function(specimen){
    library(stringr)
    library(pander)    

    cmd <- sprintf('Rscript %1$s/correct_read_assignment.R -d %2$s -c %1$s -s %3$s',
                   args$codeDir, args$outputDir, specimen)
    
    pander(sprintf("System call for processing: %1$s \n", specimen))
    pander(cmd)
    cmdOut <- system(cmd, intern = TRUE)
    pander(paste0(cmdOut, collapse = '\n'))
  })
  stopCluster(buster)
  
}else if(args$process == "serial"){
  stop("Serial processing currently in development.")
}

pander("Script completed.")
q()
