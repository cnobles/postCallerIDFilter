options(stringsAsFactors = FALSE)
#Load required dependancies
rPackages <- c("plyr", 
               "dplyr", 
               "GenomicRanges", 
               "Biostrings", 
               "igraph", 
               "argparse",
               "devtools") 

stopifnot(all(sapply(rPackages, 
                     require, 
                     character.only=TRUE, 
                     quietly=TRUE, 
                     warn.conflicts=FALSE)))

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
  parser$add_argument("-d", default = getwd(), 
                      help = "Primary analysis directory.")
  parser$add_argument("-c", "--codeDir", type="character", nargs=1,
                      default=codeDir,
                      help= "Code directory.")
  
  arguments <- parser$parse_args()
  arguments
}

arguments <- setArguments()
print(arguments)

primeDir <- arguments$d
codeDir <- arguments$codeDir

if(!"postCallerIDData" %in% list.files(primeDir)){
  system("mkdir postCallerIDData")
}

setwd(paste0(primeDir, "/postCallerIDData"))

#Load all data needed for analysis
source(paste0(codeDir, "/utilities.R"))
sampleInfo <- read.delim(paste0(primeDir, "/sampleInfo.tsv"))
sampleInfo$specimen <- sapply(strsplit(sampleInfo$alias, "-"), "[[", 1)

message("Loading the following specimens:")
print(unique(sampleInfo$specimen))

allSites <- lapply(sampleInfo$alias, 
                   load_intSiteCaller_data, 
                   dataType = "allSites", 
                   dataDir = primeDir)
primerIDs <- lapply(sampleInfo$alias,
                    load_intSiteCaller_data, 
                    dataType = "primerIDData", 
                    dataDir = primeDir)

names(allSites) <- names(primerIDs) <- sampleInfo$alias

if(exists("allSites")){message("Unique sites loaded.")}
if(exists("primerIDs")){message("PrimerIDs loaded.")}

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
  save(sites, file = paste0(names(spSites[i]), "_prefilReads.RData"))
  })

lapply(names(spSites), function(specimen){
  bsub(jobName=sprintf("BushmanPostCallerProcessing_%s", specimen),
       maxmem=64000,
       logFile=paste0("logs/", specimen, "_output.txt"),
       command=paste0("Rscript ", codeDir, "/correct_read_assignment.R ",
                      "-d ", primeDir, "/postCallerIDData ",
                      "-c ", codeDir, " ",
                      "-s ", specimen)
  )
})