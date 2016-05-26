#Utility functions
load_intSiteCaller_data <- function(sampleName, dataType, dataDir){
  content <- list.files(path = paste(dataDir, sampleName, sep = "/"))
  isThere <- any(grepl(dataType, content))
  
  if(isThere){
    file <- grep(dataType, content, value = TRUE)
    filePath <- paste(dataDir, sampleName, file, sep = "/")
    load(filePath)
  }
  
  if(dataType == "allSites" & isThere){
    allSites
  }else if(dataType == "primerIDData" & isThere){
    primerIDs
  }else if(dataType == "keys" & isThere){
    keys
  }else if(isThere){
    stop("dataType not supported by this function. Check dataType.")
  }
}

graphOverlaps <- function(sites, gap){
  overlaps <- findOverlaps(sites, maxgap = gap)
  edgelist <- matrix(c(queryHits(overlaps), subjectHits(overlaps)), ncol = 2)
  graph <- graph.edgelist(edgelist, directed = FALSE)
  graph
}

generate_posID <- function(sites=NULL, seqnames=NULL, strand=NULL, start=NULL, end=NULL, ...){
  if(length(sites) != 0){
    if(class(sites) == "GRanges"){
      chr <- as.character(seqnames(sites))
      strand <- as.vector(strand(sites))
      pos <- ifelse(strand == "+", start(sites), end(sites))
      posID <- paste0(chr, strand, pos)
    }else{
      message("Sites provided not a GRanges object, please use alternative inputs.")
      stop()
    }
  }else{
    if(length(seqnames) != 0 & length(strand) != 0 & length(start) != 0 & length(end) != 0){
      chr <- as.character(seqnames)
      strand <- as.vector(strand)
      start <- as.integer(start)
      end <- as.integer(end)
      sites.df <- data.frame(chr, strand, start, end)
      sites.df$pos <- ifelse(strand == "+", sites.df$start, sites.df$end)
      posID <- paste0(sites.df$chr, sites.df$strand, sites.df$pos)
    }else{
      message("Please supply seqnames, strand, start, and end info.")
      stop()
    }}
  return(posID)
}

bsub <- function(queue="normal", cpus=1, maxmem=NULL, wait=NULL, jobName=NULL, logFile=NULL, command=NULL){
  stopifnot(!is.null(maxmem))
  stopifnot(!is.null(command))
  
  cmd <- paste0("bsub -q ", queue, " -n ", as.character(cpus), " -M ", maxmem)
  ##cmd <- sprintf("bsub -q %s -n %s -M %s", queue, cpus, maxmem)
  
  if(!is.null(wait)){
    LSF.VERSION <- system2("bsub", "-V", stdout=TRUE, stderr=TRUE)[1]
    if( grepl("openlava", LSF.VERSION, ignore.case=TRUE) ) {
      wait <- sub("done", "ended", wait)
    }
    cmd <- paste0(cmd, " -w \"", wait, "\"")
  }
  
  if(!is.null(jobName)) cmd <- paste0(cmd, " -J \"", jobName, "\"")
  if(!is.null(logFile)) cmd <- paste0(cmd, " -o ", logFile)
  
  cmd <- paste0(cmd, " ", command)
  message(cmd)
  system(cmd)
}

assign_sampleName_by_primerID <- function(sites){
  IDs <- table(sites$primerID)
  two_or_more <- sapply(
    names(IDs[IDs >= 2]), 
    function(x) unique(sites[sites$primerID == x]$sampleName)
  )
  two_or_more <- two_or_more[sapply(two_or_more, length) >= 2]
  sites_to_correct <- sites[sites$primerID %in% names(two_or_more)]
  cluster_graph <- graphOverlaps(
    flank(sites_to_correct, -1, start = TRUE),
    gap = 5L
  )
  sites_to_correct$clusID <- paste(
    as.character(sites_to_correct$primerID),
    clusters(cluster_graph)$membership,
    sep = ":"
  )
  sites_to_correct <- split(sites_to_correct, sites_to_correct$clusID)
  sites_to_correct <- sites_to_correct[sapply(sites_to_correct, function(x){
    length(unique(x$sampleName)) >= 2
  })]
  
  corrected_reads <- lapply(
    sites_to_correct,
    function(sites_matching_primerID){
      sampleName_freq <- table(sites_matching_primerID$sampleName)
      origin <- names(sampleName_freq[sampleName_freq == max(sampleName_freq)])
      modified_reads <- sites_matching_primerID[
        !sites_matching_primerID$sampleName == origin
        ]
      modified_reads$sampleName <- origin
      modified_reads
    }
  )
  
  reassignment_frame <- data.frame(
    row.names = unlist(sapply(1:length(corrected_reads), function(i){
      names(corrected_reads[[i]])
    })),
    "reassign_sampleName" = unlist(sapply(1:length(corrected_reads), function(i){
      corrected_reads[[i]]$sampleName
    }))
  )
  reassignment_frame$prev_sampleName <- 
    sites[rownames(reassignment_frame)]$sampleName
  
  sites_reassigned <- sites
  sites_reassigned[rownames(reassignment_frame)]$sampleName <- 
    reassignment_frame$reassign_sampleName
  sites_reassigned$clusID <- NULL
  output <- list(sites_reassigned, reassignment_frame)
  names(output) <- c("sites_reassigned", "reassignment_frame")
  output
}
