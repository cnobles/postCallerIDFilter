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

assign_groups <- function(sites, start = TRUE){
  origin_order <- names(sites)
  sites_flank <- flank(sites, -1, start = start)
  sites_red <- reduce(sites_flank, min.gapwidth = 5L, with.revmap = TRUE)
  revmap <- sites_red$revmap
  groups <- as.character(Rle(
    values = sapply(revmap, "[", 1), 
    lengths = sapply(revmap, length)
  ))
  sites <- sites[unlist(revmap)]
  sites$clusID <- groups
  sites <- sites[origin_order]
  as.numeric(sites$clusID)
}

group_sites <- function(sites){
  position_grp <- assign_groups(sites, start = TRUE)
  breakpoint_grp <- assign_groups(sites, start = FALSE)
  sites$clusID <- paste(sites$primerID, 
                        position_grp, 
                        breakpoint_grp,
                        sep = ":")
  sites
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
  
  if(length(sites_to_correct) > 0){
    sites_to_correct <- group_sites(sites_to_correct)
    sites_to_correct <- split(sites_to_correct, sites_to_correct$clusID)
    sites_to_correct <- sites_to_correct[sapply(sites_to_correct, function(x){
      length(unique(x$sampleName)) >= 2
    })]
  }

  if(length(sites_to_correct) > 0){
    corrected_reads <- lapply(
      sites_to_correct,
      function(sites_matching_primerID){
        sampleName_freq <- table(sites_matching_primerID$sampleName)
        origin <- names(sampleName_freq[sampleName_freq == max(sampleName_freq)])
        modified_reads <- sites_matching_primerID[
          !sites_matching_primerID$sampleName == origin[1]
          ]
        modified_reads$sampleName <- origin[1]
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
  }else{
    message("No miss assigned reads found.")
    sites_reassigned <- sites
    reassignment_frame <- NULL
  }
  output <- list(sites_reassigned, reassignment_frame)
  names(output) <- c("sites_reassigned", "reassignment_frame")
  output
}
