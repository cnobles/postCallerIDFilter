# postCallerIDFilter
Filter output from intSiteCaller by primerIDs and resort filtered reads to their correct replicates.

## Usage:
```
Rscript path/to/postCallerIDFilter.R
Rscript path/to/postCallerIDFilter.R -d path/to/primaryAnalysisDirectory
Rscript path/to/postCallerIDFilter.R -c path/to/codeDirectory
Rscript path/to/postCallerIDFilter.R -r hg38
Rscript path/to/postCallerIDFilter.R -p r-parallel --cores 10
Rscript path/to/postCallerIDFilter.R -p bsub
```

## Arguments
**[-h, --help]** Help information regarding input format and arguments available.

**[-d, --analysisDir]** File path to the intSiteCaller primary analysis directory, default is current path.

**[-o, --outputDir]** Output directory name (to be created within the primary analysis directory, default is 'postCallerIDData'.

**[-c, --codeDir]** File path to code directory, default is taken from the system call.

**[-r, --refGenome]** Reference genome used to process samples, default is hg38.

**[-p, --process]** Parallel processing method, options include: serial **(not currently implemented)**, r-parallel, and bsub. Default is r-parallel.

**[--cores]** Specify the number of cores to use during r-parallel processing, default is 1.

## Dependencies
* argparse
* pander
* dplyr
* GenomicRanges
* Biostrings
* parallel (if using parallel processing on non-HPC machine)
