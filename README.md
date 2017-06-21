# postCallerIDFilter
Filter output from intSiteCaller by primerIDs and resort filtered reads to their correct replicates.

## Usage:
```
Rscript path/to/postCallerIDFilter.R
Rscript path/to/postCallerIDFilter.R -d path/to/primaryAnalysisDirectory
Rscript path/to/postCallerIDFilter.R -c path/to/codeDirectory
Rscript path/to/postCallerIDFilter.R -r hg38
```

## Arguments
**[-h, --help]** Help information regarding input format and arguments available.

**[-d, --analysisDir]** File path to the intSiteCaller primary analysis directory, default is current path.

**[-c, --codeDir]** File path to code directory, default is taken from the system call.

**[-r, --refGenome]** Reference genome used to process samples, default is hg38.

## Dependencies
* argparse
* pander
* dplyr
* GenomicRanges
* Biostrings
