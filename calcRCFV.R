#!/usr/bin/env Rscript

# Setup -------------------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = F)

# Load libraries ----------------------------------------------------------

suppressMessages(require(getopt))

# Constants ---------------------------------------------------------------

AA <- c("a", "c", "d", "e", "f", "g", "h", "i", "k", "l", "m", "n", "p", "q", "r", "s", "t", "v", "w", "y")
NT <- c("a", "c", "g", "t")

# Functions ---------------------------------------------------------------

read.fasta.simple <- function(path){
  con = file(path, 'r')
  seqstrings <- list()
  header <- NULL
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ){
      break
    } else if ( grepl("^>", line) ){
      header <- sub("^>", "", line)
    } else if ( is.null(header) ) {
      stop("First line of ", path, " does not start with >, are you sure it's",
           " a fasta file?")
    } else {
      seqstrings[header] <- paste0(seqstrings[[header]], tolower(line))
    }
  }
  close(con)
  do.call(rbind, strsplit(unlist(seqstrings), ""))
}

read.partition <- function(path){
  
  raxmlregex <- "^[A-Za-z ]+, *[^ =]+ *= *[0-9]+ *- *[0-9]+\\\\?[0-9]*"
  parsepartspec <- function(spec){
    s <- as.numeric(strsplit(spec, "[-/\\]")[[1]])
    return(seq(s[1], s[2], by = c(s[3], 1)[is.na(s[3])+1]))
  }
  
  lines <- readLines(path)
  
  partitions <- if( lines[1] == "#nexus" ){
    nexpart <- lapply(lines[grepl("^ *charset", lines)], function(x){
      y <- strsplit(x, "[ =;]+")[[1]]
      y <- y[! y %in% c("", "charset")]
      list(y[1],unlist(lapply(y[-1], parsepartspec)))
    })
  } else if( grepl(raxmlregex, lines[1]) ){
    lapply(lines, function(x){
      y <- strsplit(gsub(" ", "", x), "[,=]")[[1]][-1]
      list(y[1], unlist(sapply(y[-1], parsepartspec)))
    })
  }
  return(setNames(lapply(partitions, "[[", 2), lapply(partitions, "[[", 1)))

}

rcfv <- function(m, charset){
  # Calculate the frequency of each character in each sequence
  cfreq <- do.call(
    cbind, 
    apply(m, 1, function(x) tabulate(factor(x, levels = charset)), simplify = F)
  )
  rownames(cfreq) <- charset
  # Drop sequences not present in the dataset
  present <- colSums(cfreq) > 0
  missing <- colnames(cfreq)[!present]
  cfreq <- cfreq[, present]
  if ( ncol(cfreq) == 0 ){
    return(list(csRCFV = NULL,
                tsRCFV = NULL,
                RCFV = NULL,
                dropped = missing))
  }
  # Convert to the relative frequency of each character in each sequence = uij
  # I.e. standardised by sequence length
  relfreq <- t(t(cfreq)/colSums(cfreq))
  # Calculate the mean frequency for each character across each sequence = mean(uj)
  # I.e. sum of frequencies for that character divided by the number of sequences
  meanfreq <- rowSums(relfreq)/sum(relfreq)
  # Standardise the relative frequency of each character in each sequence by the
  # mean for that character across sequences, then divide by the number of sequences
  # I.e. abs(uij - mean(uj))/n
  finfreq <- abs(relfreq - meanfreq)/ncol(cfreq)
  return(list(csRCFV = rowSums(finfreq),
              #ncsRCFV = rowSums(finfreq)/(ncol(m)^(-.5) * 100),
              tsRCFV = colSums(finfreq)[rownames(m)],
              #ncsRCFV = colSums(finfreq)[rownames(m)]/
              #  (ncol(m)^(-.5) * ncol(finfreq)^(-1) * length(charset) * 100)
              RCFV = sum(finfreq),
              #nRCFV = sum(finfreq)/
              #  (ncol(m)^(-.5) * ncol(finfreq)^(0.01) * length(charset) * 100),
              dropped = missing))
}

rcfv.partitioned <- function(m, charset, partitions){
  m = seqmatrix
  out <- lapply(partitions, function(p) rcfv(m[,p], charset))
  dropped <- lapply(out, '[[', 4)
  for( pn in names(dropped) ){
    if( length(dropped[[pn]]) > 0 )
      warning(paste("Partition", pn, "dropped", length(dropped[[pn]]), 
                    "terminals with no data"))
  }
  
  return(list(csRCFV = do.call(rbind, lapply(out, "[[", 1)),
              tsRCFV = do.call(rbind, lapply(out, '[[', 2)),
              RCFV = unlist(lapply(out, '[[', 3))))
       
}


# Set up options ----------------------------------------------------------
# col1: long flag name
# col2: short flag name
# col3: 0 = no argument, 1 = required, 2 = optional
# col3: logical, integer, double, complex, character
# col5: optional, brief description

spec <- matrix(c(
  'help'       , 'h', 0, "logical"  , "show this helpful message",
  'matrix'     , 'm', 1, "character", "path to a supermatrix",
  'datatype'   , 'd', 1, "character", "the type of data, either AA or NT",
  'partitions' , 'p', 2, "character", "path to a nexus or raxml format partition file (optional)",
  'outprefix'  , 'o', 1, "character", "file name prefix for output files"
), byrow = T, ncol = 5)

# Read options and do help -----------------------------------------------

opt <- getopt(spec)

if ( is.null(opt) | !is.null(opt$help) ){
  message(getopt(spec, usage = T))
  q(status = 1)
}

rm(spec)

# Check inputs ------------------------------------------------------------

if ( ! opt$datatype %in% c("AA", "NT") ){
  stop("--datatype/-d must be one of AA or NT")
}

# Set character set -------------------------------------------------------

charset <- if ( opt$datatype == "AA" ) AA else NT

# Read sequence data ------------------------------------------------------

seqmatrix <- read.fasta.simple(opt$matrix)
seqmatrix[! seqmatrix %in% charset ] <- NA


# Calculate RCFV, partitioned or otherwise --------------------------------

results <- if( ! is.null(opt$partitions) ){
  rcfv.partitioned(seqmatrix, charset, read.partition(opt$partitions))
} else {
  rcfv(seqmatrix, charset)
}


# Output ------------------------------------------------------------------

for( d in names(results)[1:3] ){
  if(d == "RCFV" & is.null(opt$partitions)){
    message("RCFV: ", results$RCFV)  
  } else {
    write.csv(results[[d]], paste0(opt$outprefix, "_", d, ".csv"), 
              row.names = T, quote = F)
  }
}
