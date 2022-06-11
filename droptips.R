#!/usr/bin/env Rscript

# Setup -------------------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = F)

# Load libraries ----------------------------------------------------------

suppressMessages(require(getopt))
suppressMessages(require(ape))

# Set up options ----------------------------------------------------------
  # col1: long flag name
  # col2: short flag name
  # col3: 0 = no argument, 1 = required, 2 = optional
  # col3: logical, integer, double, complex, character
  # col5: optional, brief description

spec <- matrix(c(
  'help'   , 'h', 0, "logical",
  'phylogeny', 'p', 1, "character",
  'tips'   , 't', 1, "character"
), byrow = T, ncol = 4)

# Read options and do help -----------------------------------------------

opt <- getopt(spec)

if ( !is.null(opt$help) ){
  cat(getopt(spec, usage = T))
  q(status = 1)
}


# Load in data -----------------------------------------------------------

phy <- read.tree(opt$phylogeny)

dropl <- readLines(opt$tips)

message("Read ", length(dropl), " entries to drop")

intips <- dropl %in% phy$tip.label
if( !all(intips) ){
  message("Ignoring ", sum( !intips), " tips to drop that aren't in phylogeny")
  dropl <- dropl[intips]
}

phy <- drop.tip(phy, dropl)

write.tree(phy, stdout())