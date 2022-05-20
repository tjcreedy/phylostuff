#!/usr/bin/env Rscript

# Setup -------------------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = F)

# Load libraries ----------------------------------------------------------

suppressMessages(require(getopt))
suppressMessages(require(ape))


# Load global variables -----------------------------------------------------------------------

taxlevels <- readLines("https://raw.githubusercontent.com/tjcreedy/constants/main/taxlevels.txt")

# Load functions ----------------------------------------------------------


# Set up options ----------------------------------------------------------
# col1: long flag name
# col2: short flag name
# col3: 0 = no argument, 1 = required, 2 = optional
# col3: logical, integer, double, complex, character
# col5: optional, brief description

spec <- matrix(c(
  'help'       , 'h', 0, "logical"  , "show this helpful message",
  'constraint' , 'c', 1, "character", "either the path to a phylogeny, or a newick tree as a quoted character string",
  'taxonomy'   , 'y', 1, "character", "path to a taxonomy table to search for tips matching the taxa in the constraint",
  'subset'     , 's', 2, "character", "path to a text file listing tips to use",
  'output'     , 'o', 1, "character", "path to write the output constraint"),
  byrow = T, ncol = 5)

# Read options and do help -----------------------------------------------

opt <- getopt(spec)

if ( !is.null(opt$help) ){
  cat(getopt(spec, usage = T))
  q(status = 1)
}

# Read template constraint and taxonomy -------------------------------------------------------

if( grepl("^\\(.*[);]$", opt$constraint) ){
  opt$constraint <- if( !grepl(";$", opt$constraint) ) paste0(opt$constraint, ";")
  template <- read.tree(text = opt$constraint)
} else {
  template <- read.tree(opt$constraint)
}

taxonomy <- read.csv(opt$taxonomy)
taxonomy <- cbind(id=taxonomy[, 1], taxonomy[, taxlevels[taxlevels %in% names(taxonomy)]])

if( !is.null(opt$subset) ){
  exclude <- readLines(opt$subset)
  taxonomy <- taxonomy[taxonomy$id %in% exclude, ]
}

# Find set of ids for each tip on the template ------------------------------------------------

# Get info
tipinfo <- setNames(lapply(template$tip.label, function(tip){
  tiplocations <- which(taxonomy == tip, arr.ind = T)
  if( length(unique(tiplocations[,2])) > 1 ){
    message("Warning: the tip ", tip, " is found in multiple taxonomic levels, taking ids with the most frequent taxonomic grade only")
    tipfreqs <- table(tiplocations[,2])
    tiplocations <- tiplocations[tiplocations[,2] == names(tipfreqs)[which.max(tipfreqs)], ]
  }
  list(level = unique(tiplocations[,2]),
       ids = taxonomy$id[tiplocations[,1]])
}), template$tip.label)

# Separate out into tips and levels

tipids <- lapply(tipinfo, "[[", "ids")
tiplevels <- sort(sapply(tipinfo, "[[", "level"))

# Build output constraint from lowest levels up ----------------------------------------------

usedids <- c()
output <- template
for(tip in names(tiplevels)){
  # Pull ids
  ids <- tipids[[tip]]
  if( length(ids) == 0){
    message("Warning: the tip ", tip, " has no matching ids in the taxonomy. This tip will be pruned from the constraint")
    output <- drop.tip(output, tip)
    next
  }
  # Remove any already used
  ids <- ids[! ids %in% usedids]
  if( length(ids) == 0 ){
    next
  }
  # Add to used
  usedids <- c(usedids, ids)
  # Make polytomy
  tiptree <- read.tree(text = paste0("(", paste(ids, collapse=","), ");"))
  # Graft to template
  output <- bind.tree(output, tiptree, which(output$tip.label == tip))
}

# Output --------------------------------------------------------------------------------------

write.tree(output, opt$output)
