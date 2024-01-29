#!/usr/bin/env Rscript

# Setup -------------------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = F)

# Load libraries ----------------------------------------------------------

suppressMessages(require(getopt))
suppressMessages(require(ape))


# Load functions ------------------------------------------------------------------------------

suppressMessages(source(pipe(paste("wget -O -", "https://raw.githubusercontent.com/tjcreedy/phylostuff/main/phylofuncs.R"))))

# Set up options ----------------------------------------------------------
  # col1: long flag name
  # col2: short flag name
  # col3: 0 = no argument, 1 = required, 2 = optional
  # col3: logical, integer, double, complex, character
  # col5: optional, brief description

spec <- matrix(c(
  'help'     , 'h', 0, "logical",   "show this helpful message",
  'phylogeny', 'p', 1, "character", "path to a phylogeny to calculate tRI for",
  'taxonomy' , 't', 1, "character", "path to a taxonomy csv with tree tip labels in the first column and taxonomic levels as other columns",
  'taxlevel' , 'l', 1, "character", "the taxonomic level for which to calculate an index",
  'drop'     , 'd', 2, "character", "a comma-separated list of tip labels to drop from the tree before monophyly assessment - not recommended",
  'exclude'  , 'e', 2, "character", "a comma-separated list of taxa of the given level to exclude from reporting (but include in monophyly assesment), e.g. an outgroup",
  'output'   , 'o', 2, "character", "if desired, a file path to which a CSV of individual taxon indices will be written for taxa with >1 representative",
  'alloutput', 'a', 2, "character", "if desired, a file path to which a CSV of individual taxon indices will be written for all taxa"
), byrow = T, ncol = 5)

# Read options and do help -----------------------------------------------

opt <- getopt(spec)

if ( !is.null(opt$help) ){
  cat(getopt(spec, usage = T))
  q(status = 1)
}

# Load in data -----------------------------------------------------------

phy <- read.tree(opt$phylogeny)
taxonomy <- read.csv(opt$taxonomy)

# Set defaults -----------------------------------------------------------

if( ! is.null(opt$drop) ){
  opt$drop <- strsplit(opt$drop, ",")[[1]]
  missing <- opt$drop[ ! opt$drop %in% phy$tip.label]
  if( length(missing) > 0 ) { stop("Error: one or more members of the drop list is not in the phylogeny:", paste0('"', missing, '"', collapse = ', ')) }
}

if( ! is.null(opt$exclude) ){
  opt$exclude <- strsplit(opt$exclude, ",")[[1]]
  missing <- opt$exclude[ ! opt$exclude %in% phy$tip.label]
  if( length(missing) > 0 ) { stop("Error: one or more members of the exclusion list is not in the phylogeny:", paste0('"', missing, '"', collapse = ', ')) }
}

# Filter taxonomy -----------------------------------------------------------------------------

  # Check for missing tips
missingtips <- phy$tip.label[!phy$tip.label %in% taxonomy[,1]]
if( length(missingtips) > 0 ){
  if ( length(missingtips) == Ntip(phy) ){
    stop("Error: taxonomy table is missing all tip labels from the phylogeny, or column one is not tip labels") 
  } else {
    stop("taxonomy table is missing the following tips - ", paste(missingtips, collapse = ", "))
  }
}

  # Get taxonomic levels and check
taxa <- taxonomy[match(phy$tip.label, taxonomy[,1]), opt$taxlevel]
if( is.null(taxa) | is.data.frame(taxa) && ncol(taxa) == 0 ) { 
  stop("Error: supplied taxlevel is not present in the taxonomy table") 
}
if( length(taxa) < Ntip(phy) ) { stop("Error: taxonomy table is missing tip labels from the phylogeny, or column one is not tip labels") }

# Do calculation and output -------------------------------------------------------------------

indices <- calculate_taxonomic_indices(phy, taxa, exclude = opt$exclude, drop.missing = T, drop.tips = opt$drop)

write("Tree\tN taxa\tN informative taxa\tMean taxonomic CI\tMean taxonomic RI", stderr())
write(paste(c(opt$phylogeny, indices$summary), collapse = "\t"), stdout())

if( ! is.null(opt$output) ) {
  write.csv(indices$informative, opt$output, row.names = F, quote = F)
}

if( ! is.null(opt$alloutput) ) {
  write.csv(indices$all, opt$alloutput, row.names = F, quote = F)
}
