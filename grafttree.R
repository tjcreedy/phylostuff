#!/usr/bin/env Rscript

# Setup -------------------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = F)

# Load libraries ----------------------------------------------------------

suppressMessages(require(getopt))
suppressMessages(require(ape))
suppressMessages(require(taxize))


# Load functions ----------------------------------------------------------

listdescendants <- function(tree, n, nodes = T, tips = T, inc.n = F){
  if(nodes == F & tips == F){
    stop("Nothing to return!")
  }
  if(n < 1 | n > Ntip(tree) + tree$Nnode){
    stop("Argument to n is not a valid node or tip")
  }
  # if(inc.n == T){
  #   if(nodes == F & ! n <= Ntip(tree)) warning("inc.n = T ignored because n is an internal node and nodes = F")
  # }
  
  chs <- tree$edge[tree$edge[,1] == n, 2]
  
  out <- unlist(sapply(chs, function(ch){
    if(ch <= Ntip(tree)){
      if(tips == T){
        return(ch)
      } else if(tips == "labels"){
        return(tree$tip.label[ch])
      }
    } else {
      ot <- listdescendants(tree, ch, nodes, tips)
      if(nodes == T){
        ot <- c(ch,ot)
      }
      return(ot)
    }
  }))
  
  if(inc.n == T & (nodes == T | n <= Ntip(tree))){
    out <- c(n, out)
  }
  return(out)
}

# Set up options ----------------------------------------------------------
  # col1: long flag name
  # col2: short flag name
  # col3: 0 = no argument, 1 = required, 2 = optional
  # col3: logical, integer, double, complex, character
  # col5: optional, brief description

spec <- matrix(c(
  'help'   , 'h', 0, "logical",
  'host'   , 's', 1, "character",
  'graft'  , 'g', 1, "character",
  'taxon'  , 't', 2, "character",
  'ott'    , 'o', 2, "integer"
), byrow = T, ncol = 4)

# Read options and do help ------------------------------------------------

opt <- getopt(spec)

if ( !is.null(opt$help) ){
  cat(getopt(spec, usage = T))
  q(status = 1)
}


# Set defaults ------------------------------------------------------------

if( is.null(opt$taxon) & is.null(opt$ott) ){
  exit("Error: one of -o/--ott or -t/--taxon is required")
}


# Parse taxon location ----------------------------------------------------

if( is.null(opt$ott) ){
  ott <- suppressWarnings(get_tolid(opt$taxon, ask = F, messages = F))  
  
  if( is.null(ott) ){
    quit("Taxon name given to -t/--taxon does not return a single unambiguous OTT ID from Tree Of Life. 
         Please search the taxon yourself at https://tree.opentreeoflife.org/taxonomy/browse and supply the ID to -o/--ott", 1)
  }
} else {
  ott <- opt$ott
}

if( !grepl('^ott', ott) ){
  ott <- paste0('ott', ott)
}

# Load in host and check taxon is present ---------------------------------

host <- read.tree(opt$host)

if( ! ott %in% host$node.label ){
  quit("OTT ID given or derived from taxon given not found in supplied host tree internal nodes", 1)
}

message("Host tree has ", Ntip(host), " terminals")

# Load in graft -----------------------------------------------------------

graft <- read.tree(opt$graft)

message("Graft tree has ", Ntip(graft), " terminals")

# Do grafting -------------------------------------------------------------

node <- which(host$node.label == ott) + Ntip(host)
dumptips <- listdescendants(host, node, nodes = F)

message("Graft location found at ", ott, " with ", length(dumptips), " descendants, these will be removed.")

out <- suppressWarnings(bind.tree(host, graft, where = node))
out <- drop.tip(out, dumptips)

message("Final tree has ", Ntip(out), " terminals")

write.tree(out, stdout())