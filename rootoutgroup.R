#!/usr/bin/env Rscript

# Setup -------------------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = F)

# Load libraries ----------------------------------------------------------

suppressMessages(require(getopt))
suppressMessages(require(ape))


# Load functions ------------------------------------------------------------------------------

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
  
  out <- unlist(lapply(chs, function(ch){
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

find_monophyletic_subtrees <- function(tree, tips, start = Ntip(tree)+1){
  for(tip in tips){
    if(!tip %in% tree$tip.label){
      stop(paste("Tip", tip, "not found in the supplied tree"))
    }
  }
  
  currtips <- NULL
  if(start %in% 1:Ntip(tree)){
    currtips <- start
  } else {
    currtips <- listdescendants(tree = tree, n  = start, nodes = F, tip = T, inc.n = F)
  }

  intips <- tree$tip.label[currtips] %in% tips
  
  if(all(intips)){
    return(c(start))
  } else if(! any(intips)){
    return(NULL)
  } else if(sum(intips) == 1){
    return(currtips[intips])
  } else {
    chs <- tree$edge[tree$edge[,1] == start, 2]
    return(unlist(sapply(chs, function(ch) find_monophyletic_subtrees(tree, tips, start = ch))))
  }
}


find_largest_outgroup_parent <- function(tree, tips, ignore = NULL){
  if( ! is.null(ignore) ){
    fulltree <- tree 
    tree <- drop.tip(tree, tree$tip.label[grepl(ignore, tree$tip.label)])
  }
  subtreenodes <- find_monophyletic_subtrees(tree, tips)
  if( length(subtreenodes)  == length(tips)) {
    stop("no monophyletic subtrees comprising >1 of the supplied tips found. Do you have extra tips that should be ignored?")
  }
  subtreelength <- sapply(subtreenodes, function(n){
    length(listdescendants(tree, n, nodes = F))
  })
  lop <- subtreenodes[which.max(subtreelength)]
  if( ! is.null(ignore) ){
    lop <- getMRCA(fulltree, tree$tip.label[listdescendants(tree, lop, nodes = F, tips = T)])
  }
  return(lop)
}

root_outgroup_fuzzy <- function(tree, outgroup, ignore = NULL){
  return(ladderize(root(tree, 
                        node = find_largest_outgroup_parent(tree, outgroup, ignore), 
                        resolve.root = T)))
}

# Set up options ----------------------------------------------------------
# col1: long flag name
# col2: short flag name
# col3: 0 = no argument, 1 = required, 2 = optional
# col3: logical, integer, double, complex, character
# col5: optional, brief description

spec <- matrix(c(
  'help'     , 'h', 0, "logical",
  'phy'      , 'p', 1, "character",
  'outgroup' , 'o', 1, "character",
  'ignore'   , 'g', 2, "character"
), byrow = T, ncol = 4)

# Read options and do help -----------------------------------------------

opt <- getopt(spec)

if ( !is.null(opt$help) ){
  cat(getopt(spec, usage = T))
  q(status = 1)
}


# Load in data -----------------------------------------------------------

phy <- read.tree(opt$phy)
outgroup <- readLines(opt$outgroup)

# Check outgroup ------------------------------------------------------------------------------

oginphy <- outgroup[outgroup %in% phy$tip.label]
if( length(oginphy) == 0 )
  stop("cannot find any outgroup tips in phylogeny")
if( length(oginphy) < length(outgroup))
  message("Warning: not all outgroup tips in phylogeny")

write.tree(root_outgroup_fuzzy(phy, oginphy, ignore = opt$ignore))