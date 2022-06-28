#!/usr/bin/env Rscript

# Setup -------------------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = F)

# Load libraries ----------------------------------------------------------

suppressMessages(require(getopt))
suppressMessages(require(ape))

# Load functions

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

unresolve <- function(phy, nodes){
  # nodes = unnodes
  unresolve_node <- function(phy, n){
    # n = node
    # If there are node names, get for later reference, otherwise use "NA"
    nodenames <- if( "node.label" %in% names(phy) ) phy$node.label else "NA"
    # Get the branch length above the node and add it to the branches below, to preserve total 
    # tip heights
    heightabove <- phy$edge.length[phy$edge[,2] == n]
    phy$edge.length[phy$edge[,1] == n] <- phy$edge.length[phy$edge[,1] == n] + heightabove
    # Extract the clade from the node, leaving behind one or more branches to tips names with 
    # the previous node label
    graft <- extract.clade(phy, n)
    phy <- drop.tip(phy, graft$tip.label, trim.internal = F)
    # Find the new number of this node in the pruned tree, and that of its parent
    newnode <- if( sum(phy$tip.label %in% nodenames) == 1 ){
      which(phy$tip.label %in% nodenames)
    } else {
      getMRCA(phy, phy$tip.label[phy$tip.label %in% nodenames])
    }
    parent <- phy$edge[,1][phy$edge[,2] == newnode]
    # Bind the extracted clade to the parent, and remove all the leftover branches from prior drop
    phy <- bind.tree(phy, graft, parent)
    phy <- drop.tip(phy, phy$tip.label[phy$tip.label %in% nodenames])
    return(phy)
  }
  
  unresolve_nodes <- function(phy, nodes){
    # Because unresolving nodes will keep changing the node numbers on the tree, define nodes 
    # instead by their children
    unresolve.tips <- lapply(nodes, function(n){
      phy$tip.label[listdescendants(phy, n, nodes = F, tips = T, inc.n = F)]})
    
    # Rename the nodes with unique ids, storing the original labels
    # Because if there are duplicate node labels, when tips are getting trimmed with 
    # drop.tip(trim.internal = F), node labels will end up as tip labels and then there'll be errors
    # in identifying the right place to graft
    nlreplace <- data.frame(inp = phy$node.label,
                            rep = paste0("n", sprintf("%03d", Ntip(phy) + (1:phy$Nnode))))
    phy$node.label <- nlreplace$rep
    for(tips in unresolve.tips){
      phy <- unresolve_node(phy, getMRCA(phy, tips))
    }
    # Return the input labels of the remaining nodes
    phy$node.label <- setNames(nlreplace$inp, nlreplace$rep)[phy$node.label]
    return(phy)
  }
  
  if( length(nodes) > 1 ){
    return(unresolve_nodes(phy, nodes))
  } else {
    return(unresolve_node(phy, nodes[1]))
  }
}


unresolve_by_support <- function(phy, threshold, 
                                 support = NULL, supporti = NULL, splitchar = "/", na.keep = T){
  # threshold = opt$threshold; supporti = opt$supporti; splitchar = opt$supportsep; support = NULL; na.keep = T
  # If reading support values from node labels, separate the values out if multiple supports, and in
  # either case check that they are coercable into numerics
  if ( is.null(support) ){
    if ( ! is.null(supporti) ){
      support <- sapply(strsplit(phy$node.label, splitchar), '[', supporti)
      supportcheck <- as.numeric(support)
    } else {
      supportcheck <- as.numeric(phy$node.label)
    }
    if ( sum(is.na(suppressWarnings(supportcheck))) == length(phy$node.label) ) {
      stop("Cannot coerce node labels to support values, are there multiple values? If so, supply an index to supporti")
    }
    support <- supportcheck 
  }
  
  if( max(support, na.rm = T) > 1  & threshold <= 1 ){
    warning("It looks like the support values are 0-100, but you've supplied a threshold <=1")
  } else if ( max(support, na.rm = T) <= 1 & threshold > 1 ){
    stop("It looks like the support values are 0-1, but you've supplied a threshold > 1 and so all nodes would be unresolved")
  } else if ( max(support, na.rm = T) < threshold ){
    stop("The threshold supplied is greater than the maximum support, all nodes would be unresolved")
  }
  
  # Get the list of nodes to unresolve (always excluding the root node, of course)
  unnodes <- support[-1] < threshold | (is.na(support[-1]) & !na.keep)
  unnodes <- (Ntip(phy)+2:phy$Nnode)[unnodes]
  
  return(unresolve(phy, unnodes))
} 

# Set up options ----------------------------------------------------------
# col1: long flag name
# col2: short flag name
# col3: 0 = no argument, 1 = required, 2 = optional
# col3: logical, integer, double, complex, character
# col5: optional, brief description

spec <- matrix(c(
  'help'      , 'h', 0, "logical"  , "show this helpful message",
  'phylogeny' , 'p', 1, "character", "path to a phylogeny to unresolve",
  'threshold' , 't', 1, "numeric"  , "support threshold below which nodes should be unresolved",
  'output'    , 'o', 1, "character", "path to write the unresolved tree",
  'supportsep', 's', 2, "character", "if the nodes contain multiple support values, how are they separated",
  'supporti'  , 'i', 2, "numeric"  , "if the nodes contain multiple support values, which one should be used?",
  'keepbl'    , 'k', 2, "logical"  , "keep branch lengths, total tip height will be maintained but patristic distances will not, not recommended"
  ), byrow = T, ncol = 5)

# Read options and do help -----------------------------------------------

opt <- getopt(spec)

if ( !is.null(opt$help) ){
  cat(getopt(spec, usage = T))
  q(status = 1)
}

# Do unresolving ------------------------------------------------------------------------------

phy <- read.tree(opt$phylogeny)

phy <- unresolve_by_support(phy, opt$threshold, supporti = opt$supporti, splitchar = opt$supportsep)

if ( is.null(opt$keepbl) ){
  phy$edge.length <- NULL
}

write.tree(phy, opt$output)

