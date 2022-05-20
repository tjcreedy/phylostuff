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
  
  currtips
  
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

count_monophyletic_subtrees_by_group <- function(tree, group){
  if( any(is.na(group)) ){
    message("Warning: grouping variable contains NAs, no information will be returned for these")
  }
  gr <- na.omit(sort(unique(group)))
  out <- cbind(data.frame(taxon = gr), 
               t(sapply(gr, function(g){
                 tips <- tree$tip.label[group == g & !is.na(group)]
                 if( length(tips) > 1){
                   subtree <- extract.clade(tree, getMRCA(tree, tips))
                   intips <- subtree$tip.label[!subtree$tip.label %in% tips]
                 } else {
                   intips <- NULL
                 }
                 c(ntips = length(tips),
                   nmono = length(find_monophyletic_subtrees(tree, tips)),
                   ninsert = length(find_monophyletic_subtrees(tree, intips)))
               })))
  row.names(out) <- NULL
  return(out)
}

consistency_index <- function(min, obs) min/obs

retention_index <- function(min, max, obs) (max - obs)/(max - min)

calculate_taxonomic_indices <- function(tree, taxonomy, exclude = NULL, drop.missing = F){
  remove <- NULL
  if( drop.missing ){
    remove <- tree$tip.label[ taxonomy == "" | is.na(taxonomy) ]
  }
  if( !is.null(exclude) ){
    remove <- unique(c(remove, exclude))
  }
  if( length(remove) > 0 ){
    taxonomy <- taxonomy[! tree$tip.label %in% remove]
    tree <- drop.tip(tree, remove)
  }
  bt <- count_monophyletic_subtrees_by_group(tree, taxonomy)
  bt$transitions <- with(bt, ifelse(nmono == 1, 1, ifelse(nmono < ninsert, nmono, ninsert + 1)))
  bt$TCI = consistency_index(1, bt$transitions)
  bt$TRI = retention_index(1, bt$ntips, bt$transitions)
  bti <- bt[bt$ntips > 1, ]
  sm <- c(n_taxa = nrow(bt), n_informative_taxa = nrow(bti), 
          CTCI = mean(bti$TCI), CTRI = mean(bti$TRI))
  return(list(summary = round(sm, 3),
              informative = bti,
              all = bt))
}


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
  'exclude'  , 'x', 2, "character", "a comma-separated list of tip labels to exclude from index calculation, e.g. an outgroup",
  'output'   , 'o', 2, "character", "if desired, a file path to which a CSV of individual taxon indices will be written for taxa with >1 representative",
  'alloutput', 'a', 2, "character", "if desired, a file path to which a CSV of individual taxon indices will be written for all taxa"
), byrow = T, ncol = 5)

# Read options and do help -----------------------------------------------

opt <- getopt(spec)

if ( !is.null(opt$help) ){
  cat(getopt(spec, usage = T))
  q(status = 1)
}

# Set defaults -----------------------------------------------------------

if( ! is.null(opt$exclude) ){
  opt$exclude <- strsplit(opt$exclude, ",")[[1]]
  missing <- opt$exclude[ ! opt$exclude %in% phy$tip.label]
  if( length(missing) > 0 ) { stop("Error: one or more members of the exclusion list is not in the phylogeny:", paste0('"', missing, '"', collapse = ', ')) }
}

# Load in data -----------------------------------------------------------

phy <- read.tree(opt$phylogeny)
taxonomy <- read.csv(opt$taxonomy)

taxonomy <- taxonomy[match(phy$tip.label, taxonomy[,1]), opt$taxlevel]
if( is.null(taxonomy) ) { stop("Error: supplied taxlevel is not present in the taxonomy table") }
if( length(taxonomy) < Ntip(phy) ) { stop("Error: taxonomy table is missing tip labels from the phylogeny, or column one is not tip labels") }

# Do calculation and output -------------------------------------------------------------------

indices <- calculate_taxonomic_indices(phy, taxonomy, exclude = opt$exclude, drop.missing = T)

write("Tree\tN taxa\tN informative taxa\tMean taxonomic CI\tMean taxonomic RI", stderr())
write(paste(c(opt$phylogeny, indices$summary), collapse = "\t"), stdout())

if( ! is.null(opt$output) ) {
  write.csv(indices$informative, opt$output, row.names = F, quote = F)
}

if( ! is.null(opt$alloutput) ) {
  write.csv(indices$all, opt$alloutput, row.names = F, quote = F)
}
