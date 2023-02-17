#!/usr/bin/env Rscript

# Setup -------------------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = F)

# Load libraries ----------------------------------------------------------

suppressMessages(require(getopt))
suppressMessages(require(geiger))
suppressMessages(require(taxize))
suppressMessages(require(plyr))

# Global variables ----------------------------------------------------------------------------

taxlevels <- readLines("https://raw.githubusercontent.com/tjcreedy/constants/main/taxlevels.txt")

# Load functions ----------------------------------------------------------

listancestors <- function(tree, n, inc.n = F){
  root <- Ntip(tree)+1
  
  if(n %in% tree$tip.label){
    n <- which(tree$tip.label == n)
  }
  if(n == root){
    return(0)
  } else {
    p <- tree$edge[tree$edge[,2] == n,1]
    
    if(p == root){
      return(p)
    } else {
      out <- c(p,listancestors(tree, p, inc.n = F))
      if(inc.n == T){
        return(c(n, out))
      } else {
        return(out)
      }
    }
  }
}

get_uid_local_and_remote <- function(terms, taxcache = NULL, auth){
  if( !is.null(taxcache) & length(taxcache) > 0 ){
    taxtable <- do.call("rbind", taxcache)
    taxtable <- unique(taxtable)
    
    termcounts <- sapply(terms, function(t) sum(taxtable$name == t))
    localterms <- taxtable$name %in% names(termcounts)[termcounts == 1]
    localoutput <- setNames(taxtable$id[localterms], taxtable$name[localterms])
  } else {
    localoutput <- NULL
  }
  remoteterms <- terms[ ! terms %in% names(localoutput) ]
  suppressWarnings(suppressMessages(
    remotesearch <- get_uid(remoteterms, ask = F, messages = F)
  ))
  if ( is.null(auth) ) Sys.sleep(0.5)
  attributes(remotesearch) <- NULL
  remoteoutput <- setNames(remotesearch, remoteterms)
  return(c(localoutput, remoteoutput)[terms])
}

get_taxonomy_from_taxids <- function(taxids, taxcache, auth){
  uuids <- unique(taxids)
  
  # Extract any from taxcache
  inlocal <- uuids[uuids %in% names(taxcache)]
  taxlocal <- list()
  if ( !is.null(taxcache) & length(inlocal) > 0 ) {
    taxlocal <- taxcache[inlocal]
    message(paste("Taxonomy retrieved from local NCBI cache for", length(taxlocal), "unique NCBI taxids,"))
  }
  
  # Get from NCBI
  uuids <- uuids[! uuids %in% inlocal]
  taxncbi <- list()
  if ( length(uuids) > 0 ){
    suppressMessages(suppressWarnings(
      taxncbi <- classification(uuids, db = "ncbi")
    ))
    if ( is.null(auth) ) Sys.sleep(0.5)
    taxncbi <- taxncbi[!is.na(taxncbi)]
    message(paste("Taxonomy retrieved from remote NCBI search for", length(taxncbi), "unique NCBI taxids,"))
  }
  
  # Concatenate
  taxall <- c(taxlocal, taxncbi)
  taxall <- setNames(taxall[taxids], names(taxids))
  
  # Convert to data frame
  taxonomy <- do.call('rbind.fill', lapply(taxall, function(cls){
    data.frame(matrix(rev(cls$name[cls$rank != 'no rank']), nrow = 1, 
                      dimnames = list(cls$id[nrow(cls)], 
                                      rev(cls$rank[cls$rank != 'no rank']))))
  }))
  
  # Order data frame
  taxonomy <- taxonomy[, taxlevels[taxlevels %in% colnames(taxonomy)]]
  
  return(list(taxonomy, c(taxcache, taxncbi)))
}

# Set up options ----------------------------------------------------------
# col1: long flag name
# col2: short flag name
# col3: 0 = no argument, 1 = required, 2 = optional
# col3: logical, integer, double, complex, character
# col5: optional, brief description

spec <- matrix(c(
  'help'     , 'h', 0, "logical"  , "show this helpful message",
  'phylo'    , 'p', 1, "character", "the taxonomised phylogeny",
  'tip'      , 'e', 1, "character", "regular expression denoting tips to infer taxonomy for ( e.g. \'^otu\')",
  'output'   , 'o', 1, "character", "path to write taxonomy table for classified tips",
  'usencbi'  , 'u', 2, "logical",   "search taxonomy information against NCBI to fill in missing levels",
  'taxcache' , 'c', 2, 'character', "path to a .RDS cache of taxonomy data to read from and/or write to",
  'auth'     , 'a', 2, "character", "an ncbi_authentication text file with your API key as the second line not beginning with #"
), byrow = T, ncol = 5)

# Read options and do help -----------------------------------------------

opt <- getopt(spec)

if ( is.null(opt) | !is.null(opt$help) ){
  message(getopt(spec, usage = T))
  q(status = 1)
}

rm(spec)

# Set defaults -----------------------------------------------------------

if( is.null(opt$usencbi) ){
  if( !is.null(opt$taxcache) | !is.null(opt$auth) ){
    stop("Specifying --taxcache and/or --auth is redundant if not specifying --usencbi")
  }
  opt$usencbi <- FALSE
} else {
  if( is.null(opt$auth) ){
    warning("Warning: trying to access NCBI without API authentication. This might not work. You might want to use --auth")
  }
}

# Load in tree and get tips ----------------------------------------------

message(paste("Reading tree", opt$phylo))

phy <- read.tree(opt$phylo)
tips <- phy$tip

noveltips <- tips[grepl(opt$tip, tips)]
rm(tips)

message(paste("Read tree with", length(phy$tip.label), "tips, found", length(noveltips), "to classify"))


# Check the node labels -----------------------------------------------------------------------

nodelabels <- c("Coleoptera", "Buprestidae")

nodelabels <- unlist(strsplit(phy$node.label, "/"))
nodelabels <- suppressWarnings(na.omit(as.numeric(nodelabels)))

if( length(nodelabels) > 1 ){
  if( all(0 <= nodelabels) & ( all(nodelabels <= 1) | all(nodelabels <= 100) )){
    stop("the node labels in ", opt$phylo, " appear to be support values, not taxonomic levels. Have you run phylabel.R with the --taxonomise/-t option?")
  } else {
    message("Warning: there appear to be numeric nodelables of ", opt$phylo, ". These might cause issues with taxonomisation")
  }
}

# Find inferred taxonomies ------------------------------------------------

noveltaxonomy <- lapply(noveltips, function(tip){
  # tip <- noveltips[4543]
  # Get ancestor nodes for this tip
  ancnods <- listancestors(phy, tip)
  ancnods <- rev(sort(ancnods))
  # Retrieve any names
  nodes <- phy$node.label[ancnods - length(phy$tip.label)]
  nodes <- nodes[nodes != '']
  gsub("_", " ", nodes)
})
names(noveltaxonomy) <- noveltips

# Output if not using NCBI ------------------------------------------------

if ( !opt$usencbi ) {
  taxonomy <- sapply(noveltaxonomy, function(tax) paste(tax, collapse = ','))
  estimated <- grepl('\"', taxonomy)
  write.csv(cbind(estimated, gsub('\"', '', taxonomy)), opt$output, quote = F)
  message(paste("Successfully assigned taxonomy to", sum(taxonomy != ""), "out of", length(noveltips), "tips."))
  q()
}

# Get cache if present ---------------------------------------------------

taxcache <- list()
if ( !is.null(opt$taxcache) && file.exists(opt$taxcache) ){
  taxcache <- readRDS(opt$taxcache)
}

# Load API key if present ------------------------------------------------

if ( !is.null(opt$auth) ){
  authlines <- readLines(opt$auth)
  authlines <- authlines[! grepl("^#", authlines)]
  options(ENTREZ_KEY = authlines[2])
}


# Parse taxonomies --------------------------------------------------------

taxids <- NULL
estimated <- NULL
i <- 1

while( length(taxids) < length(noveltaxonomy) ){
  taxa <- sapply(noveltaxonomy[ ! names(noveltaxonomy) %in% names(taxids) ], function(tax){
    if ( length(tax) >= i ){
      return(tax[i])
    } else {
      return(NA)
    }
  })
  taxa <- taxa[!is.na(taxa)]
  if ( length(taxa) == 0 ){
    break
  }
  est <- setNames(grepl('\"', taxa), names(taxa))
  taxa <- gsub('\"', '', taxa)
  uniqtaxa <- unique(taxa)
  uids <- get_uid_local_and_remote(unique(taxa), taxcache, opt$auth)
  uids <- uids[!is.na(uids)]
  foundtaxa <- taxa[taxa %in% names(uids)]
  taxids <- append(taxids, setNames(uids[foundtaxa], names(foundtaxa)))
  estimated <- append(estimated, est[names(foundtaxa)])
  i <- i + 1
}

gtftreturn <- get_taxonomy_from_taxids(taxids, taxcache, opt$auth)
taxidtaxonomy <- gtftreturn[[1]]
rownames(taxidtaxonomy) <- names(taxids)
taxcache <- gtftreturn[[2]]

output <- cbind("estimated" = estimated[rownames(taxidtaxonomy)], taxidtaxonomy)

message(paste("Successfully assigned taxonomy to", nrow(output), "out of", length(noveltips), "tips."))

missing <- names(noveltaxonomy[sapply(noveltaxonomy, length) == 0])
if ( length(missing) > 0 ){
  output <- rbind(output, matrix('', nrow = length(missing), ncol = ncol(output), dimnames = list(missing, colnames(output))))
}

write.csv(output, opt$output, row.names = T, quote = F)
