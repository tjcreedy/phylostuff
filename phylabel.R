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

read.tabular <- function(path){
  if( grepl("\\.csv$", tolower(path)) ){
    return(read.csv(path, row.names = 1))
  } else {
    tab <- read.table(path)
    if(ncol(tab) == 2){
      rownames(tab) <- tab[,1]
      if( ! is.numeric(tab[1, 2]) ){
        colnames(tab) <- tab[1, ]
        tab <- tab[-1, ]
      }
      return(tab[, 2, drop = F])
    } else {
      colnames(tab) <- tab[1, ]
      return(tab[-1, ])
    }
  }
}

generate_uid_table <- function(cache){
  utab <- do.call("rbind", cache)
  utab <- unique(utab)
  row.names(utab) <- NULL
  return(utab)
}

get_uid_local.uidtable <- function(terms, localcache, rank = NULL){
  utab <- generate_uid_table(localcache)
  if( !is.null(rank) ){ # This option is more elegant but substantially slower than the sapply version in the parent function
    utab <- utab[utab$rank == rank, ]  
  }
  alluids <- setNames(utab$id, utab$name)
  return(unname(alluids[terms]))
}

get_uid_local <- function(terms, localcache, rank = NULL, method = if( is.null(rank) ) "table" else "search"){
  method <- match.arg(method)
  # terms = uniqtaxa
  # localcache = taxcache
  # rank = preslevels[level]
  if( !is.null(rank) ){
    if(method == "search"){
      taxnames <- sapply(localcache, function(td){
        if ( rank %in% td$rank ){
          return(td$name[td$rank == rank])
        } else {
          return(NA)
        }
      })
      taxnames <- taxnames[!is.na(taxnames)]
      taxnames <- setNames(names(taxnames), taxnames)
      return(unname(taxnames[terms]))
    } else if(method == "table"){
      return(get_uid_local.uidtable(terms, localcache, rank = rank))
    }
  } else {
    if(method == "table"){
      return(get_uid_local.uidtable(terms, localcache, rank = NULL))
    } else {
      stop("Error: when rank is not supplied, only table method is available")
    }
  }
}

validate_taxcache <- function(taxcache, path = NULL, exit.text = NULL){
  baseuids <- sapply(taxcache, function(tc) rev(tc$id)[1])
  check <- names(taxcache) == baseuids
  if( ! all(check) ){
    taxcache <- taxcache[check]
    if( ! is.null(path) ){
      saveRDS(taxcache, path)
    }
    if( ! is.null(exit.text) ){
      stop(exit.text)
    }
  }
  return(taxcache)
}

update_taxcache <- function(taxids, taxcache, auth, indent = ""){
  taxids <- taxids[! taxids %in% names(taxcache) ]
  if( length(taxids) == 0 ){
    return(list())
  }
  message(indent, "Updating taxonomy cache for new ", length(taxids), " taxids")
  newtaxcache <- get_taxonomy_from_taxids(taxids, taxcache, auth, indent = paste0(indent, "\t"))[[2]][taxids]
  return(newtaxcache)
}

generate_taxonlist <- function(df, taxcolumns,  onlyifonly = NULL, ntaxa = NULL){
  x <- df[, taxcolumns[taxcolumns %in% colnames(df)], drop = F]
  if(is.null(ntaxa)){
    ntaxa <- ncol(x)
  }
  x <- apply(x, 1, function(y){
    y <- y[! is.na(y) & ! y == '']
    if(!is.null(onlyifonly)){
      othernames <- names(y)[!names(y) %in% onlyifonly]
      if(length(othernames) > 0){
        y <- y[othernames]
      }
      y <- y[]
    }
    if(length(y) > ntaxa){
      y <- y[1:ntaxa]
    }
    return(paste(y, collapse = ';'))
  })
  return(x)
}

generate_tobycodes <- function(df, taxcolumns){
  get_uniq_substrings <- function(x){
    
    split_to_matrix <- function(y){
      ysplit <- strsplit(y, '')
      mxlen <- max(sapply(ysplit, length))
      ysplit <- sapply(ysplit, function(ys) c(ys, rep('', mxlen - length(ys)))) 
      return(t(ysplit))
    }  
    
    uniq_substr <- function(y, n = nrow(y), nc = 1){
      s <- apply(y[,1:nc, drop = F], 1, paste0, collapse = '')
      nodups <- !s %in% s[duplicated(s)]
      y[nodups, 1] <- s[nodups]
      y[nodups, -1] <- ''
      if(ncol(y) > 1 && any(y[, 2] != '')){
        y <- uniq_substr(y = y, n = n, nc = nc + 1)
      } else {
        y <- y[, 1]
      }
      y
    }
    
    x <- as.character(x)
    uniq <- unique(x)
    uniqmat <- split_to_matrix(uniq)
    renamer <- setNames(uniq_substr(uniqmat), uniq)
    return(renamer[x])
  }
  x <- df[, taxcolumns[taxcolumns %in% colnames(df)]]
  x[x == ''] <- NA
  if(all(c('morphospecies', 'species') %in% taxcolumns)){
    x$species <- ifelse(is.na(x$species), x$morphospecies, x$species)
    x <- subset(x, select = -morphospecies)
  }
  x <- apply(x, 2, get_uniq_substrings)
  xcolnames <- colnames(x)
  x <- generate_taxonlist(x, taxcolumns)
  return(sapply(x, function(y){
    paste0(which(xcolnames %in% names(y)), y, collapse = '')
  }))
}

parse_presence <- function(path, genes = c('ND1', 'CYTB', 'ND6', 'ND4L', 'ND4', 'ND5', 'ND3', 'COX3', 'ATP6', 'ATP8', 'COX2', 'COX1', 'ND2')){
  presence <- readLines(path)
  presence <- strsplit(presence, ',|\t') 
  return(sapply(presence, function(l){
    setNames(paste0(ifelse(genes %in% l[-1], 'X', '-'), collapse = ''), l[1])
  }))
}

parse_metadata <- function(path, names){
  metadata <- read.tabular(path)
  metadata <- metadata[, names[names %in% colnames(metadata)], drop = F]
  metadata <- setNames(apply(metadata, 1, function(x) paste(x, collapse = '-')), rownames(metadata))
  return(metadata)
}

get_taxids_from_file <- function(path, tips){
  # Load data
  taxiddata <- read.tabular(path)
  
  # Find taxids
  if( ncol(taxiddata) == 1 ){
    taxids <- setNames(taxiddata[, 1], rownames(taxiddata))
  } else {
    taxidnames <- c('taxon_id', 'taxid', 'ncbi_taxid')
    taxidname <- taxid_names[taxidnames %in% colnames(taxiddata)]
    if ( length(taxid_name) > 0 ){
      taxids <- setNames(taxiddata[, taxidname[1]], rownames(taxiddata))
    } else {
      stop(paste("Error: cannot find column named taxon_id, taxid or ncbi_taxid in", path))
    }
  }
  
  # Clean taxids
  taxids <- taxids[! is.na(taxids) & taxids != '']
  
  # Filter to matching tips
  matchtips <- intersect(tips, names(taxids))
  if( length(matchtips) > 0 ){
    taxids <- taxids[matchtips]
  } else {
    exit(paste("Error: cannot find any rows in", path, "matching to tips on the tree"))
  }
  message("\tRead ", length(taxids), " NCBI taxids matching tips on the tree in ", path)
  return(taxids)
}

get_taxonomy_from_file <- function(path, tips){
  # Load data
  taxondata <- read.tabular(path)
  colnames(taxondata) <- tolower(colnames(taxondata))
  
  # Filter to matching tips
  matchtips <- intersect(tips, rownames(taxondata))
  if( length(matchtips) > 0 ){
    taxondata <- taxondata[matchtips, ]
  } else {
    exit(paste("Error: cannot find any rows in", path, "matching to tips on the tree"))
  }
  
  # Find taxonomy information
  matchtax <- intersect(c('morphospecies', taxlevels), colnames(taxondata))
  if( length(matchtax) > 0 ){
    taxondata <- taxondata[,  c('morphospecies', taxlevels)[ c('morphospecies', taxlevels) %in% colnames(taxondata) ] ]
  } else {
    exit(paste("Error: cannot find any columns in", path, "that match to known taxonomy levels"))
  }
  message("\tRead ", nrow(taxondata), " rows of ", ncol(taxondata), " taxonomy levels ",
          if('morphospecies' %in% colnames(taxondata)) "(including morphospecies) ",
          "matching tips on the tree in ", path)
  return(taxondata)
}

get_uid_remote <- function(names, auth, indent = ""){
  message(indent, "Running remote NCBI search on ", length(names), " taxon terms...", appendLF = F)
  start <- Sys.time()
  suppressWarnings(suppressMessages(
    uids <- get_uid(names, messages = F, ask = F)
  ))
  stop <- Sys.time()
  searchtime <- as.numeric(difftime(stop, start, units = "secs"))
  message("completed in ", round(searchtime, 1), " seconds")
  if ( is.null(auth) ) Sys.sleep(0.5)
  return(uids)
}

get_taxids_from_taxonomy <- function(taxonomytable, taxcache, auth){
  require(taxize)
  # Extract the present taxonomy levels in order
  preslevels <- taxlevels[taxlevels %in% colnames(taxonomytable)]
  
  # Setup taxid container
  taxids <- NULL
  
  # Setup unmatched name recording
  unmatchednames <- lapply(preslevels, function(x) NULL)
  setNames(unmatchednames, preslevels)
  
  # Find taxon ids using taxize by working up the available taxonomic levels
  taxonomytable <- taxonomytable[, preslevels]
  message("\tSearching a taxonomy table with ", 
          nrow(taxonomytable), 
          " entries and ", 
          ncol(taxonomytable), 
          " taxonomy levels to find NCBI taxids:")
  level <- 1
  while( level <= length(preslevels) & nrow(taxonomytable) > 0 ){
    taxa <- taxonomytable[,level]
    uniqtaxa <- unique(taxa)
    uniqtaxa <- na.omit(uniqtaxa[uniqtaxa != ''])
    message("\t\tLevel ", level, " (", preslevels[level], ") has data for ", sum(taxa != ''), " entries, ", length(uniqtaxa), " unique names")
    if( length(uniqtaxa) > 0 ){
      locuids <- c()
      if ( length(taxcache) > 0 ){
        locuids <- get_uid_local(uniqtaxa, taxcache, preslevels[level])
        locuids <- setNames(locuids, uniqtaxa)
        uniqtaxa <- uniqtaxa[is.na(locuids)]
        message("\t\t\tMatched ", sum(!is.na(locuids)), " names to taxids from local taxonomy cache")
      }
      remuids <- c()
      if ( length(uniqtaxa) > 0 ){
        remuids <- get_uid_remote(uniqtaxa, auth, indent = "\t\t\t")
        remuids <- setNames(remuids, uniqtaxa)
        uniqtaxa <- uniqtaxa[is.na(remuids)]
        message("\t\t\tMatched ", sum(!is.na(remuids)), " names to taxids from remote taxonomy")
      }
      uids <- c(locuids[!is.na(locuids)], remuids[!is.na(remuids)])
      foundtaxids <- setNames(uids[taxa], rownames(taxonomytable))
      foundtaxids <- foundtaxids[!is.na(foundtaxids)]
      taxids <- append(taxids, foundtaxids)
      taxonomytable <- taxonomytable[! rownames(taxonomytable) %in% names(taxids), ]
      message("\t\t\tFound taxids for ", length(uids), " unique names, ", length(uniqtaxa), " could not be matched")
      unmatchednames[[preslevels[level]]] <- uniqtaxa
    }
    level <- level + 1
    message("\t\t\tTaxonomy table has ", nrow(taxonomytable), " unknown rows remaining")
  }
  
  if( nrow(taxonomytable) > 0 ){
    taxids <- append(taxids, setNames(rep(NA, nrow(taxonomytable)), rownames(taxonomytable)))
  }
  message("\tCompleted taxid search, successfully found taxids for ", sum(!is.na(taxids)), "/", length(taxids), " input taxonomies")
  # Return
  return(taxids)
}

get_rank_from_uid <- function(uid, knownranks, localcache, auth, indent){
  # Set up accessory function and vector of unranked ranks
  convert_ranks <- function(ranks, knownranks) unname(setNames(1:length(knownranks), knownranks)[ranks])
  shitranks <- c("no rank", "clade")
  # Get all taxonomy into table to find ranks
  cache <- list()
  tries <- 0
  while(length(cache) < length(uid) & tries < 3){
    tries <- tries + 1
    wrap <- if(tries > 1) suppressMessages else identity
    cache <- wrap(get_taxonomy_from_taxids(uid, localcache, auth, "\t\t\t\t")[[2]][uid])
    cache <- cache[!is.na(cache)]
    
    cache <- validate_taxcache(cache)
  }
  
  
  
  uid[!uid %in% unique(tab$id)]
  uidm <- "1343364"
  cache[uidm]
  classification("1343364", db = "ncbi")
  # Filter down to a single entry for each taxid and extract ranks for uids
  tab <- tab[!duplicated(tab$id), ]
  ranks <- setNames(tab$rank, tab$id)[uid]
  
  # Identify any ranks that are unranked
  shitranki <- ranks %in% shitranks
  uids.norank <- names(ranks)[shitranki]
  ranks <- setNames(convert_ranks(ranks[!shitranki], knownranks), names(ranks))
  if( length(uids.norank) > 0 ){
    ranks.norank <- sapply(uids.norank, function(nru){
      # Get ranks for a species that contains this uid
      ranktable <- localcache[[which(sapply(cache, function(ut) nru %in% ut$id))[1]]]
      # Convert the ranks to numbers
      rankns <- convert_ranks(ranktable$rank, knownranks)
      # Interpolate the missing shitranks
      ranknap <- approx(x = which(!is.na(rankns)), y = rankns[!is.na(rankns)], xout = which(is.na(rankns)))
      rankns[ranknap$x] <- ranknap$y
      # Set cellular organisms and any other higher shitranks to Infinite
      rankns[1:(which(!is.na(rankns))[1]-1)] <- Inf
      # Set any remaining (i.e low end) shitranks to low values
      rankns[is.na(rankns)] <- seq(min(rankns, na.rm = T), 0, length.out = sum(is.na(rankns)) + 1)[-1]
      return(rankns[ranktable$id == nru])
    })
    ranks <- c(ranks.norank, ranks)[uid]
  }
  return(ranks)
}

get_taxids_from_tiptaxonomy <- function(tips, localcache, auth){
  localcache <- if( length(localcache) > 0 ) localcache else NULL
  message("\tParsing ", length(tips), " tip labels for taxonomy and searching against NCBI")
  
  # Split strings by common separators
  tipsplit <- setNames(strsplit(tips, '[-;,~_ ]'), tips)
  # Retain only strings that are words
  tipsplit <- lapply(tipsplit, function(ts)  ts[grepl('^[A-Za-z][a-z]+$', ts)] )
  # Rejoin any potential species names, retaining separate words in case incorrect
  tipsplit.binom <- lapply(tipsplit, function(ts){
    potspec <- paste(ts[-length(ts)], ts[-1])
    return(potspec[grepl('^[A-Z][a-z]+ [a-z]+$', potspec)])
  })
  
  # Generate unique corpora
  taxcorpus <- unique(unlist(tipsplit))
  taxcorpus.binom <- unique(unlist(tipsplit.binom))
  message("\t\tFound potential taxonomy in ", sum(sapply(tipsplit, length) > 0), " tip labels, with ", 
          length(taxcorpus), " unique taxon words and ", length(taxcorpus.binom), " unique potential bionomials")
  
  # Set up accessory function
  list_assign_uid_by_corpus <- function(corpus, splitlist, localcache = NULL, rank = NULL){
    name <- if( !is.null(rank) && rank == "species") " potential binomials " else " unique taxon words "
    if( !is.null(localcache) ){
      locuids <- setNames(get_uid_local(corpus, localcache, rank), corpus)
      corpleft <- corpus[is.na(locuids)]
      locuids <- locuids[!is.na(locuids)]
      message("\t\t\tMatched ", length(locuids), "/", length(corpus), name, "to taxids from local taxonomy cache")
    } else {
      corpleft <- corpus
    }
    remuids <- setNames(get_uid_remote(corpleft, auth, "\t\t\t"), corpleft)
    message("\t\t\tMatched ", sum(!is.na(remuids)), "/", length(corpleft), name, "to taxids from remote taxonomy")
    uids <- c(locuids, remuids)[corpus]
    return(lapply(splitlist, 
                  function(ts) return(na.omit(uids[ts])))
    )
  }
  
  # Search binomials and filter main corpus
  if( length(taxcorpus.binom) > 0 ){
    uidlist.binom <- list_assign_uid_by_corpus(taxcorpus.binom, tipsplit.binom, localcache, "species")
    hasuid <- sapply(uidlist.binom, length) > 0
    tipsplit.sub <- tipsplit[!hasuid]
    taxcorpus.sub <- unique(unlist(tipsplit.sub))
    message("\t\t\tFound taxids for ", sum(hasuid), " tips, ", length(tipsplit.sub), " remain with ", length(taxcorpus.sub), " unique taxon words")
  } else {
    uidlist.binom <- NULL
    tipsplit.sub <- tipsplit
    taxcorpus.sub <- taxcorpus
  }
  
  # Search main corpus
  uidlist.main <- list_assign_uid_by_corpus(taxcorpus.sub, tipsplit.sub, localcache)
  uidlist <- c(uidlist.binom, uidlist.main)
  message("\t\t\tFound taxids for ", sum(sapply(uidlist, length) > 0), " further tips")
  
  # Filter hits by rank
  message("\t\t\tFiltering taxids by rank for each tip")
  uidlist <- lapply(uidlist, unname)
  uidlist[sapply(uidlist, function(x) "1343364" %in% x)]
  
  
  
  uidranks <- get_rank_from_uid(unique(unlist(uidlist)), taxlevels, taxcache, auth, "\t\t\t\t")
  uidlist.multi <- sapply(uidlist, length) > 1
  
  # Finalise and return
  out <- setNames(c(unlist(uidlist[!uidlist.multi]), 
                    sapply(uidlist[uidlist.multi], function(pu) {
                      ranki <- uidranks[pu]
                      return(pu[ranki == min(ranki)][1])
                    }))[tips], tips)
  message("\tCompleted search, successfully found taxids for ", sum(!is.na(out)), "/", length(out), " tips labels from taxonomy")
  return(out[!is.na(out)])
}

get_taxids_from_tipGBaccessions <- function(tips, auth){
  # Extract genbank accessions from the tips
  gbaccessions <- sapply(tips, function(tip){
    regmatches(tip, regexpr("(?:[^A-Za-z0-9]|^)[A-Z]{1,2}_{0,1}[0-9]{5,8}(?:[^A-Za-z0-9]|$)", tip))[1]
  })
  gbout <- NULL
  if ( sum(!is.na(gbaccessions)) > 0 ) {
    # Create dataframe and drop any rows without putative genbank accessions
    gbdf <- data.frame(gbacc = gbaccessions, id = unknowntips)[gbaccessions != '' & !is.na(gbaccessions),]
    
    # Retrieve the taxon ids from genbank and drop rows without taxon ids
    gbdf$taxon_id <- suppressMessages(suppressWarnings(
      unlist(genbank2uid(gbdf$gbacc))
    ))
    if ( is.null(auth) ) Sys.sleep(0.5)
    gbout <- gbdf[! is.na(gbdf$taxon_id), 'taxon_id']
    names(gbout) <- row.names(gbdf)
  }
  
  # Return
  message(paste("Retrieved taxon ids for", length(gbout), "tips from GB accession numbers in tip names"))
  return(gbout)
}

get_taxonomy_from_taxids <- function(taxids, taxcache, auth, indent = ""){
  taxids <- as.character(taxids)
  uuids <- unique(taxids)
  
  # Extract any from taxcache
  taxlocal <- list()
  inlocal <- uuids[uuids %in% names(taxcache)]
  if ( length(inlocal) > 0){
    taxlocal <- taxcache[inlocal]
    message(indent, "Found taxonomy from local NCBI cache for ", length(taxlocal), " unique taxids")
  }
  
  # Get from NCBI
  uuids <- uuids[! uuids %in% names(taxlocal)]
  taxncbi <- list()
  if ( length(uuids) > 0 ){
    message(indent, "Running remote NCBI search on ", length(uuids), " taxids...", appendLF = F)
    start <- Sys.time()
    suppressWarnings(suppressMessages(
      taxncbi <- classification(uuids, db = "ncbi")
    ))
    stop <- Sys.time()
    taxncbi <- taxncbi[!is.na(taxncbi)]
    searchtime <- as.numeric(difftime(stop, start, units = "secs"))
    message("completed in ", round(searchtime, 1), " seconds, found taxonomy for ", length(taxncbi), " taxids")
    if ( is.null(auth) ) Sys.sleep(0.5)
  }
  
  # Concatenate
  taxall <- c(taxlocal, taxncbi)
  taxall <- setNames(taxall[as.character(taxids)], names(taxids))
  
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
  'help'        , 'h', 0, "logical"  , "show this helpful message",
  'phylo'       , 'p', 1, "character", "the phylogeny to label",
  'output'      , 'o', 1, "character", "path to write labelled phylogeny", 
  'taxonomise'  , 'x', 0, "logical"  , "assign taxonomic levels to nodes on the tree",
  'excludetip'  , 'e', 2, "character", "regular expression denoting tips that should not be used for taxonomisation ( e.g. \'^otu\')",
  'rename'      , 'r', 0, "logical"  , "rename tips using taxonomy, presence and/or metadata information given",
  'nameorder'   , 'd', 2, "character", "comma-separated list giving the order of \'name\', \'genepresence\', \'taxonomy\', \'species\' and/or \'metadata\' in new tip names (default uses any supplied in that order)",
  'taxonomy'    , 'y', 2, "character", "path to a table containing taxonomies to add to tip labels and/or taxonomise the tree with",
  'taxlevels'   , 'v', 2, "character", "comma-separated list of taxonomic levels to use in renaming tips (default: order,family)",
  'usencbi'     , 'u', 2, "logical"  , "search taxonomy information against NCBI to fill in missing levels",
  'ncbitaxids'  , 'i', 2, "character", "path to a table containing ncbi taxids to retrieve taxonomy for adding to tip labels and/or taxonomisation",
  'taxcache'    , 'c', 2, 'character', "path to a .RDS cache of taxonomy data to read from and/or write to",
  'auth'        , 'a', 2, "character", "an ncbi_authentication text file with your API key as the second line not beginning with #",
  'metadata'    , 'm', 2, "character", "path to a table containing metadata to add to tip labels",
  'metanames'   , 'l', 2, "character", "comma separated list of columns to use from metadata table",
  'tobycodes'   , 'b', 0, "logical"  , "use abbreviated taxonomy names (\"tobycodes\") instead of full taxonomy",
  'genepresence', 'g', 2, "character", "path to a gene presence file to parse and add to tip labels",
  'threads'     , 't', 2, 'integer'  , "number of threads to run on (only used for taxonomisation)",
  'strict'      , 's', 0, 'logical'  , "strictly taxonomise the tree"
), byrow = T, ncol = 5)

# Read options and do help -----------------------------------------------


opt <- getopt(spec)

#setwd("/home/thomas/work/iBioGen_postdoc/MMGdatabase/phylogeny/reftree498_project/")
#opt$phylo <- "6_trees/6_nt_ultboot_nopart.contree"
#opt$taxonomy <- "6_renamedata.csv"
#opt$metadata <- opt$taxonomy
#opt$output <- "testout.tre"
#opt$metanames <- "locality,site,morphospecies"
#opt$rename <- T
#opt$taxlevels <- "order,family,subfamily,genus"
#opt$taxonomise <- T
#opt$taxcache <- "NCBI_taxonomy.RDS"
#opt$auth <- "~/passwords/api_tokens/ncbi_authentication.txt"
#opt$threads <- 6
#opt$usencbi <- T

if ( is.null(opt) | !is.null(opt$help) ){
  message(getopt(spec, usage = T))
  q(status = 1)
}

rm(spec)

# Set defaults and error check -----------------------------------------

if ( is.null(opt$phylo)     )  { stop("Phylogeny is required")                }
if ( is.null(opt$threads)   )  { opt$threads    <- 1                          }
if ( is.null(opt$strict)    )  { opt$strict  <- FALSE                      }
if ( is.null(opt$taxonomise) ) { opt$taxonomise <- FALSE       }
if ( is.null(opt$rename) ) { opt$rename <- FALSE   }
if ( is.null(opt$tobycodes) ) { opt$tobycodes <- FALSE }
if ( is.null(opt$usencbi) ) { opt$usencbi <- FALSE }
if (  !opt$usencbi & (!is.null(opt$ncbitaxids) | !is.null(opt$taxcache) | !is.null(opt$auth)) ){
  stop("Specified --ncbitaxids, --taxcache and/or --auth, but not specified --usencbi")
}
if ( !opt$rename & !opt$taxonomise ){
  stop("Specify to --rename and/or --taxonomise")
}
if ( opt$usencbi & is.null(opt$auth) ){
  warning("Warning: trying to access NCBI without API authentication. This might not work. You might want to use --auth")
}
if ( !opt$rename & 
     ( !is.null(opt$metadata) | !is.null(opt$genepresence) | opt$tobycodes | !is.null(opt$metanames) | !is.null(opt$nameorder) | !is.null(opt$taxlevels)) ){
  stop("When only taxonomising, --metadata, --genepresence, --tobycodes, --metanames, --nameorder, --taxlevels are redundant")
}
if ( ! is.null(opt$taxlevels) ) {
  opt$taxlevels <- strsplit(opt$taxlevels, ",")[[1]]
  if ( ! all(opt$taxlevels %in% c('morphospecies', taxlevels)) ){
    stop("Error: unrecognised taxonomic level passed to --taxlevels")
  } else {
    opt$taxlevels <- c('morphospecies', taxlevels)[match(opt$taxlevels, c('morphospecies', taxlevels))]
  }
} else {
  opt$taxlevels <- c('order', 'family')
} 
if ( !is.null(opt$nameorder) ){
  opt$nameorder <- strsplit(opt$nameorder, ',', fixed = T)
  if ( any( c('species', 'taxonomy') %in% opt$nameorder) & is.null(opt$taxonomy) & is.null(opt$ncbitaxids) ) {
    stop("To rename with taxonomy and/or species, supply table(s) to --taxonomy and/or --taxids")
  }
  if ( 'species' %in% opt$nameorder & 'species' %in% opt$taxlevels ){
    stop("\'species\' is in both --nameorder and --taxlevels, pick one otherwise you'll have redundant names")
  }
  if ( 'metadata' %in% opt$nameorder & is.null(opt$metadata) ){
    stop("To rename with metadata, supply table to --metadata")
  }
  if ( 'genepresence' %in% opt$nameorder & is.null(opt$genepresence) ){
    stop("To rename with genepresence, supply table to --genepresence")
  }
} else {
  if ( 'species' %in% opt$taxlevels ){
    stop("--nameorder is not specified and \'species\' is in --taxlevels. If not specified, nameorder includes \'species\' by default, so this would produce redundant names - either remove \'species\' from --taxlevels or explicitly specify --nameorder")
  }
  opt$nameorder <- 'name'
  if ( !is.null(opt$genepresence) ) { opt$nameorder <- append(opt$nameorder, 'genepresence') }
  if ( !is.null(opt$taxonomy) | !is.null(opt$ncbitaxids) ) { opt$nameorder <- append(opt$nameorder, c('taxonomy', 'species')) }
  if ( !is.null(opt$metadata) ) { opt$nameorder <- append(opt$nameorder, 'metadata') }
}
if ( opt$taxonomise & (!opt$usencbi & is.null(opt$taxonomy))){
  stop("To taxonomise without NCBI, supply taxonomy information to --taxonomy")
}

# Load in tree and get tips ----------------------------------------------

message(paste("Reading tree", opt$phylo))

phy <- read.tree(opt$phylo)
unknowntips <- phy$tip.label
if ( ! is.null(opt$excludetip) ){
  unknowntips <- unknowntips[ ! grepl(opt$excludetip, unknowntips) ]
}

message(paste("Read tree with", length(phy$tip.label), "tips"))


# Get any gene presence and/or metadata -----------------------------------

renamelist <- list()

if ( !is.null(opt$genepresence) ){
  renamelist[['genepresence']] <- parse_presence(opt$genepresence)
  message("Successfully parsed gene presence data")
}

if ( !is.null(opt$metadata) ){
  names <- strsplit(opt$metanames, ',')[[1]]
  renamelist[['metadata']] <- parse_metadata(opt$metadata, names)
  message("Successfully parsed metadata")
}

# Get taxonomy ------------------------------------------------------------

message("Starting taxonomy operations")

if ( opt$usencbi ){
  
  # Get cache if present 
  taxcache <- list()
  if ( !is.null(opt$taxcache) && file.exists(opt$taxcache) ){
    taxcache <- readRDS(opt$taxcache)
    taxcache <- validate_taxcache(taxcache, opt$taxcache)
    message("Read supplied taxonomy cache")
  }
  
  # Load API key if present 
  if ( !is.null(opt$auth) ){
    authlines <- readLines(opt$auth)
    authlines <- authlines[! grepl("^#", authlines)]
    options(ENTREZ_KEY = authlines[2])
  }
  
  taxids <- NULL
  tipswotaxonomy <- unknowntips
  
  message("Starting search for taxonomy for ", length(tipswotaxonomy), " tips")
  
  # Retrieve taxids from file if present
  if ( ! is.null(opt$ncbitaxids) ){
    taxids <- get_taxids_from_file(opt$ncbitaxids, tipswotaxonomy)
    message("Retrieved taxonomy ids for ", sum(names(taxids) %in% tipswotaxonomy), " tips from taxids file")
    tipswotaxonomy <- tipswotaxonomy[ !tipswotaxonomy %in% names(taxids) ]
    newtaxcache <- update_taxcache(taxids, taxcache, opt$auth, "\t")
    taxcache <- validate_taxcache(c(taxcache, newtaxcache), exit.text = 
                                  "Error: discrepancy found in local taxonomy cache. This is likely due to your taxids file having an outdate taxid. Please check your taxids and start again")
  }
  
  # Retrieve taxonomy details from file if present
  if ( ! is.null(opt$taxonomy) ){
    taxonomy <- get_taxonomy_from_file(opt$taxonomy, tipswotaxonomy)
    newtaxids <- na.omit(get_taxids_from_taxonomy(taxonomy, taxcache, opt$auth))
    taxids <- append(taxids, newtaxids)
    message("Retrieved taxonomy ids for ", sum(names(newtaxids) %in% tipswotaxonomy), " tips from taxonomy table")
    tipswotaxonomy <- tipswotaxonomy[ !tipswotaxonomy %in% names(taxids) ]
    newtaxcache <- update_taxcache(unique(taxids), taxcache, opt$auth, "\t")
    taxcache <- validate_taxcache(c(taxcache, newtaxcache), opt$taxcache, 
                                  "Error: discrepancy found in local taxonomy cache. This is likely due to a taxonomy id merge in the NCBI data. The offending entries have been removed, please start again")
  }
  
  # If taxonomising nodes, search tips for taxonomy info if requested
  if ( opt$taxonomise & length(tipswotaxonomy) > 0 ){
    taxids <- append(taxids, get_taxids_from_tiptaxonomy(tipswotaxonomy, taxcache, opt$auth))
    tipswotaxonomy <- tipswotaxonomy[ !tipswotaxonomy %in% names(taxids) ]
    
    taxids <- append(taxids, get_taxids_from_tipGBaccessions(tipswotaxonomy))
    tipswotaxonomy <- tipswotaxonomy[ !tipswotaxonomy %in% names(taxids) ]
  }
  
  # Generate final complete taxonomy table
  gtftreturn <- get_taxonomy_from_taxids(taxids, taxcache, opt$auth)
  taxidtaxonomy <- gtftreturn[[1]]
  rownames(taxidtaxonomy) <- names(taxids)
  
  if ( exists("taxonomy") && 'morphospecies' %in% colnames(taxonomy) ){
    taxidtaxonomy <- cbind( taxonomy[rownames(taxidtaxonomy), 'morphospecies'], taxidtaxonomy )
  }
  if ( exists("taxonomy") && any( !rownames(taxonomy) %in% rownames(taxidtaxonomy) ) ){
    taxonomy <- rbind.fill(taxidtaxonomy, taxonomy[! rownames(taxonomy) %in% rownames(taxidtaxonomy), ])
  } else {
    taxonomy <- taxidtaxonomy
  }
  if ( ! is.null(opt$taxcache) ) {
    newcache <- gtftreturn[[2]]
    newcache <- newcache[!names(newcache) %in% names(taxcache)]
    taxcache <- c(taxcache, newcache)
    saveRDS(taxcache, opt$taxcache)
  }

} else {
  
  # Retrieve taxonomy details from file if present
  if ( ! is.null(opt$taxonomy) ){
    taxonomy <- get_taxonomy_from_file(opt$taxonomy, unknowntips)
  }
  
}
 
# Parse taxonomy for renaming, if doing
if ( opt$rename ){
  if ( 'taxonomy' %in% opt$nameorder ) {
    if ( opt$tobycodes ) {
      renamelist[['taxonomy']] <- setNames(generate_tobycodes(taxonomy, opt$taxlevels), rownames(taxonomy))
    } else {
      renamelist[['taxonomy']] <- setNames(generate_taxonlist(taxonomy, opt$taxlevels), rownames(taxonomy))
    }
  }
  if ( 'species' %in% opt$nameorder ) {
    if ( 'species' %in% colnames(taxonomy) ) {
      species <- taxonomy$species
    } else {
      species <- rep(NA, nrow(taxonomy))
    }
    if ( 'morphospecies' %in% colnames(taxonomy) ){
      morphospecies <- taxonomy$morphospecies
    } else {
      morphospecies <- rep(NA, nrow(taxonomy))
    }
    renamelist[['species']] <- setNames(ifelse(!is.na(species) & !species != '', morphospecies, species), rownames(taxonomy))
  }
}

# Do taxonomisation, if doing ---------------------------------------------
if ( opt$taxonomise ){
  message(paste("Applying taxonomy. The next message will probably be warnings from geiger::nodelabel.phylo - don't be alarmed:",
                "> The message 'redundant labels encountered' means that these levels are fully nested, no big deal.",
                "> The message 'redundant labels encountered at root' means that the entire tree falls inside these taxa, also no big deal unless the first value is a very low taxonomic level.",
                "> The message 'labels missing from phy' means that it couldn't place these taxonomic levels. We don't expect it to be able to place everything, but if this seems like a lot you could try running again without -strict (if you aren't already).",
                sep = '\n'))
  phy <- nodelabel.phylo(phy, taxonomy, strict = opt$strict, ncores = opt$threads)
}

# Do renaming, if doing ---------------------------------------------------

if ( opt$rename ){
  renamelist <- lapply(renamelist, function(l){
    l <- l[phy$tip.label]
    return(l)
  }) 
  renamelist <- do.call('cbind', renamelist)
  renamelist <- cbind(phy$tip.label, renamelist)
  phy$tip.label <- apply(renamelist, 1, function(x) paste(x[!is.na(x) & x != ''], collapse = '~'))
}

# Output ------------------------------------------------------------------
write.tree(phy, opt$output)
message(paste("Completed successfully, tree written to",  opt$output))