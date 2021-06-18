#!/usr/bin/env Rscript

# Setup -------------------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = F)

# Global variables --------------------------------------------------------

taxlevels <- c("subspecies","species","superspecies",
               "subgenus","genus",
               "infratribe","subtribe","tribe","supertribe",
               "infrafamily","subfamily","family","superfamily",
               "parvorder","infraorder","suborder","order","superorder","magnorder",
               "cohort","legion",
               "parvclass","subteclass","infraclass","subclass","class","superclass",
               "microphylum","infraphylum","subphylum","phylum","superphylum",
               "infrakingdom","subkingdom","kingdom","superkingdom") # From lowest to highest

# Load libraries ----------------------------------------------------------

suppressMessages(require(getopt))
suppressMessages(require(geiger))
suppressMessages(require(taxize))
suppressMessages(require(plyr))

# Load functions ----------------------------------------------------------

read.tabular <- function(path){
  if( grepl("\\.csv$", path) ){
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

get_uid_local <- function(terms, localcache, rank){
  # terms = uniqtaxa
  # localcache = taxcache
  # rank = preslevels[level]
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
  message("Found ", length(taxids), " NCBI taxids matching tips on the tree in ", path)
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
  message("Found ", nrow(taxondata), " rows of ", ncol(taxondata), " taxonomy levels matching tips on the tree in ", path)
  return(taxondata)
}

get_taxids_from_taxonomy <- function(taxonomytable, taxcache, auth){
  # Extract the present taxonomy levels in order
  preslevels <- taxlevels[taxlevels %in% colnames(taxonomytable)]
  
  # Setup taxid container
  taxids <- NULL
  
  # Find taxon ids using taxize by working up the available taxonomic levels
  taxonomytable <- taxonomytable[, preslevels]
  level <- 1
  while( level <= length(preslevels) & nrow(taxonomytable) > 0 ){
    taxa <- taxonomytable[,level]
    uniqtaxa <- unique(taxa)
    uniqtaxa <- na.omit(uniqtaxa[uniqtaxa != ''])
    if( length(uniqtaxa) > 0 ){
      locuids <- c()
      if ( length(taxcache) > 0 ){
        locuids <- get_uid_local(uniqtaxa, taxcache, preslevels[level])
        locuids <- setNames(locuids, uniqtaxa)
        uniqtaxa <- uniqtaxa[is.na(locuids)]
      }
      remuids <- c()
      if ( length(uniqtaxa) > 0 ){
        suppressWarnings(suppressMessages(
          remuids <- get_uid(uniqtaxa, messages = F, ask = F)
        ))
        if ( is.null(auth) ) Sys.sleep(0.5)
        remuids <- setNames(remuids, uniqtaxa)
      }
      uids <- c(locuids[!is.na(locuids)], remuids[!is.na(remuids)])
      foundtaxids <- setNames(uids[taxa], rownames(taxonomytable))
      foundtaxids <- foundtaxids[!is.na(foundtaxids)]
      taxids <- append(taxids, foundtaxids)
      taxonomytable <- taxonomytable[! rownames(taxonomytable) %in% names(taxids), ]
    }
    level <- level + 1
  }
  
  if( nrow(taxonomytable) > 0 ){
    taxids <- append(taxids, setNames(rep(NA, nrow(taxonomytable)), rownames(taxonomytable)))
  }
  
  # Return
  return(taxids)
}

get_taxids_from_tiptaxonomy <- function(tips, auth){
  tipuids <- sapply(tips, function(tip){
    out <- NA
    
    # Remove any underscores
    tiprn <- gsub('_+', ' ', tip)
    
    # Split up by common separators
    tipsplit <- strsplit(tiprn, ';,')[[1]]
    
    # If only one name still, try getting the taxon_id just by searching the whole name
    if ( length(tipsplit) == 1 ){
      suppressWarnings(suppressMessages(
        out <- get_uid(tiprn, ask = F, messages = F)[1]
      ))
      if ( is.null(auth) ) Sys.sleep(0.5)
      # If no luck, try breaking it up
      if ( is.na(out) ) {
        tipsplit <- strsplit(tiprn, ' ')[[1]]
      }
    }
    
    if ( is.na(out) & length(tipsplit) > 1) {
      # Try each subset from the end
      tipsplit <- rev(tipsplit)
      for ( t in tipsplit ) {
        suppressWarnings(suppressMessages(
          out <- get_uid(t, messages = F, ask = F)[1]
        ))
        if ( is.null(auth) ) Sys.sleep(0.5)
        if ( !is.na(out) ) break
      }
    }
    
    return(out)
  })
  
  # Drop any that failed to return
  tipout <- tipuids[tipuids != '' & ! is.na(tipuids)]
  
  # Return
  message(paste("Retrieved taxon ids for", length(tipout), "tips from taxonomy in tip names"))
  return(tipout)
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

get_taxonomy_from_taxids <- function(taxids, taxcache, auth){
  taxids <- as.character(taxids)
  uuids <- unique(taxids)
  
  # Extract any from taxcache
  inlocal <- uuids[uuids %in% names(taxcache)]
  taxlocal <- list()
  if ( !is.null(taxcache) & length(inlocal) > 0  ){
    taxlocal <- taxcache[inlocal]
    message(paste("Taxonomy retrieved from local NCBI cache for", length(taxlocal), "unique NCBI taxids,"))
  }
  
  # Get from NCBI
  uuids <- uuids[! uuids %in% inlocal]
  taxncbi <- list()
  if ( length(uuids) > 0 ){
    suppressWarnings(suppressMessages(
      taxncbi <- classification(uuids, db = "ncbi")
    ))
    if ( is.null(auth) ) Sys.sleep(0.5)
    taxncbi <- taxncbi[!is.na(taxncbi)]
    message(paste("Taxonomy retrieved from remote NCBI search for", length(taxncbi), "unique NCBI taxids,"))
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
}

if ( !is.null(opt$metadata) ){
  names <- strsplit(opt$metanames, ',')[[1]]
  renamelist[['metadata']] <- parse_metadata(opt$metadata, names)
}

# Get taxonomy ------------------------------------------------------------

if ( opt$usencbi ){
  
  # Get cache if present 
  taxcache <- list()
  if ( !is.null(opt$taxcache) && file.exists(opt$taxcache) ){
    taxcache <- readRDS(opt$taxcache)
  }
  
  # Load API key if present 
  if ( !is.null(opt$auth) ){
    authlines <- readLines(opt$auth)
    authlines <- authlines[! grepl("^#", authlines)]
    options(ENTREZ_KEY = authlines[2])
  }
  
  taxids <- NULL
  tipswotaxonomy <- unknowntips
  
  # Retrieve taxids from file if present
  if ( ! is.null(opt$ncbitaxids) ){
    taxids <- get_taxids_from_file(opt$ncbitaxids, tipswotaxonomy)
    tipswotaxonomy <- tipswotaxonomy[ !tipswotaxonomy %in% names(taxids) ]
  }
  
  # Retrieve taxonomy details from file if present
  if ( ! is.null(opt$taxonomy) ){
    taxonomy <- get_taxonomy_from_file(opt$taxonomy, tipswotaxonomy)
    taxids <- append(taxids, get_taxids_from_taxonomy(taxonomy, taxcache, opt$auth))
    tipswotaxonomy <- tipswotaxonomy[ !tipswotaxonomy %in% rownames(taxonomy) ]
  }
  
  # If taxonomising nodes, search tips for taxonomy info if requested
  if ( opt$taxonomise & length(tipswotaxonomy) > 0){
    taxids <- append(taxids, get_taxids_from_tiptaxonomy(tipswotaxonomy, opt$auth))
    tipswotaxonomy <- tipswotaxonomy[ !tipswotaxonomy %in% names(taxids) ]
    
    taxids <- append(taxids, get_taxids_from_tipGBaccessions(tipswotaxonomy))
    tipswotaxonomy <- tipswotaxonomy[ !tipswotaxonomy %in% names(taxids) ]
  }
  
  # Generate final complete taxonomy table
  gtftreturn <- get_taxonomy_from_taxids(taxids, taxcache, opt$auth)
  taxidtaxonomy <- gtftreturn[[1]]
  rownames(taxidtaxonomy) <- names(taxids)
  taxcache <- gtftreturn[[2]]
  
  if ( exists("taxonomy") && 'morphospecies' %in% colnames(taxonomy) ){
    taxidtaxonomy <- cbind( taxonomy[rownames(taxidtaxonomy), 'morphospecies'], taxidtaxonomy )
  }
  if ( exists("taxonomy") && any( !rownames(taxonomy) %in% rownames(taxidtaxonomy) ) ){
    taxonomy <- rbind.fill(taxidtaxonomy, taxonomy[! rownames(taxonomy) %in% rownames(taxidtaxonomy), ])
  } else {
    taxonomy <- taxidtaxonomy
  }
  if ( ! is.null(opt$taxcache) ) {
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