suppressMessages(require(ape))

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

rootdist <- function(tree, n){
  sum(tree$edge.length[tree$edge[,2] %in% listancestors(tree, n, inc.n = T)])
}

find_subtrees_by_group <- function(tree, groups, p, start = Ntip(tree)+1, skip.p = 0, skip.n = 0){
  if(length(unique(groups)) < 2){
    stop("Argument 'groups' must have at least 2 different values")
  }
  if(length(groups) != length(tree$tip.label)){
    stop("Argument 'groups' must be the same length as the number of tips in the tree")
  }
  if(p < 0 | p > 1){
    stop("Argument 'p' must be greater than or equal to 0 and less than or equal to 1")
  }
  if(skip.p < 0 | skip.p > 1){
    stop("Argument 'skip.p' must be greater than or equal to 0 and less than or equal to 1")
  }
  if(skip.p > p){
    stop("Argument 'skip.p' cannot be greater than argument 'p'")
  }
  if(skip.n < 0){
    stop("Argument 'skip.n' must be greater than or equal to 0")
  }
  if(skip.n >= Ntip(tree)){
    stop("Argument 'skip.n' is greater than number of tips in the tree")
  }
  
  chs <- tree$edge[tree$edge[,1] == start, 2]
  
  chdata <- sapply(chs, function(ch){
    tips <- NULL
    if(ch %in% 1:Ntip(tree)){
      tips <- ch
    } else {
      tips <- listdescendants(tree, ch, nodes = F)
    }
    c(ch = ch,
      ntips = length(tips),
      prop_repd = length(unique(groups[tips]))/length(unique(groups)))
  }) %>% t()
  
  chdata

  if(all(chdata[, "prop_repd"] < p)){
    return(start)
  } else if(all(chdata[, "prop_repd"] >= p)){
    return(unlist(sapply(chs, function(ch){
      find_subtrees_by_group(tree, groups, p, ch, skip.p, skip.n)
    })))
  } else if(all(chdata[,"prop_repd"] > skip.p) & all(chdata[, "ntips"] > skip.n)){
    return(start)
  } else {
    drop.ch <- ! chdata[,"prop_repd"] >= p
    message(paste("Skipped subtree", chdata[drop.ch,"ch"]))
    return(find_subtrees_by_group(tree, groups, p, unname(chdata[!drop.ch, "ch"]), skip.p, skip.n))
  }
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
    return(c(find_monophyletic_subtrees(tree, tips, start = chs[1]), find_monophyletic_subtrees(tree, tips, start = chs[2])))
  }
}

mean_patristic_distance <- function(tree){
  m <- sapply(tree$edge[,2], function(x) length( listdescendants(tree, x , nodes = F, inc.n = T) ) )
  return( sum( m * (Ntip(tree) - m) * tree$edge.length )/( ( Ntip(tree)^2 - Ntip(tree) ) / 2 ) ) 
}

mean_sd_cophenetic <- function(tree){
  td <- cophenetic.phylo(tree)
  td <- td[lower.tri(td)]
  out <- c(mean(td), sd(td))
  gc()
  return(out)
}

sym_diff <- function(a, b) unique(c(setdiff(a, b), setdiff(b, a)))

patristic_distance <- function(tree, n1, n2){
  return(sum(tree$edge.length[tree$edge[,2] %in% sym_diff(listancestors(tree, n1, inc.n = T), listancestors(tree, n2, inc.n = T))]))
}

