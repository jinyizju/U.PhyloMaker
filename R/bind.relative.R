bind.relative <- function (sp.list, tree, gen.list, nodes.type = 1, output.sp.list = TRUE)
{
  if (!all(sub("_.*", "", tree$tip.label) %in% gen.list$genus)) {
    tree <- drop.tip(tree, which(sub("_.*", "", tree$tip.label) %in% setdiff(sub("_.*", "", tree$tip.label), gen.list$genus)))
  }
  treeX <- tree
  d1 <- data.frame(species = tree$tip.label, genus = sub("_.*", "", tree$tip.label))
  d1 <- merge(d1, gen.list, all.x = T)[, union(colnames(d1), colnames(gen.list))]
  if (nodes.type == 1)   nodes <- build.nodes.1(tree, d1)
  if (nodes.type == 2)   nodes <- build.nodes.2(tree, d1)
  nodesN <- nodes
  if (is.null(tree$node.label))   tree$node.label <- rep("", tree$Nnode)
  dimnames(sp.list)[[2]][1:3] <- c("species", "genus", "family")
  sp.list[sapply(sp.list, is.factor)] <- lapply(sp.list[sapply(sp.list, is.factor)], as.character)
  sp.list <- sp.list[sp.list$species != "", ]
  sp.list$species <- gsub(" ", "_", sp.list$species)
  sp.list$species <- gsub("(^[[:alpha:]])", "\\U\\1", sp.list$species, perl = TRUE)
  sp.list$genus <- gsub("(^[[:alpha:]])", "\\U\\1", sp.list$genus, perl = TRUE)
  sp.list$family <- gsub("(^[[:alpha:]])", "\\U\\1", sp.list$family, perl = TRUE)
  sp.list$species <- sub("_+$", "", sp.list$species)
  if (any(duplicated(sp.list$species))) {
    print("Duplicated species detected and removed.")
    print(data.frame(species=unique(sp.list$species[duplicated(sp.list$species)])))
  }
  sp.list <- sp.list[!duplicated(sp.list$species), ]
  sp.list$genus <- sub(" +$", "", sp.list$genus)
  sp.list$family <- sub(" +$", "", sp.list$family)
  sp.list$species.relative <- sub(" +$", "", sp.list$species.relative)
  sp.list$genus.relative <- sub(" +$", "", sp.list$genus.relative)
  sp.list[sp.list == ""] <- NA
  if (any(is.na(sp.list$genus))|any(is.na(sp.list$family)))
  {
    y <- which(is.na(sp.list$genus))
    z <- which(is.na(sp.list$family))
    sp.list$genus[y] <- sub("_.*", "", sp.list$species[y])
    sp.list$family[z] <- gen.list$family[match(sp.list$genus[z], gen.list$genus)]
  }
  sp.list$species.relative <- gsub(" ", "_", sp.list$species.relative)
  sp.list$species.relative <- gsub("(^[[:alpha:]])", "\\U\\1", sp.list$species.relative, perl = TRUE)
  sp.list$genus.relative <- gsub("(^[[:alpha:]])", "\\U\\1", sp.list$genus.relative, perl = TRUE)
  sp.list.original <- sp.list
  oriN <- tree$node.label
  tree$node.label <- paste("N", 1 : tree$Nnode, sep = "")
  add.tip <- sp.list[!(sp.list$species %in% tree$tip.label), ]
  status <- rep("prune", dim(sp.list)[1])
  if (dim(add.tip)[1] == 0 & length(which(sp.list$species %in% tree$tip.label)) == 0)
    stop("Incorrect format of species list.")
  if (all(sp.list$species %in% treeX$tip.label) & !is.null(sp.list$species)) {
    print("All species in sp.list are present on tree, return the phylogeny pruned from the tree.")
    splis <- NULL # sp.list.original
    treeX <- drop.tip(treeX, setdiff(treeX$tip.label, sp.list$species))
    return(treeX)
    stop()
  }
  f <- which(sp.list$species %in% tree$tip.label)
  x <- which(sp.list$species.relative %in% tree$tip.label)
  h.sp <- setdiff(x, f)
  fG <- which(sp.list$genus %in% nodes[nodes$level == "G", ]$genus)
  xG <- which(sp.list$genus.relative %in% nodes[nodes$level == "G", ]$genus)
  h.gen <- setdiff(xG, fG)
  h.gen <- setdiff(h.gen, union(x, f))
  h <- union(h.sp, h.gen)
  sel.sp <- sp.list[h.sp, ]
  sel.gen <- sp.list[h.gen, ]
  sel.gen <- sel.gen[!duplicated(sel.gen$genus), ]
  if (dim(sel.gen)[1] > 0) {
    for (i in 1:dim(sel.gen)[1]) {
      n <- match(sel.gen$genus.relative[i], nodes$genus)
      x <- length(tree$tip.label) + which(tree$node.label ==
                                            nodes$bn[n])
      m <- data.frame(level = "G", family = nodes$family[n],
                      genus = sel.gen$genus[i], rn = nodes$bn[n], rn.bl = nodes$rn.bl[n],
                      bn = nodes$bn[n], bn.bl = nodes$bn.bl[n], gen.n = 1,
                      sp.n = 1, taxa = sel.gen$species[i], stringsAsFactors = FALSE)
      nodesN <- rbind(nodesN, m)
      tree <- at.node(tree, x, sel.gen$species[i])
    }
  }
  if (dim(sel.sp)[1] > 0) {
    for (i in 1:dim(sel.sp)[1]) {
      n <- which(tree$edge[, 2] == match(sel.sp$species.relative[i],
                                         tree$tip.label))
      tree <- at.node(tree, tree$edge[n, 1], sel.sp$species[i])
    }
  }
  tree$edge.length <- as.numeric(tree$edge.length)
  tree$node.label <- oriN
  status <- rep("", dim(sp.list)[1])
  status[h] <- "add to relative"
  sp.list.original$status.relative <- status
  if (output.sp.list == FALSE)   sp.list.original <- NULL
  phylo <- list(tree = tree, species.list = sp.list.original,
                nodes.info = nodesN)
  phylo[sapply(phylo, is.null)] <- NULL
  return(phylo)
}
