phylo.maker <- function (sp.list, tree, gen.list, nodes.type = 1, output.sp.list = TRUE,
                         output.tree = FALSE, scenarios = 3, r = 1)
{
  if (!all(sub("_.*", "", tree$tip.label) %in% gen.list$genus)&isFALSE(output.tree)) {
    tree <- drop.tip(tree, which(sub("_.*", "", tree$tip.label) %in% setdiff(sub("_.*", "", tree$tip.label), gen.list$genus)))
  }
  treeX <- tree
  d1 <- data.frame(species = tree$tip.label, genus = sub("_.*", "", tree$tip.label))
  d1 <- merge(d1, gen.list, all.x = T)[, union(colnames(d1), colnames(gen.list))]
  if (nodes.type == 1)   nodes <- build.nodes.1(tree, d1)
  if (nodes.type == 2)   nodes <- build.nodes.2(tree, d1)
  nodes[sapply(nodes, is.factor)] <- lapply(nodes[sapply(nodes, is.factor)], as.character)

  if (is.null(tree$node.label))  tree$node.label <- rep("", tree$Nnode)
  rnN <- data.frame(node.label = paste("N", 1 : tree$Nnode,
                                       sep = ""), oriN = tree$node.label, stringsAsFactors = FALSE)
  tree$node.label <- paste("N", 1 : tree$Nnode, sep = "")  #

  if (is.factor(sp.list)|(is.character(sp.list)))  sp.list <- data.frame(species = sp.list, genus = NA, family = NA)
  if (is.data.frame(sp.list) & dim(sp.list)[2] < 3)   sp.list[, 2 : 3] <- NA

  dimnames(sp.list)[[2]][1:3] <- c("species", "genus", "family")
  sp.list[sapply(sp.list, is.factor)] <- lapply(sp.list[sapply(sp.list, is.factor)], as.character)
  sp.list$species <- gsub("(^[[:alpha:]])", "\\U\\1", sp.list$species, perl = TRUE)
  sp.list$genus <- gsub("(^[[:alpha:]])", "\\U\\1", sp.list$genus, perl = TRUE)
  sp.list$family <- gsub("(^[[:alpha:]])", "\\U\\1", sp.list$family, perl = TRUE)

  sp.list$species <- sub("_+$", "", sp.list$species)
  sp.list <- sp.list[sp.list$species != "", ]
  sp.list$species <- gsub(" ", "_", sp.list$species)

  if (any(duplicated(sp.list$species))) {
    print("Duplicated species detected and removed.")
    print(data.frame(species=unique(sp.list$species[duplicated(sp.list$species)])))
    sp.list <- sp.list[!duplicated(sp.list$species), ]
  }

  sp.list$genus <- sub(" +$", "", sp.list$genus)
  sp.list$family <- sub(" +$", "", sp.list$family)
  sp.list[sp.list == ""] <- NA

  sp.list.original <- sp.list   #
  sp.list.original$species <- gsub("_"," ", sp.list.original$species)  #

  sp.list$genus <- sub("_.*", "", sp.list$species) #

  m <- gen.list
  dimnames(m)[[2]] <- gsub("family", "family_in_gen.list", dimnames(m)[[2]])
  m0 <- sp.list[, c("species", "genus", "family")]
  dimnames(m0)[[2]][3] <- "family_in_sp.list"
  mm <- merge(m0, m, all.x = T)        ##
  g <- mm[-which(mm$family_in_sp.list == mm$family_in_gen.list),]
  sp.list.original$family.in.genus.list <- NA
  if (dim(g)[1] > 0) {
    print("Taxonomic classification not consistent between sp.list and gen.list.")
    print(g[, c(1,3,4)])
    g$species <- gsub("_"," ", g$species)
    sp.list.original$family.in.genus.list[match(g$species, sp.list.original$species)] <- as.character(g$family_in_gen.list)
  }
  if (any(is.na(sp.list$family)))    #
  {                                                         #
    z <- which(is.na(sp.list$family))   #
    sp.list$family[z] <- sp.list.original$family.in.genus.list[z] <- as.character(gen.list$family[match(sp.list$genus[z], gen.list$genus)])
  }   #


  add.tip <- sp.list[!(sp.list$species %in% tree$tip.label), ]
  status <- rep("present in megatree", dim(sp.list)[1])
  if (dim(add.tip)[1] == 0 & length(which(sp.list$species %in% tree$tip.label)) == 0)
    stop("Incorrect format of species list.")
  if (all(sp.list$species %in% treeX$tip.label) & !is.null(sp.list$species)) {
    print("All species in sp.list are present on the tree, return the phylogeny pruned from the tree.")
   # splis <- NULL # sp.list.original
    treeX <- drop.tip(treeX, setdiff(treeX$tip.label, sp.list$species))
    return(treeX)
    stop()
  }
  add.tip$sort <- ""
  add.tip$sort[add.tip$genus %in% nodes[nodes$level == "G", ]$genus] <- "G1"
  add.tip$sort[(!add.tip$genus %in% nodes[nodes$level == "G", ]$genus) & (add.tip$family %in% nodes[nodes$level == "F", ]$family)] <- "F1"
  add.tip$sort[add.tip$sort == "F1"][duplicated(add.tip[add.tip$sort == "F1", ]$genus)] <- "F2"
  status[match(add.tip$species, sp.list$species)] <- "insertion based on genus"  #
  status[match(add.tip[add.tip$sort != "G1", ]$species, sp.list$species)] <- "insertion based on family"  #
  a <- which(add.tip$sort == "")
  if (length(a) > 0) {
    print(paste("Note:", length(a), "species fail to be binded to the tree.", sep = " "))
    print(add.tip$species[a])
    status[match(add.tip$species[a], sp.list$species)] <- "no insertion"
  }
  sp.list.original$output.note <- status
  if (1 %in% scenarios) {
    t1 <- tree
    rnN1 <- rnN
    nG <- nodes[nodes$level == "G", ]
    nF <- nodes[nodes$level == "F", ]
    data <- add.tip[add.tip$sort == "F1" | add.tip$sort == "F2", ]
    if (dim(data)[1] > 0) {
      for (i in 1:dim(data)[1]) {
        num <- nF$bn[match(data$family[i], nF$family)]
        t1 <- at.node(t1, location.node = num, tip.label = data$species[i])
      }
    }
    data <- add.tip[add.tip$sort == "G1", ]
    if (dim(data)[1] > 0) {
      for (i in 1:dim(data)[1]) {
        num <- nG$bn[match(data$genus[i], nG$genus)]
        t1 <- at.node(t1, location.node = num, tip.label = data$species[i])
      }
    }
    t1$edge.length <- as.numeric(t1$edge.length)
    tree1 <- t1
    tree1$node.label <- rnN1$oriN[match(tree1$node.label, rnN1$node.label)]
    toDrop <- setdiff(1:length(t1$tip.label), which(t1$tip.label %in% sp.list$species))
    t1 <- drop.tip(t1, tip = toDrop)
    Re <- which(t1$node.label %in% rnN1$node.label)
    noRe <- which(!t1$node.label %in% rnN1$node.label)
    t1$node.label[Re] <- rnN1$oriN[match(t1$node.label, rnN1$node.label)[Re]]
    t1$node.label[noRe] <- ""
  }
  else {
    t1 <- NULL
    tree1 <- NULL
  }
  if (2 %in% scenarios) {
    t2r <- replicate(r, list())
    names(t2r) <- paste("run", 1:r, sep = ".")
    tree2r <- replicate(r, list())
    names(tree2r) <- paste("run", 1:r, sep = ".")
    for (o in 1:r) {
      t2 <- tree
      rnN2 <- rnN
      nG <- nodes[nodes$level == "G", ]
      nF <- nodes[nodes$level == "F", ]
      data <- add.tip[add.tip$sort == "F1", ]
      if (dim(data)[1] > 0) {
        for (i in 1:dim(data)[1]) {
          n <- match(data$family[i], nF$family)
          g <- nF$gen.n[n]
          s <- nF$sp.n[n]
          if (g == 1 & s == 1) {
            num <- match(nF$taxa[n], t2$tip.label)
            nlabel <- paste("N", t2$Nnode + 1, sep = "")
            t2 <- ext.node(t2, location.tip = num, tip.label = data$species[i],
                           node.label = nlabel, position = 2/3)
            nF$gen.n[n] <- g + 1
            nF$sp.n[n] <- s + 1
            x <- which(t2$node.label == nlabel)
            xx <- rbind(rnN2[1:(x - 1), ], data.frame(node.label = nlabel, oriN = ""))
            xx <- rbind(xx, rnN2[x:dim(rnN2)[1], ])
            rnN2 <- xx
            nF$bn[n] <- nlabel
          }
          else
          {
            num <- sample(nG$bn[nG$family %in% data$family[i]], 1)
            t2 <- at.node(t2, location.node = num, tip.label = data$species[i])
            nF$gen.n[n] <- g + 1
            nF$sp.n[n] <- s + 1
          }
        }
      }
      data <- add.tip[add.tip$sort == "F2", ]
      if (dim(data)[1] > 0) {
        for (i in 1 : dim(data)[1]) {
          n <- which(sub("_.*", "", t2$tip.label) == data$genus[i])
          nlabel <- paste("N", t2$Nnode + 1, sep = "")
          if (length(n) == 1) {
            num <- n
            t2 <- ext.node(t2, location.tip = num, tip.label = data$species[i],
                           node.label = nlabel, position = 1/2)
            x <- which(t2$node.label == nlabel)
            xx <- rbind(rnN2[1 : (x - 1), ], data.frame(node.label = nlabel, oriN = ""))
            xx <- rbind(xx, rnN2[x : dim(rnN2)[1], ])
            rnN2 <- xx
          }
          if (length(n) > 1) {
            num <- t2$edge[which(t2$edge[, 2] == n[1]),
                           1]
            t2 <- at.node(t2, location.node = num, tip.label = data$species[i])
          }
        }
      }
      data <- add.tip[add.tip$sort == "G1", ]
      if (dim(data)[1] > 0) {
        for (i in 1 : dim(data)[1]) {
          n0 <- match(data$genus[i], nG$genus)
          n <- nG$sp.n[n0]
          nlabel <- paste("N", t2$Nnode + 1, sep = "")
          if (n == 1) {
            num <- t2$tip.label[match(nG$taxa[n0], t2$tip.label)]
            t2 <- ext.node(t2, location.tip = num, tip.label = data$species[i],
                           node.label = nlabel, position = 1/2)
            x <- which(t2$node.label == nlabel)
            xx <- rbind(rnN2[1 : (x - 1), ], data.frame(node.label = nlabel,
                                                        oriN = "", stringsAsFactors = FALSE))
            xx <- rbind(xx, rnN2[x : dim(rnN2)[1], ])
            rnN2 <- xx
            nG$sp.n[n0] <- nG$sp.n[n0] + 1
          }
          if (n > 1) {
            num <- which(t2$node.label == nG$bn[n0]) +
              length(t2$tip.label)
            num1 <- which(t2$edge[, 1] %in% num)
            part1 <- t2$edge[1 : min(num1), ]
            n1 <- max(which(part1[, 1] < num), 0) + 1
            part2 <- t2$edge[max(num1) : dim(t2$edge)[1],
            ]
            n2 <- min(which(part2[, 1] < num), dim(part2)[1] +
                        1) + max(num1) - 2
            sect <- t2$edge[n1 : n2, ]
            sect <- sort(unique(c(sect[, 1], sect[, 2])))
            sect <- sect[which(sect > length(t2$tip.label))]
            num2 <- sect[sample(1 : length(sect), 1)]
            t2 <- at.node(t2, location.node = num2, tip.label = data$species[i])
          }
        }
      }
      t2$edge.length <- as.numeric(t2$edge.length)
      tree2 <- t2
      tree2$node.label <- rnN2$oriN[match(tree2$node.label,
                                          rnN2$node.label)]
      toDrop <- setdiff(1:length(t2$tip.label), which(t2$tip.label %in% sp.list$species))
      t2 <- drop.tip(t2, tip = toDrop)
      Re <- which(t2$node.label %in% rnN2$node.label)
      noRe <- which(!t2$node.label %in% rnN2$node.label)
      t2$node.label[Re] <- rnN2$oriN[match(t2$node.label, rnN2$node.label)[Re]]
      t2$node.label[noRe] <- ""
      t2r[[o]] <- t2
      tree2r[[o]] <- tree2
    }
  }
  else {
    t2r <- NULL
    tree2r <- NULL
  }
  if (3 %in% scenarios) {
    t3 <- tree
    rnN3 <- rnN
    nG <- nodes[nodes$level == "G", ]
    nF <- nodes[nodes$level == "F", ]
    data <- add.tip[add.tip$sort == "F1", ]
    if (dim(data)[1] > 0) {
      for (i in 1 : dim(data)[1]) {
        n <- match(data$family[i], nF$family)
        g <- nF$gen.n[n]
        s <- nF$sp.n[n]
        if (g == 1 & s == 1) {
          num <- match(nF$taxa[n], t3$tip.label)
          nlabel <- paste("N", t3$Nnode + 1, sep = "")
          t3 <- ext.node(t3, location.tip = num, tip.label = data$species[i],
                         node.label = nlabel, position = 2/3)
          nF$gen.n[n] <- g + 1
          nF$sp.n[n] <- s + 1
          x <- which(t3$node.label == nlabel)
          xx <- rbind(rnN3[1 : (x - 1), ], data.frame(node.label = nlabel,
                                                      oriN = ""))
          xx <- rbind(xx, rnN3[x : dim(rnN3)[1], ])
          rnN3 <- xx
          nF$bn[n] <- nlabel
        }
        if (g == 1 & s > 1) {
          nlabel <- paste("N", t3$Nnode + 1, sep = "")
          if ((2/3) * nF$rn.bl[n] <= 3 * nF$bn.bl[n]) {
            len <- (nF$rn.bl[n] - nF$bn.bl[n])/2
          }
          if ((2/3) * nF$rn.bl[n] > 3 * nF$bn.bl[n]) {
            len <- nF$rn.bl[n] * 2/3 - nF$bn.bl[n]
          }
          port <- len / (nF$rn.bl[n] - nF$bn.bl[n])
          t3 <- int.node(t3, location.node = nF$bn[n],
                         tip.label = data$species[i], node.label = nlabel,
                         position = port)
          nF$gen.n[n] <- g + 1
          nF$sp.n[n] <- s + 1
          x <- which(t3$node.label == nlabel)
          xx <- rbind(rnN3[1 : (x - 1), ], data.frame(node.label = nlabel,
                                                      oriN = "", stringsAsFactors = FALSE))
          xx <- rbind(xx, rnN3[x : dim(rnN3)[1], ])
          rnN3 <- xx
          nF$bn[n] <- nlabel
        }
        if (g > 1) {
          t3 <- at.node(t3, location.node = nF$bn[n],
                        tip.label = data$species[i])
        }
      }
    }
    data <- add.tip[add.tip$sort == "F2", ]
    if (dim(data)[1] > 0) {
      for (i in 1 : dim(data)[1]) {
        n <- which(sub("_.*", "", t3$tip.label) == data$genus[i])
        nlabel <- paste("N", t3$Nnode + 1, sep = "")
        if (length(n) == 1) {
          t3 <- ext.node(t3, location.tip = n, tip.label = data$species[i],
                         node.label = nlabel, position = 1/2)
          nG$sp.n[match(data$genus[i], nG$genus)] <- length(n) + 1
          x <- which(t3$node.label == nlabel)
          xx <- rbind(rnN3[1 : (x - 1), ], data.frame(node.label = nlabel,
                                                      oriN = "", stringsAsFactors = FALSE))
          xx <- rbind(xx, rnN3[x : dim(rnN3)[1], ])
          rnN3 <- xx
          nG$bn[match(data$genus[i], nG$genus)] <- nlabel
        }
        if (length(n) > 1) {
          num <- ancestor(t3, min(n), max(n))
          t3 <- at.node(t3, location.node = num, tip.label = data$species[i])
        }
      }
    }
    data <- add.tip[add.tip$sort == "G1", ]
    if (dim(data)[1] > 0) {
      for (i in 1 : dim(data)[1]) {
        n0 <- match(data$genus[i], nG$genus)
        s <- nG$sp.n[n0]
        nlabel <- paste("N", t3$Nnode + 1, sep = "")
        if (s == 1) {
          num <- t3$tip.label[match(nG$taxa[n0], t3$tip.label)]
          t3 <- ext.node(t3, location.tip = num, tip.label = data$species[i],
                         node.label = nlabel, position = 1/2)
          nG$sp.n[n0] <- nG$sp.n[n0] + 1
          x <- which(t3$node.label == nlabel)
          xx <- rbind(rnN3[1 : (x - 1), ], data.frame(node.label = nlabel,
                                                      oriN = "", stringsAsFactors = FALSE))
          xx <- rbind(xx, rnN3[x : dim(rnN3)[1], ])
          rnN3 <- xx
          nG$bn[n0] <- nlabel
        }
        if (s > 1) {
          t3 <- at.node(t3, location.node = nG$bn[n0],
                        tip.label = data$species[i])
        }
      }
    }
    t3$edge.length <- as.numeric(t3$edge.length)
    tree3 <- t3
    tree3$node.label <- rnN3$oriN[match(tree3$node.label, rnN3$node.label)]
    toDrop <- setdiff(1 : length(t3$tip.label), which(t3$tip.label %in% sp.list$species))
    t3 <- drop.tip(t3, tip = toDrop)
    Re <- which(t3$node.label %in% rnN3$node.label)
    noRe <- which(!t3$node.label %in% rnN3$node.label)
    t3$node.label[Re] <- rnN3$oriN[match(t3$node.label, rnN3$node.label)[Re]]
    t3$node.label[noRe] <- ""
  }
  else {
    t3 <- NULL
    tree3 <- t3
  }
  if (output.sp.list == FALSE) {
    sp.list.original <- NULL
  }
  if (output.tree == FALSE) {
    tree1 <- tree2r <- tree3 <- NULL
  }
  if (r == 1) {
    t2r <- t2r$run.1
    tree2r <- tree2r$run.1
  }
  phylo <- list(s1 = t1, s2 = t2r, s3 = t3,
                sp.list = sp.list.original, tree.s1 = tree1,
                tree.s2 = tree2r, tree.s3 = tree3)
  phylo[sapply(phylo, is.null)] <- NULL
  if (length(phylo)==1)  phylo <- phylo[1]
  if (length(scenarios) == 1)  names(phylo) [1] <- "phylo"
  return(phylo)
}
