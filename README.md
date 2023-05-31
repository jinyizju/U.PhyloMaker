# A brief introduction to the U.PhyloMaker package


# 1 Introduction
This package aims to generate phylogenetic trees for a list of species of plants or animals, based on an existent backbone phylogeny. This package contains seven functions adopted from the V.PhyloMaker and V.PhyloMaker2 packages (Jin & Qian, 2019, 2022). Two of the functions, including the main function "phylo.maker" and "bind.relative", were modified to better deal with various situations of manipulating different backbone phylogenies. Below is the introduction of the two modified functions, introduction of the other functions can be found in the V.PhyloMaker and V.PhyloMaker2 packages.

## Install the "U.PhyloMaker" package
> devtools::install_github("jinyizju/U.PhyloMaker")

# 2. phylo.maker
This function is named the same as the "phylo.maker" functions in V.PhyloMaker and V.PhyloMaker2 packages, but the usage of this function differs from those of the latter two. To run "phylo.maker" of this package, you need to input three files, the first is the species list of yours, the second is the backbone phylogeny, the third is a genus-family relationship file that includes at least all the genera that appear on the backbone phylogeny. Below is an example,

## Load the "U.PhyloMaker" package
> library("U.PhyloMaker")

## Make up a species list
> c1 <- c("Genus1 add1", "Genus0 add2", "Genus0 add3", "Genus1 t2", "Genus1 add5", "Genus2 t6", "Genus2 add6", "Genus3 add7")\
> c2 <- c("Genus1", "Genus0", "Genus0", "Genus1", "Genus1", "Genus2","Genus2","Genus3")\
> c3 <- c("Family1", "Family1", "Family1", "Family1", "Family1", "Family1", "Family1", "Family3")\
> c4 <- c("Genus1_t1", rep("", 7))\
> sample1 <- data.frame(species = c1, genus = c2, family = c3, species.relative = c4, genus.relative = NA)\
> sample1       

## Make up a backbone phylogeny
> x <- "((((((t1:1,t2:1):1,t3:2):3,t4:5):1,(t5:3,t6:3):3):1,(t7:4,t8:4):3):1,t9:8);"
cat(x, file = "tree.tre", sep = "\n")\
> megatree <- read.tree("tree.tre")\
> g <- rep("Genus2", 9); g[c(1:4, 7)] <- "Genus1"\
> f <- rep("Family2", 9); f[c(1:4, 7)] <- "Family1"\
> megatree$tip.label <- paste(g, sort(megatree$tip.label), sep = "_")

## Make up the genus-family relationship file
> genus.list <- data.frame(genus = c("Genus0", "Genus1", "Genus2", "Genus3"), family = c(rep("Family1", 2), "Family2","Family3"))

## Run the phylo.maker function with the three files, and output the results
> result1 <- U.PhyloMaker::phylo.maker(sp.list = sample1, tree = megatree, gen.list = genus.list, nodes.type = 1, scenario = 3)\
> result1

By default, the output of "phylo.maker" includes two files, one is the phylogeny of "sample1", the other is "sample1" with two more columns, one includes the assigned family name(s) in the "genus.list" for the species when different from the family name(s) specified in "sample1", and the other includes the information about if the species was directly pruned from the backbone phylogeny or was bound to the backbone phylogeny and pruned later.  
The "phylo.maker" fucntion can deal with species list with incorrect or missing information. For example,

## Make up a species list with incorrect and missing information
> c1 <- c("Genus1 add1", "Genus0 add2", "Genus0 add3", "Genus1 t2", "Genus1 add5", "Genus2 t6", "Genus2 add6", "Genus3 add7")\
> c2 <- c("Genus2", "", "", "", "Genus1", "Genus2","Genus2","Genus3")\
> c3 <- c("Family1", "Family1", "", "", "Family1", "Family2", "Family1", "")\
> c4 <- c("Genus1 t1", rep("", 7)) \
> sample2 <- data.frame(species = c1, genus = c2, family = c3, species.relative = c4, genus.relative = NA)\
> sample2       

## Use "sample2" to replace "sample1", and re-run the example above, and output the results
> result2 <- U.PhyloMaker::phylo.maker(sp.list = sample2, tree = megatree, gen.list = genus.list, nodes.type = 1, scenario = 3)\
> result2

In the extreme case, the "phylo.maker" function can also deal with a species list that includes only the first column, namely the names of the species.

## Make up a species list with only species names
> sample3 <- c("Genus1 add1", "Genus0 add2", "Genus0 add3", "Genus1 t2", "Genus1 add5", "Genus2 t6", "Genus2 add6", "Genus3 add7")

## Use "sample3" to replace "sample2", and re-run the example above, and output the results.
> result3 <- U.PhyloMaker::phylo.maker(sp.list = sample2, tree = megatree, gen.list = genus.list, nodes.type = 1, scenario = 3)\
> result3

## Compare the three phylogenies respectively generated based on the three sample specie lists above.
> par(mfrow = c(1, 3), mar = c(1, 0, 1, 1))\
> plot.phylo(result1$phylo, cex =1.5, main = "sample1")\
> plot.phylo(result2$phylo, cex =1.5, main = "sample2")\
> plot.phylo(result3$phylo, cex =1.5, main = "sample3")

# 3 bind.relative
This function is named the same as the "bind.relative" functions in V.PhyloMaker and V.PhyloMaker2 packages, however, the usage of this function is different from the previous ones. Similar to the "phylo.maker" function in V.PhyloMaker and V.PhyloMaker2, to run "bind.relative" of this package, you need to input three files, the first is the species list of yours, the second is the backbone phylogeny, the third is a genus-family relationship file that includes at least all the genera that appear on the backbone phylogeny.
You can use the species list "sample1", the backbone phylogeny "megatree" and the genus-family relationship file "genus.list" that made up above to test "bind.relative", as shown:

## Run the function "bind.relative", and output the results.
> result4 <- U.PhyloMaker::bind.relative(sp.list = sample1, tree = megatree, gen.list = genus.list, nodes.type = 1, output.sp.list = TRUE)\
> result4

Running "bind.relative" will result in two files, one is the input species list with a column indicating the species that have been bound to the relative on the megatree, the other is the megatree with the species bound to the relative on the megatree. 

## Plot the megatree with species (in red color) binded to its designated relative on the megatree.
> plot.phylo(result4$tree, cex = 1.5, tip.color = c("red", rep("black",9)), main = "bind.relative")

Then, you can run the "phylo.maker" function with the output files of "bind.relative" and the "genus.list", to bind/prune the rest of species in your species list to the megatree and generate the more resolved phylogeny compared with just using "phylo.maker" to generate phylogeny on all the species in your species list.

## run the function bind.relative, and output the results
> result5 <- U.PhyloMaker::phylo.maker(sp.list = result4$species.list, tree = result4$tree, gen.list = genus.list, nodes.type = 1, output.tree = TRUE, scenario = 3)\
> result5

## Compare the generated phylogeny of your species list, between using only "phylo.maker" and using "bind.relative" in combination with "phylo.maker". 
> par(mfrow = c(1, 2), mar = c(1, 0, 1, 1))\
> plot.phylo(result1$phylo, cex =1.5, main = "phylo.maker", tip.color = c(rep("black", 4),"red","black","black"))\
> plot.phylo(result5$phylo, cex =1.5, tip.color = c(rep("black", 4),"red","black","black"), main = "bind.relative + phylo.maker")


# 4 References
Jin Y, & Qian H. (2019) V.PhyloMaker: an R package that can generate very large phylogenies for vascular plants. Ecography, 42, 1353–1359.

Jin Y, & Qian H. (2022) V. PhyloMaker2: An updated and enlarged R package that can generate very large phylogenies for vascular plants. Plant Diversity, 44, 335–339.

