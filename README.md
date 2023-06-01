# A brief introduction to the U.PhyloMaker package


# 1 Introduction
This package aims to generate phylogenetic trees for a list of species of plants or animals, based on an existent backbone phylogeny. This package contains seven functions adopted from the V.PhyloMaker and V.PhyloMaker2 packages (Jin & Qian, 2019, 2022). Two of the functions, including the main function "phylo.maker" and "bind.relative", were modified to better deal with various situations of manipulating different backbone phylogenies. Below is the introduction of the two modified functions, introduction of the other functions can be found in the V.PhyloMaker and V.PhyloMaker2 packages.

## Install the "U.PhyloMaker" package
> devtools::install_github("jinyizju/U.PhyloMaker")

# 2. phylo.maker
This function is named the same as the "phylo.maker" functions in V.PhyloMaker and V.PhyloMaker2 packages, but the usage of this function differs from those of the latter two. To run "phylo.maker" of this package, you need to input three files, the first is the species list of yours, the second is the backbone phylogeny, the third is a genus-family relationship file that includes at least all the genera that appear on the backbone phylogeny. Below is an example,

## Load the "U.PhyloMaker" package
> library("U.PhyloMaker")

## Download a sample species of amphibians from the online database of megatrees (https://github.com/megatrees)
> sp.list <- read.csv('https://raw.githubusercontent.com/megatrees/amphibian_20221117/main/amphibian_sample_species_list.csv', sep=",")\
> sp.list       

## Download the megatree of amphibians from the online database of megatrees
> megatree <- read.tree('https://raw.githubusercontent.com/megatrees/amphibian_20221117/main/amphibian_megatree.tre')

## Download the the genus-family relationship file of amphibians from the online database of megatrees
> gen.list <- read.csv('https://raw.githubusercontent.com/megatrees/amphibian_20221117/main/amphibian_genus_list.csv', sep=",")

## Run the phylo.maker function with the three files, and output the results
> result <- phylo.maker(sp.list, megatree, gen.list, nodes.type = 1, scenario = 3)\
> result

By default, the output of "phylo.maker" includes two files, one is the phylogeny of "sp.list", the other is "sp.list" with two more columns, one includes the assigned family name(s) in the "genus.list" for the species when different from the family name(s) specified in "sp.list", and the other includes the information about if the species was directly pruned from the backbone phylogeny or was bound to the backbone phylogeny and pruned later.  
### [Note] When importing data from the megatree database, make sure it is the URL of the 'Raw' option on Github page where the data is stored.

# 3 References
Jin Y, & Qian H. (2019) V.PhyloMaker: an R package that can generate very large phylogenies for vascular plants. Ecography, 42, 1353–1359.

Jin Y, & Qian H. (2022) V. PhyloMaker2: An updated and enlarged R package that can generate very large phylogenies for vascular plants. Plant Diversity, 44, 335–339.

