## gene tree visualization - STAMPS 2022 
## Sarah Teichman, Amy Willis 
## Friday 29 July 2022 

# ------------------------------ Introduction ---------------------------------

# This tutorial will work through visualizing a phylogenomic tree and gene trees. 
# The R package to do this analysis is still under development.

# In this tutorial we will be consider a set of Prevotella genomes. We generated
# this genome set in order to have broad representation across the Prevotella genus 
# and work with high quality genomes. We started with the 383 species-level representatives
# in the Genome Taxonomy Database (GTDB). The GTDB is a database described in 
# Parks et al. "A standardized bacterial taxonomy based on genome phylogeny 
# substantially revises the tree of life". Nature biotechnology, 36(10):996–1004, 
# 2018. It contains microbial genomes, and was created in an effort to build a 
# microbial taxonomy based on phylogeny. To develop a set of genes to analyze, we 
# used HMMER to identify all protein families in the Pfam database with hits in at
# least 90% of our 383 Prevotella genomes. We identified 63 genes, and found
# that 73 of our genomes contained all 63 genes of interest. This gave us a set
# of genomes that each contain all of our genes of interest, covering much of the 
# breadth of the Prevotella genus as well as a variety of genes.

# We curated this set of genes and genomes using functionality implemented in
# GToTree and another one of Mike's collections of useful tools, bit:
# https://github.com/AstrobioMike/bit 
# https://f1000research.com/articles/11-122/v1

# ------------------------------ Loading packages -------------------------------

# To analyze trees, we'll use the packages `ape` and `ggtree`.
# only run this command once in the console if `ape` is not yet installed
#install.packages("ape")
library(ape)
# only run this command once in the console if `ggtree` is not yet installed
#if (!require("remotes", quietly = TRUE)) {
#library(BiocManager)
#BiocManager::install("ggtree")
library(ggtree)

# Because we're working on the RStudio Server, we already have `grove` installed. 
# However, if you're working at home, uncomment the follow lines and run them in 
# your console (not in this script).
# check that remotes is installed
#  install.packages("remotes") 
#}
#remotes::install_github("statdivlab/grove")

# We'll start with the output from GToTree, which includes a 
# phylogenomic tree and gene level alignments. To ask GToTree to 
# save your gene level alignments, run it with flag -k. 

# ------------------------- Loading in data -------------------------------------

# Before visualizing a set of gene trees and a phylogenomic tree, there are 
# several data processing steps. These are omitted here for clarity. However,
# a tutorial on processing your data from GToTree output to a set of gene trees
# and a phylogenomic tree is under development. If you want to run this tool on
# your own genomes and genes, send me an email at teichs@uw.edu and I'll send you
# recommendations for data processing. The general steps are the following: 

# (1) Look at gene level alignments and identify which genes are present in which
#      genomes. 

# (2) Remove genes or genomes until you have a subset of genes that are present in
#      a subset of genomes. You can choose this subset algorithmically or based on
#      certain genes or genomes that you want to retain for scientific reasons. 

# (3) Build a new phylogenomic tree using only the genes and genomes that you've 
#      selected for your analysis. 

# (4) Estimate gene level trees. I typically do this with IQTree with the following
#      command: IQTree -s $gene_file_name -mset WAG,LG -bb 1000
# Citation: 
# B.Q. Minhet al. (2020) "IQ-TREE 2: New models and efficient methods for 
# phylogenetic inference in the genomic era". Mol. Biol. Evol., 37:1530-1534.

# For our purposes, I've already done all of these steps. We can load in the
# phylogenomic tree and the set of gene trees. 

# Phylogenomic trees

phy_genom <- read.tree("https://raw.githubusercontent.com/statdivlab/stamps2022/main/Friday-morning/labs/gene_tree_lab/data/trees_txt/phy_genom.txt")

# Let's look at our tree object! What information does it have?
# You can use the $ operator to see what information is saved within the object
# phy_genom

phy_genom

# A phylo object contains the following: 
# edge: a list of branches that connect one node of the tree to another
# edge.length: a list of branch lengths for each branch 
# Nnode: the number of internal nodes on the tree, these represent ancestor organisms
# tip.label: a vector of names of tips (in our case, genomes)

# We can plot our tree using ggtree (this is a very cool plotting package that 
# extends ggplot2 for phylogenetic trees). 

ggtree(phy_genom) + 
  geom_tiplab(size = 2)

# Gene trees 

gene_trees <- read.tree("https://raw.githubusercontent.com/statdivlab/stamps2022/main/Friday-morning/labs/gene_tree_lab/data/trees_txt/gene_trees.txt")

# Let's take a look at our gene trees

gene_trees

# Here we can see that this is a set of 63 phylogenetic trees. If we want to look
# at a specific gene tree, we can use the `[[]]` subsetting operator.

gene_trees[[1]]

# We can also plot any of our gene trees.

ggtree(gene_trees[[1]]) + 
  geom_tiplab(size = 2)

# What do we notice? 
# This tree looks different than the phylogenomic tree! Actually, it turns
# out that the topology (branching structure) is slightly different between each
# gene tree and between the gene trees and phylogenomic tree. 

# The last thing we need to load in is a list of our gene names. We can do this 
# with: 

gene_names <- as.character(read.delim("https://raw.githubusercontent.com/statdivlab/stamps2022/main/Friday-morning/labs/gene_tree_lab/data/gene_names.txt", 
                                      header = FALSE)$V1)

# -------------------- Building ordination plot ---------------------------------

# First we need a vector that contains the pathways to each of our individual
# gene trees. We'll use the function paste0, which combines two strings together
# into one string.

pathways <- paste0("https://raw.githubusercontent.com/statdivlab/stamps2022/main/Friday-morning/labs/gene_tree_lab/data/trees_txt/",
                   gene_names, ".txt")

#-------
vectors <- compute_logmap(cons_path = "https://raw.githubusercontent.com/statdivlab/stamps2022/main/Friday-morning/labs/gene_tree_lab/data/trees_txt/phy_genom.txt", 
                          tree_paths = pathways,
                          add_pendant_branches = TRUE,
                          cons_tree = phy_genom,
                          trees_complete = c(phy_genom, trees))

other_paths <- paste0("Friday-morning/labs/gene_tree_lab/data/trees_txt/", 
                      gene_names, ".txt")
# this one works! 
vectors <- compute_logmap(cons_path = "Friday-morning/labs/gene_tree_lab/data/trees_txt/phy_genom.txt", 
                          tree_paths = other_paths,
                          add_pendant_branches = TRUE,
                          cons_tree = phy_genom,
                          trees_complete = c(phy_genom, gene_trees))

vectors <- compute_logmap(cons_path = "https://raw.githubusercontent.com/statdivlab/stamps2022/main/Friday-morning/labs/gene_tree_lab/data/trees_txt/phy_genom.txt", 
                          tree_paths = other_paths,
                          add_pendant_branches = TRUE,
                          cons_tree = phy_genom,
                          trees_complete = c(phy_genom, gene_trees))
#-------------
lm_res <- plot_logmap(vectors = vectors, 
                      base_name = base_name, 
                      gene_names = gene_names,
                      cons_exists = FALSE)
lm_res$plot
















# First, we're going to check which genes were found in which genomes 
# We'll read in gene_names from the file "gene_names.txt"

gene_names <- as.character(read.delim("https://raw.githubusercontent.com/statdivlab/stamps2022/main/Friday-morning/labs/gene_tree_lab/data/gene_names.txt", 
                                      header = FALSE)$V1)

# Note, in order to run `grove`, you need a complete set of genes and genomes.
# This means that each gene that you want to investigate needs to have been
# observed each of your genomes. In practice, this is rarely the case due to 
# variation in gene presence/absence across taxa and errors in the sequencing 
# and processing of data. To deal with this, you will need to subsample either 
# the set of genes or genomes or both to consider a set of genomes that each
# include all genes. Stay tuned for a tutorial with recommendations on how to 
# do this! 

# Because we chose the genomes for this analysis specifically so that they would
# include a set of 63 genes that we would like to examine in Prevotella, we are
# starting with a set of genomes that include all of our genes of interest. 

# Next, read in the trees. 
# Start with the phylogenomic tree.

phy_genom <- read.tree("https://raw.githubusercontent.com/statdivlab/stamps2022/main/Friday-morning/labs/gene_tree_lab/data/tree_files/Prevotella-GToTree-res.tre")

# let's take a look at the phylogenomic tree 
# with a tree this large, we won't add tip labels for now 

ggtree(phy_genom, size = 0.5) 

# how many tips does this tree have? 

length(phy_genom$tip.label)

# It has 380. This is because the concatenated tree accounts 

# Next, we will load each gene tree into our environment.
# Gene trees were each generated with the following command line code:
# IQTree -s $gene_file_name -mset WAG,LG -bb 1000
# Citation: 
# B.Q. Minhet al. (2020) "IQ-TREE 2: New models and efficient methods for 
# phylogenetic inference in the genomic era". Mol. Biol. Evol., 37:1530-1534.

# We want to load each gene tree into our environment. 
curr_gene <- gene_names[1]
# read the estimated gene tree (output from IQTree)
tree <- read.tree(paste0("https://raw.githubusercontent.com/statdivlab/stamps2022/main/Friday-morning/labs/gene_tree_lab/data/tree_files/",
                         curr_gene, "_aln.faa.treefile"))
# remove the node labels from our trees 
tree$node.label <- NULL
# save tree in a set of trees
trees <- tree
for (i in 2:length(gene_names)) {
  curr_gene <- gene_names[i]
  tree <- read.tree(paste0("https://raw.githubusercontent.com/statdivlab/stamps2022/main/Friday-morning/labs/gene_tree_lab/data/tree_files/",
                           curr_gene, "_aln.faa.treefile"))
  tree$node.label <- NULL
  trees <- c(trees, tree)
}

# we now have a multiPhylo object containing all of our gene trees
trees

# look at presence/absence 

# drop tips for phylogenomic tree (if necessary)

# plot phylogenomic tree

# make pca plot 

# identify interesting trees 

# plot trees 



