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

# ------------------------------ Loading packages -------------------------------

# To analyze trees, we'll use the packages `ape` and `ggtree`.
# only run this command once in the console if `ape` is not yet installed
#install.packages("ape")
library(ape)
# only run this command once in the console if `ggtree` is not yet installed
#if (!require("remotes", quietly = TRUE)) {
# check that remotes is installed
#  install.packages("BiocManager") 
#}
#library(BiocManager)
#BiocManager::install("ggtree")
library(ggtree)

# Because we're working on the RStudio Server, we already have `grove` installed. 
# However, if you're working at home, uncomment the follow lines and run them in 
# your console (not in this script).
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



