## gene tree visualization - STAMPS 2022 
## Sarah Teichman, Amy Willis 
## Friday 29 July 2022 

# --------------------- Introduction -------------------------

# This tutorial will work through visualizing a phylogenomic 
# tree and gene trees. 
# The R package to do this analysis is still under construction.

# ------------------ Loading packages and data ---------------

# To analyze trees, we'll use the package `ape` 
install.packages("ape") # only run this command once 
library(ape)

# To make our gene tree plots, we'll use the package `TreeVizPackage`
# First make sure that we have `remotes`
if (!require("remotes", quietly = TRUE)) {
  install.packages("remotes") # check that remotes is installed
}
remotes::install_github("svteichman/TreeVizPackage")

# We'll start with the output from GToTree, which includes a 
# phylogenomic tree and gene level alignments. To ask GToTree to 
# save your gene level alignments, run it with flag -k. 

# First, read in the trees. 
# Start with the phylogenomic tree.
phy_genom <- read.tree("https://raw.githubusercontent.com/statdivlab/stamps2022/main/Friday-afternoon/labs/gene_tree_lab/data/tree_files/AI-2E_transport_aln.faa.treefile")

# read in gene trees 

# look at presence/absence 

# drop tips for phylogenomic tree (if necessary)

# plot phylogenomic tree

# make pca plot 

# identify interesting trees 

# plot trees 



