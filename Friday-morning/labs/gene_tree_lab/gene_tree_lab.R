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

# To analyze trees, we'll use the packages `ape`, `plotly` and `ggtree`.
# only run this command once in the console if `ape` is not yet installed
#install.packages("ape")
library(ape)

# only run this command in the console if `plotly` is not yet installed
#install.packages("plotly")
library(plotly)

# only run this command once in the console if `ggtree` is not yet installed
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
library(grove)

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

# For our purposes, these steps have already been done for you. We can load in the
# phylogenomic tree and the set of gene trees. 

# Uncomment the following lines and run in your console. This will create a new 
# folder in your current working directly that has txt files will all of our 
# trees. 

#download.file("https://github.com/statdivlab/stamps2022/raw/main/Friday-morning/labs/gene_tree_lab/data/trees_txt.zip", 
#"trees_txt.zip")
#unzip("trees_txt.zip")

# Phylogenomic trees

phy_genom <- read.tree("trees_txt/phy_genom.txt")

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

gene_trees <- read.tree("trees_txt/gene_trees.txt")

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

pathways <- paste0("trees_txt/", gene_names, ".txt")

# Now we will run the function `compute_logmap` to approximate all of our trees
# as vectors in Euclidean space

# the path to the phylogenomic tree, to do our approximation
# we need a base tree to compare all other trees to, 
# and we use the phylogenomic tree as the base.
vectors <- compute_logmap(cons_path = "trees_txt/phy_genom.txt", 
                          # a vector of paths to all gene trees we want to study
                          tree_paths = pathways,
                          # this lets us consider the lengths of branches to the 
                          # tips, by default these are not considered but we like
                          # to consider them because the amount of evolution from
                          # the most recent common ancestor to the current organism
                          # may be of biological interest. 
                          add_pendant_branches = TRUE,
                          # the phylogenomic tree
                          cons_tree = phy_genom,
                          # a set of trees including the phylogenomic tree and 
                          # all of the gene trees
                          trees_complete = c(phy_genom, gene_trees))

# Let's take a quick look at this output

dim(vectors)

# This output has 64 rows (for the phylogenomic tree and 63 gene trees) and 153
# columns. These 153 columns approximate the information encoded in the topology
# (branching order) and branch lengths from the set of trees in Euclidean space.

# Next, let's build up our ordination plot. The next function will run PCA (our
# preferred ordination procedure) and plot the first two principal components.
# The idea behind this is that we want to visualize our set of trees in two 
# dimensions in a way that will retain as much variation between the trees as 
# possible. This let's us get a sense of which trees are more similar, which are
# more different, and which might be outliers that could be interesting to 
# investigate.

# The approximation of our trees in Euclidean space
lm_res <- plot_logmap(vectors = vectors,
                      # The name of our base tree (we use the phylogenomic
                      # tree for this)
                      base_name = "Phylogenomic",
                      # A vector of gene names
                      gene_names = gene_names)

# Let's check out our ordination plot. Each point represents a single tree.
# The phylogenomic tree is shown in red and all other trees are shown in black.

lm_res$plot

# What do you notice from this plot? 
# I notice a few things here. 
# (1) The phylogenomic tree is in the middle of the majority of the gene trees.
# (2) There are two potential outliers in the first principal component. We should
#     check this out to see if we can tease out why!
# (3) There is one potential outlier in the second principal component. We should
#     also check this out! 

# It would be awesome if this plot were interactive and we could mouse over the 
# points to see which points represent which genes. It turns out we can! `plotly` 
# is a great visualization package that lets us do this. We can wrap up our ggplot
# with the command ggplotly and specify the variable that we want to use to label
# each point.

ggplotly(lm_res$plot, tooltip = "name")

# By mousing over our plot, we can see that the two potential outliers in the 
# first principal component are GTP_cyclohydroI and DMRL_synthase, and the potential
# outlier in the second principal component is BacA. 

# ------------------------- Visualizing individual trees -----------------------

# Now that we've visualized all of our trees, we have a few trees that we'd like
# to investigate. Let's start with GTP_cyclohydroI and DMRL_synthase. 

# First we need to figure out which trees in our `gene_trees` object correspond
# with these genes. We'll use the function which to give us the indices. 

which(gene_names == "GTP_cyclohydroI")
which(gene_names == "DMRL_synthase")

# We can start by visualizing GTP_cyclohydroI. 

ggtree(gene_trees[[11]]) + 
  geom_tiplab(size = 2)

# Here we can see that there is a very long branch connecting the genome with 
# accession GCA_00239436. More material about this genome could be found here:
# https://gtdb.ecogenomic.org/genome?gid=GCA_002394365.1 

ggtree(gene_trees[[6]]) + 
  geom_tiplab(size = 2)

# In DMRL_synthase, we also see a long branch to the same genome. This could help
# generate hypotheses about the connection between these two genes and this genome.
# However, this could also be a sign that something unexpected is happening this 
# genome alignment for these two genes. 

# Next, let's consider the second principal component and BacA. 

which(gene_names == "BacA")

ggtree(gene_trees[[42]]) + geom_tiplab(size = 2)

# Here, we see another long branch, but one that separates one large clade of the 
# tree from the rest. This represents a large amount of evolution between these two
# sets of genomes. We aren't the biologists (that's all of you!), but it could be
# worth investigating this gene further to see what is happening between these two
# group of genomes. 

# If you have any questions or want to talk about using this tool for your own 
# metagenomics data, please reach out and send me an email at teichs@uw.edu. 
