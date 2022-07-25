## happi lab - STAMPS 2022 
## Pauline Trinh, David Clausen, Sarah Teichman, Amy Willis 
## Friday 29 July 2022 

# --------------------- Introduction -------------------------

# In this lab, we are going to learn about using `happi` for  modeling gene presence 
# while accounting for differences in genome quality! 

# `happi` uses a hierarchical modeling approach to incorporate genome quality information  
# when conducting gene enrichment testing in pangenomics. 
# For example, if you were interested in understanding whether a particular gene X 
# is more present in *E. coli* genomes from sea otters compared to 
# *E. coli* genomes from narwhals then `happi` is well-suited to help you answer this question! 

# `happi` models the association between 
# a covariate (e.g. narwhal E.coli genomes vs. sea otter E.coli genomes) and gene presence (gene X) where the 
# **covariate** is the **primary predictor** of interest and 
# **gene presence** is the **outcome** of interest. 

# The key difference between `happi` and existing methods for modeling gene presence in pangenomics 
# is that `happi` incorporates information about the quality of each genome in its modeling approach. 


# --------------------- Installing happi -------------------------

# To get started let's install `happi` if you haven't already done so! 
# To install: 
if (!require("remotes", quietly = TRUE))
  install.packages("remotes") # check that remotes is installed

remotes::install_github("statdivlab/happi", build_vignettes = TRUE, dependencies = TRUE) # install happi using remotes and build vignettes

# if you want to check that `happi` is installed on your system you can run: 
"happi" %in% rownames(installed.packages()) # This should return [1] TRUE if happi is installed


# --------------------- Load required packages for this lab -------------------------

# Let's make sure we have all the packages and data that we need to follow along in this lab 
library(happi)
library(tidyverse)
library(ggplot2)

prausnitzii_data <- read_csv("https://raw.githubusercontent.com/statdivlab/stamps2022/main/Thursday-afternoon/labs/happi_lab/stamps_prausnitzii.csv")

# --------------------- Data description -------------------------

# The data that we will be using for this lab focuses on 
# 21 Faecalibacterium prausnitzii metagenome-assembled genomes (MAGs). 
# These 21 MAGs were recovered from shotgun metagenomic sequencing data of stool samples 
# collected from a cohort of farm workers and community controls. 

# We are interested in understanding whether there are functions that are associated with either of these groups 
# Put another way, are there functions that are more enriched in F. prausnitzii genomes 
# belonging to farm workers compared to community controls? 

# Let's look at the dimensions of our data 
dim(prausnitzii_data) # dim 21 x 1356 
# Our data has 21 rows and 1356 columns. 
# The 21 rows correspond to the 21 F. prausnitzii metagenome-assembled genomes (MAGs)
# These genomes were recovered from 11 farm workers and 10 community controls 
# Columns 3-1355 correspond to Clusters of Orthologous Groups (COG) functions that have been 
# identified as present (=1) or absent (=0) in each genome 

head(prausnitzii_data)

# Our covariates are as follows:

#  - `farm_worker`: Indicator of whether the genome was recovered from a farm worker or community control
#  - `group`: Categorical variable designating association with farm worker or community control 
#  - `mean_coverage`: the mean coverage* of each MAG/genome 
# *mean coverage refers to the average number of reads that align to or "cover" a reference area (e.g. genome size or locus size)

# We are going to use the mean coverage of each genome as the genome quality variable we want to use with `happi`. 


# --------------------- Why account for genome quality? -------------------------
# Metagenome-assembled genomes (MAGs) are frequently incomplete or can contain errors 
# (i.e., contamination of fragments from other, missing genes due to assembly issues or shallow sequencing depth). 
# These errors can impact the probability of detecting a gene in our genomes so we want to be able to account 
# for differential genome quality factors and improve our statistical inference! 


# --------------------- Let's use happi! -------------------------
# We are interested in understanding whether the presence of a gene, 
# we'll choose the gene that encodes for L-rhamnose isomerase,
# is more enriched (aka prevalent) in F. prausnitzii MAGs in one group compared to another. 

# Let's subset our data to look at this particular gene for L-rhamnose isomerase
prausnitzii_data %>% 
  dplyr::select(group, farm_worker, mean_coverage,
                `L-rhamnose isomerase`) -> isomerase

# And let's plot this data to see what we're working with! 
isomerase %>% ggplot() +
  geom_jitter(aes(x = mean_coverage, y = `L-rhamnose isomerase`, col = group, pch = group), height=0.08, width=0.00) +
  xlab("Mean coverage") + ylab("") +
  theme_bw() + 
  scale_colour_manual(values= c("mediumseagreen", "dodgerblue")) + 
  theme(legend.position="right") +
  scale_y_continuous(breaks = c(0,1),
                     label = c("Not detected", "Detected"), limits=c(-0.32, 1.1)) 

# So we see that there are more community control-associated MAGs that do not have 
# `L-rhamnose isomerase` detected in their genomes and these appear to be MAGs that have lower mean coverage. 

# Looking at this plot, there may be a potential difference between the presence/absence of 
# `L-rhamnose isomerase` by  site where this gene appears to be more prevalent in 
# farm worker-associated MAGs than in community-associated MAGs. 
# However, are we conflating this difference in gene detection with 
# differences in genome quality (aka mean coverage)? 

# Let's compare an existing method for gene enrichment testing with `happi`. 
# One existing method to test the hypothesis of whether there is a difference in gene presence 
# by some covariate of interest is to use a generalized linear model (GLM) with Rao score test. 

# So let's use that and see what we get. 
ha <- glm(`L-rhamnose isomerase` ~ farm_worker, family="binomial", data = isomerase)
h0 <- glm(`L-rhamnose isomerase` ~ 1, family="binomial", data = isomerase)
anova(ha, h0, test = "Rao")[2, "Pr(>Chi)"]

# Using a GLM + Rao score test we see that the p-value when testing for differences between 
# gene presence in farm worker compared to community control associated F. prausnitzii MAGs is 
# [1] 0.07246124

# Recall though, we were concerned about conflating of this difference in gene detection with 
# differences in mean coverage. We would expect that if this were the case, a method 
# that does not account for genome quality would produce smaller p-values than a method 
# that does account for genome quality. 

# `happi` accounts for genome quality in its modeling of gene presence and allows user flexibility 
# to specify which genome quality variable is relevant to their experimental condition. 
# Let's see how `happi` does! 

x_matrix <- model.matrix(~farm_worker, data = isomerase) # create design matrix using your main predictor of interest

# the main function we'll be using is happi() which has various options that you can specify 
# We'll describe a few here but to see the full list of available options you can run: 

?happi 

# The main options we'll use to run happi() include: 
#  - outcome: Your gene presence/absence variable that needs to be coded as 0 or 1
#  - covariate: the main predictor/covariate of interest formatted as a design matrix (as we did up above with x_matrix)
#  - quality_var: the quality variable of interest. Currently happi only accepts 1 quality variable but we will expand this capability to multiple quality variables 
#  - max_iterations: the maximum number of E-M steps the algorithm will run for. A larger number of steps may mean a longer run time if the algorithm doesn't reach the change_threshold.  
#  - change_threshold: the maximum percent change of the likelihood for 5 iterations in a row for both the alternative and null models that would signal termination of the algorithm 
#  - epsilon: the fixed probability of observing a gene given that it shouldn't be present (aka contamination of the genome from other genomes)
#  - method: method for estimating f. Defaults to "splines" which fits a monotone spline with df determined by argument spline_df; "isotone" for isotonic regression fit
#  - firth: uses Firth penalty (default is TRUE)
#  - spline_df: degrees of freedom (in addition to intercept) to use in monotone spline fit (default is 3)

happi_results <- happi(outcome=isomerase$`L-rhamnose isomerase`, 
                       covariate=x_matrix, 
                       quality_var=isomerase$mean_coverage,
                       max_iterations=100, 
                       change_threshold=0.001, 
                       epsilon=0, 
                       method = "splines", 
                       firth = T, 
                       spline_df = 4)
# Let's look  at our results! 
# ** note that the p-value results are stored in a vector (happi_results$loglik$pvalue) that is the length of 
# the number of max_iterations you specified. If the algorithm reaches the change_threshold requirements 
# before the max_iterations, this will result in a lot of NAs at the end because the algorithm completes earlier. 
# The final p-value result is the last element of happi_results$loglik$pvalue before the trailing NAs. 
# In this particular example that would be element 42 as can be seen using happi_results$loglik$pvalue[42] or...

happi_results$loglik$pvalue %>% tibble() %>% filter(!is.na(.)) %>% tail(1)

# We see that the p-value we get from `happi` is p = 0.180 and that this p-value is 
# larger than the p-value (p = 0.07) we got using GLM + Rao when we didn't account for mean coverage. 

# We think this is great! In  situations where  the pattern of detection or non-detection 
# could be attributable to genome quality, we think statements about significance should be moderated. 
# In this case we want to see larger p-values! 

# If we want to get the beta estimates from happi we can run the following: 

happi_results$beta %>% tibble() %>% drop_na() %>% tail(1)

#.[,1]  [,2]
# <dbl> <dbl>
# -0.699  2.98
# and we see that our estimates are -0.699 for our intercept beta_0
# and our estimate for beta_1 which corresponds to our main predictor of interest is 2.98 
# Based on our results, we see that the gene encoding for L-rhamnose isomerase is more enriched in 
# F. prausnitzii MAGs found in farm workers but that this difference is not significant at the 5% significance level (p = 0.180). 


# --------------------- But what if I have thousands of genes I want to look at? -------------------------

# You can parallelize your analyses! To do that you need to make sure you have the package `parallel` installed. 

if (!require("parallel", quietly = TRUE))
  install.packages("parallel") # check that parallel is installed. If not then install. 

library(parallel) # load parallel package 

# to run your analyses in parallel let's first make our design matrix x_matrix_prausnitzii
# and also a function run_happi_prausnitzii() that contains all the specifications of happi() we want to use 

x_matrix_prausnitzii <- model.matrix(~farm_worker, data = prausnitzii_data)

run_happi_prausnitzii <- function(colnum) {
  happi(outcome=unlist(prausnitzii_data[,colnum]), 
        covariate=x_matrix_prausnitzii, 
        quality_var=prausnitzii_data$mean_coverage,
        method="splines", 
        firth=T, 
        spline_df=4,
        max_iterations=100, 
        change_threshold=0.01, 
        epsilon=0)
}

# For the purposes of this lab we will not be running all the COG functions and will focus on only 50. 
# Additionally for quicker run-time we're going to keep the max_iterations low at 100 
# and a higher change threshold of 0.01 

# To run happi in parallel: 
prausnitzii_results <- mclapply(3:52, run_happi_prausnitzii, mc.cores=6) # this should only take about a minute or two to finish running... 
# 3:52 denotes the COG functions I want to run from columns 3 - 52 
# run_happi_prausnitzii is the function we created above 
# and mc.cores allows me to specify how many cores I want to allocate to this computational task 

# we store our results in prausnitzii_results and can consolidate our beta estimates along with p-values using

pvalue_prausnitzii <- lapply(prausnitzii_results, function(x) tail(x$loglik$pvalue[!is.na(x$loglik$pvalue)], 1)) %>% unlist
beta_prausnitzii <- lapply(prausnitzii_results, function(x) tail(x$beta[!is.na(x$beta[,1]),],1)) %>% do.call("rbind",.)

prausnitzii_hyp_results <- tibble("gene" = colnames(prausnitzii_data)[3:52], # grab the names of each gene from prausnitzii_data
                                  pvalue_prausnitzii, # Combine with our pvalues and betas
                                  beta_prausnitzii[,1],
                                  beta_prausnitzii[,2]) %>% 
  arrange(pvalue_prausnitzii)

head(prausnitzii_hyp_results)

# And here we see the variables
#  - gene: contains the COG function names
#  - pvalue_prausnitzii: the happi p-values 
#  - beta_prausnitzii[, 1]: the estimate for beta0 that corresponds to our intercept 
#  - beta_prausnitzii[, 2]: the estimate for beta1 that corresponds to our primary covariate (farm_worker)

### Remember to choose your favorite FDR or FWER method to address multiple comparisons! 


# --------------------- Final Notes -------------------------
# - `happi` can be extended for use when testing for gene presence between metagenomes of varying quality (e.g., sequencing depth)
# - `happi` currently can only take in 1 quality variable but we are prioritizing incorporation of multiple quality variables
# - If you have additional questions or issues with using `happi` please open an issue on our GitHub https://github.com/statdivlab/happi

# Please cite our work if you use `happi`! 
# Pauline Trinh, David S. Clausen, and Amy D. Willis. happi: a hierarchical approach to pangenomics inference. bioRxiv e-prints, April 2022.











