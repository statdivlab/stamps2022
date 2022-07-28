### This lab examines a sequencing experiment performed by Costea et al.
### (2017) [1] and subsequently analyzed by McLaren et al. (2019) [2]. Briefly,
### the data we'll look at in this lab consists of output
### from shotgun sequencing of fecal specimens collected from 10 unique
### study participants. Each specimen was split into three samples
### which were sequenced according to one of three protocols (labeled
### H,Q, and W).

### Before each sample was sequenced, a mock
### community containing 10 bacterial species was spiked into it.
### In addition to sequencing data, Costea et al. (2017) published
### flow cytometry measurements taken on each of the bacterial isolates
### combined to create this mock community. Our focus in this lab is how
### shotgun-sequencing-based estimates of relative abundance within this
### mock community compare to (presumably more accurate) estimates based on
### flow cytometry.

# [1] Costea, Paul I., et al.
# "Towards standards for human fecal sample processing in metagenomic studies."
# Nature biotechnology 35.11 (2017): 1069-1076.

# [2] McLaren, Michael R., Amy D. Willis, and Benjamin J. Callahan.
# "Consistent and correctable bias in metagenomic sequencing experiments."
# Elife 8 (2019): e46923.

# Alright let's get started!
### Let's load libraries we'll need

library(tidyverse)

# Because we're working on the RStudio Server, we already have `tinyvamp` installed. However,
# if you're working at home, uncomment the follow lines and run them in your console (not in
# this script).

# To install:
#if (!require("remotes", quietly = TRUE))
#  install.packages("remotes") # check that remotes is installed

# Install tinyvamp using remotes and build vignettes:
#remotes::install_github("https://github.com/statdivlab/tinyvamp")

library(tinyvamp)

# ### And load data as well
# read in metaphlan2 profiles from McLaren et al.
costea2017_metaphlan2_profiles <-
  read_csv("https://raw.githubusercontent.com/statdivlab/stamps2022/main/Thursday-afternoon/labs/tinyvamp/costea2017_metaphlan2_profiles.csv")

# read in flow cytometry data from Costea et al.
costea2017_mock_composition <-
  read_csv("https://raw.githubusercontent.com/statdivlab/stamps2022/main/Thursday-afternoon/labs/tinyvamp/costea2017_mock_composition.csv")
# read in sample metadata
costea2017_sample_data <-
  read_csv("https://raw.githubusercontent.com/statdivlab/stamps2022/main/Thursday-afternoon/labs/tinyvamp/costea2017_sample_data.csv")



### Let's take a look at 'costea2017_metaphlan2_profiles'
head(costea2017_metaphlan2_profiles)
head(costea2017_metaphlan2_profiles$Clade)

### It looks like the 'Clade' column includes metaPhlan2 output for
### all levels of taxonomy between kingdom and species
### so let's...

################# Filter data (species level, mock taxa only) ##################
# first let's store the names of taxa in the spike-in mock
mock_taxa <- costea2017_mock_composition$Taxon

# now we'll filter out observations that 1) aren't for those taxa or
# 2) aren't at species level
costea2017_metaphlan2_profiles_species <- costea2017_metaphlan2_profiles %>%
  filter(!str_detect(Clade, "t__")) %>%
  filter(str_detect(Clade, paste(mock_taxa, collapse="|"))) %>%
  mutate(Clade = str_remove(Clade, ".*s__"))

### now let's take a look at "costea2017_mock_composition"
### which contains flow cytometry on isolates used to
### create the mock communities we have shotgun sequencing output for
head(costea2017_mock_composition)

### let's pull out "Taxon" and create a new matrix containing
### only flow cytometry measurements
mock_cols <- costea2017_mock_composition$Taxon

mock_mat <- t(as.matrix(costea2017_mock_composition[,-c(1:2)]))
colnames(mock_mat) <- mock_cols

### tinyvamp is expecting an observation matrix whose rows are samples
### and whose columns are taxa... so let's make that:

W <- costea2017_mock_composition %>%
  as.data.frame() %>%
  dplyr::rename("cells_per_ml" = "bacterial cells/ml OD 1") %>% ## rename cell conc
                                                            ## column
  group_by(Taxon) %>%
  summarize(cells_per_ml = mean(cells_per_ml)) %>% ## we have two flow cytometry
                                                   ## measurements for *most*
                                                   ## isolates; combine these
                                                   ## by taking an average
  (function(x){ m <- matrix(x$cells_per_ml,nrow = 1)
  colnames(m) <- x$Taxon
  return(m)})

# let's figure out how to order our MetaPhlan2 output
# so that its columns are in the same order as W's
W_reorder <- sapply(colnames(W), # if you're wondering if there's a simpler way
                    function(x) which( # to do this, yes - there probably is
                      sapply(costea2017_metaphlan2_profiles_species$Clade,
                             function(d) grepl(x,d,fixed = TRUE)))) %>%
  as.numeric()

# and now we'll add that output to W
W <- rbind(W,
           costea2017_metaphlan2_profiles_species[W_reorder,-c(1:2)] %>%
             as.matrix() %>%
             t)


# now we organize sequencing protocol information
protocol_df <- W %>%
  as_tibble(rownames="Run_accession") %>%
  full_join(costea2017_sample_data) %>%
  select(Run_accession, Protocol)

protocol_df

#give "Run_accession" column for the flow cytometry measurements value
# "Flow Cytometry" and also call  the corresponding protocol "Flow Cytometry"

protocol_df$Run_accession[1] <- "Flow Cytometry"
protocol_df$Protocol[1] <- "Flow Cytometry"
#make sure first row of W (contains flow cytometry measurements)
#is named in the same way
rownames(W)[1] <- "Flow Cytometry"

# how many protocols are there?
unique(protocol_df$Protocol)

# we're about to inner_join with costea2017_sample_data,
# but we'll lose our flow cytometry measurements if we don't
# add a row to this tibble corresponding to these measurements
costea2017_sample_data <-
  rbind(tibble("...1" = 0,
               "Sample" = "Flow Cytometry",
               "Protocol" = "Flow Cytometry",
               "Individual" = "none",
               "SI_sample" = "irrelevant",
               "Run_accession" = "Flow Cytometry"),
        costea2017_sample_data)

#and inner_join to relate samples to the study participants ("Individual" column)
#who provided them
protocol_df <-
  protocol_df %>%
  inner_join(costea2017_sample_data)

# and put these things together in a dataframe called "measurements with metadata"
measurements_with_metadata <- W %>%
  as.data.frame() %>%
  rownames_to_column(var = "Run_accession") %>%
  pivot_longer(-Run_accession) %>%
  group_by(Run_accession) %>%
  inner_join(x = ., protocol_df)

### Let's see how empirical proportions from metaPhlan2 compare to
### flow cytometry measurements
measurements_with_metadata %>%
  group_by(Run_accession) %>%
  mutate(empirical_proportion = value/sum(value)) %>%
  ungroup %>%
  ggplot() +
  geom_bar(aes(x = Run_accession, y = empirical_proportion,
                group = Run_accession, fill = name),
           position="stack", stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(~Protocol, scales = "free_x") +
  scale_fill_viridis_d() +
  ylab("Empirical Proportion") +
  xlab("Sample") +
  labs(fill = "Species")

### What seems to be over-represented? Under-represented?

### Another view of this
measurements_with_metadata %>%
  mutate(flow_cytometry = measurements_with_metadata$value[
    measurements_with_metadata$Protocol=="Flow Cytometry"
  ]) %>%
  mutate(flow_cytometry = flow_cytometry/sum(flow_cytometry)) %>%
  filter(Protocol != "Flow Cytometry") %>%
  group_by(Run_accession) %>%
  mutate(empirical_proportion = value/sum(value)) %>%
  ungroup %>%
  ggplot() +
  geom_line(aes(x = name, y = empirical_proportion,
                group = Run_accession, color = Protocol)) +
  geom_line(aes(x = name, y = flow_cytometry,
                group = Run_accession),linetype = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(Protocol~.) +
  ylab("Empirical Relative Abundances \n(Flow Cytometry Shown as Dotted Line)") +
  xlab("Species") +
  scale_y_log10()
#Flow-cytometry-based relative abundance estimates are shown as black dotted lines
#in the plot above.
#How well do you think the empirical relative abundances within the mock community
#that we measured under each protocol reflect the relative abundances as measured
#by flow cytometry? Are we close to being in agreement with flow cytometry?

### Note that each protocol gives fairly consistent results on
### this mock community -- within-protocol technical variation
### is low! But none of them accurately reflect our flow cytometry
### data (which we are taking as a standard)

### A brief digression... let's talk a bit about bias and variance.
### This is on the relative abundance scale -- scale *really* matters when we're
### talking about bias and variance. We would have different results if we
### log or log-ratio transformed our observations (in part because we would be
### estimating a different parameter).

### Calculate and store empirical proportions
measurements_for_bias_figure <-
  measurements_with_metadata %>%
  mutate(flow_cytometry = measurements_with_metadata$value[
    measurements_with_metadata$Protocol=="Flow Cytometry"
  ]) %>%
  group_by(Run_accession) %>%
  mutate(empirical_proportion = value/sum(value)) %>%
  ungroup

### Transform flow_cytometry column to relative abundance scale
measurements_for_bias_figure <-
  measurements_for_bias_figure %>%
  group_by(Run_accession) %>%
  mutate(flow_cytometry = flow_cytometry/sum(flow_cytometry))


### Now let's just look at within-protocol technical variation at
### the relative abundance / empirical proportion scale:
measurements_for_bias_figure %>%
  filter(Protocol !="Flow Cytometry") %>% #we're comparing shotgun sequencing
                                           #protocols to flow cytometry,
                                           #so cut out flow cytometry rows
                                           #(we added flow cytometry data as
                                           #a column in this tibble)
  group_by(Protocol,name) %>%
  summarize(within_protocol_variance = var(empirical_proportion),
            bias_in_rel_abundance_estimator = mean(empirical_proportion) -
              flow_cytometry,
            flow_cytometry = flow_cytometry) %>%
  mutate(relative_bias = bias_in_rel_abundance_estimator/flow_cytometry) %>%
  ggplot() +
  geom_line(aes(x = name, y = within_protocol_variance,
                group = Protocol,
                color= Protocol)) +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Within-protocol Variance of Estimated Relative Abundances") +
  xlab("Species")

# ~~~~****Technical variation is pretty low!****~~~~

#How does bias look for empirical proportions (relative to flow cytometry)?

measurements_for_bias_figure %>%
  filter(Protocol !="Flow Cytometry") %>% #we're comparing shotgun sequencing
  #protocols to flow cytometry,
  #so cut out flow cytometry rows
  #(we added flow cytometry data as
  #a column in this tibble)
  group_by(Protocol,name) %>%
  summarize(within_protocol_variance = var(empirical_proportion),
            bias_in_rel_abundance_estimator = mean(empirical_proportion) -
              flow_cytometry,
            flow_cytometry = flow_cytometry) %>%
  mutate(relative_bias = bias_in_rel_abundance_estimator/flow_cytometry) %>%
  ggplot() +
  geom_line(aes(x = name, y =  bias_in_rel_abundance_estimator,
                group = Protocol,
                color= Protocol)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Empirical Bias")


#not looking quite so hot (esp. Prevotella), but...
#what if we look at *relative* bias (i.e., the bias
#divided by the flow cytometry measurement -- so
#how large is bias as a proportion of the "true"/flow-cytometry value

measurements_for_bias_figure %>%
  filter(Protocol !="Flow Cytometry") %>% #we're comparing shotgun sequencing
  #protocols to flow cytometry,
  #so cut out flow cytometry rows
  #(we added flow cytometry data as
  #a column in this tibble)
  group_by(Protocol,name) %>%
  summarize(within_protocol_variance = var(empirical_proportion),
            bias_in_rel_abundance_estimator = mean(empirical_proportion) -
              flow_cytometry,
            flow_cytometry = flow_cytometry) %>%
  mutate(relative_bias = bias_in_rel_abundance_estimator/flow_cytometry) %>%
  ggplot() +
  geom_line(aes(x = name, y =  relative_bias,
                group = Protocol,
                color= Protocol)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Empirical Relative Bias")

# So our bias is sometimes ~ 6x larger than the true value of the proportion
# we're trying to estimate

#We can also just directly look at the ratio of empirical proportions
#to flow cytometry measurements to get another quantification of how far
#we're off by
measurements_for_bias_figure %>%
  filter(Protocol !="Flow Cytometry") %>%
  ggplot() +
  geom_line(aes(x = name, y =  empirical_proportion/flow_cytometry,
                group = Run_accession,
                color= Protocol)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Ratio: metaPhlan2 proportions over flow cytometry proportions") +
  geom_abline(aes(slope = 0, intercept = 0),
              linetype = 2) +
  scale_y_log10()

# so one thing this is telling us is that it's not unusual for our
# empirical relative abundances (using these shotgun sequencing protocols)
# to be about ***an order of magnitude*** off (e.g., B hansenii,
# Y. pseudotuberculosis), and sometimes we are even further off (e.g.,
# L plantarum under protocol H, which gives us an estimated relative abundance
# that is off by about ***100 fold***).


### Using tinyvamp to estimate "catchability" across taxa and protocols
# I'm going to assert without a huge amount of proof that the phenomenon
# we are observing in the Costea data above is differing levels of "catchability"
# in differing taxa -- i.e., that under each protocol, some microbes are easier
# to detect than others. This results in multiplicative over- and under-detection
# of differing taxa under each protocol. We'll refer to this over/under-detection
# as a "detection effect."
#
# The R package tinyvamp allows us to estimate detection effects and estimate
# relative abundances after 'removing' these effects.

# ... but tinyvamp does not currently accept formulas as input (I know, boo)
# so we'll need to use the model.matrix() function to create our design
# matrices

# First we'll create a matrix telling tinyvamp what detection effects to estimate

X <- model.matrix(~-1 + Protocol,
             data = protocol_df)
head(X)
# We're going to use flow cytometry data as a standard, so we *don't* want
# to estimate an effect for it. To accomplish this, we'll remove the first
# column of X:
X <- X[,-1]

# tinyvamp also needs us to provide a matrix Z explaining how the metagenomic
# measurements we feed to it (the matrix W in this case) are related to unique
# specimens whose relative abundance we want to estimate
# For the moment, we're treating all of our observations as coming from the
# same specimen (since we're focusing on contents of a spike-in control),
# so Z is just a single column of ones
Z <- matrix(1,nrow = 30, ncol = 1)
# This Z will tell tinyvamp that samples 1, ..., 30 came from a single specimen
# To model measurements on multiple specimens, we would add a column in Z for
# each additional specimen

# Z has a cousin named Z_tilde -- this is relevant to estimating contamination,
# which we're not doing here since we're looking just at taxa truly present
# in a mock community, so we'll set all elements of this matrix equal to zero
Z_tilde <- matrix(0,nrow= 30, ncol = 1)
Z_tilde_gamma_cols <- 1 #this is another argument we need to specify -- it
                        #doesn't matter at all to our model fit here though

# we also need to give tinyvamp an initial estimate of sampling
# "intensities" gamma (similar to log read depths -- here our shotgun
# observations are on the relative abundance scale, so read depth isn't really
# meaningful, but we still need to estimate this parameter to formally fit
# our model)
gammas <- apply(W,1,function(x) log(sum(x)))
# we also tell tinyvamp to estimate all gammas and not hold any of them fixed
gammas_fixed_indices <- rep(FALSE, length(gammas))

# We also want to estimate a relative abundance profile of our mock community
# this will be represented by P, so we give this P initial values:
P <- matrix(W[1,]/sum(W[1,]),nrow =1, ncol = 10)
# and tell tinyvamp not to hold any elements of P fixed (i.e., estimate P)
P_fixed_indices <- matrix(FALSE,nrow = 1, ncol = 10)

# this is another input related to contamination -- we can ignore it here
X_tilde <- matrix(0,ncol = 3, nrow= 1)

# Now we tell tinyvamp what detection effects (over/under-representation
# of different taxa in different protocols) we want it to estimate
# We do this with matrix B, which we'll give starting value 0
B <- matrix(0,
            ncol = 10,
            nrow = 3)
# We actually *do* need to set some elements of B to be fixed
# This is essentially because we can only estimate how
# under- or over-represented each taxon is relative to
# some other taxon. We set the last (10th) column of B
# to be fixed at zero, which means we are estimating
# degree of over/under-representation of each other taxon
# relative to our 10th taxon, Y. pseudotuberculosis
B_fixed_indices <- matrix(FALSE,ncol = 10,nrow = 3)
B_fixed_indices[,10] <- TRUE

# It doesn't actually matter what taxon we set to
# be fixed -- we will get the same between-taxon
# relationships in detection regardless of what
# choice we make

# More contamination parameters to ignore!
P_tilde <- P*0
P_tilde_fixed_indices <- !P_fixed_indices
gamma_tilde <- matrix(0,ncol = 1, nrow = 1)
gamma_tilde_fixed_indices <- TRUE

# (For the record, we are working on the tinyvamp interface to
# make it slightly less... like this)

# Ok we can fit a model now!
full_model  <- estimate_parameters(W = W,
                      X = X,
                      Z = Z,
                      Z_tilde = Z_tilde,
                      Z_tilde_gamma_cols = 1,
                      gammas = gammas,
                      gammas_fixed_indices = gammas_fixed_indices,
                      P = P,
                      P_fixed_indices = P_fixed_indices,
                      B = B,
                      B_fixed_indices = B_fixed_indices,
                      X_tilde = X_tilde,
                      P_tilde = P_tilde,
                      P_tilde_fixed_indices = P_tilde_fixed_indices,
                      gamma_tilde = gamma_tilde,
                      gamma_tilde_fixed_indices = gamma_tilde_fixed_indices,
                      alpha_tilde = NULL,
                      Z_tilde_list = NULL,
                      barrier_t = 1, #starting value of reciprocal barrier penalty coef.
                      barrier_scale = 10, #increments for value of barrier penalty
                      max_barrier = 1e12, #maximum value of barrier_t
                      initial_conv_tol = 1000,
                      final_conv_tol = 0.1,
                      # final_f = 1e-6,
                      constraint_tolerance = 1e-10,
                      hessian_regularization = 0.01,
                      criterion = "Poisson",
                      # subproblem_method = "Newton",
                      profile_P = FALSE,
                      profiling_maxit = 25,
                      wts = NULL,
                      verbose = FALSE)

# Ok we fit our model -- so what does it say?
# Let's look at the estimated detection effects B
full_model$B

# Cool there are numbers. Do they mean something?
# To clarify, let's give B informative column and row names
rownames(full_model$B) <- colnames(X) #label effects by protocol
colnames(full_model$B) <- colnames(W) #label effects by taxon

# Now we can exponentiate B to get more interpretable results:
round(exp(full_model$B),2)
# We interpret the estimate for B. hansenii to indicate that
# under protocol H, if we sequence a specimen consisting of equal
# parts Y. pseudotuberculosis and B. hansenii (as measured by flow cytometry),
# we will on average
# observe 0.2 B. hansenii reads for each Y. pseudotuberculosis read
# (or, more or less equivalently, we can say that the ratio of 1) the MetaPhlan2
# estimate of relative abundance in B. hansenii to 2) the estimate of
# relative abundance in Y. pseudotuberculosis will on average be
# 1/5 as large as the ratio of true relative abundances)

# On the other hand, we estimate that under protocol H, we will observe
# about 65 P melanogenica reads for each Y. pseudotuberculosis read
# if we sequence an even mixture of these two species
# Does this align with what you see in the plots we produced earlier in
# this lab?


# Does estimating these detection effects help us estimate relative abundances?
# Briefly, yes, but we can examine this with a cross-validation. Specifically,
# we conduct a 10-fold cross-validation, and for each fold we hold out all
# samples from a single participant's specimen.* We fit our model on the
# remainder of the samples and predict on the held-out samples.

#*strictly, this is for the first 8 folds. Then we have a fold consisting of
#measurements on a pure mock sample as the 9th fold and a 10th fold consisting of
#measurements on samples "A" and "B" (human fecal samples that were sequenced
#only under protocol Q).
#)

# We'll do the cross-validation now:
# construct folds
folds <- vector(10, mode = "list")
available <- 2:30
unique_individuals <- c(1:8,"M")

for(i in 1:9){
  folds[[i]] <- which(protocol_df$Individual == unique_individuals[i]) +1
}

folds[[10]] <- which(protocol_df$Individual %in% c("A","B")) +1

full_cv <- vector(10, mode = "list")
for(whichfoldout in 1:10){
  print(paste("Fitting fold", whichfoldout))
  heldout <- folds[[whichfoldout]]
  nheldout <- length(heldout)
  Z_cv <- Z
  Z_cv <- cbind(Z_cv - Z_cv*(1:30 %in% heldout))
  P_cv <- P
  P_fixed_indices_cv <- P_fixed_indices
  for(k in 1:nheldout){
    Z_cv <- cbind(Z_cv,as.numeric(1:30 == heldout[k]))
    P_cv <- rbind(P_cv,P)
    P_fixed_indices_cv <- rbind(P_fixed_indices_cv,
                                P_fixed_indices)
  }

  full_cv[[whichfoldout]]  <-
    estimate_parameters(W = W,
                        X = X,
                        Z = Z_cv,
                        Z_tilde = Z_tilde,
                        Z_tilde_gamma_cols = 1,
                        gammas = gammas,
                        gammas_fixed_indices = gammas_fixed_indices,
                        P = P_cv,
                        P_fixed_indices = P_fixed_indices_cv,
                        B = B,
                        B_fixed_indices = B_fixed_indices,
                        X_tilde = X_tilde,
                        P_tilde = P_tilde,
                        P_tilde_fixed_indices = P_tilde_fixed_indices,
                        gamma_tilde = gamma_tilde,
                        gamma_tilde_fixed_indices = gamma_tilde_fixed_indices,
                        alpha_tilde = NULL,
                        Z_tilde_list = NULL,
                        barrier_t = 1, #starting value of reciprocal barrier penalty coef.
                        barrier_scale = 10, #increments for value of barrier penalty
                        max_barrier = 1e12, #maximum value of barrier_t
                        initial_conv_tol = 1000,
                        final_conv_tol = 0.1,
                        # final_f = 1e-6,
                        constraint_tolerance = 1e-10,
                        hessian_regularization = 0.01,
                        criterion = "Poisson",
                        # subproblem_method = "Newton",
                        profile_P = FALSE,
                        profiling_maxit = 25,
                        wts = NULL,
                        verbose = FALSE)


}
# Extract predictions on held-out folds
full_cv_predictions <- lapply(1:10,
                              function(x)
                                full_cv[[x]]$varying[
                                  full_cv[[x]]$varying$param == "P"&
                                    full_cv[[x]]$varying$k>1,])

for(i in 1:10){
  full_cv_predictions[[i]]$k <- sapply(full_cv_predictions[[i]]$k,
                                       function(x) folds[[i]][x - 1])

}

full_cv_predictions <- do.call(rbind,full_cv_predictions)
# generate flow cytometry estimates of relative abundance
fc_values <- W[1,]/sum(W[1,])
full_cv_predictions$fc_value <- sapply(full_cv_predictions$j,
                                       function(d) fc_values[d])
full_cv_predictions$protocol <-
  sapply(full_cv_predictions$k,
         function(d) protocol_df$Protocol[d-1]) #  d - 1 bc k starts at 2
#  (k = 1 is fc data)

# label predictions by specimen predicted on
full_cv_predictions$specimen <-
  sapply(full_cv_predictions$k,
         function(d) protocol_df$Individual[d - 1]) #  d - 1 bc k starts at 2
#  (k = 1 is fc data)

# label these predictions as being from tinyvamp
full_cv_predictions$model <- "tinyvamp"

# generate naive / plug-in estimates of mock community relative abundances
# on the basis of MetaPhlan2 output under each of the sequencing protocols
W_prop <- W[-1,] #don't use flow cytometry measurements (1st row of W)
for(i in 1:nrow(W_prop)){
  W_prop[i,] <- W_prop[i,]/sum(W_prop[i,])
}

#collect naive predictions
naive_predictions <- full_cv_predictions[numeric(0),]

for(i in 1:nrow(W_prop)){
  protocol <- protocol_df$Protocol[i]
  specimen <- protocol_df$Individual[i]
  for(j in 1:ncol(W)){
    naive_predictions <- rbind(naive_predictions,
                               data.frame(value = W_prop[i,j],
                                          param = "P",
                                          k = i,
                                          j  = j,
                                          fc_value = fc_values[j],
                                          protocol = protocol,
                                          specimen = specimen,
                                          model = "Plug-in"))
  }
}


#plot!

rbind(full_cv_predictions,
      naive_predictions) %>%
  filter(protocol != "Flow Cytometry") %>%
  mutate(protocol = sapply(protocol, function(x) paste("Protocol ",x,
                                                       sep = "",
                                                       collapse = ""))) %>%
  # filter(!is.na(protocol)) %>% #why would protocol be NA? Check!
  ggplot() +
  geom_point(aes(x = fc_value, y = value, color= specimen),
             size = .5) +
  geom_line(aes(x = fc_value, y = value, color = specimen,
                group = as.factor(k)), size = .5) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 2) +
  facet_grid(model~protocol,scales = "free_y") +
  scale_color_viridis_d() +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() +
  guides(color=guide_legend(title="Specimen")) +
  xlab("Relative Abundance Measured by Flow Cytometry") +
  ylab("Cross-Validated Estimated Relative Abundance")

# The plot above provides flow cytometry estimates of relative abundance
# on the x-axis and either plug-in (top row) or cross-validated tinyvamp
# (bottom row) estimates of relative abundance on the y-axis. Sequencing
# protocols are split across columns, and all relative abundance estimates
# are color-coded according to the specimen on which they were taken.
# The line y = x is shown as a dotted line above.

# What are your takeaways? How do the protocols compare in terms of bias
# and variance of plug-in relative abundance based on measurements taken
# according to them? How reliable are the plug-in estimates of relative
# abundance?

##### Thoughts about experimental design, analysis, and interpretation #####

# In general, it's risky territory to interpret read proportions / metaphlan
# output / etc. as if they were true relative abundances (or even close to true)
# in a sequenced sample. We only looked at one community (a mock) in this lab,
# but it's worth considering as you look at read proportions from an experiment
# that
      # ***detection effects impact plug-in relative abundance
      # estimates in different communities differently!***
# (For example, consider we sequence two samples, an even mixture of
# species A and B and an even mixture of species A and C. What happens to
# plug-in relative abundances estimates of species A if species B is much
# easier to detect than species C?)

# As a result, any analysis that relies on observed read proportions to be
# similar to true relative abundance profiles may be misleading... but it's
# hard to predict exactly when this will happen. Very large differences between
# communities will likely still be visible... but detection effects can also
# produce spurious differences between communities, so again care is
# warranted.

# If you are able to sequence positive controls with your biological samples,
# this is always a good idea. Better yet if you have the wherewithal to
# include mocks containing taxa you are particularly interested in
# characterizing in your samples of interest.

# In a similar vein, since detection effects can vary from batch to batch,
# if you have to sequence samples in multiple batches,
# it's generally best to make sure each batch contains a similar
# array of samples (so don't, for example, sequence all samples from
# one environment in batch A and all samples from another in batch B).
# But!! This doesn't guarantee that comparisons of observed relative abundances
# won't be misleading.

# As well, some analytical methods are more robust to detection effects
# than others. In general, "compositional data" (a misnomer imho, but I
# digress) methods are fairly robust to detection effects. At the same time,
# these methods generally make fairly restrictive assumptions on taxon
# presence (usually all taxa are assumed to be truly present in every sample,
# which may not be realistic in a variety of settings --
# in other cases, such as ancomII, it is assumed that we can identify which
# taxa are present in which sample... which also strikes me as a dubious
# proposition). radEmu was specifically designed to be robust to detection
# effects and does not assume all taxa are present in every sample... so yeah.
# With all of this said, none of these methods have been demonstrated to be
# robust to measurement error that includes detection effects *and* other
# types of error (like contamination or misclassification of taxa), so
# it is in general important to consider how measurement error might impact
# results of a microbiome data analysis.


