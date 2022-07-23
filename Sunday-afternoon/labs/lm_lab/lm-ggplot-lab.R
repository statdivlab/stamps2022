# Linear regression and ggplot lab 
## prepared by Amy Willis & Sarah Teichman

# In this lab, we'll practice running linear regression with the function `lm`. 
# We'll also manipulate our data with the `tidyverse` plot our data with `ggplot`. 

# Let's begin by loading in all the packages we'll need for this tutorial.

library(tidyverse)
library(lme4)

# ------------------------ loading the data -----------------------

# Next, let's load our data. We'll start with two data sets. This data comes from
# a study of Cystic Fibrosis. It consists of sputum and saliva samples from ten
# people with Cystic Fibrosis. More information can be found here: 
# https://journals.asm.org/doi/full/10.1128/mSystems.00296-20?rfr_dat=cr_pub++0pubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&

ddpcr <- read_csv("https://raw.githubusercontent.com/statdivlab/stamps2022/main/Sunday-afternoon/labs/lm_lab/data/ddPCR.csv")
meta <- read_csv("https://raw.githubusercontent.com/statdivlab/stamps2022/main/Sunday-afternoon/labs/lm_lab/data/meta_data.csv")

# ------------------------ exploring the data ---------------------

# We can take a look at the beginning of each dataset with the command head

head(ddpcr)
head(meta)

# We'd like to combine these two data sets together. They both contain the variable 
# `sample_name`, so we will match observations based on this variable. 

# Note that when we use `tidyverse` syntax, the pipe symbol %>% is our friend. The
# pipe lets us chain commands together, using the output of the previous command as
# the input to the next command. 

both <- meta %>% 
  inner_join(ddpcr, by = "sample_name") 

# Next, we will manipulate our data to make it easier to work with. We'll rename the
# variable `Average Replicates, stock DNA (copies /µl)` (remember, we like variable
# names to be concise and avoid spaces) and tell R that this variable is a number. 
# We'll do the same for FEV1.

both <- both %>%
  rename(ddpcr = `Average Replicates, stock DNA (copies /µl)`) %>%
  mutate(ddpcr = gsub(",", "", ddpcr)) %>%
  mutate(ddpcr = as.numeric(ddpcr)) %>%
  mutate(FEV1 = as.numeric(FEV1))

# Earlier we looked at the first few observations of each dataset. To see a full list
# of variables in our dataset, we can use the `names()` function.

names(both)

# We'll be focusing on a few variables from this dataset: ddpcr, FEV1, Treatment 
# Group, Sample Type, and Subject ID. 

# ddprc: digital PCR, the total bacterial load in the sample
# FEV1: Forced expiratory volume (FEV) measures how much air a person can exhale 
# during a forced breath
# Treatment Group: "ON" if the sample is from a subject treated with antibiotics, 
# "OFF" otherwise 
# Sample Type: the location of the sample, either Saliva or Sputum
# Subject ID: an ID identifying the participant 

# Before we do inference, we like to plot our data. This let's us get a sense of 
# general trends and build intuition about our data. 

# We'll use ggplot for plotting. Remember that in ggplot we start with the base
# of our plot, where we specify our dataset and variables, and then we can build
# up our plot by adding layers with the `+` operator. 

ggplot(data = both, aes(x = `Sample Type`, y = ddpcr, col = `Subject ID`)) +
  geom_jitter(height = 0, width = 0.3) + 
  labs(y = "Average observed ddPCR\n(copies/ul)",
       title = "Digital PCR by Sample Type") + 
  scale_y_log10() + 
  theme(plot.title = element_text(hjust = 0.5))

# What do we learn from this? 
# It seems that there are sputum samples generally have higher pcr values
# than saliva samples. 

ggplot(data = both, aes(x = `Treatment Group`, y = ddpcr, col = `Subject ID`)) +
  geom_jitter(height = 0, width = 0.3) + 
  labs(y = "Average observed ddPCR\n(copies/ul)",
       title = "Digital PCR by Treatment") + 
  scale_y_log10() + 
  theme(plot.title = element_text(hjust = 0.5))

# What do we learn from this? 
# There are a lot more non-treated samples than treated samples. The general trend
# seems to be higher pcr values from non-treated samples than treated samples.

ggplot(data = both, aes(x = as.numeric(FEV1), y = ddpcr, col = `Subject ID`)) +
  geom_point() + 
  labs(x = "FEV Value",
       y = "Average observed ddPCR\n(copies/ul)",
       title = "Digital PCR by FEV Value") + 
  scale_y_log10()
  theme(plot.title = element_text(hjust = 0.5))

# What do we learn from this? 

# There might be a slight negative correlation between PCR values and FEV values
# There seem to be two observations for each value on the x-axis. Why might that be?

# It turns out, saliva and sputum samples are paired, and share a value of FEV1!
# This insight from our plot might lead us to analyze our data differently to account
# for these paired samples. In our case, we will ignore this. 

# ----------------- Fitting regression models ----------------

# Now, let's look at the association between digital PCR measurements and whether
# the sample is from saliva or sputum. We'll use the lm command. 
  
# Note that we are deciding to model our outcome on the linear scale. Some people
# might model the log of this variable. Both options are valid, as long as 
# you interpret the coefficients in terms any transformations you've done.

mod_type <- lm(ddpcr ~ `Sample Type`, data = both)

# Notice the syntax we use here. The first argument is called a formula. It takes
# the form "dependent variable ~ independent variable". The second argument specifies
# the dataset that contains the variables you are using. 

# We use `` around the variable name `Sample Type` because it has a space in it. 
# This is another good reason not to have spaces in your variable names! 

# To look at the results of our model, we can use the function `summary`

summary(mod_type)

# How do we interpret this output? 

# We have two rows. The first row corresponds with the intercept of our model.
# Following rows correspond with covariates. 

# We have 4 columns. The first column gives the estimate of the coefficient for 
# the intercept or covariate on the row. The second column gives the standard error
# (remember, this is how variable our estimate is). 

# `lm` will by default run a hypothesis test with the hypothesis that the 
# coefficient for each covariate is equal to 0. If this hypothesis is rejected, 
# this suggests a statistically significant relationship between our covariate and 
# outcome. The third column of our output gives the test statistic for that hypothesis
# test, and the fourth column gives a p-value. 

# What does the second row of the output tell us? 

# Our estimate of 1,006,955 for the coefficient for `Sample Type`Sputum tells us 
# that the estimated mean observed ddPCR for a sample of sputum is 1,000,955 
# units higher than the estimated mean obseved ddPCR for a sample from saliva. 

# We have a p-value of 0.01355. If we use a alpha level of 0.05, this gives us enough
# evidence to reject the null hypothesis of no relationship between ddpcr and sample
# type. 

# Now, let's fit a model with multiple covariates. Let's consider Treatment and 
# Sample Type. 

mod_treat_type <- lm(ddpcr ~ `Sample Type` + `Treatment Group`, data = both)

# We can again use `summary` to check out the results. 

summary(mod_treat_type)

# How would you interpret these results? 
# Here, we have a significant relationship between Sample Type and ddpcr, but not
# Treatment and ddpcr. 

# What if we think there is an interaction between the effects of Treatment and 
# Sample Type on ddpcr? 

# We can fit an interaction model by replacing the `+` above with `*`. 

mod_interact <- lm(ddpcr ~ `Sample Type` * `Treatment Group`, data = both)
summary(mod_interact)

# Now we have an extra row. Not only do we have rows for `Sample Type`Sputum and 
# Treatment, but we also have an interaction row. The estimate for the interaction
# can be interpreted as the following: 

# -409,238 is the estimated difference in mean bacterial concentration in sputum for
# people on antibiotics compared to people not on antibiotics minus the estimated 
# difference in mean bacterial concentration in saliva for people on antibiotics 
# compared to people not on antibiotics.

# Note that our estimates and standard errors for `Sample Type`Sputum and 
# `Treatment Group`ON are different for our two models. This is because whenever we
# change something (add a covariate, add an interaction, etc) we are fitting a 
# different model. There are methods available to help choose between models, but
# ONE TAKEAWAY - ONLY ADD VARIABLES THAT YOU CARE ABOUT? SOMETHING ELSE?