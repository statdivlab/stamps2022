# Linear regression and ggplot lab 
## prepared by Amy Willis & Sarah Teichman

# In this lab, we'll practice running linear regression with the function `lm`. 
# We'll also manipulate our data with the `tidyverse` plot our data with `ggplot`.
# There will be a few intentional errors here for you to fix. If a line of coding
# isn't working, think about how you can fix it to make it work. 

# Let's begin by loading in the tidyverse, which we'll need for this tutorial.

library(tidyverse)

# ------------------------ loading the data -----------------------

# Next, let's load our data. We'll start with two data sets. This data comes from
# a study of people living with cystic fibrosis. It consists of 16S and droplet 
# digital PCR (ddPCR) data from both saliva and sputum (coughed up mucus)
# from study participants. More information can be found here: 
# https://journals.asm.org/doi/full/10.1128/mSystems.00296-20?rfr_dat=cr_pub++0pubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&

ddpcr <- read_csv("https://raw.githubusercontent.com/statdivlab/stamps2022/main/Sunday-afternoon/labs/lm_lab/data/ddPCR.csv")
meta <- read_csv("https://raw.githubusercontent.com/statdivlab/stamps2022/main/Sunday-afternoon/labs/lm_lab/data/meta_data.csv")

# ------------------------ exploring the data ---------------------

# We can take a look at the beginning of each dataset with the function head

head(dpcr)
head(meta)

# We'd like to combine these two data sets together. They both contain the variable 
# `sample_name`, so we will match observations based on this variable. `inner_join` 
# combines these two datasets into a single one, cleverly using the `by` argument 
# to figure out how to rows match up. 

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

# You can run each of these lines individually to further explore what each command
# is doing 

# Earlier we looked at the first few observations of each dataset. To see a full list
# of variables in our dataset, we can use the `names()` function.

names(both)

# We'll be focusing on a few variables from this dataset: ddpcr, FEV1, Treatment 
# Group, Sample Type, and Subject ID. 

# ddpcr: droplet digital PCR, giving the observed total bacterial load in the sample in copies/uL
# FEV1: Forced expiratory volume (FEV) measures how much air a person can exhale
# during a forced breath (L/s)
# Treatment Group: "ON" if the sample is from a subject receiving antibiotics, 
# "OFF" otherwise 
# Sample Type: the body site of origin of the sample, either Saliva or Sputum
# Subject ID: an ID identifying the participant 

# Before we do inference, we should plot our data. This lets us get a sense of 
# general trends and build intuition about our data. 

# We'll use ggplot for plotting. Remember that in ggplot we start with the base
# of our plot, where we specify our dataset and variables, and then we can build
# up our plot by adding layers with the `+` operator. We'll build this first plot
# step by step.

ggplot(data = both, aes(x = `Sample Type`, y = ddpcr, col = `Subject ID`))

# This sets up our axis labels. We are now ready to add a plotting layer. 

ggplot(data = both, aes(x = `Sample Type`, y = ddpcr, col = `Subject ID`)) +
  geom_jitter(height = 0, width = 0.3)

# geom_jitter is a great way to make a scatterplot with a categorical variable on 
# the x-axis. It spreads the points out so that they do not all line up exactly 
# for the same values of the x variable.

ggplot(data = both, aes(x = `Sample Type`, y = ddpcr, col = `Subject ID`)) +
  geom_jitter(height = 0, width = 0.3) + 
  labs(y = "Average observed ddPCR\n(copies/ul)",
       title = "Digital PCR by Sample Type")

# We like to add informative labels with the `labs` layer to improve on our variable
# names. 

ggplot(data = both, aes(x = `Sample Type`, y = ddpcr, col = `Subject ID`)) +
  geom_jitter(height = 0, width = 0.3) + 
  labs(y = "Average observed ddPCR\n(copies/ul)",
       title = "Digital PCR by Sample Type") + 
  scale_y_log10()

# The ddpcr data is skewed, with a few larger outliers taking up a lot of space on
# our plot. We can transform ddpcr to the log of ddpcr to visualize this in a 
# different way. 

ggplot(data = both, aes(x = `Sample Type`, y = ddpcr, col = `Subject ID`)) +
  geom_jitter(height = 0, width = 0.3) + 
  labs(y = "Average observed ddPCR\n(copies/ul)",
       title = "Digital PCR by Sample Type") + 
  scale_y_log10() + 
  theme(plot.title = element_text(hjust = 0.5))

# ggplot by default left-justifies titles. The theme layer centers the title. 

# What do we learn from this plot? 
# It seems that there are sputum samples generally have higher concentration (ddPCR) values
# than saliva samples. 

ggplot(data = both, aes(x = `Treatment Group`, y = ddpcr, col = `Subject ID`)) 
  geom_jitter(height = 0, width = 0.3) +
  labs(y = "Average observed ddPCR\n(copies/ul)",
       title = "Digital PCR by Treatment") + 
  scale_y_log10() + 
  theme(plot.title = element_text(hjust = 0.5))

# What do we learn from this plot? 
# There are a lot more non-treated samples than treated samples. The general trend
# seems to be higher observed concentration values from non-treated
# samples than treated samples.

ggplot(data = both, aes(x = as.numeric(FEV1), y = ddpcr, col = `Subject ID`)) +
  geom_point() + 
  labs(x = "FEV Value",
       y = "Average observed ddPCR\n(copies/ul)",
       title = "Digital PCR by FEV Value") + 
  scale_y_log10() +
  theme(plot.title = element_text(hjust = 0.5))

# What do we learn from this? 

# There might be a slight negative correlation between concentration values and FEV values
# There seem to be two observations for each value on the x-axis. Why might that be?

# It turns out, saliva and sputum samples are paired (one each for each subject), 
# and so share a value of FEV1! Whoops - we forgot about this - or maybe the postdoc
# left and we weren't told. 

# Insights like this from our plot might lead us to analyze our data. 

# ----------------- Fitting regression models ----------------

# Now, let's look at the association between digital PCR measurements and whether
# the sample is from saliva or sputum. We'll use the lm command. 
  
# Note that we are deciding to model the mean of our outcome on the linear scale. We could also
# choose to model the mean of the log of this variable. Both options are valid, as long as 
# you interpret the coefficients in terms any transformations you've done.

mod_type <- lm(ddpcr ~ `Sample Type`, data = both)

# Notice the syntax we use here. The first argument is called a formula. It takes
# the form "response variable ~ predictor variable". The second argument specifies
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
# (remember, this is an estimate of the standard deviation of our estimate). 

# `lm` will by default run a hypothesis test with the hypothesis that the 
# coefficient for each covariate is equal to 0. If this hypothesis is rejected, 
# we determine that there is statistically significant relationship between our covariate and 
# outcome. The third column of our output gives the test statistic for that hypothesis
# test, and the fourth column gives a p-value. 

# What does the second row of the output tell us? 

# Our estimate of 1,006,955 for the coefficient for `Sample Type`Sputum tells us 
# that the estimated mean observed ddPCR reading for a sample of sputum is 1,000,955 
# counts/uL higher than the estimated mean observed ddPCR for a sample from saliva. 

# We have a p-value of 0.01355. If we use a alpha level of 0.05, we
# reject the null hypothesis that there is no difference in the true mean concentration 
# in Sputum and Saliva samples. 

# Now, let's fit a model with multiple covariates. Let's consider Treatment and 
# Sample Type. 

mod_treat_type <- lm(ddpcr ~ Sample Type + Treatment Group, data = both)

# We can again use `summary` to check out the results. 

summary(mod_treat_type)

# How would we interpret these results? 
# Here, we have a significant difference in mean ddPCR across 
# Sample Types (at the alpha=0.05 level), but not across Treatment groups. 

# What if we think there is an interaction between the effects of Treatment and 
# Sample Type on ddpcr? i.e., a different change between Treatment ON and OFF
# for Sputum vs Saliva. 

# We can fit an interaction model by replacing the `+` above with `*`. 

mod_interact <- lm(both, ddpcr ~ `Sample Type` * `Treatment Group`)
summary(mod_interact)

# Now we have an extra row. Not only do we have rows for `Sample Type`Sputum and 
# Treatment, but we also have an interaction row. The estimate for the interaction
# can be interpreted as the following: 

# -349,375 is the estimated value of the difference in mean bacterial concentration in sputum for
# people on antibiotics compared to people not on antibiotics minus the 
# difference in mean bacterial concentration in saliva for people on antibiotics 
# compared to people not on antibiotics.

# Note that our estimates and standard errors for `Sample Type`Sputum and 
# `Treatment Group`ON are different for our two models. This is because whenever we
# change something (add a covariate, add an interaction, etc) we are fitting a 
# different model. Our strongest advice on the subject of "which model to choose?"
# is to be motivated by your scientific question, and consider 
# the classification of variables presented in the lecture: predictor of interest (always include),
# precision variables (good to include most important ones), 
# confounders (always include, but use the actual definition, not the way most people
# talk about it) and relevant effect modifiers (include if scientifically relevant)). 