### radEmu lab
### David Clausen
### July 24th, 2022


### In this lab we'll explore a dataset published by Wirbel et al. (2019).
### (https://www.nature.com/articles/s41591-019-0406-6)
### This is a meta-analysis of case-control studies, meaning that Wirbel
### et al. collected raw sequencing data from studies other researchers
### conducted and re-analyzed it (in this case, they also collected some
### new data of their own).

### Wirbel et al. published two pieces of data we'll focus on today:

# metadata giving demographics and other information about participants

# an mOTU table

### We'll look at a subset of all 849 mOTUs Wirbel et al. published

### We're most interested to compare microbial abundance in cases diagnosed
### with colorectal cancer to abundances in controls (without this diagnosis)

### First let's load libraries we'll need
library(tidyverse)
library(Matrix)
library(remotes)
library(monotone)
install_github("https://github.com/statdivlab/radEmu")
library(radEmu)



metadata <-
  read_csv("https://raw.githubusercontent.com/statdivlab/stamps2022/main/Sunday-afternoon/labs/rademu_lab/wirbel_et_al_metadata.csv")
head(metadata)

### Let's see how many observations we have among cases ("CRC") and
### controls ("CTR")
metadata %>%
  group_by(Group) %>%
  summarize(n = length(Group))

### metadata$Group tells us which participants are cases and which are controls

### We have data from studies in 5 different countries
### How much from each study, you ask? Let's find out!
metadata %>%
  group_by(Country) %>%
  summarize(n = length(Group))

### Let's see how many cases and controls were enrolled in each study as well
metadata %>%
  with(table(Country, Group))

### Now let's load the mOTU table
mOTU_table <-
  read_csv("https://raw.githubusercontent.com/statdivlab/stamps2022/main/Sunday-afternoon/labs/rademu_lab/wirbel_et_al_mOTUs.csv")


# let's take a peek at the mOTU table
head(mOTU_table)
# how are the columns in this object named?
# (Looking at the metadata again may be informative)

### column "X1" of our mOTU table contains names of our mOTUs
### let's store them in a separate object
mOTU_names <- mOTU_table$mOTU


### now we'll clean up the mOTU table a bit
mOTU_table <- mOTU_table %>%
  dplyr::select(-1) %>% #remove first column (mOTU names)
  as.matrix() %>% # make this whole deal a matrix
  t() %>% # transpose so that rows are samples and columns are mOTU names
  (Matrix::Matrix) #now make our transposed matrix into a ~fancy~ Matrix
# using the Matrix package

### we need to keep track of which columns of our transposed mOTU_table
### correspond to which mOTUs, so let's name them!
colnames(mOTU_table) <- mOTU_names

### we're not *quite* done with data cleaning and organization
### first, if we compare the dimensions of mOTU_table and metadata...
dim(mOTU_table)
dim(metadata)

### yike -- we have more measurements than metadata!
### in general, this would require further investigation...
### but you can take my word that we don't need the extra observations
### let's get rid of them!

# fun fact: I (David) use sapply a lot
# apparently very few other people do
# in any case, let me talk you through what this sapply is doing:
rows_in_metadata <- sapply( #sapply basically says "take X and do Y to each
  #element of X
  rownames(mOTU_table), # <- this is X, so
  #so far we've told sapply that we want it to do
  #something using the row names of mOTU_table
  #(the row names here tell us what sample each row
  #of the mOTU table contains observations on)
  #this next bit specifies what to do with the row names
  #namely, we're asking it to tell us if each row name
  #exists inside metadata$Sample_ID -- in other words,
  #to tell us which samples we have mOTU data for we
  #also have metadata for
  function(x) x %in% metadata$Sample_ID # <- this is Y
)


#so now that we've used everyone's favorite R function, sapply,
#we have an object "rows_in_metadata" that tells us which rows
#of mOTU_table correspond to rows in our metadata
#we'll use this to subset mOTU_table so we have data for the same
#participants in metadata and mOTU_table
mOTU_table <- mOTU_table[rows_in_metadata,]

#note: I'm using base R to do this -- you could also use tidyverse
#... which some of my colleagues might prefer


### in any case, let's check our dimensions again now
dim(mOTU_table)
dim(metadata)

### success!
### let's check that rows of metadata and mOTU_table are in the same order
order_for_metadata <- sapply(rownames(mOTU_table),
                             function(x) which(x == metadata$Sample_ID))

# if rows are in same order, this plot should look like y = x
plot(order_for_metadata)

# ok nope - not in same order, so let's reorder:
metadata <- metadata[order_for_metadata,]

# now let's check our work and look at the ordering again
order_for_metadata <- sapply(rownames(mOTU_table),
                             function(x) which(x == metadata$Sample_ID))

plot(order_for_metadata)

# yay
# again -- I'm using base R here. You could also use inner_join() to do the
# same thing (I am told)

### if we look at the first few rownames of mOTU_table and sample_IDs
### in metadata, we see the order also matches

metadata$Sample_ID %>% head
rownames(mOTU_table) %>% head

### one last data organizing / cleaning / beautifying step...
### Wirbel et al. published their mOTU table...
### after dividing each row (i.e., count data for each sample) by its total
### and then rounding anything below 1e-6 to zero
### ... please don't do this! Publish your count data as is -- other
### researchers can easily transform counts, but un-transforming is
### harder!

### in any case, Wirbel et al. also published library size (i.e., count totals)
### by sample, so we can recover most of the count data by multiplying each
### row by the corresponding library size (we still don't get back counts that
### were rounded to zero)

for(i in 1:nrow(mOTU_table)){
  mOTU_table[i,] <- mOTU_table[i,]*metadata$Library_Size[i]
}

mOTU_table[1:5, 1:5]

### now we'll pull out some taxa we might be interested in to fit models on
### (we can also fit a model to all mOTUs, but this takes longer)
#yay another sapply
fuso <- sapply(mOTU_names, #go through the names of our mOTUs
               #and for each name:
               
               function(x) grepl("Fusobac", #see if the name contains the string
                                 #"Fusobac"
                                 #if so, return TRUE; otherwise FALSE
                                 x,
                                 fixed = TRUE))
#now do the same thing for some other genera:
prevo <- sapply(mOTU_names,function(x) grepl("Prevotella ",x,fixed = TRUE))
porph <-  sapply(mOTU_names,function(x) grepl("Porphyromonas ",x,fixed = TRUE))
clostr <- sapply(mOTU_names, function(x) grepl("Clostridium",x, fixed = TRUE))
faecali <- sapply(mOTU_names, function(x) grepl("Faecalibact",x, fixed = TRUE))
eubact <- sapply(mOTU_names, function(x) grepl("Eubact",x, fixed = TRUE))

# take a look at "eubact" -- what does this look like? What information
# does it contain?

head(eubact)

### "Group" in metadata indicates whether each participant is case (diagnosed
### with colorectal cancer) or a control (no such diagnosis)
### Let's make this a factor
unique(metadata$Group)

### make "CTR" the reference category (i.e., the first category)
metadata$Group <- factor(metadata$Group,
                         levels = c("CTR","CRC"))

### we'll stick to the genera Eubacterium, Porphyromonas,
### and Fusobacterium for now
# store names of the mOTUs in these genera
restricted_mOTU_names <- mOTU_names[eubact|porph|fuso]

# Figure out which columns of mOTU_table contain observations in
# Eubacterium, Porphyromonas, or Fusobacterium
which_mOTU_names <- which(eubact|porph|fuso)

# Among just the mOTUs in Eubacterium, Porphyromonas, or Fusobacterium,
# flag those that are in Eubacterium
eubact_restr <- sapply(restricted_mOTU_names,
                       function(x) grepl("Eubact",x, fixed = TRUE))

### Ooh ok now we're going to define the constraint we'll use with radEmu
### For now, we'll make this a median over Eubacterium
### (so we require the median of effects we estimate for Eubacterium to be
### zero

# fun question time: how does using this constraint affect the interpretation of
# our estimates?

constraint_fn <- function(x){median(x[eubact_restr])}

### Great! We're ready to start fitting models!

# We'll begin with data from a Chinese study Wirbel et al. analyzed
ch_study_obs <- which(metadata$Country %in% c("CHI"))

# let's fit a model!
# ... emuFit will talk at you a bit, but don't worry about it
ch_fit <-
  emuFit(~ Group, # this is a formula telling radEmu what predictors to
         # use in fitting a model
         # we are using Group -- i.e., an indicator for which
         # participants have CRC diagnoses and which do not
         # as our predictor
         covariate_data = metadata[ch_study_obs, #ch_study obs = we're
                                   # only looking at rows
                                   # containing observations
                                   # from the Chinese study
         ], # covariate_data
         # contains our predictor
         # data
         Y = mOTU_table[ch_study_obs,
                        which_mOTU_names], #which_mOTU_names = we're limiting the taxa
         # we're looking at to Eubacterium,
         # Porphyromonas, and Fusobacterium
         # (which_mOTU_names was constructed above)
         
         constraint_fn = constraint_fn, #make sure radEmu is using our
         #constraint function
         optim_only = TRUE, # you can ignore the rest of these arguments
         # for now
         tolerance = 1,
         method = "FL",
         verbose= TRUE,
         reweight = TRUE)

### this will take a moment and emuFit will not talk to you (computation
### is happening in parallel in environments we don't have direct access to)

ch_fit_cis <- emuCI(emuMod = ch_fit, #we have to give emuCI a fitted
                    #radEmu model
                    nboot = 100, # radEmu uses a bootstrap to form confidence
                    # intervals -- in a publication, we recommend
                    # more bootstrap iterations (1000 is great)
                    parallel = TRUE, # speed up by
                    # doing bootstrap fits in parallel?
                    seed = 0, #set seed for reproducibility
                    ncore = 6) #how many cores to parallelize over?

### Ok we have estimates and confidence intervals for the group effect
### Let's take a look:

ch_fit_cis%>%
  filter(row ==2) %>%
  mutate(taxon = factor(restricted_mOTU_names,
                        levels = mOTU_names[c(
                          which(eubact),which(porph),which(fuso)
                        )]
  )) %>%
  mutate(Genus = sapply(taxon, function(x) ifelse(
    grepl("Eubact",x,fixed = TRUE),"Eubacterium",
    ifelse(grepl("Fuso",x,fixed = TRUE),
           "Fusobacterium",
           "Porphyromonas")
  ))) %>%
  ggplot() + geom_point(aes(x = taxon, y = estimate,
                            color = Genus),
                        size = .5) +
  geom_errorbar(aes(x = taxon, ymin = lower, ymax = upper,
                    color = Genus),
                width = .25)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

### The y-axis is a bit wild because of large uncertainty in
### the estimate for one mOTU
### We can fix this with coord_cartesian()

### here's the help information for this function:
?coord_cartesian

### thoughts on how to make the y scale more reasonable
### in the plot below?

ch_fit_cis%>%
  filter(row ==2) %>%
  mutate(taxon = factor(restricted_mOTU_names, #define "taxon" to be a factor
                        levels = mOTU_names[c( #put the levels of this factor
                          #in a particular order:
                          which(eubact), #Eubacterium mOTUs first
                          which(porph),  #then Porphyromonas
                          which(fuso)    #then Fusobacterium
                        )] #(specifying the order of the levels puts the
                        #x-axis of the plot below in prettier order)
  )) %>%
  mutate(Genus = sapply(taxon, #sapply! to define a new variable: Genus
                        function(x) ifelse(
                          grepl("Eubact",x,fixed = TRUE), #if the mOTU contains
                          #the string "Eubact",
                          "Eubacterium", #return genus Eubacterium
                          ifelse(grepl("Fuso",x,fixed = TRUE), #otherwise if it contains string "Fuso",
                                 "Fusobacterium", #return genus Fusobacterium
                                 "Porphyromonas") #otherwise return genus Porphyromonas
                        ))) %>%
  ggplot() + geom_point(aes(x = taxon, y = estimate,
                            color = Genus),
                        size = .5) +
  geom_errorbar(aes(x = taxon, ymin = lower, ymax = upper,
                    color = Genus),
                width = .25)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian( # here for you to fill in woooo
    
  )

# ok now that (hopefully) the plot is easier to read...
# what patterns do you see across mOTUs?

# many of the estimates for Eubacterium are close to zero -- why?


### Let's look at a French study Wirbel et al. analyized
# -- do we see similar patterns?
fr_study_obs <- which(metadata$Country %in% c("FRA"))

fr_fit <-
  emuFit(~ Group,
         covariate_data = metadata[fr_study_obs,],
         Y = mOTU_table[fr_study_obs,which_mOTU_names],
         optim_only = TRUE,
         tolerance = 1,
         constraint_fn = constraint_fn,
         method = "FL", # use Firth-penalized likelihood
         verbose= TRUE,
         reweight = TRUE) # tell us what's happening during optimization

### this will take a moment

fr_fit_cis <- emuCI(emuMod = fr_fit,
                    nboot = 100,
                    parallel = TRUE,
                    seed = 0,
                    ncore = 6)

fr_fit_cis%>%
  filter(row ==2) %>%
  mutate(taxon = factor(restricted_mOTU_names,
                        levels = mOTU_names[c(
                          which(eubact),which(porph),which(fuso)
                        )]
  )) %>%
  mutate(Genus = sapply(taxon, function(x) ifelse(
    grepl("Eubact",x,fixed = TRUE),"Eubacterium",
    ifelse(grepl("Fuso",x,fixed = TRUE),
           "Fusobacterium",
           "Porphyromonas")
  ))) %>%
  ggplot() + geom_point(aes(x = taxon, y = estimate,
                            color = Genus),
                        size = .5) +
  geom_errorbar(aes(x = taxon, ymin = lower, ymax = upper,
                    color = Genus),
                width = .25)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim = c(-5,14))


### Let's look at those together
rbind(
  ch_fit_cis %>%
    mutate(study = "Chinese"),
  fr_fit_cis %>%
    mutate(study = "French")
) %>%
  filter(row ==2) %>%
  mutate(taxon = factor(rep(restricted_mOTU_names,2),
                        levels = mOTU_names[c(
                          which(eubact),which(porph),which(fuso)
                        )]
  )) %>%
  mutate(Genus = sapply(taxon, function(x) ifelse(
    grepl("Eubact",x,fixed = TRUE),"Eubacterium",
    ifelse(grepl("Fuso",x,fixed = TRUE),
           "Fusobacterium",
           "Porphyromonas")
  ))) %>%
  ggplot() + geom_point(aes(x = taxon, y = estimate,
                            group = study,
                            color = Genus),
                        size = .5) +
  geom_errorbar(aes(x = taxon, ymin = lower, ymax = upper,
                    group = study,
                    linetype = study,
                    color = Genus),
                width = .25)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim = c(-5,14))

### Hard to read overlapping CIs -- let's fix

rbind(
  ch_fit_cis %>%
    mutate(study = "Chinese"),
  fr_fit_cis %>%
    mutate(study = "French")
) %>%
  filter(row ==2) %>%
  mutate(taxon = factor(rep(restricted_mOTU_names,2),
                        levels = mOTU_names[c(
                          which(eubact),which(porph),which(fuso)
                        )]
  )) %>%
  mutate(Genus = sapply(taxon, function(x) ifelse(
    grepl("Eubact",x,fixed = TRUE),"Eubacterium",
    ifelse(grepl("Fuso",x,fixed = TRUE),
           "Fusobacterium",
           "Porphyromonas")
  ))) %>%
  ggplot() + geom_point(aes(x = taxon, y = estimate,
                            group = study,
                            color = Genus),
                        position = position_dodge(.25), #the magic fix
                        #for overlapping CIs
                        size = .5) +
  geom_errorbar(aes(x = taxon, ymin = lower, ymax = upper,
                    group = study,
                    linetype = study,
                    color = Genus),
                position = position_dodge(.25), #use position_dodge() here too
                #so error bars are adjusted too
                width = .25)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim = c(-5,14))


### In the Chinese study, samples were taken from patients undergoing
### colonoscopy, with some patients providing samples before
### colonoscopy and some after -- maybe we should try to adjust for this
### by including timing relative to colonoscopy as a predictor in our model
### (Why?)

ch_fit_timing <-
  emuFit(~ Group + Sampling_rel_to_colonoscopy,
         covariate_data = metadata[ch_study_obs,],
         Y = mOTU_table[ch_study_obs,which_mOTU_names],
         optim_only = TRUE,
         tolerance = 1,
         constraint_fn = constraint_fn,
         method = "FL", # use Firth-penalized likelihood
         verbose= TRUE,
         reweight = TRUE)

### Error :(
### "Matrices must have same number of rows for arithmetic" -- seems like
### we might not be giving emuFit inputs that have the correct dimensions

### On that theory, let's look at "Sampling_rel_to_colonoscopy"
metadata$Sampling_rel_to_colonoscopy[ch_study_obs]
### What do you notice about the values of this variable?








### ... spoiler below









### that's right! there is an NA in Sampling_rel_to_colonoscopy
### let's record which observations we have for which
#### Sampling_rel_to_colonoscopy is not NA
timing_available <- !is.na(metadata$Sampling_rel_to_colonoscopy[ch_study_obs])

###... and fit a model using only those observations
ch_fit_timing <-
  emuFit(~ Group + Sampling_rel_to_colonoscopy,
         covariate_data = metadata[ch_study_obs,][timing_available, #only rows
                                                  # with non-NA values of
                                                  # Sampling_rel_to_colonoscopy
         ],
         Y = mOTU_table[ch_study_obs,which_mOTU_names][timing_available, #same
                                                       # thing except now we're
                                                       # subsetting an already-
                                                       # subsetted version of
                                                       # mOTU_table...
                                                       # I write the most elegant code
         ],
         optim_only = TRUE,
         tolerance = 1,
         constraint_fn = constraint_fn,
         method = "FL",
         verbose= TRUE,
         reweight = TRUE)

chi_fit_timing_cis <- emuCI(emuMod = ch_fit_timing,
                            nboot = 100,
                            parallel = TRUE,
                            seed = 0,
                            ncore = 6)

### Let's plot the timing-adjusted and unadjusted fits together
rbind(
  ch_fit_cis %>%
    mutate(adjustment = "No adjustment"),
  chi_fit_timing_cis %>%
    mutate(adjustment = "Adjustment for timing")
) %>%
  filter(row ==2) %>%
  mutate(taxon = factor(rep(restricted_mOTU_names,2),
                        levels = mOTU_names[c(
                          which(eubact),which(porph),which(fuso)
                        )]
  )) %>%
  mutate(Genus = sapply(taxon, function(x) ifelse(
    grepl("Eubact",x,fixed = TRUE),"Eubacterium",
    ifelse(grepl("Fuso",x,fixed = TRUE),
           "Fusobacterium",
           "Porphyromonas")
  ))) %>%
  ggplot() + geom_point(aes(x = taxon, y = estimate,
                            group = adjustment,
                            color = Genus),
                        position = position_dodge(.25),
                        size = .5) +
  geom_errorbar(aes(x = taxon, ymin = lower, ymax = upper,
                    group = adjustment,
                    linetype = adjustment,
                    color = Genus),
                position = position_dodge(.25),
                width = .25)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim = c(-5,14))

### seems like point estimates and CIs are pretty similar for adjusted
### and unadjusted fits -- does this mean we don't need to worry about
### the impact of timing of sample collection? Why or why not?

### We can also fit a single model to data from multiple studies

fr_ch_study_obs <- which(metadata$Country %in% c("CHI","FRA"))

fr_ch_fit <-
  emuFit(~ Group + Country,
         covariate_data = metadata[fr_ch_study_obs,],
         Y = mOTU_table[fr_ch_study_obs,which_mOTU_names],
         optim_only = TRUE,
         tolerance = 1,
         constraint_fn = constraint_fn,
         method = "FL", # use Firth-penalized likelihood
         verbose= TRUE,
         reweight = TRUE) # tell us what's happening during optimization

### this will take a moment (longer than previous fits):

fr_ch_fit_cis <- emuCI(emuMod = fr_ch_fit,
                       nboot = 100,
                       parallel = TRUE,
                       seed = 0,
                       ncore = 6)

## Let's look at results
fr_ch_fit_cis%>%
  filter(row ==2) %>%
  mutate(taxon = factor(restricted_mOTU_names,
                        levels = mOTU_names[c(
                          which(eubact),which(porph),which(fuso)
                        )]
  )) %>%
  mutate(Genus = sapply(taxon, function(x) ifelse(
    grepl("Eubact",x,fixed = TRUE),"Eubacterium",
    ifelse(grepl("Fuso",x,fixed = TRUE),
           "Fusobacterium",
           "Porphyromonas")
  ))) %>%
  ggplot() + geom_point(aes(x = taxon, y = estimate,
                            color = Genus),
                        size = .5) +
  geom_errorbar(aes(x = taxon, ymin = lower, ymax = upper,
                    color = Genus),
                width = .25)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


### Notice that we didn't include Sampling_rel_to_colonoscopy in this regression
### (Optional exercise: fit the model above with Sampling_rel_to_colonoscopy
### added to the model)
### (Well, everything is optional -- you can fit a model adjusting for age or
### BMI as well. Up to you!)

