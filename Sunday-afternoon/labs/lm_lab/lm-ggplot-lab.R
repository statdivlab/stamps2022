## Amy Willis & Sarah Teichman

library(tidyverse)
library(readxl)
ddpcr <- read_csv("lm_lab/data/ddPCR.csv")
meta <- read_xlsx("lm_lab/data/Sample Metadata.xlsx")

both <- meta %>% 
  inner_join(ddpcr, by = "sample_name") 

both <- both %>%
  rename(ddpcr = `Average Replicates, stock DNA (copies /µl)`) %>%
  mutate(ddpcr = gsub(",", "", ddpcr)) %>%
  mutate(ddpcr = as.numeric(ddpcr)) 

both %>%
  ggplot(aes(x = `Sample Type`, y = ddpcr, col = `Subject ID`)) +
  geom_jitter(height = 0, width = 0.3)

both %>%
  ggplot(aes(x = `Sample Type`, y = ddpcr, col = `Subject ID`)) +
  geom_jitter(height = 0, width = 0.3) +
  scale_y_log10() +
  ylab("Average observed ddPCR\n(copies/ul)")

both %>%
  filter(`Sample Type` == "Sputum") %>%
  ggplot(aes(x = `Treatment Group`, y = ddpcr, col = `Subject ID`)) +
  geom_jitter(height = 0, width = 0.3) +
  scale_y_log10() +
  ylab("Average observed ddPCR\n(copies/ul)")

both %>%
  filter(`Sample Type` == "Sputum") %>%
  lm(ddpcr ~ `Treatment Group`, .) %>%
  summary



mod1 <- lm(ddpcr ~ `Treatment Group` + `Sample Type`, data = both)

mod2 <- lm(ddpcr ~ `Treatment Group`*`Sample Type`, data = both)



