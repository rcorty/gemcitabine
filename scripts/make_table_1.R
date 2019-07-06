spsm <- suppressPackageStartupMessages
spsm(library(haven))
spsm(library(tidyverse))

source('scripts/utils_topdown.R')

demog <- read_demog()

demog %>% glimpse()

demog$age %>% median()

demog$age %>% range()

demog$sex %>% table()
demog$sex %>% table() %>% `/`(430) %>% round(2)

demog$race %>% table()
demog$race %>% table() %>% `/`(430) %>% round(2)


performance <- read_baseline_performance()

performance %>% glimpse()

performance %>%
  pull(kps) %>%
  table(useNA = 'ifany')

performance %>%
  pull(kps) %>%
  table() %>%
  `/`(429) %>%
  round(2)


# trusting NEJM article for lesions



