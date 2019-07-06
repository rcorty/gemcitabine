spsm <- suppressPackageStartupMessages
spsm(library(haven))
spsm(library(tidyverse))
spsm(library(furrr))
plan(multiprocess)
spsm(library(IRanges))

source('scripts/utils_topdown.R')
source('scripts/constants_topdown.R')

total_days <- sum(read_disc()$dsc_day)

grade_at_least <- 3L
lag_days <- 0

list(ae_any = read_aes(mingrade = grade_at_least)) -> events_ae_any

future_map(.x = common_ae_codes_by_class,
           .f = read_aes,
           mingrade = grade_at_least) -> events_ae_byclass

future_map(.x = names_of_common_aes,
           .f = read_aes,
           mingrade = grade_at_least) -> events_ae_bycode

common_ae_table

# number of pts with each AE,
# number of times each AE was experienced,
# and duration each AE was experienced
map(.x = c(events_ae_any, events_ae_byclass, events_ae_bycode),
    .f = ~ c(num_pts = sum(unlist(map(.x = .$event,
                          .f = ~ length(.) > 0))),
             num_exps = sum(unlist(map(.x = .$event,
                            .f = ~ length(.)))),
             num_days = sum(unlist(map(.x = .$event,
                            .f = width))))) %>%
  do.call(what = bind_rows, args = .) %>%
  mutate(ae = names(c(events_ae_any, events_ae_byclass, events_ae_bycode)),
         frac_pts = 100*round(num_pts/430, 2)) %>%
  # filter(num_exps > 13) %>%
  unite(col = num_pts, num_pts, frac_pts, sep = ' (') %>%
  mutate(num_pts = paste0(num_pts, '%)'),
         incid_rate = round(1000*num_exps/total_days, 1)) %>%
  arrange(desc(num_exps)) %>%
  select(`Adverse Event` = ae,
         `Number of Patients` = num_pts,
         `Number of Events` = num_exps,
         `Total Duration (days)` = num_days,
         `Incidence rate (events per 1000 patient-days)` = incid_rate) ->
  table2

# write_csv(x = table2, path = 'results/table2.csv')
write_csv(x = table2, path = 'results/tableS1.csv')

print(x = table2, n = 30)
