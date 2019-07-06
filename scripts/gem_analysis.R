library(furrr)
plan(multiprocess)
library(survival)
library(IRanges)
library(tidyverse)

source('scripts/gem_utils.R')
source('scripts/gem_constants.R')

mingrades <- c(2L, 3L)
lags <- c(0, 7, 14)

params <- bind_rows(crossing(mingrade = mingrades, lag = lags, permanent = FALSE),
                    crossing(mingrade = mingrades, lag = 0, permanent = TRUE))

outfile_name <- paste0(params_to_string(mingrades, lags), '.RDS')
print(outfile_name)

if (get_os() == 'unix') {

  future_pmap(.l = list(mingrade = params$mingrade,
                        lag_days = params$lag,
                        permanent_exposure = params$permanent),
              .f = run_analysis) -> analyses

  saveRDS(object = analyses, file = outfile_name)
}

if (get_os() == 'mac') {

  analyses <- readRDS(outfile_name)

  # map(.x = analyses, .f = ~ paste(names(attributes(.)[4:7]), attributes(.)[4:7]))

  map(.x = analyses,
      .f = make_results_table) -> results_tables

  walk(.x = results_tables,
       .f = make_plots)

}
