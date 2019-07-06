make_plots <- function(results_table) {

  mingrade <- attr(x = results_table, which = 'mingrade')
  perm_exp <- attr(x = results_table, which = 'permanent_exposure')
  lag_days <- ifelse(test = perm_exp,
                     yes = 'inf',
                     no = attr(x = results_table, which = 'lag_days'))

  results_table %>%
    separate(col = outcome, into = c('abx', 'event'), sep = '__') %>%
    filter(!is.na(events_on_exp)) %>%
    mutate(event_fac = factor(x = event,
                              levels = levels_of_events_for_plotting,
                              labels = names_of_events_for_plotting)) %>%
    filter(!duplicated(event_fac)) ->
    to_plot

  to_plot %>%
    filter(grepl('any', abx),
           grepl('ae_any|class', event)) %>%
    mutate(event_fac = factor(paste0(event_fac, ' (', events_on_exp, '/', total_events, ')')),
           event_fac = factor(event_fac, levels = levels(event_fac)[order(total_events)])) ->
    toplot_fig1

  toplot_fig1 %>%
    ggplot(mapping = aes(x = HR, y = event_fac)) +
    geom_vline(xintercept = 1, color = 'gray') +
    geom_point(size = 2) +
    geom_segment(aes(x = LL, xend = UL, yend = event_fac)) +
    geom_text(data = toplot_fig1 %>% filter(p_val < 1),
              mapping = aes(x = HR, y = event_fac,
                            label = paste0('HR = ', round(HR, 2), ', ',
                                           'CI = (', round(LL, 2), ', ', round(UL, 2), '), ',
                                           'p = ', signif(p_val, 2))), #, ', q = ', signif(qvalues, 2))),
              vjust = -1, hjust = 0.5, size = 3) +
    scale_x_continuous(trans = 'log', breaks = c(0.5, 1, 2, 4)) +
    # scale_x_continuous(limits = c(-0.7, 2.5), breaks = -1:3) +
    # scale_x_exp() +
    # theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          strip.background = element_rect(fill = 'lightgray'),
          strip.placement = 'outside',
          # axis.text.y = element_text(margin = margin(0, -10, 0, 10)),
          plot.title = element_text(hjust = 0),
          axis.ticks = element_blank()) +
    labs(title = paste0('Hazard Ratio of Antibiotic Exposure ',
                        ifelse(test = perm_exp,
                               yes = 'with Permanent Effect',
                               no = paste0('with lag ', lag_days, ' days')),
                        '\non Composite Adverse Events Grade ',
                        mingrade, '+'))

  ggsave(filename = paste0('plots/fig1_',
                           params_to_string(mingrade, lag_days),
                           '.pdf'),
         height = nrow(toplot_fig1)/2 + 1,
         width = 7)

  to_plot %>%
    filter(grepl('any', abx),
           grepl('aecode', event)) %>%
    mutate(categ = common_ae_table$class[match(x = tolower(event_fac), table = tolower(common_ae_table$ae_code_american))],
           categ = factor(x = categ,
                          levels = paste0('aeclass_', c('hem', 'gi', 'const', 'hep', 'pulm')),
                          labels = ae_class_names_for_plotting_short)) %>%
    mutate(event_fac = paste0(event_fac, ' (', events_on_exp, '/', total_events, ')'),
           event_fac = factor(x = event_fac, levels = unique(event_fac)[order(total_events)])) ->
    toplot_fig2


  toplot_fig2 %>%
    ggplot(mapping = aes(x = HR, y = event_fac)) +
    geom_vline(xintercept = 1, color = 'gray') +
    geom_segment(mapping = aes(x = LL, xend = UL, yend = event_fac)) +
    geom_point() +
    facet_grid(rows = vars(categ), drop = TRUE, scales = 'free_y', space = 'free_y', switch = 'y') +
    geom_text(data = toplot_fig2 %>% filter(p_val < 1),
              mapping = aes(x = HR, y = event_fac,
                            label = paste0('HR = ', round(HR, 2), ', ',
                                           'CI = (', round(LL, 2), ', ', round(UL, 2), '), ',
                                           'p = ', signif(p_val, 2))), #, ', q = ', signif(qvalues, 2))),
              vjust = -1, hjust = 0.5, size = 3) +
    scale_x_continuous(trans = 'log', breaks = c(0.5, 1, 2, 4, 8)) +
    # theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          strip.background = element_rect(fill = 'lightgray'),
          strip.placement = 'outside',
          # axis.text.y = element_text(margin = margin(0, -10, 0, 10)),
          plot.title = element_text(hjust = 0),
          axis.ticks = element_blank()) +
    labs(title = paste0('Hazard Ratio of Antibiotic Exposure ',
                        ifelse(test = perm_exp,
                               yes = 'with Permanent Effect',
                               no = paste0('with lag ', lag_days, ' days')),
                        '\non Specific Adverse Events Grade ',
                        mingrade, '+'))

  ggsave(filename = paste0('plots/fig2_',
                           params_to_string(mingrade, lag_days),
                           '.pdf'),
         height = nrow(toplot_fig2)/2 + 1,
         width = 7)
}

make_results_table <- function(analysis) {

  analysis$result %>%
    keep(is_tibble) %>%
    bind_rows(.id = 'outcome') ->
    result_table

  attributes(result_table) <- c(attributes(result_table),
                                mingrade = attr(x = analysis, which = 'mingrade'),
                                lag_days = attr(x = analysis, which = 'lag_days'),
                                permanent_exposure = attr(x = analysis, which = 'permanent_exposure'))

  return(result_table)
}

run_analysis <- function(mingrade = 3L,
                         lag_days = 0,
                         permanent_exposure = TRUE,
                         do_cat3_abx = TRUE,
                         do_composite_aes = TRUE,
                         do_individual_aes = TRUE) {

  #### read in exposure data (abx) ####
  message('Reading antibacterial exposure data...')
  exposures <- list(abx_any = read_abx(lag_days = lag_days,
                                       permanent_exposure = permanent_exposure))

  if (do_cat3_abx) {
    exposures <- c(exposures,
                   list(abx_bycat3 = read_abx(
                     classify_by = quos(
                       abx_cat3 == 'QUINOLONE ANTIBACTERIALS' ~ 'abx_quin',
                       abx_cat3 == 'BETA-LACTAM ANTIBACTERIALS, PENICILLINS' ~ 'abx_pen',
                       abx_cat3 == 'OTHER BETA-LACTAM ANTIBACTERIALS' ~ 'abx_betalac',
                       TRUE ~ 'abx_other'
                     ),
                     lag_days = lag_days,
                     permanent_exposure = permanent_exposure
                   )
                   )
    )
  }

  message('Done reading antibacterial exposure data.')

  #### read in event data ####
  message('Reading adverse event data...')
  events <- list(ae_any = read_aes(mingrade = mingrade))
  if (do_composite_aes) {
    events <- c(events,
                future_map(.x = common_ae_codes_by_class,
                           .f = read_aes,
                           mingrade = mingrade))
  }

  if (do_individual_aes) {
    events <- c(events,
                future_map(.x = names_of_common_aes,
                           .f = read_aes,
                           mingrade = mingrade))
  }
  message('Done reading adverse event data.')

  # maybe filter 'out' AE's with <15 events here?
  # to save the time of making a surv_df for them

  #### set up names ####
  result_names <- as.vector(outer(X = names(exposures),
                                  Y = names(events),
                                  FUN = paste, sep = '__'))

  #### run analysis ####
  message('Combining antibacterial and adverse event data...')
  cross_df(list(exposures = exposures,
                events = events,
                baseline = list(select(read_pt_data(), -dsc_day)))) %>%
    mutate(surv_df = pmap(list(exposures, events, baseline), make_surv_df),
           result = map(surv_df, fit_ag_model),
           result = set_names(result, result_names)) ->
    results
  message('Done fitting survival model')

  attr(x = results, which = 'mingrade') <- mingrade
  attr(x = results, which = 'lag_days') <- lag_days
  attr(x = results, which = 'permanent_exposure') <- permanent_exposure
  attr(x = results, which = 'do_cat3_abx') <- do_cat3_abx

  return(results)

}

#' @title fit anderson gill model
#'
#' @param exposure A tibble where each row is a patient. It must have two columns.
#' The first must be called 'pt_id" and be a factor of patient ID's.
#' The second must be called 'exposure' and be an IRanges object.
#' This IRanges object must have columns o
#' @param events a tibble yada yada...TODO
#' @param baseline baseline data
#'
#' @return
#' @export
#'
#' @examples
fit_ag_model <- function(agdata, ag_formula = NULL) {

  if (is.character(agdata))
    return('fewer than 15 events total')

  # drop races that had <5 events, estimates are unstable
  # TODO: maybe combine rare races rather than drop
  agdata %>%
    filter(event) %>%
    count(race) %>%
    filter(n >= 5) %>%
    dplyr::select(-n) ->
    races_prevalences_table

  # if only 1 race present, no use modeling it
  # otherwise, filter out races w unstable estimates
  if (nrow(races_prevalences_table) == 1) {
    agdata %>%
      dplyr::select(-race) ->
      agdata
  } else {
    right_join(agdata, races_prevalences_table, by = 'race') %>%
      mutate(race = droplevels(race)) ->
      agdata
  }

  # only estimate effects of exposures w 5+ events
  # return null if no such exposure
  if (any(grepl('on_abx', names(agdata)))) {

    agdata %>%
      filter(event) %>%
      gather(key = 'expo', value = 'on_off', matches('on_abx')) %>%
      group_by(expo) %>%
      summarise(n = sum(on_off)) %>%
      filter(n >= 5) ->
      events_by_exposure

    if (nrow(events_by_exposure) == 0)
      return('no exposure with 5+ events')

    an <- names(agdata)
    good_names <- an[!grepl('on_abx', an) | an %in% events_by_exposure$expo]

    agdata <- dplyr::select(agdata, good_names)
  }

  if (is.null(ag_formula)) {
    start_idx <- which(names(agdata) == 'event') + 1
    covar_idxs <- start_idx:ncol(agdata)
    covar_names <- names(agdata)[covar_idxs]
    covar_string <- paste(covar_names, collapse = ' + ')
    ag_formula <- as.formula(
      paste0('Surv(start, end, event) ~ ', covar_string)
    )
  }

  # fit the model
  agfit <- coxph(formula = ag_formula,
                 data = agdata,
                 control = coxph.control(iter.max = 50))

  # set up contrasts
  abx_effects <- grep(pattern = 'on_abx', x = names(coef(agfit)), value = TRUE)

  # if we're working on an 'empty' exposure df
  if (length(abx_effects) == 0) {
    return(agfit)
  }

  compare_w_zero <- paste(abx_effects, '== 0')
  if (length(abx_effects) > 1) {
    compare_w_other <- paste(unlist(combn(x = abx_effects,
                                          m = 2,
                                          FUN = paste,
                                          simplify = FALSE,
                                          collapse =' - ')), '== 0')
  } else {
    compare_w_other <- NULL
  }

  # browser()
  # good_terms <- grepl(pattern = 'abx', x = tidy(agfit)$term)
  #
  # tibble(
  #   hypothesis = tidy(agfit)$term[good_terms]
  #   est = tidy(agfit)$estimate[good_terms],
  #   std_err = tidy(agfit)$std.error[good_terms],
  #   p_val = tidy(agfit)$p.value[good_terms]
  # )

  hts <- multcomp::glht(model = agfit,
                        linfct = c(compare_w_zero, compare_w_other))

  tibble(
    hypothesis = c(compare_w_zero, compare_w_other),
    est = get_coef(ht = hts),
    std_err = get_se(ht = hts),
    p_val = get_pval(ht = hts)
  ) %>%
    mutate(total_events = sum(agdata$event),
           events_on_exp = c(events_by_exposure$n,
                             rep(NA, length(compare_w_other))),
           HR = exp(est),
           LL = exp(est - 2*std_err),
           UL = exp(est + 2*std_err),
           num_events = c(events_by_exposure$n,
                          rep(NA, length(compare_w_other))),
           signif = case_when(p_val < 0.001 ~ '***',
                              p_val < 0.01 ~ '**',
                              p_val < 0.05 ~ '*',
                              TRUE ~ 'n.s.')) %>%
    dplyr::select(hypothesis, total_events, events_on_exp,
                  HR, LL, UL, num_events, p_val, signif)
}

make_surv_df <- function(exposures = read_abx(),
                         events = read_aes(),
                         baseline = read_pt_data()) {

  num_events <- sum(unlist(map(.x = events$event, .f = length)))
  if (num_events < 15)
    return('fewer than 15 events total')

  inner_join(x = exposures,
             y = events,
             by = 'pt_id') %>%
    mutate(surv_df = future_map2(.x = exposure,
                                 .y = event,
                                 .f = make_surv_df_)) %>%
    dplyr::select(pt_id, surv_df) %>%
    unnest(cols = surv_df) %>%
    left_join(baseline, by = 'pt_id')
}

make_surv_df_ <- function(exps, evs) {

  dj <- IRanges::disjoin(x = c(exps, evs), with.revmap = TRUE)
  at_risk <- dj[sapply(X = mcols(dj)$revmap, FUN = length) == 1, ]

  end(at_risk) <- end(at_risk) + 1

  if (is.null(mcols(exps))) {
    out_mcols <- NULL
  } else {
    out_mcols <- as.data.frame(
      mcols(exps)[unlist(mcols(at_risk)$revmap),,drop=FALSE]
    )
  }

  return(bind_cols(as.data.frame(at_risk),
                   list(event = end(at_risk) %in% start(evs)),
                   out_mcols))
}

get_coef <- function(ht) {
  summary(ht)$test$coef
}

get_se <- function(ht) {
  summary(ht)$test$sigma
}

get_pval <- function(ht) {
  summary(ht)$test$pvalue
}


#' @title Read in and process antibiotic exposure data
#'
#' @param filter_by inclusion criteria
#' @param classify_by list of quo()'s for case_when to classify abx's
#' @param pt_df patient data
#'
#' @return a tibble with columns pt_id, abx_class, on_abx, and off_abx where
#' pt_id is a factor encoding patient id,
#' abx_class is a factor encoding the class of abx,
#' on_abx is an IRanges of when the pt was on that abx_class,
#' and off_abx is the IRanges of when the pt was not on that abx_class
#'
read_abx <- function(
  filter_by = TRUE,
  classify_by = list(TRUE ~ 'abx'),
  pt_df = read_pt_data(),
  permanent_exposure = FALSE,
  lag_days = 0
) {

  filter_by <- enquo(filter_by)

  pt_id_and_dsc_day <- pt_df %>% dplyr::select(pt_id, dsc_day)

  # read in data
  read_abx_basic(pt_id_and_dsc_day = pt_id_and_dsc_day) %>%
    # filter and classify according to input
    filter(!!filter_by) %>%
    mutate(class = case_when(!!!classify_by)) %>%
    # tidy up and ensure all necessary rows present (even if no data)
    dplyr::select(-matches('abx_cat'), -abx_name) %>%
    complete(pt_id, class) %>%
    nest(abx_data = c(abx_start, abx_end)) %>%
    right_join(pt_id_and_dsc_day, by = 'pt_id') %>%
    # for each pt and exposure class, make an iranges
    mutate(one_exposure = future_pmap(.l = list(d = abx_data,
                                                dsc_day = dsc_day,
                                                out_prefix = class),
                                      .f = one_exposure_iranges_from_df,
                                      in_prefix = 'abx',
                                      permanent = permanent_exposure,
                                      lag_days = lag_days)) %>%
    # for each pt, combine exposure classes
    dplyr::select(pt_id, class, one_exposure) %>%
    nest(data = c(class, one_exposure)) %>%
    mutate(exposure = future_map(.x = data,
                                 .f = ~ combine_exposures(.$one_exposure))) %>%
    dplyr::select(-data)
}

read_abx_basic <- function(pt_id_and_dsc_day = read_pt_data()) {

  celgene_read_sas(fn = 'cm') %>%
    select(pt_id = RUSUBJID,
           abx_cat2 = CMATC2,
           abx_cat3 = CMATC3,
           abx_cat4 = CMATC4,
           abx_name = CMDECOD,
           abx_start = CMSTDY,
           abx_end = CMENDY) %>%
    # keep levels of ID column congruent across data tables for easy joins
    mutate(pt_id = factor(x = pt_id, levels = levels(pt_id_and_dsc_day$pt_id))) %>%
    # filter to abx with valid start and end
    filter(
      abx_cat2 == 'ANTIBACTERIALS FOR SYSTEMIC USE',
      abx_start <= abx_end,
      0 <= abx_end
    ) %>%
    # filter out abx that start after pt left study
    right_join(pt_id_and_dsc_day, by = 'pt_id') %>%
    filter(abx_start <= dsc_day) %>%
    # and immediately rj again to keep all pts in the df...or don't?
    # dplyr::select(-dsc_day) %>%
    # right_join(pt_id_and_dsc_day, by = 'pt_id') %>%
    # trim any abx periods that go too low or high
    mutate(
      abx_start = pmax(abx_start, 1),
      abx_end = pmin(abx_end, dsc_day)
    ) %>%
    select(-dsc_day)
}

one_exposure_iranges_from_df <- function(
  d,
  dsc_day,
  in_prefix = '',
  out_prefix = in_prefix,
  permanent = FALSE,
  lag_days = 0
) {

  if (all(is.na(d))) {

    ir <- IRanges(start = 0, end = dsc_day)
    mcols(ir)$placeholder <- FALSE
    names(mcols(ir)) <- paste('on', out_prefix, sep = '_')
    return(ir)

  } else {

    start_name <- paste(in_prefix, 'start', sep = '_')
    end_name <- paste(in_prefix, 'end', sep = '_')

    if (permanent) {
      on_ir <- IRanges(
        start = d[[start_name]],
        end = dsc_day
      )
    } else {
      on_ir <- IRanges(
        start = d[[start_name]],
        end = pmin(d[[end_name]] + lag_days, dsc_day)
      )
    }
    on_ir <- IRanges::reduce(on_ir)
    mcols(on_ir)$placeholder <- TRUE
    names(mcols(on_ir)) <- paste('on', out_prefix, sep = '_')

    off_ir <- IRanges::gaps(x = on_ir, start = 0, end = dsc_day)
    mcols(off_ir)$placeholder <- FALSE
    names(mcols(off_ir)) <- paste('on', out_prefix, sep = '_')

    return(sort(c(on_ir, off_ir)))
  }
}

combine_exposures <- function(exposures) {

  cd <- do.call(what = c, args = exposures)
  mc <- mcols(cd)
  dj <- disjoin(x = cd, with.revmap = TRUE)
  r <- mcols(dj)$revmap

  # idxs <- rep(0:(ncol(mc) - 1), times = nrow(mc)) * nrow(mc) + unlist(r)
  # as.matrix(mc)[idxs],

  idxs <- cbind(unlist(r), rep(x = 1:ncol(mc), times = length(dj)))

  mcols(dj) <- matrix(data = as.matrix(mc)[idxs],
                      nrow = length(dj),
                      byrow = TRUE,
                      dimnames = list(NULL, names(mc)))
  return(dj)
}

# ae_body_system = AEBODSYS,
# ae_high_level_term = AEHLTTX,
# ae_low_level_term = AELLTTX,
# ae_high_level_group_term = AEHLGTTX,
## ae_outcome = AEOUT,
# TODO: consider filtering on other things to classify groups of AE's,
# like body_system, etc
read_aes <- function(
  ae_codes = 'ae_any',
  mingrade = 3L,
  pt_df = read_pt_data()
) {

  stopifnot(mingrade %in% 1L:5L)

  pt_id_and_dsc_day <- pt_df %>% select(pt_id, dsc_day)

  # read in data, filter per arguments, and select useful bits
  read_aes_basic(pt_id_and_dsc_day = pt_id_and_dsc_day) %>%
    filter(ae_grade >= mingrade) %>%
    cond_filter(ae_code %in% ae_codes,
                execute = !identical(ae_codes, 'ae_any')) %>%
    # cond_filter(ae_code == filter_ae_class,
    #             execute = !is.null(filter_ae_class)) %>%
    select(pt_id, ae_start, ae_end, ae_code) ->
    aes_raw

  aes_raw %>%
    filter(!is.na(ae_start) & !is.na(ae_end)) ->
    aes_no_na

  # trim AEs to study period
  aes_no_na %>%
    complete(pt_id) %>%
    right_join(pt_id_and_dsc_day, by = 'pt_id') %>%
    mutate(ae_end = pmin(ae_end, dsc_day),
           ae_start = pmax(ae_start, 1)) %>%
    select(pt_id, ae_start, ae_end) %>%
    nest(data = c(ae_start, ae_end)) %>%
    mutate(event = future_map(.x = data,
                              .f = one_ae_iranges_from_df,
                              in_prefix = 'ae')) %>%
    select(pt_id, event)
}


read_aes_basic <- function(pt_id_and_dsc_day = read_pt_data()) {

  celgene_read_sas(fn = 'ae') %>%
    select(pt_id = RUSUBJID,
           ae_code = AEDECOD,
           ae_grade = AETOXGR,
           ae_body_system = AEBODSYS,
           ae_high_level_term = AEHLTTX,
           ae_low_level_term = AELLTTX,
           ae_high_level_group_term = AEHLGTTX,
           ae_start = AESTDY,
           ae_end = AEENDY) %>%
    mutate_if(is.numeric, as.integer) %>%
    mutate(pt_id = factor(x = pt_id,
                          levels = levels(pt_id_and_dsc_day$pt_id))) %>%
    filter(!(is.na(ae_start) & is.na(ae_end))) %>%
    right_join(pt_id_and_dsc_day, by = 'pt_id') %>%
    filter(is.na(ae_start) | ae_start <= dsc_day) %>%
    mutate(ae_end = pmin(ae_end, dsc_day))
}

one_ae_iranges_from_df <- function(
  d,
  in_prefix = '',
  out_prefix = in_prefix
) {

  if (all(is.na(d))) {

    return(IRanges())

  } else {

    start_name <- paste(in_prefix, 'start', sep = '_')
    end_name <- paste(in_prefix, 'end', sep = '_')

    ir <- IRanges(
      start = d[[start_name]],
      end = d[[end_name]]
    )

    return(IRanges::reduce(ir))
  }
}


read_pt_data <- function() {

  purrr::reduce(
    list(
      read_demog(),
      read_disc(),
      read_baseline_performance()
    ),
    left_join,
    by = 'pt_id'
  )
}

read_demog <- function() {
  celgene_read_sas(fn = 'dm') %>%
    select(pt_id = RUSUBJID,
           age = AGE,
           sex = SEX,
           race = RACEGEN) %>%
    mutate(pt_id = factor(x = pt_id),
           age = as.integer(age),
           sex = factor(x = sex, levels = c('F', 'M')),
           race = factor(x = race, levels = c('WHITE', 'BLACK', 'ASIAN', 'OTHER')))
}

read_disc <- function(pt_ids) {
  celgene_read_sas(fn = 'ds') %>%
    select(pt_id = RUSUBJID,
           dsc_day = DSDY) %>% #,
    # dsc_reason = DSTERM) %>%
    mutate(pt_id = factor(x = pt_id),
           dsc_day = as.integer(dsc_day))
}

read_baseline_performance <- function() {
  celgene_read_sas_factor(fn = 'kp') %>%
    filter(VISIT == 'BASE') %>%
    select(pt_id = RUSUBJID,
           kps = KPORRESN) %>%
    mutate(kps = as.integer(kps))
}

celgene_folder <- 'PDS_data/Celgene_132/CA046 data'

celgene_read_sas <- function(fn) {
  fn %>%
    paste0(., '.sas7bdat') %>%
    file.path(celgene_folder, .) %>%
    haven::read_sas()
}

celgene_read_sas_factor <- function(fn) {
  celgene_read_sas(fn = fn) %>%
    mutate_if(is.character, factor)
}

conditionally <- function(fun){
  function(first_arg, ..., execute){
    if(execute) return(fun(first_arg, ...))
    else return(first_arg)
  }
}

cond_filter <- conditionally(filter)
# cond_select <- conditionally(select)

get_os <- function() {
  if (.Platform$OS.type == "windows") {
    "win"
  } else if (Sys.info()["sysname"] == "Darwin") {
    "mac"
  } else if (.Platform$OS.type == "unix") {
    "unix"
  } else {
    stop("Unknown OS")
  }
}

params_to_string <- function(...) {
  dots <- substitute(list(...))[-1]
  paste(sapply(dots, deparse),
        sapply(list(...), paste0, collapse = ','),
        sep = '=',
        collapse = '_')
}
