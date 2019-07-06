aes_basic <- read_aes_basic()

dplyr::filter(aes_basic, !is.na(ae_start) & !is.na(ae_end)) %>%
  count(ae_code) %>%
  arrange(desc(n)) %>%
  dplyr::slice(1:20) ->
  common_ae_table

names_of_common_aes <- as.character(common_ae_table$ae_code)
names_of_common_aes <- set_names(x = names_of_common_aes,
                                 nm = paste('aecode', names_of_common_aes, sep = '_'))

classes_by_common_ae <- paste('aeclass',
                              c('hem', 'hem', 'gi', 'hem', 'const',
                                'const', 'gi', 'gi', 'hep', 'gi',
                                'gi', 'hem', 'hem', 'hep', 'gi',
                                'const', 'hep', 'hem', 'pulm', 'const'),
                              sep = '_')

classes_of_common_aes <- sort(unique(classes_by_common_ae))
classes_of_common_aes <- set_names(x = classes_of_common_aes,
                                   nm = classes_of_common_aes)

common_ae_table %>%
  mutate(class = classes_by_common_ae,
         ae_code_american = recode(ae_code,
                                   Diarrhoea = 'Diarrhea',
                                   Anaemia = 'Anemia',
                                   Pyrexia = 'Fever',
                                   `Oedema peripheral` = 'Edema',
                                   `Alanine aminotransferase increased` = 'ALT elevation',
                                   `Aspartate aminotransferase increased` = 'AST elevation',
                                   Asthenia = 'weakness',
                                   Dyspnoea = 'dyspnea')) ->
  common_ae_table

common_ae_table %>%
  group_by(class) %>%
  summarise(codes = list(ae_code)) %>%
  pull(codes) %>%
  set_names(nm = sort(classes_of_common_aes)) ->
  common_ae_codes_by_class

# TODO: make color scheme for abx exposures
# eg black for abx_any
# blue for this class and variants of blue for all specifcis in that...

ae_class_levels <- c("aeclass_hem", "aeclass_gi", "aeclass_const", "aeclass_hep", 'aeclass_pulm')
ae_code_levels <- c("aecode_Neutropenia", "aecode_Thrombocytopenia",
                    "aecode_Nausea", "aecode_Anaemia", "aecode_Pyrexia", "aecode_Fatigue",
                    "aecode_Vomiting", "aecode_Diarrhoea", "aecode_Oedema peripheral",
                    "aecode_Constipation", "aecode_Abdominal pain", "aecode_Leukopenia",
                    "aecode_Platelet count decreased", "aecode_Alanine aminotransferase increased",
                    "aecode_Decreased appetite", "aecode_Asthenia", "aecode_Aspartate aminotransferase increased",
                    "aecode_Neutrophil count decreased", "aecode_Dyspnoea")

ae_class_names_for_plotting <- c("hematologic", "gastrointestinal",
                                 "constitutional", "hepatologic",
                                 'pulmonary')
ae_class_names_for_plotting_short <- c("hematologic", "gastro.", "const.", "hep.", 'pulm.')
ae_code_names_for_plotting <- c("neutropenia", "thrombocytopenia", "nausea", "anemia",
                                'fever', 'fatigue', 'vomiting', 'diarrhea', 'edema',
                                'constipation', 'abdominal pain', 'leukopenia', 'thrombocytopenia',
                                'ALT elevation', 'loss of appetite', 'weakness', 'AST elevation',
                                'neutropenia', 'dyspnea')

levels_of_events_for_plotting <-
  c("ae_any",
    ae_class_levels,
    ae_code_levels)



names_of_events_for_plotting <-
  c("any adverse event",
    ae_class_names_for_plotting,
    ae_code_names_for_plotting)

names_of_hypothesis_for_plotting <-
  c('any abx', 'penicillin abx', 'quinolone abx', 'cipro', 'pen with BLI',
    'other abx', 'non-pen. beta lactam', '1st gen cephalospin')
