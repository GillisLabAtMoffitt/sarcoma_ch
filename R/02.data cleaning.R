################################################################################# II ### Data cleaning
sarcoma_info <- 
  sarcoma_info %>% 
  distinct() %>% 
  arrange(deidentified_patient_id) %>% 
  left_join(., 
            Demographic, 
            by = "deidentified_patient_id") %>% 
  mutate(across(where(is.character), ~str_to_lower(.)))
write_rds(sarcoma_info, "sarcoma_info.rds")

# Chemot
Chemot <- Chemot %>% 
  filter(chemotherapy_type == "CHEMO NOS" |
           chemotherapy_type == "MULTI-AGENT CHEMO" |
           chemotherapy_type == "NONE, NOT PLANNED" | # clean more
           chemotherapy_type == "RECOMMENDED,UNKN IF GIVEN" | # clean more
           chemotherapy_type == "SINGLE-AGENT CHEMO" |
           chemotherapy_type == "UNKNOWN; DC ONLY"  # clean more
  ) %>% 
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  # Make it easier to not have na for future filtering, rescue the ones with a drug name
  mutate(chemotherapy_completion_status_first_course = case_when(
    !is.na(chemotherapy_drug) |
      !is.na(chemotherapy_start_date)                        ~ coalesce(chemotherapy_completion_status_first_course, "chemo given")
  )) %>% 
  
  # filtering
  # Remove NA in both drug name and date
  filter(!(is.na(chemotherapy_drug) & is.na(chemotherapy_start_date) & is.na(chemotherapy_end_date))) %>%  # remove 341
  # Remove no chemo given in chemotherapy_drug
  filter(chemotherapy_drug != "no chemo given" | # remove 137
           is.na(chemotherapy_drug)) %>% 
  # Remove no chemotherapy when 18.. or 2300 dates only but keep if real date
  filter(!(str_detect(chemotherapy_start_date, "12:00:00 am") & # remove 271
             chemotherapy_completion_status_first_course == "no chemotherapy")) %>% 
  # To help after removing the unknown date
  # mutate(drugs_unk = case_when(
  #   str_detect(chemotherapy_start_date, "12:00:00 am") &
  #     is.na(chemotherapy_drug)                              ~ "unk drug, 1800 date",
  #   !is.na(chemotherapy_start_date) &
  #     is.na(chemotherapy_drug)                              ~ "unk drug",
  #   is.na(chemotherapy_start_date) &
  #     is.na(chemotherapy_drug)                              ~ "unk drug, NA date"
  # )) %>% 
  mutate(data_unk = case_when(
    str_detect(chemotherapy_start_date, "12:00:00 am") &
      !is.na(chemotherapy_drug)                             ~ "known drug, 1800 date",
    !is.na(chemotherapy_start_date) &
      !is.na(chemotherapy_drug)                             ~ "known drug",
    is.na(chemotherapy_start_date) &
      !is.na(chemotherapy_drug)                             ~ "unk drug, NA date",
    str_detect(chemotherapy_start_date, "12:00:00 am") &
      is.na(chemotherapy_drug)                              ~ "unk drug, 1800 date",
    !is.na(chemotherapy_start_date) &
      is.na(chemotherapy_drug)                              ~ "unk drug, known date",
  )) %>% 
  # Dates
  # remove the data 18.., 23..
  mutate(chemotherapy_start_date = case_when(
    str_detect(chemotherapy_start_date, "12:00:00 am")       ~ NA_character_,
    TRUE                                                     ~ chemotherapy_start_date
  ), 
  chemotherapy_start_date = as.Date(as.numeric(chemotherapy_start_date), 
                                    origin = "1899-12-30")
  ) %>% 
  mutate(chemotherapy_end_date = case_when(
    str_detect(chemotherapy_end_date, "12:00:00 am")         ~ NA_character_,
    TRUE                                                     ~ chemotherapy_end_date
  ), 
  chemotherapy_end_date = as.Date(as.numeric(chemotherapy_end_date), 
                                  origin = "1899-12-30")
  ) %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  mutate(chemotherapy_start_date = 
           coalesce(chemotherapy_start_date, 
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>%  # 4,293 dates created
  # Fix the 2300 dates
  mutate(treatment_unk = case_when(
    str_detect(chemotherapy_start_date, "2300|2301")          ~ "Unknown when and if given 2300",
    TRUE                                                      ~ NA_character_
  )) %>% 
  mutate(chemotherapy_drug = coalesce(treatment_unk, chemotherapy_drug, data_unk)) %>% 
  mutate(chemotherapy_start_date = case_when(
    str_detect(chemotherapy_start_date, "2300|2301")          ~ as.Date("1700-01-01", origin = "1899-12-30"),
    TRUE                                                      ~ chemotherapy_start_date
  )) %>% 
  
  # Start pivot
  select(deidentified_patient_id, chemotherapy_drug, chemotherapy_start_date, chemotherapy_end_date) %>% 
  distinct() %>% 
  group_by(deidentified_patient_id, chemotherapy_start_date, chemotherapy_end_date) %>%
  summarise_at(vars(chemotherapy_drug), str_c, collapse = "; ") %>% 
  ungroup()

# pivot wider, use dcast bc better to keep date class
chemot <- dcast(setDT(Chemot), deidentified_patient_id ~ rowid(deidentified_patient_id),
                value.var = c(
                  "chemotherapy_drug",
                  "chemotherapy_start_date",
                  "chemotherapy_end_date"
                )
)

Chemot <- Chemot %>% 
  select(deidentified_patient_id, treatment_start_date = chemotherapy_start_date, 
         treatment_end_date = chemotherapy_end_date, treatment = chemotherapy_drug) %>% 
  mutate(treatment_type = "chemo")
# write_rds(Chemot, "Chemot.rds")

# Immnunot
Immnunot <- Immnunot %>% 
  filter(immunotherapy_type == "IMMUNO ADMINISTERED" |
           immunotherapy_type == "NONE, NOT PLANNED" |
           immunotherapy_type == "RECOMMENDED,UNKN IF GIVEN" | # clean more
           immunotherapy_type == "UNKNOWN; DC ONLY"  # clean more
  ) %>% 
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  mutate(immunotherapy_drug = case_when(
    is.na(immunotherapy_start_date) &
      !is.na(immunotherapy_drug)                             ~ "unk drug, NA date",
    str_detect(immunotherapy_start_date, "12:00:00 am") &
      is.na(immunotherapy_drug)                              ~ "unk drug, 1800 date",
    !is.na(immunotherapy_start_date) &
      is.na(immunotherapy_drug)                              ~ "unk drug, known date",
  )) %>% 
  filter(!is.na(immunotherapy_drug)) %>% 
  # Dates
  # remove the data 18.., 23..
  mutate(immunotherapy_start_date = case_when(
    str_detect(immunotherapy_start_date, "12:00:00 am")       ~ NA_character_,
    TRUE                                                     ~ immunotherapy_start_date
  ), 
  immunotherapy_start_date = as.Date(as.numeric(immunotherapy_start_date), 
                                     origin = "1899-12-30")
  ) %>% 
  mutate(immunotherapy_end_date = case_when(
    str_detect(immunotherapy_end_date, "12:00:00 am")         ~ NA_character_,
    TRUE                                                     ~ immunotherapy_end_date
  ), 
  immunotherapy_end_date = as.Date(as.numeric(immunotherapy_end_date),
                                   origin = "1899-12-30")
  ) %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  mutate(immunotherapy_start_date = 
           coalesce(immunotherapy_start_date, 
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>% # 3,647 dates created
  # No 2300 date
  
  # Start pivot
  select(deidentified_patient_id, immunotherapy_drug, immunotherapy_start_date, immunotherapy_end_date) %>% 
  distinct() %>% 
  group_by(deidentified_patient_id, immunotherapy_start_date, immunotherapy_end_date) %>%
  summarise_at(vars(immunotherapy_drug), str_c, collapse = "; ") %>% 
  ungroup()

immnunot <- dcast(setDT(Immnunot), deidentified_patient_id ~ rowid(deidentified_patient_id),
                  value.var = c(
                    "immunotherapy_drug",
                    "immunotherapy_start_date",
                    "immunotherapy_end_date"
                  )
)

Immnunot <- Immnunot %>% 
  select(deidentified_patient_id, treatment_start_date = immunotherapy_start_date, 
         treatment_end_date = immunotherapy_end_date, treatment = immunotherapy_drug) %>% 
  mutate(treatment_type = "immunoT")



# Radiot
Radiot <- Radiot %>% 
  filter(reason_for_no_radiation == "RAD THERAPY PERFORMED"| 
           reason_for_no_radiation == "RECOMMENDED,UNKN IF GIVEN" | 
           reason_for_no_radiation == "UNKNOWN; DC ONLY" |
           is.na(reason_for_no_radiation)) %>% 
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  
  mutate(data_unk = case_when(
    str_detect(radiation_start_date, "12:00:00 am") &
      !is.na(boost_dose_c_gy)                             ~ "known dose, 1800 date",
    !is.na(radiation_start_date) &
      !is.na(boost_dose_c_gy)                             ~ "known dose",
    is.na(radiation_start_date) &
      !is.na(boost_dose_c_gy)                             ~ "unk dose, NA date",
    str_detect(radiation_start_date, "12:00:00 am") &
      is.na(boost_dose_c_gy)                              ~ "unk dose, 1800 date",
    !is.na(radiation_start_date) &
      is.na(boost_dose_c_gy)                              ~ "unk dose, known date",
  )) %>% 
  filter(!is.na(data_unk)) %>% 
  # Dates
  # remove the data 18.., 23..
  mutate(radiation_start_date = case_when(
    str_detect(radiation_start_date, "12:00:00 am")       ~ NA_character_,
    TRUE                                                     ~ radiation_start_date
  ), 
  radiation_start_date = as.Date(as.numeric(radiation_start_date), 
                                 origin = "1899-12-30")
  ) %>% 
  mutate(radiation_end_date = case_when(
    str_detect(radiation_end_date, "12:00:00 am")         ~ NA_character_,
    TRUE                                                     ~ radiation_end_date
  ), 
  radiation_end_date = as.Date(as.numeric(radiation_end_date), ################################################ Fix bug create a 2300 date, lubridate::origin
                               origin = "1899-12-30")
  ) %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  mutate(radiation_start_date = 
           coalesce(radiation_start_date, 
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>% # 3,647 dates created
  # Fix the 2300 dates
  mutate(treatment_unk = case_when(
    str_detect(radiation_start_date, "2300|2301")             ~ "Unknown when and if given 2300",
    TRUE                                                      ~ NA_character_
  )) %>%
  mutate(boost_dose_c_gy = coalesce(treatment_unk, as.character(boost_dose_c_gy), data_unk)) %>% 
  mutate(radiation_start_date = case_when(
    str_detect(radiation_start_date, "2300|2301")             ~ as.Date("1700-01-01", origin = "1899-12-30"),
    TRUE                                                      ~ radiation_start_date
  )) %>%
  
  # Start pivot
  select(deidentified_patient_id, boost_dose_c_gy, radiation_start_date, radiation_end_date) %>% 
  distinct() %>% 
  group_by(deidentified_patient_id, radiation_start_date, radiation_end_date) %>%
  summarise_at(vars(boost_dose_c_gy), str_c, collapse = "; ") %>% 
  ungroup()


radiot <- dcast(setDT(Radiot), deidentified_patient_id ~ rowid(deidentified_patient_id),
                value.var = c(
                  "boost_dose_c_gy",
                  "radiation_start_date",
                  "radiation_end_date"
                )
)


Radiot <- Radiot %>% 
  select(deidentified_patient_id, treatment_start_date = radiation_start_date, 
         treatment_end_date = radiation_end_date, treatment = boost_dose_c_gy) %>% 
  mutate(treatment_type = "radioT")




# Combine
treatment <- bind_rows(Chemot, Immnunot, Radiot) %>% 
  arrange(deidentified_patient_id, treatment_start_date) %>% 
  group_by(deidentified_patient_id, treatment_type) %>% 
  mutate(treatment_line = row_number(deidentified_patient_id)) %>% 
  unite(treatment_line, c(treatment_type, treatment_line), sep = "_", remove = FALSE)

Treatment <- full_join(chemot, immnunot, by = "deidentified_patient_id") %>% 
  full_join(., radiot, by = "deidentified_patient_id") %>% 
  mutate(had_chemo = case_when(
    !is.na(chemotherapy_start_date_1) ~ "Yes"
  )) %>% 
  mutate(had_immuno = case_when(
    !is.na(immunotherapy_start_date_1) ~ "Yes"
  )) %>% 
  mutate(had_rad = case_when(
    !is.na(radiation_start_date_1) ~ "Yes"
  )) %>% 
  mutate(had_chemo_rad = case_when(
    !is.na(chemotherapy_start_date_1) &
      !is.na(radiation_start_date_1) ~ "Yes"
  )) %>% 
  mutate(had_treatment = case_when(
    !is.na(chemotherapy_start_date_1) &
      !is.na(immunotherapy_start_date_1) &
      !is.na(radiation_start_date_1)             ~ "Yes"
  ))

write_rds(Treatment, "Treatment.rds")



# DNA, clean and add same sample/same date on the same row
sarcoma_dna <- sarcoma_DNA %>% 
  filter(derived_tissue_type == "Blood",sample_type != "WBC/RBC") %>% 
  mutate(deidentified_patient_id = str_to_lower(deidentified_patient_id)) %>% 
  select(deidentified_patient_id, sample_family_id_sf, sample_id,
         specimen_collection_date) %>%
  mutate(specimen_collection_date = as.Date(specimen_collection_date, format = "%Y-%M-%d")) %>%
  arrange(deidentified_patient_id, specimen_collection_date) %>% 
  group_by(deidentified_patient_id, sample_family_id_sf, specimen_collection_date) %>% 
  summarise_at(vars(sample_id), str_c, collapse = "; ") %>%
  ungroup()

# write_rds(sarcoma_dna, "sarcoma_dna.rds")

sarcoma_dna %>% 
  distinct(deidentified_patient_id, specimen_collection_date, .keep_all = TRUE) %>% 
  group_by(deidentified_patient_id) %>% 
  mutate(sample_count = factor(row_number(deidentified_patient_id))) %>% 
  ungroup() %>% 
  arrange(desc(sample_count)) %>% 
  distinct(deidentified_patient_id, .keep_all = TRUE) %>% 
  select(sample_count) %>% 
  tbl_summary()




# sarcoma_dna1 <- sarcoma_dna %>% left_join(., Treatment, by = "deidentified_patient_id") %>% 
#   mutate(blood_bf_chemo = case_when(
#     specimen_collection_date <= chemotherapy_start_date_1                ~ "Yes",
#     specimen_collection_date > chemotherapy_start_date_1                 ~ "No",
#     is.na(chemotherapy_start_date_1)                                     ~ "not administred",
#     TRUE                                                            ~ NA_character_
#   )) %>% 
#   mutate(blood_bf_hormone = case_when(
#     specimen_collection_date <= hormone_therapy_start_date_1                ~ "Yes",
#     specimen_collection_date > hormone_therapy_start_date_1                 ~ "No",
#     is.na(hormone_therapy_start_date_1)                                     ~ "not administred",
#     TRUE                                                            ~ NA_character_
#   )) %>% 
#   mutate(blood_bf_immuno = case_when(
#     specimen_collection_date <= immunotherapy_start_date_1                ~ "Yes",
#     specimen_collection_date > immunotherapy_start_date_1                 ~ "No",
#     is.na(immunotherapy_start_date_1)                                     ~ "not administred",
#     TRUE                                                            ~ NA_character_
#   )) %>% 
#   mutate(blood_bf_rad = case_when(
#     specimen_collection_date <= radiation_start_date_1                ~ "Yes",
#     specimen_collection_date > radiation_start_date_1                 ~ "No",
#     is.na(radiation_start_date_1)                                     ~ "not administred",
#     TRUE                                                            ~ NA_character_
#   )) %>% 
#   mutate(blood_bf_chemo_rad = case_when(
#     specimen_collection_date <= chemotherapy_start_date_1 &
#       specimen_collection_date <= radiation_start_date_1                ~ "Yes",
#     specimen_collection_date > chemotherapy_start_date_1 |
#       specimen_collection_date > radiation_start_date_1                 ~ "No",
#     is.na(chemotherapy_start_date_1) |
#       is.na(radiation_start_date_1)                                 ~ "not administred",
#     TRUE                                                            ~ NA_character_
#   )) %>% 
#   mutate(blood_bf_treatment = case_when(
#     specimen_collection_date <= chemotherapy_start_date_1 &
#       specimen_collection_date <= hormone_therapy_start_date_1 &
#       specimen_collection_date <= immunotherapy_start_date_1 &
#       specimen_collection_date <= radiation_start_date_1                ~ "Yes",
#     if_any(contains("start_date_1"), ~ . < specimen_collection_date)    ~ "No",
#     # if_all(contains("start_date_1"), ~ . < specimen_collection_date)    ~ "Nope",
#     is.na(chemotherapy_start_date_1) |
#       is.na(hormone_therapy_start_date_1) |
#       is.na(immunotherapy_start_date_1) |
#       is.na(radiation_start_date_1)                                     ~ "not administred",
#     TRUE                                                                ~ NA_character_
#   )) %>% 
#   mutate(blood_bf_30_days_chemo = case_when(
#     specimen_collection_date >= (chemotherapy_start_date_1 - days(30)) &
#       specimen_collection_date <= (chemotherapy_start_date_1 + days(30))              ~ "Yes",
#     specimen_collection_date > (chemotherapy_start_date_1 + days(30))                 ~ "No",
#     is.na(chemotherapy_start_date_1)                                     ~ "not administred",
#     TRUE                                                            ~ NA_character_
#   )) %>%
#   mutate(blood_bf_30_days_hormone = case_when(
#     specimen_collection_date <= (hormone_therapy_start_date_1 + days(30))                ~ "Yes",
#     specimen_collection_date > (hormone_therapy_start_date_1 + days(30))                 ~ "No",
#     is.na(hormone_therapy_start_date_1)                                     ~ "not administred",
#     TRUE                                                            ~ NA_character_
#   )) %>%
#   mutate(blood_bf_30_days_immuno = case_when(
#     specimen_collection_date <= (immunotherapy_start_date_1 + days(30))                ~ "Yes",
#     specimen_collection_date > (immunotherapy_start_date_1 + days(30))                 ~ "No",
#     is.na(immunotherapy_start_date_1)                                     ~ "not administred",
#     TRUE                                                            ~ NA_character_
#   )) %>%
#   mutate(blood_bf_30_days_rad = case_when(
#     specimen_collection_date <= (radiation_start_date_1 + days(30))                ~ "Yes",
#     specimen_collection_date > (radiation_start_date_1 + days(30))                 ~ "No",
#     is.na(radiation_start_date_1)                                     ~ "not administred",
#     TRUE                                                            ~ NA_character_
#   )) %>%
#   mutate(blood_bf_30_days_chemo_rad = case_when(
#     specimen_collection_date <= (chemotherapy_start_date_1 + days(30)) &
#       specimen_collection_date <= (radiation_start_date_1 + days(30))                ~ "Yes",
#     specimen_collection_date > (chemotherapy_start_date_1 + days(30)) &
#       specimen_collection_date > (radiation_start_date_1 + days(30))                 ~ "No",
#     is.na(chemotherapy_start_date_1) |
#       is.na(radiation_start_date_1)                                 ~ "not administred",
#     TRUE                                                            ~ NA_character_
#   )) %>%
#   mutate(blood_bf_30_days_treatment = case_when(
#     specimen_collection_date <= (chemotherapy_start_date_1 + days(30)) &
#       specimen_collection_date <= (hormone_therapy_start_date_1 + days(30)) &
#       specimen_collection_date <= (immunotherapy_start_date_1 + days(30)) &
#       specimen_collection_date <= (radiation_start_date_1 + days(30))                ~ "Yes",
#     if_any(contains("start_date_1"), ~ (.  + days(30)) < specimen_collection_date)    ~ "No",
#     # if_all(contains("start_date_1"), ~ . < specimen_collection_date)    ~ "Nope",
#     is.na(chemotherapy_start_date_1) |
#       is.na(hormone_therapy_start_date_1) |
#       is.na(immunotherapy_start_date_1) |
#       is.na(radiation_start_date_1)                                     ~ "not administred",
#     TRUE                                                                ~ NA_character_
#   )) %>%
#   mutate(across(contains("blood_bf_"), ~ factor(., levels = c("Yes", "No", "not administred")))) %>% 
# arrange(deidentified_patient_id, specimen_collection_date, blood_bf_treatment) %>% 
#   group_by(deidentified_patient_id) %>% 
# 
# select(deidentified_patient_id, specimen_collection_date, #sample_lag, has_a_good_sample, has_a_good_seq_sample,
#        "blood_bf_chemo", "blood_bf_hormone", 
#        "blood_bf_immuno", "blood_bf_rad",
#        "blood_bf_chemo_rad",
#        "blood_bf_treatment", everything())
# 
# 
# sarcoma_dna1 %>% filter(chemotherapy_start_date_1 == "1700-01-01") %>% nrow()
# sarcoma_dna1 %>% filter(hormone_therapy_start_date_1 == "1700-01-01") %>% nrow()
# sarcoma_dna1 %>% filter(immunotherapy_start_date_1 == "1700-01-01") %>% nrow()
# sarcoma_dna1 %>% filter(radiation_start_date_1 == "1700-01-01") %>% nrow()
# sarcoma_dna1 %>% filter(chemotherapy_start_date_1 == "1700-01-01" & 
#                          radiation_start_date_1 == "1700-01-01") %>% nrow()
# sarcoma_dna1 %>% filter(chemotherapy_start_date_1 == "1700-01-01" & 
#                          hormone_therapy_start_date_1 == "1700-01-01" &
#                          immunotherapy_start_date_1 == "1700-01-01" &
#                          radiation_start_date_1 == "1700-01-01") %>% nrow()






# sarcoma_dna2 <- sarcoma_dna1 %>% 
#   mutate(had_good_sample_chemo = case_when(
#     blood_bf_chemo == "Yes"              ~ "Yes"
#   )) %>% 
#   group_by(deidentified_patient_id) %>% 
#   fill(had_good_sample_chemo, .direction = "updown") %>% 
#   mutate(seq_sample_chemo = case_when(
#     had_good_sample_chemo == "Yes" &
#       blood_bf_chemo == "No" ~ "Yes",
#     TRUE ~ "No"
#   )) %>% 
#   mutate(had_good_sample_hormone = case_when(
#     blood_bf_hormone == "Yes"              ~ "Yes"
#   )) %>% 
#   group_by(deidentified_patient_id) %>% 
#   fill(had_good_sample_hormone, .direction = "updown") %>% 
#   mutate(seq_sample_hormone = case_when(
#     had_good_sample_hormone == "Yes" &
#       blood_bf_hormone == "No" ~ "Yes",
#     TRUE ~ "No"
#   )) %>% 
#   mutate(had_good_sample_immuno = case_when(
#     blood_bf_immuno == "Yes"              ~ "Yes"
#   )) %>% 
#   group_by(deidentified_patient_id) %>% 
#   fill(had_good_sample_immuno, .direction = "updown") %>% 
#   mutate(seq_sample_immuno = case_when(
#     had_good_sample_immuno == "Yes" &
#       blood_bf_immuno == "No" ~ "Yes",
#     TRUE ~ "No"
#   )) %>% 
#   mutate(had_good_sample_rad = case_when(
#     blood_bf_rad == "Yes"              ~ "Yes"
#   )) %>% 
#   group_by(deidentified_patient_id) %>% 
#   fill(had_good_sample_rad, .direction = "updown") %>% 
#   mutate(seq_sample_rad = case_when(
#     had_good_sample_rad == "Yes" &
#       blood_bf_rad == "No" ~ "Yes",
#     TRUE ~ "No"
#   )) %>% 
#   mutate(had_good_sample_chemo_rad = case_when(
#     blood_bf_chemo_rad == "Yes"              ~ "Yes"
#   )) %>% 
#   group_by(deidentified_patient_id) %>% 
#   fill(had_good_sample_chemo_rad, .direction = "updown") %>% 
#   mutate(seq_sample_chemo_rad = case_when(
#     had_good_sample_chemo_rad == "Yes" &
#       blood_bf_chemo_rad == "No" ~ "Yes",
#     TRUE ~ "No"
#   )) %>% 
#   mutate(had_good_sample_treatment = case_when(
#     blood_bf_treatment == "Yes"              ~ "Yes"
#   )) %>% 
#   group_by(deidentified_patient_id) %>% 
#   fill(had_good_sample_treatment, .direction = "updown") %>% 
#   mutate(seq_sample_treatment = case_when(
#     had_good_sample_treatment == "Yes" &
#       blood_bf_treatment == "No" ~ "Yes",
#     TRUE ~ "No"
#   )) %>% 
#   ungroup()








################################################################################# III ### Merge data
Global_data <- full_join(sarcoma_info, sarcoma_dna, by = "deidentified_patient_id") %>% 
  # full_join(., sarcoma_info, by = "deidentified_patient_id") %>% 
  full_join(., Chemot, by = "deidentified_patient_id") %>% 
  full_join(., Immnunot, by = "deidentified_patient_id") %>% 
  full_join(., Radiot, by = "deidentified_patient_id")

write_rds(Global_data, "Global_data.rds")

blood_patients <- Global_data %>% 
  # filter to patients who have blood samples
  filter(!is.na(specimen_collection_date))

# End cleaning
