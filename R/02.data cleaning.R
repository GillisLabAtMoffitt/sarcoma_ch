################################################################################# II ### Data cleaning
# DNA, clean and add same sample/same date on the same row
sarcoma_dna <- sarcoma_DNA %>% 
  mutate(mrn = as.character(mrn),
         mrn = coalesce(mrn, party_id)) %>% 
  mutate_at(c("mrn"), ~str_to_lower(.)) %>% 
  # Select sample type needed for CH detection
  filter(collection_site_tissue_type == "Blood", 
         str_detect(sample_type, "Buffy|Genomic|Unprocessed|CD138|MNC$")) %>% 
  mutate(across(where(is.character), ~str_to_sentence(.))) %>% 
  select(mrn, party_id, sample_family_id, sample_id,
         specimen_collection_date, gender_cancer_registry) %>%
  # add same sample/same date on the same row
  arrange(mrn, specimen_collection_date) %>% 
  # Summarize to have 1 sample/day per row and not 1 row for each aliquot of the same sample 
  group_by(mrn, party_id, sample_family_id, specimen_collection_date, gender_cancer_registry) %>% 
  summarise_at(vars(sample_id), str_c, collapse = "; ") %>%
  # separate(col = sample_id, paste("sample_id", 1:3, sep="_"), sep = "; ", extra = "drop", fill = "right")
  ungroup()

write_rds(sarcoma_dna, "sarcoma_dna.rds")


# Cancer Characteristics----
sarcoma_info <- sarcoma_info %>%
  select(-c(group_name, approached_tcc_consent_oncore, is_active_tcc)) %>% 
  
  filter(!is.na(mrn)) %>% 
  mutate(mrn = as.character(mrn)#,
         # mrn = coalesce(mrn, party_id)
  ) %>% 
  mutate(across(where(is.character), ~str_to_sentence(.))) %>% 
  filter(!str_detect(primary_site_group, "Breast")) %>%
  
  mutate(date_of_diagnosis = as.Date(as.numeric(date_of_diagnosis), 
                                         origin = "1899-12-30")) %>% 
  filter(!is.na(date_of_diagnosis)) %>% 
  
  distinct(mrn, date_of_diagnosis, .keep_all = TRUE) %>% 
  arrange(mrn, date_of_diagnosis) %>% 
  # Summarize to have 1 row per patients
  group_by(mrn, #party_id, 
           # gender_cancer_registry, 
           date_of_birth#, 
           # gender_derived, gender_cerner
           ) %>%
  summarise_at(vars(date_of_diagnosis, tumor_sequence_number, 
                    clinical_tnm_group_stage, histology, 
                    tnm_cs_mixed_group_stage, primary_site_group), 
               str_c, collapse = "; ") %>% 
  ungroup() %>% 
  separate(date_of_diagnosis, paste("date_of_diagnosis", 1:10, sep = ""), 
           sep = "; ", remove = TRUE, 
           extra = "warn", fill = "right") %>% 
  separate(tnm_cs_mixed_group_stage, paste("tnm_cs_mixed_group_stage", 1:10, sep = ""), 
           sep = "; ", remove = TRUE, 
           extra = "warn", fill = "right") %>% 
  separate(primary_site_group, paste("primary_site_group", 1:10, sep = ""), 
           sep = "; ", remove = TRUE, 
           extra = "warn", fill = "right") %>% 
  purrr::keep(~!all(is.na(.)))

# sarcoma_info <- 
#   sarcoma_info %>% 
#   distinct() %>% 
#   arrange(deidentified_patient_id) %>% 
#   left_join(., 
#             Demographic, 
#             by = "deidentified_patient_id") %>% 
#   mutate(across(where(is.character), ~str_to_lower(.)))
# write_rds(sarcoma_info, "sarcoma_info.rds")


# Demographic----
Demographic <- Demographic %>%
  select(-c(group_name, approached_tcc_consent_oncore, is_active_tcc)) %>% 
  mutate(mrn = as.character(mrn),
         mrn = coalesce(mrn, party_id)) %>% 
  mutate(across(where(is.character), ~str_to_sentence(.)))


################################################################################# III ### Merge data
sarcoma_patients <- sarcoma_dna %>% 
  # Merge with Cancer Char for patients with samples available
  left_join(., sarcoma_info, by = c("mrn"#, "party_id"
                                    )) %>% 
  # Merge with Demographic for patients with samples available
  left_join(., Demographic, 
            c("mrn", "party_id", "date_of_birth"))

sarcoma_patients_id <- paste0(sarcoma_patients$mrn, collapse = "|")


################################################################################# IV ### Clean treatments

# Chemotherapy----
Chemot <- Chemot %>% 
  mutate(mrn = as.character(mrn),
         mrn = coalesce(mrn, party_id)) %>% 
  # Limit to Breast cancer patients
  filter(str_detect(mrn, sarcoma_patients_id)) %>% 
  # Transform number to date and 12:00:00 am character as NA
  mutate(chemotherapy_end_date = as.Date(as.numeric(chemotherapy_end_date), 
                                         origin = "1899-12-30")) %>% 
  mutate(chemotherapy_start_date = as.Date(as.numeric(chemotherapy_start_date),
                                           origin = "1899-12-30")) %>% 
  # Fix the 2300 dates
  mutate(chemotherapy_start_date = case_when(
    str_detect(chemotherapy_start_date, "2300")                   ~ NA_Date_,
    TRUE                                                          ~ chemotherapy_start_date
  )) %>% 
  mutate(chemotherapy_end_date = case_when(
    str_detect(chemotherapy_end_date, "2300")                     ~ NA_Date_,
    TRUE                                                          ~ chemotherapy_end_date
  )) %>% 
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  # Remove no chemo given in chemotherapy_drug
  filter(chemotherapy_drug != "no chemo given" | is.na(chemotherapy_drug)) %>% 
  # Remove drugs not for breast cancer
  # filter(!str_detect(chemotherapy_drug, "AG-013736|BRENTUXUMAB|DEPOCYT") |
  #        is.na(chemotherapy_drug)) %>% 
  mutate(chemotherapy_drug = case_when(
    chemotherapy_drug == "cisplatin - c1"                         ~ "cisplatin",
    chemotherapy_drug == "abraxane (form of taxol)"               ~ "paclitaxel",
    chemotherapy_drug %in% 
      c("niraparib", "ag-013736 (investigationa", 
        "brentuxumab", "depocyt")                                 ~ "non-breast",
    TRUE                                                          ~ chemotherapy_drug
  )) %>% 
  filter(chemotherapy_drug != "non-breast" | is.na(chemotherapy_drug)) %>% 
  # Remove NA in both drug name and date
  filter_at(vars(chemotherapy_drug, chemotherapy_start_date,
                 chemotherapy_end_date), any_vars(!is.na(.)))

Chemot1 <- Chemot %>% 
  # remove the chemo_type CONTRAINDICATED, PT DIED, RECOMMENDED, NOT GIVEN, "REFUSED"
  filter(chemotherapy_type == "chemo nos" |
           chemotherapy_type == "multi-agent chemo" |
           chemotherapy_type == "contraindicated" |
           chemotherapy_type == "none, not planned" | # clean more
           chemotherapy_type == "recommended,unkn if given" |
           chemotherapy_type == "single-agent chemo" |
           chemotherapy_type == "unknown; dc only"
  ) %>% 
  # Make it easier to not have na for future filtering, rescue the ones with a drug name
  mutate(chemotherapy_completion_status_first_course = case_when(
    !is.na(chemotherapy_drug) |
      !is.na(chemotherapy_start_date)         
    ~ coalesce(chemotherapy_completion_status_first_course, "chemo given")
  ))

Chemot1 <- Chemot1 %>% 
  # Clean contraindicated and none, not planned
  mutate(remove = case_when(
    chemotherapy_type == "contraindicated" &
      is.na(chemotherapy_drug)                                    ~ 1,
    chemotherapy_type == "none, not planned" &
      (is.na(chemotherapy_drug) & 
         is.na(chemotherapy_start_date))                          ~ 1,
    TRUE                                                          ~ 0
  )) %>% 
  filter(remove == 0) %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  mutate(chemotherapy_start_date =
           coalesce(chemotherapy_start_date,
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>%
  distinct(mrn, chemotherapy_drug, chemotherapy_start_date, chemotherapy_end_date, .keep_all = TRUE) 

Chemot <- Chemot1 %>%
  arrange(mrn, chemotherapy_start_date, chemotherapy_drug) %>% 
  # Combine drugs into regimen
  group_by(mrn, chemotherapy_start_date) %>%
  summarise_at(vars(chemotherapy_drug, chemotherapy_end_date, chemotherapy_type), 
               str_c, collapse = "; ") %>% 
  ungroup() %>% 
  separate(chemotherapy_end_date, paste("chemotherapy_end_date", 10:1, sep = ""), 
           sep = "; ", remove = TRUE, 
           extra = "warn", fill = "left") %>% 
  mutate(across(starts_with("chemotherapy_end_date"), ~ as.Date(.))) %>% 
  purrr::keep(~!all(is.na(.))) %>% 
  
  # Remove drugs that are before breast diagnosis
  # left_join(., 
  #           sarcoma_patients %>% 
  #             distinct(mrn, .keep_all = TRUE) %>% 
  #             select(mrn, date_of_diagnosis1), by= "mrn") %>% 
  # mutate(remove = case_when(
  #   chemotherapy_start_date < date_of_diagnosis1             ~ "remove",
  #   TRUE                                                     ~ NA_character_
  # )) %>% 
  # filter(is.na(remove)) %>% 
  
  distinct(mrn, chemotherapy_drug, chemotherapy_start_date, chemotherapy_end_date1, 
           .keep_all = TRUE) %>%  # 705976
  # select(-c(remove, date_of_diagnosis1, AC, AC_start_date, paclitaxel, pac_end_date)) %>% 
  mutate(linenumber = row_number())


# Chemot <- Chemot %>% 
#   filter(chemotherapy_type == "CHEMO NOS" |
#            chemotherapy_type == "MULTI-AGENT CHEMO" |
#            chemotherapy_type == "NONE, NOT PLANNED" | # clean more
#            chemotherapy_type == "RECOMMENDED,UNKN IF GIVEN" | # clean more
#            chemotherapy_type == "SINGLE-AGENT CHEMO" |
#            chemotherapy_type == "UNKNOWN; DC ONLY"  # clean more
#   ) %>% 
#   mutate(across(where(is.character), ~str_to_lower(.))) %>% 
#   # Make it easier to not have na for future filtering, rescue the ones with a drug name
#   mutate(chemotherapy_completion_status_first_course = case_when(
#     !is.na(chemotherapy_drug) |
#       !is.na(chemotherapy_start_date)                        ~ coalesce(chemotherapy_completion_status_first_course, "chemo given")
#   )) %>% 
#   
#   # filtering
#   # Remove NA in both drug name and date
#   filter(!(is.na(chemotherapy_drug) & is.na(chemotherapy_start_date) & is.na(chemotherapy_end_date))) %>%  # remove 341
#   # Remove no chemo given in chemotherapy_drug
#   filter(chemotherapy_drug != "no chemo given" | # remove 137
#            is.na(chemotherapy_drug)) %>% 
#   # Remove no chemotherapy when 18.. or 2300 dates only but keep if real date
#   filter(!(str_detect(chemotherapy_start_date, "12:00:00 am") & # remove 271
#              chemotherapy_completion_status_first_course == "no chemotherapy")) %>% 
#   # To help after removing the unknown date
#   # mutate(drugs_unk = case_when(
#   #   str_detect(chemotherapy_start_date, "12:00:00 am") &
#   #     is.na(chemotherapy_drug)                              ~ "unk drug, 1800 date",
#   #   !is.na(chemotherapy_start_date) &
#   #     is.na(chemotherapy_drug)                              ~ "unk drug",
#   #   is.na(chemotherapy_start_date) &
#   #     is.na(chemotherapy_drug)                              ~ "unk drug, NA date"
#   # )) %>% 
#   mutate(data_unk = case_when(
#     str_detect(chemotherapy_start_date, "12:00:00 am") &
#       !is.na(chemotherapy_drug)                             ~ "known drug, 1800 date",
#     !is.na(chemotherapy_start_date) &
#       !is.na(chemotherapy_drug)                             ~ "known drug",
#     is.na(chemotherapy_start_date) &
#       !is.na(chemotherapy_drug)                             ~ "unk drug, NA date",
#     str_detect(chemotherapy_start_date, "12:00:00 am") &
#       is.na(chemotherapy_drug)                              ~ "unk drug, 1800 date",
#     !is.na(chemotherapy_start_date) &
#       is.na(chemotherapy_drug)                              ~ "unk drug, known date",
#   )) %>% 
#   # Dates
#   # remove the data 18.., 23..
#   mutate(chemotherapy_start_date = case_when(
#     str_detect(chemotherapy_start_date, "12:00:00 am")       ~ NA_character_,
#     TRUE                                                     ~ chemotherapy_start_date
#   ), 
#   chemotherapy_start_date = as.Date(as.numeric(chemotherapy_start_date), 
#                                     origin = "1899-12-30")
#   ) %>% 
#   mutate(chemotherapy_end_date = case_when(
#     str_detect(chemotherapy_end_date, "12:00:00 am")         ~ NA_character_,
#     TRUE                                                     ~ chemotherapy_end_date
#   ), 
#   chemotherapy_end_date = as.Date(as.numeric(chemotherapy_end_date), 
#                                   origin = "1899-12-30")
#   ) %>% 
#   # Create a really early date to add to NAs date to make sure to exclude these patients
#   mutate(chemotherapy_start_date = 
#            coalesce(chemotherapy_start_date, 
#                     as.Date("1700-01-01", origin = "1899-12-30"))) %>%  # 4,293 dates created
#   # Fix the 2300 dates
#   mutate(treatment_unk = case_when(
#     str_detect(chemotherapy_start_date, "2300|2301")          ~ "Unknown when and if given 2300",
#     TRUE                                                      ~ NA_character_
#   )) %>% 
#   mutate(chemotherapy_drug = coalesce(treatment_unk, chemotherapy_drug, data_unk)) %>% 
#   mutate(chemotherapy_start_date = case_when(
#     str_detect(chemotherapy_start_date, "2300|2301")          ~ as.Date("1700-01-01", origin = "1899-12-30"),
#     TRUE                                                      ~ chemotherapy_start_date
#   )) %>% 
#   
#   # Start pivot
#   select(deidentified_patient_id, chemotherapy_drug, chemotherapy_start_date, chemotherapy_end_date) %>% 
#   distinct() %>% 
#   group_by(deidentified_patient_id, chemotherapy_start_date, chemotherapy_end_date) %>%
#   summarise_at(vars(chemotherapy_drug), str_c, collapse = "; ") %>% 
#   ungroup()

# pivot wider, use dcast bc better to keep date class
chemot <- dcast(setDT(Chemot), mrn ~ rowid(mrn),
                value.var = c(
                  "chemotherapy_drug",
                  "chemotherapy_start_date",
                  "chemotherapy_end_date1",
                  "chemotherapy_type"
                )) %>% 
  select(mrn, starts_with("chemotherapy_drug"),
         starts_with("chemotherapy_start"),
         starts_with("chemotherapy_end"),
         chemotherapy_type_1,
         chemotherapy_type_2)

# chemot <- dcast(setDT(Chemot), deidentified_patient_id ~ rowid(deidentified_patient_id),
#                 value.var = c(
#                   "chemotherapy_drug",
#                   "chemotherapy_start_date",
#                   "chemotherapy_end_date"
#                 )
# )
Chemot <- Chemot %>%
  select(mrn, treatment_start_date = chemotherapy_start_date,
         treatment_end_date = chemotherapy_end_date1, treatment = chemotherapy_drug) %>%
  mutate(treatment_type = "chemo")

# write_rds(Chemot, "Chemot.rds")

# Immnunotherapy----
Immnunot <- Immnunot %>% 
  mutate(mrn = as.character(mrn),
         mrn = coalesce(mrn, party_id)) %>% 
  # Limit to Breast cancer patients
  filter(str_detect(mrn, sarcoma_patients_id)) %>% 
  # Transform number to date and 12:00:00 am character as NA
  mutate(across(ends_with("date"), ~ as.Date(as.numeric(.), 
                                             origin = "1899-12-30"))) %>% 
  # No 2300 date
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  # Remove no hormone given in immunotherapy_drug
  # filter(immunotherapy_drug != "no  given" | is.na(immunotherapy_drug)) %>% 
  # remove rows with no drugs info
  filter_at(vars(immunotherapy_drug, immunotherapy_start_date,
                 immunotherapy_end_date), any_vars(!is.na(.)))

Immnunot1 <- Immnunot %>% 
  # remove the chemo_type CONTRAINDICATED, PT DIED, RECOMMENDED, NOT GIVEN, "REFUSED"
  filter(immunotherapy_type == "immuno administered" |
           immunotherapy_type == "none, not planned" |
           immunotherapy_type == "contraindicated" |
           immunotherapy_type == "unknown; dc only" |
           immunotherapy_type == "recommended, unk if given"
  )

Immnunot1 <- Immnunot1 %>% 
  # Clean contraindicated and none, not planned
  mutate(remove = case_when(
    immunotherapy_type == "contraindicated" &
      (is.na(immunotherapy_drug) & 
         is.na(immunotherapy_start_date))                        ~ 1,
    immunotherapy_type == "none, not planned" &
      (is.na(immunotherapy_drug) & 
         is.na(immunotherapy_start_date))                        ~ 1,
    TRUE                                                         ~ 0
  )) %>% 
  filter(remove == 0) %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  mutate(immunotherapy_start_date =
           coalesce(immunotherapy_start_date,
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>%
  distinct(mrn, immunotherapy_drug, immunotherapy_start_date, immunotherapy_end_date, .keep_all = TRUE)

Immnunot <- Immnunot1 %>%
  arrange(mrn, immunotherapy_start_date, immunotherapy_drug) %>% 
  group_by(mrn, immunotherapy_start_date) %>%
  summarise_at(vars(immunotherapy_drug, immunotherapy_end_date), str_c, collapse = "; ") %>% 
  ungroup() %>% 
  separate(immunotherapy_end_date, paste("immunotherapy_end_date", 10:1, sep = ""), 
           sep = ";", remove = TRUE, 
           extra = "warn", fill = "left") %>% 
  mutate(across(starts_with("immunotherapy_end_date"), ~ as.Date(.))) %>% 
  purrr::keep(~!all(is.na(.))) #%>% 
  # Remove drugs that are before breast diagnosis
  # left_join(., 
  #           breast_patients %>% 
  #             distinct(mrn, .keep_all = TRUE) %>% 
  #             select(mrn, date_of_diagnosis1), by= "mrn") %>% 
  # mutate(remove = case_when(
  #   immunotherapy_start_date < date_of_diagnosis1             ~ "remove",
  #   TRUE                                                     ~ NA_character_
  # )) %>% 
  # filter(is.na(remove)) %>% 
  # select(-remove, -date_of_diagnosis1)


# Immnunot <- Immnunot %>% 
#   filter(immunotherapy_type == "IMMUNO ADMINISTERED" |
#            immunotherapy_type == "NONE, NOT PLANNED" |
#            immunotherapy_type == "RECOMMENDED,UNKN IF GIVEN" | # clean more
#            immunotherapy_type == "UNKNOWN; DC ONLY"  # clean more
#   ) %>% 
#   mutate(across(where(is.character), ~str_to_lower(.))) %>% 
#   mutate(immunotherapy_drug = case_when(
#     is.na(immunotherapy_start_date) &
#       !is.na(immunotherapy_drug)                             ~ "unk drug, NA date",
#     str_detect(immunotherapy_start_date, "12:00:00 am") &
#       is.na(immunotherapy_drug)                              ~ "unk drug, 1800 date",
#     !is.na(immunotherapy_start_date) &
#       is.na(immunotherapy_drug)                              ~ "unk drug, known date",
#   )) %>% 
#   filter(!is.na(immunotherapy_drug)) %>% 
#   # Dates
#   # remove the data 18.., 23..
#   mutate(immunotherapy_start_date = case_when(
#     str_detect(immunotherapy_start_date, "12:00:00 am")       ~ NA_character_,
#     TRUE                                                     ~ immunotherapy_start_date
#   ), 
#   immunotherapy_start_date = as.Date(as.numeric(immunotherapy_start_date), 
#                                      origin = "1899-12-30")
#   ) %>% 
#   mutate(immunotherapy_end_date = case_when(
#     str_detect(immunotherapy_end_date, "12:00:00 am")         ~ NA_character_,
#     TRUE                                                     ~ immunotherapy_end_date
#   ), 
#   immunotherapy_end_date = as.Date(as.numeric(immunotherapy_end_date),
#                                    origin = "1899-12-30")
#   ) %>% 
#   # Create a really early date to add to NAs date to make sure to exclude these patients
#   mutate(immunotherapy_start_date = 
#            coalesce(immunotherapy_start_date, 
#                     as.Date("1700-01-01", origin = "1899-12-30"))) %>% # 3,647 dates created
#   # No 2300 date
#   
#   # Start pivot
#   select(deidentified_patient_id, immunotherapy_drug, immunotherapy_start_date, immunotherapy_end_date) %>% 
#   distinct() %>% 
#   group_by(deidentified_patient_id, immunotherapy_start_date, immunotherapy_end_date) %>%
#   summarise_at(vars(immunotherapy_drug), str_c, collapse = "; ") %>% 
#   ungroup()


immnunot <- dcast(setDT(Immnunot), mrn ~ rowid(mrn),
                  value.var = c(
                    "immunotherapy_drug",
                    "immunotherapy_start_date",
                    "immunotherapy_end_date1"
                  ))

Immnunot <- Immnunot %>%
  select(mrn, treatment_start_date = immunotherapy_start_date,
         treatment_end_date = immunotherapy_end_date1, treatment = immunotherapy_drug) %>%
  mutate(treatment_type = "immunoT")


# Radiotion----
Radiot <- Radiot %>% 
  mutate(mrn = as.character(mrn),
         mrn = coalesce(mrn, party_id)) %>% 
  # Limit to Breast cancer patients
  filter(str_detect(mrn, sarcoma_patients_id)) %>% 
  # Transform number to date and 12:00:00 am character as NA
  mutate(across(ends_with("date"), ~ as.Date(as.numeric(.), 
                                             origin = "1899-12-30"))) %>% 
  # Fix the 2300 dates
  mutate(radiation_start_date = case_when(
    str_detect(radiation_start_date, "2300")                      ~ NA_Date_,
    TRUE                                                          ~ radiation_start_date
  )) %>% 
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  # remove rows with no rad info
  filter(!is.na(boost_dose_c_gy) | !is.na(radiation_start_date),
         !is.na(radiation_end_date) | 
           radiation_location_of_rx != "no radiation therapy")

Radiot1 <- Radiot %>% 
  # remove the PT DIED, RECOMMENDED, NOT GIVEN, "REFUSED"
  filter(reason_for_no_radiation == "rad therapy performed" |
           # reason_for_no_radiation == "not recommended/autopsy" |
           reason_for_no_radiation == "recommended, unk if given" |
           reason_for_no_radiation == "radiation contraindicated" |
           reason_for_no_radiation == "unknown/dco"
  )

Radiot1 <- Radiot1 %>% 
  # Clean contraindicated and none, not planned
  mutate(remove = case_when(
    reason_for_no_radiation == "radiation contraindicated" &
      is.na(radiation_start_date)                                 ~ 1,
    TRUE                                                          ~ 0
  )) %>% 
  filter(remove == 0) %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  mutate(radiation_start_date =
           coalesce(radiation_start_date,
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>%
  distinct(mrn, boost_dose_c_gy, radiation_start_date, radiation_end_date, .keep_all = TRUE)

Radiot <- Radiot1 %>%
  arrange(mrn, radiation_start_date) %>% 
  group_by(mrn, radiation_start_date) %>%
  summarise_at(vars(boost_dose_c_gy, radiation_end_date), str_c, collapse = "; ") %>% 
  ungroup() %>% 
  separate(radiation_end_date, paste("radiation_end_date", 10:1, sep = ""), 
           sep = ";", remove = TRUE, 
           extra = "warn", fill = "left") %>% 
  mutate(across(starts_with("radiation_end_date"), ~ as.Date(.))) %>% 
  purrr::keep(~!all(is.na(.))) #%>% 
  # Remove rad that are before breast diagnosis
  # left_join(., 
  #           breast_patients %>% 
  #             distinct(mrn, .keep_all = TRUE) %>% 
  #             select(mrn, date_of_diagnosis1), by= "mrn") %>% 
  # mutate(remove = case_when(
  #   radiation_start_date < date_of_diagnosis1                 ~ "remove",
  #   TRUE                                                     ~ NA_character_
  # )) %>% 
  # filter(is.na(remove)) %>% 
  # select(-remove, -date_of_diagnosis1)


# Radiot <- Radiot %>% 
#   filter(reason_for_no_radiation == "RAD THERAPY PERFORMED"| 
#            reason_for_no_radiation == "RECOMMENDED,UNKN IF GIVEN" | 
#            reason_for_no_radiation == "UNKNOWN; DC ONLY" |
#            is.na(reason_for_no_radiation)) %>% 
#   mutate(across(where(is.character), ~str_to_lower(.))) %>% 
#   
#   mutate(data_unk = case_when(
#     str_detect(radiation_start_date, "12:00:00 am") &
#       !is.na(boost_dose_c_gy)                             ~ "known dose, 1800 date",
#     !is.na(radiation_start_date) &
#       !is.na(boost_dose_c_gy)                             ~ "known dose",
#     is.na(radiation_start_date) &
#       !is.na(boost_dose_c_gy)                             ~ "unk dose, NA date",
#     str_detect(radiation_start_date, "12:00:00 am") &
#       is.na(boost_dose_c_gy)                              ~ "unk dose, 1800 date",
#     !is.na(radiation_start_date) &
#       is.na(boost_dose_c_gy)                              ~ "unk dose, known date",
#   )) %>% 
#   filter(!is.na(data_unk)) %>% 
#   # Dates
#   # remove the data 18.., 23..
#   mutate(radiation_start_date = case_when(
#     str_detect(radiation_start_date, "12:00:00 am")       ~ NA_character_,
#     TRUE                                                     ~ radiation_start_date
#   ), 
#   radiation_start_date = as.Date(as.numeric(radiation_start_date), 
#                                  origin = "1899-12-30")
#   ) %>% 
#   mutate(radiation_end_date = case_when(
#     str_detect(radiation_end_date, "12:00:00 am")         ~ NA_character_,
#     TRUE                                                     ~ radiation_end_date
#   ), 
#   radiation_end_date = as.Date(as.numeric(radiation_end_date), ################################################ Fix bug create a 2300 date, lubridate::origin
#                                origin = "1899-12-30")
#   ) %>% 
#   # Create a really early date to add to NAs date to make sure to exclude these patients
#   mutate(radiation_start_date = 
#            coalesce(radiation_start_date, 
#                     as.Date("1700-01-01", origin = "1899-12-30"))) %>% # 3,647 dates created
#   # Fix the 2300 dates
#   mutate(treatment_unk = case_when(
#     str_detect(radiation_start_date, "2300|2301")             ~ "Unknown when and if given 2300",
#     TRUE                                                      ~ NA_character_
#   )) %>%
#   mutate(boost_dose_c_gy = coalesce(treatment_unk, as.character(boost_dose_c_gy), data_unk)) %>% 
#   mutate(radiation_start_date = case_when(
#     str_detect(radiation_start_date, "2300|2301")             ~ as.Date("1700-01-01", origin = "1899-12-30"),
#     TRUE                                                      ~ radiation_start_date
#   )) %>%
#   
#   # Start pivot
#   select(deidentified_patient_id, boost_dose_c_gy, radiation_start_date, radiation_end_date) %>% 
#   distinct() %>% 
#   group_by(deidentified_patient_id, radiation_start_date, radiation_end_date) %>%
#   summarise_at(vars(boost_dose_c_gy), str_c, collapse = "; ") %>% 
#   ungroup()


radiot <- dcast(setDT(Radiot), mrn ~ rowid(mrn),
                value.var = c(
                  "boost_dose_c_gy",
                  "radiation_start_date",
                  "radiation_end_date1"
                )
)

Radiot <- Radiot %>%
  select(mrn, treatment_start_date = radiation_start_date,
         treatment_end_date = radiation_end_date1, treatment = boost_dose_c_gy) %>%
  mutate(treatment_type = "radioT")



# Combine
treatment <- bind_rows(Chemot, Immnunot, Radiot) %>% 
  arrange(mrn, treatment_start_date) %>% 
  group_by(mrn, treatment_type) %>% 
  mutate(treatment_line = row_number(mrn)) %>% 
  unite(treatment_line, c(treatment_type, treatment_line), sep = "_", remove = FALSE)

write_rds(treatment, "treatment_long.rds")

Treatment <- full_join(chemot, immnunot, by = "mrn") %>% 
  full_join(., radiot, by = "mrn") %>% 
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


################################################################################# III ### Merge data
Global_data <- 
  left_join(sarcoma_patients, Treatment, by = "mrn") %>% 
  # Create deidentified IDs
  mutate(rad = "sarcoma_study_") %>%
  group_by(mrn) %>% 
  mutate(id = cur_group_id()) %>%
  ungroup() %>%
  mutate(zero = 6 - nchar(id)) %>%
  mutate(ii = stringi::stri_dup("0", zero)) %>%
  select(c(rad, ii, id, mrn, everything())) %>%
  unite(deidentified_patient_id, rad:id, sep = "") %>% 
  select(c(deidentified_patient_id, mrn, everything(), -zero))

write_rds(Global_data, "Global_data.rds")

blood_patients <- Global_data %>% 
  # filter to patients who have blood samples
  filter(!is.na(specimen_collection_date))

write_rds(blood_patients, "blood_patients.rds")


# End cleaning
