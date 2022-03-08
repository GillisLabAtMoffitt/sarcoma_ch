# Table Sarcoma / Breast patients

blood_patients <- read_rds(paste0(here::here(), "/blood_patients.rds"))

breast_patients <- read_rds(paste0(here::here(), "/breast_patients.rds"))


preliminary_table <-
  bind_rows(blood_patients, breast_patients) %>% 
    mutate(primary_site_group1 = case_when(
      primary_site_group1 == "Breast"            ~ "Breast",
      primary_site_group1 != "Breast"            ~ "Sarcoma",
      is.na(primary_site_group1)                 ~ "Sarcoma"
    )) %>% 
    mutate(ethnicity_derived = case_when(
      ethnicity_derived == "Default code"        ~ NA_character_,
      ethnicity_derived == "Unknown"             ~ NA_character_,
      TRUE                                       ~ ethnicity_derived
    )#,
    # ethnicity_derived = coalesce(ethnicity_derived, ethnicity_cancer_registry, ethnicity_cerner)
    ) %>% 
    
  select(age_at_diagnosis, age_at_sample, year_at_sample,
         gender_cancer_registry,
         race_cancer_registry_1, 
         ethnicity_derived,
         primary_site_group1) %>% 
    tbl_summary(by= primary_site_group1,
                type = list(year_at_sample ~ "categorical")) %>% 
    bold_labels() %>% 
  as_gt()

gt::gtsave(preliminary_table, "preliminary table.pdf")
