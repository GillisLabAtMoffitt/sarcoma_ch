library(tidyverse)
library(gtsummary)

# Table Sarcoma / Breast patients

blood_patients <- read_rds(paste0(here::here(), "/blood_patients.rds"))

breast_patients <- read_rds(paste0(here::here(), "/breast_patients.rds"))


preliminary_table <-
  bind_rows(blood_patients, breast_patients) %>% 
  distinct(mrn, .keep_all = TRUE) %>% 
    mutate(primary_site_group1 = case_when(
      primary_site_group1 == "Breast"            ~ "Breast",
      primary_site_group1 != "Breast"            ~ "Sarcoma",
      is.na(primary_site_group1)                 ~ "Sarcoma"
    )) %>% 
  mutate(gender_at_birth = case_when(
    gender_cancer_registry == "Transsexual, natal male"    ~ "Male",
    TRUE                                                   ~ gender_cancer_registry
  )) %>% 
  mutate(race_cancer_registry_1 = case_when(
    str_detect(race_cancer_registry_1, "Asian")          ~ "Asian",
    str_detect(race_cancer_registry_1, "indian|Indian")          ~ "Am Indian",
    race_cancer_registry_1 == "Black\vother"                 ~ "Black",
    race_cancer_registry_1 == "Black\vunknown"                 ~ "Black",
    race_cancer_registry_1 == "Black\vwhite"                 ~ NA_character_,
    race_cancer_registry_1 == "Chamorran"                 ~ "Others",
    race_cancer_registry_1 == "Chinese"                 ~ "Asian",
    race_cancer_registry_1 == "Filipino"                 ~ "Asian",
    race_cancer_registry_1 == "Guamanian nos"                 ~ "Others",
    race_cancer_registry_1 == "Japanese"                 ~ "Asian",
    race_cancer_registry_1 == "Korean"                 ~ "Asian",
    str_detect(race_cancer_registry_1, "Laotian")                 ~ "Asian",
    race_cancer_registry_1 == "Other"                 ~ "Others",
    str_detect(race_cancer_registry_1, "islander")                 ~ "Others",
    race_cancer_registry_1 == "Other\vunknown"                 ~ NA_character_,
    race_cancer_registry_1 == "Other\vwhite"                 ~ "Others",
    race_cancer_registry_1 == "Pakistani"                 ~ "Asian",
    race_cancer_registry_1 == "Samoan"                 ~ "Others",
    race_cancer_registry_1 == "Thai"                 ~ "Asian",
    race_cancer_registry_1 == "Unknown"                 ~ NA_character_,
    race_cancer_registry_1 == "Unknown\vwhite"                 ~ NA_character_,
    race_cancer_registry_1 == "Vietnamese"                 ~ "Asian",
    TRUE            ~ race_cancer_registry_1
    
  )) %>% 
  mutate(ethnicity_derived = case_when(
    ethnicity_derived == "Default code"        ~ NA_character_,
    ethnicity_derived == "Unknown"             ~ NA_character_,
    TRUE                                       ~ ethnicity_derived
  )) %>% 
  mutate(ethnicity = coalesce(ethnicity_derived, ethnicity_cancer_registry, ethnicity_cerner
  ),
  ethnicity = case_when(
    ethnicity == "Cuban"                             ~ "Hispanic",
    ethnicity == "Dominican republic"                             ~ "Hispanic",
    ethnicity == "Hispanic/latino"                             ~ "Hispanic",
    ethnicity == "Mexican"                             ~ "Hispanic",
    ethnicity == "Puerto rican"                             ~ "Hispanic",
    ethnicity == "South/central american"                             ~ "Hispanic",
    ethnicity == "Spanish nos"                             ~ "Hispanic",
    ethnicity == "Spanish; hispanic"                             ~ "Hispanic",
    ethnicity == "Non-hispanic/non-latino"                             ~ "Non-Hispanic",
    str_detect(ethnicity, "Non-spanish")                             ~ "Non-Hispanic",
    ethnicity == "Missing/invalid answer"                             ~ NA_character_,
    ethnicity == "Multiple"                             ~ NA_character_,
    ethnicity == "Prefer not to answer"                             ~ NA_character_,
    ethnicity == "Pt not present"                             ~ NA_character_,
    ethnicity == "Unknown"                             ~ NA_character_,
    TRUE                             ~ ethnicity,
    
  )) %>% 
    
  select("Age at diagnosis" = age_at_diagnosis, 
         "Age at sample collection" = age_at_sample, 
         "Year of sample collection" = year_at_sample,
         "Gender at birth" = gender_at_birth, 
         "Race" = race_cancer_registry_1,
         "Ethnicity" = ethnicity,
         primary_site_group1) %>% 
    tbl_summary(by= primary_site_group1,
                sort = list(
                  c("Gender at birth", "Race", "Ethnicity")
                  ~ "frequency"),
                type = list(`Year of sample collection` ~ "categorical")) %>%
    as_hux_table()

write_csv(preliminary_table, "preliminary table breast and sarcoma patients.csv")
gt::gtsave(preliminary_table, "preliminary table breast and sarcoma patients.pdf", 
           zoom = 1)
