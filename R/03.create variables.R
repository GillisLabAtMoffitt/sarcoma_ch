## WILL USE WHEN GET THE SEQUENCING DATA

blood_patients <- read_rds(paste0(here::here(), "/blood_patients.rds"))


blood_patients <- blood_patients %>%
  # create age
  mutate(age_at_diagnosis = round(interval(start = date_of_birth, end = date_of_diagnosis1)/
           duration(n = 1, units = "years"), 1)
         ) %>%
  mutate(age_at_sample = round(interval(start = date_of_birth, end = specimen_collection_date)/
                                 duration(n = 1, units = "years"), 1)
         ) %>% 
  mutate(year_at_sample = year(specimen_collection_date))
  
write_rds(blood_patients, "blood_patients.rds")


# End create variables
