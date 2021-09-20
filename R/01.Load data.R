# Import Library
library(tidyverse)
library(data.table)
library(VennDiagram)
library(lubridate)
library(gtsummary)


################################################################################# I ### Load data
path <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", "Sarcoma CH")

Demographic <- 
  readxl::read_xlsx(paste0(path, "/raw data/Copy of Gillis_TransMed query_sarcoma pts_sent 08.16.21.xlsx"),
                    sheet = "PTE Demographics") %>% 
  janitor::clean_names()

sarcoma_DNA <- 
  readxl::read_xlsx(paste0(path, "/raw data/Copy of Gillis_TransMed query_sarcoma pts_sent 08.16.21.xlsx"),
                    sheet = "Biospecimens") %>% 
  janitor::clean_names()

sarcoma_info <- 
  readxl::read_xlsx(paste0(path, "/raw data/Copy of Gillis_TransMed query_sarcoma pts_sent 08.16.21.xlsx"),
                    sheet = "PTE CR") %>% 
  janitor::clean_names()

Chemot <- 
  readxl::read_xlsx(paste0(path, "/raw data/Copy of Gillis_TransMed query_sarcoma pts_sent 08.16.21.xlsx"),
                    sheet = "Treatment Chemotherapy") %>% 
  janitor::clean_names()

Immnunot <- 
  readxl::read_xlsx(paste0(path, "/raw data/Copy of Gillis_TransMed query_sarcoma pts_sent 08.16.21.xlsx"),
                    sheet = "Treatment Immunotherapy") %>% 
  janitor::clean_names()

Radiot <- 
  readxl::read_xlsx(paste0(path, "/raw data/Copy of Gillis_TransMed query_sarcoma pts_sent 08.16.21.xlsx"),
                    sheet = "Treatment Radiation") %>% 
  janitor::clean_names()


# End loading