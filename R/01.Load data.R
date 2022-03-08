# Import Library
library(tidyverse)
library(data.table)
library(VennDiagram)
library(lubridate)
library(gtsummary)


################################################################################# I ### Load data
path <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", "Sarcoma CH")

Demographic <- 
  readxl::read_xlsx(paste0(path, "/raw data/Breast+Sarcoma pts with DNA samples_PHI_01072020.xlsx"),
                    sheet = "PTE Demographics") %>% 
  janitor::clean_names()

sarcoma_DNA <- 
  readxl::read_xlsx(paste0(path, "/raw data/Sarcoma pts with all blood samples_02.01.22.xlsx"),
                    sheet = "CSV_Data_export_for_Biospec (2)") %>% 
  janitor::clean_names()

sarcoma_info <- 
  readxl::read_xlsx(paste0(path, "/raw data/Breast+Sarcoma pts with DNA samples_PHI_01072020.xlsx"),
                    sheet = "PTE_Cancer Characteristics") %>% 
  janitor::clean_names()

Chemot <- 
  readxl::read_xlsx(paste0(path, "/raw data/Breast+Sarcoma pts with DNA samples_PHI_01072020.xlsx"),
                    sheet = "Treatment_Chemotherapy") %>% 
  janitor::clean_names()

Immnunot <- 
  readxl::read_xlsx(paste0(path, "/raw data/Breast+Sarcoma pts with DNA samples_PHI_01072020.xlsx"),
                    sheet = "Treatment_Immuno") %>% 
  janitor::clean_names()

Radiot <- 
  readxl::read_xlsx(paste0(path, "/raw data/Breast+Sarcoma pts with DNA samples_PHI_01072020.xlsx"),
                    sheet = "Treatment_Radiation") %>% 
  janitor::clean_names()


# End loading