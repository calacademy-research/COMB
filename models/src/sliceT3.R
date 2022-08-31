#top the hack
#slice top three
# View(dataML)
library(tidyverse)
library(here)
library(data.table)
library(googledrive)
library(lubridate)

############# OPTIONS FOR SCRIPT ###################
here <- here()

# dataMLt <- fread(paste0(here, "/acoustic/data_ingest/output/dataML_tall.csv"))
# 
# sliceT3 <- function(data, foc_spp, headtt = 10^4){
#   data %>%
#   head(headtt) %>%
#   group_by(Date_Time, point, Start_Time) %>%
#   slice_max(order_by = c("species","logit"), n = 3) %>%
#   filter(species == foc_spp) -> tste #assign(x = paste0("dataMLslicedT3_foc_spp_",headtt,"_",Sys.Date(),"_out"), value = ., envir = .GlobalEnv)
#   }
# 
# sliceT3(dataMLt,"HAWO")

# drive_sync(here("acoustic/data_ingest/output/"), "https://drive.google.com/drive/folders/1eOrXsDmiIW9YqJWrlUWR9-Cgc7hHKD_5")
dataML_top3 <- fread(paste0(here, "/acoustic/data_ingest/output/dataML_top3.csv"))

?system.time

system.time(dataML_top1 <- dataML_top3 %>% 
  group_by(Date_Time, point, Start_Time) %>% 
  slice_max(order_by = logit, n = 1) %>%
  ungroup())

#about 18 minutes

dataML_top1 <- dataML_top1 %>%
  ungroup() %>%
  filter(!is.na(species))

dataML_top3 <- dataML_top3 %>%
  ungroup() %>%
  filter(!is.na(species))

dataML_top1 %>%
  filter(lubridate::year(Date_Time) == 2021) %>% 
  filter(species == "NOGO") %>% 
  group_by(point) %>% 
  summarize(n(), sum(logit_to_p_f_us(logit)), min(logit_to_p_f_us(logit)), max(logit_to_p_f_us(logit)), min(logit), max(logit)) %>% View() #HAWObyPOINT

View(HAWObyPOINT)

fwrite(dataML_top1, paste0(here, "/acoustic/data_ingest/output/dataML_top1.csv"))
