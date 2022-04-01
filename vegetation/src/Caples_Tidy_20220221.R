library(dplyr)
library(tidyverse)
library(tidyr)
library(writexl)
library(vegan)


#Import tables from Access database.  Includes data for 158 plot visits associated with avian plots (SMC + RFR).
Caples_PlotVisits_20220207 <- read_excel("~/FY2022_Tahoe/Caples/Caples_2021/Caples_Tidyverse/Caples_PlotVisits_20220207.xlsx")
TreeData_forExport_20220207c <- read_excel("~/FY2022_Tahoe/Caples/Caples_2021/Caples_Tidyverse/TreeData_forExport_20220207c.xlsx")
VegetationData_forExport_20220207c <- read_excel("~/FY2022_Tahoe/Caples/Caples_2021/Caples_Tidyverse/VegetationData_forExport_20220207c.xlsx")
#Seedling and Sapling data does not crosswalk from RFR plots
RegenData_forExport_20220128 <- read_excel("~/FY2022_Tahoe/Caples/Caples_2021/Caples_Tidyverse/RegenData_forExport_20220128.xlsx")
SpeciesComp_forExport_20220201 <- read_excel("~/FY2022_Tahoe/Caples/Caples_2021/Caples_Tidyverse/SpeciesComp_forExport_20220201.xlsx")
#CWD in RFR?
CWD_Data_forExport_20220216 <- read_excel("~/FY2022_Tahoe/Caples/Caples_2021/Caples_Tidyverse/CWD_Data_forExport_20220216.xlsx")


#Let's make some nicknames
Trees<-TreeData_forExport_20220207c
Veg<-VegetationData_forExport_20220207c
Regen<-RegenData_forExport_20220128
PlotVisits<-Caples_PlotVisits_20220207
SpComp<-SpeciesComp_forExport_20220201
CWD<-CWD_Data_forExport_20220217
#Plot_Visits for SMC plots only
SMCPlotVisits<-PlotVisits%>%
  filter(Veg_Type=="SMC")

##USE ID FOR PLOT VISITS (PV_ID) AS UNIQUE IDENTIFIER LINKING 1-TO-MANY TREE/SEEDLING/SAPLING RECORDS WITH PLOT VISITS

###VEGETATIVE COVER DATA
#################################################

#retain all variables, but detach Notes and other variables not collected in all years
Veg2<-select(Veg,-c(6,7,8:11,14,15,26,27))
#Retain NA values.  Note that heights are"0" when cover of veg classes are 0. 

###CWD
#################################################

CWD$Veg_Type<-factor(CWD$Veg_Type)
CWD$PV_ID<-as.numeric(CWD$PV_ID)
#Retain NA values. #Binned in inches, although recorded in cm (3"-12", 12"-24", 24"+) 
CWD2= CWD%>%
  filter(Diameter_Intersect>7.61999)%>%#exclude a few records that do not meet cWD diam critera(3" at transect intersect)
  filter(Length>0.99999)%>%#exclude a few records that do not meet cWD length critera(1 m at transect intersect)
  mutate("inverse_length"=(1/Length))%>%
  mutate("diam_squared"=Diameter_Intersect * Diameter_Intersect)%>%
  mutate(diam_class=cut(Diameter_Intersect, breaks=c(7.61,30.47,60.9,900),labels=c("3to12","12to24","24more")))%>%
  mutate(decay_class=cut(Decay_Class, breaks=c(0,3,6),labels=c("sound","rotten")))%>%#Decay Class 1-3 = sound, 4-5=rotten
  mutate(CWD_cat=paste0(diam_class, "_", decay_class))
 
#For CWD calcs, total transect length sampled by plot differes between SMC(45.2 m) and RFR (22.6 m) 
#Utizilizing formulas from DeVries 1973

#SMC CWD - volume per hectare
CWD_vol_ha_SMC<-CWD2%>%
  filter(Veg_Type=="SMC")%>%
  group_by(PV_ID, CWD_cat)%>%
  summarise(sum_diam_squared = sum(diam_squared))%>%
  mutate(CWD_vol_group=sum_diam_squared*0.027294)%>%
  group_by(PV_ID, CWD_cat)%>%
  summarise(CWD_vol=sum(CWD_vol_group))%>%
  spread(CWD_cat,CWD_vol)%>%
  replace(is.na(.), 0)
CWD_vol_ha_SMC2<-CWD_vol_ha_SMC%>%
  rowwise()%>%
  mutate(CWD_vol_All_ha = sum(c_across("12to24_rotten":"3to12_sound")))

CWD_vol_ha_RFR<-CWD2%>%
  filter(Veg_Type=="RFR")%>%
  group_by(PV_ID, CWD_cat)%>%
  summarise(sum_diam_squared = sum(diam_squared))%>%
  mutate(CWD_vol_group=sum_diam_squared*0.054589)%>%
  group_by(PV_ID, CWD_cat)%>%
  summarise(CWD_vol=sum(CWD_vol_group))%>%
  spread(CWD_cat,CWD_vol)%>%
  replace(is.na(.), 0)
CWD_vol_ha_RFR2<-CWD_vol_ha_RFR%>%
  rowwise()%>%
  mutate(CWD_vol_All_ha = sum(c_across("12to24_rotten":"3to12_sound")))

#SMC and RFR - CWD volume (m3) per ha
CWD_vol_ha<-bind_rows(CWD_vol_ha_SMC2,CWD_vol_ha_RFR2)%>%
  replace(is.na(.), 0)%>%
  rename("CWDvol_3to12_R"="3to12_rotten")%>%
  rename("CWDvol_3to12_S"="3to12_sound")%>%
  rename("CWDvol_12to24_R"="12to24_rotten")%>%
  rename("CWDvol_12to24_S"="12to24_sound")%>%
  rename("CWDvol_24plus_R"="24more_rotten")%>%
  rename("CWDvol_24plus_S"="24more_sound")%>%
  rename("CWDvol_ha_all"="CWD_vol_All_ha" )%>%
  select(-c("3to12_NA"))#one missing decay class value -- captured in CWDvol_ha_all
  
#SMC CWD - pieces per hectare
CWD_density_SMC<-CWD2%>%
  filter(Veg_Type=="SMC")%>%
  group_by(PV_ID, CWD_cat)%>%
  summarise(sum_inverse_lengths = sum(inverse_length))%>%
  mutate(CWD_ha_group=sum_inverse_lengths*347.5221)%>%
  group_by(PV_ID, CWD_cat)%>%
  summarise(CWD_ha=sum(CWD_ha_group))%>%
  spread(CWD_cat,CWD_ha)%>%
  replace(is.na(.), 0)
CWD_density_SMC2<-CWD_density_SMC%>%
  rowwise()%>%
  mutate(CWD_All_ha = sum(c_across("12to24_rotten":"3to12_sound")))
  
#RFR CWD - pieces per hectare
CWD_density_RFR<-CWD2%>%
  filter(Veg_Type=="RFR")%>%
  group_by(PV_ID, CWD_cat)%>%
  summarise(sum_inverse_lengths = sum(inverse_length))%>%
  mutate(CWD_ha_group=sum_inverse_lengths*695.0221239)%>%
  group_by(PV_ID, CWD_cat)%>%
  summarise(CWD_ha=sum(CWD_ha_group))%>%
  spread(CWD_cat,CWD_ha)%>%
  replace(is.na(.), 0)
CWD_density_RFR2<-CWD_density_RFR%>%
  rowwise()%>%
  mutate(CWD_All_ha = sum(c_across("12to24_rotten":"3to12_sound")))

#SMC and RFR - CWD pieces per ha
CWD_density_ha<-bind_rows(CWD_density_SMC2,CWD_density_RFR2)%>%
  replace(is.na(.), 0)%>%
  rename("CWDdensity_3to12_R"="3to12_rotten")%>%
  rename("CWDdensity_3to12_S"="3to12_sound")%>%
  rename("CWDdensity_12to24_R"="12to24_rotten")%>%
  rename("CWDdensity_12to24_S"="12to24_sound")%>%
  rename("CWDdensity_24plus_R"="24more_rotten")%>%
  rename("CWDdensity_24plus_S"="24more_sound")%>%
  rename("CWDdensity_ha_all"="CWD_All_ha" )%>%
  select(-c("3to12_NA"))#one missing decay class value -- captured in CWDdensity_ha_all

#SMC - CWD cover
CWD_cover_SMC<-CWD2%>%
  filter(Veg_Type=="SMC")%>%
  group_by(PV_ID)%>%
  summarise(sum_diam = sum(Diameter_Intersect))%>%
  mutate(CWD_cover=sum_diam*0.034752212)

#RFR - CWD cover 
CWD_cover_RFR<-CWD2%>%
  filter(Veg_Type=="RFR")%>%
  group_by(PV_ID)%>%
  summarise(sum_diam = sum(Diameter_Intersect))%>%
  mutate(CWD_cover=sum_diam*0.069504425)

CWD_cover<-bind_rows(CWD_cover_SMC,CWD_cover_RFR)%>%
  replace(is.na(.), 0)%>%
  select(-c(2))

###TREE DENSITIES BY DBH
############################################################

#Live Tree Densities by DBH and plot visit. Used size class bins from Caples HRV analysis, accounted for different plot size in RFR plots
LiveTreeTPH<-Trees %>%
  filter(Status %in% c("I","L","M"))%>%
  group_by(PV_ID, Veg_Type) %>%
  summarize(CountAllTrees = n(),
            CountLess20.3 = sum(DBH<20.3),
            Count20.3to40.6 = sum(DBH>=20.3 & DBH<40.6),
            Count40.6to61.0=sum(DBH>=40.6 & DBH<61.0),
            Count61.0to78.7=sum(DBH>=61.0 & DBH<78.7),
            CountMore78.7=sum(DBH>=78.7))%>%
  mutate(TPA = if_else(Veg_Type == "SMC", CountAllTrees*10, CountAllTrees*8.097))%>%
  mutate(TPH = if_else(Veg_Type == "SMC", CountAllTrees*24.71, CountAllTrees*20.01))%>%
  mutate(TPHLess20.3 = if_else(Veg_Type == "SMC", CountLess20.3*24.71, CountLess20.3*20.01))%>%
  mutate(TPH20.3to40.6 = if_else(Veg_Type == "SMC", Count20.3to40.6*24.71, Count20.3to40.6*20.01))%>%
  mutate(TPH40.6to61.0 = if_else(Veg_Type == "SMC", Count40.6to61.0*24.71, Count40.6to61.0*20.01))%>%
  mutate(TPH61.0to78.7 = if_else(Veg_Type == "SMC", Count61.0to78.7*24.71, Count61.0to78.7*20.01))%>%
  mutate(TPHMore78.7 = if_else(Veg_Type == "SMC", CountMore78.7*24.71, CountMore78.7*20.01))%>%
  select(-c(Veg_Type,CountAllTrees,CountLess20.3, Count20.3to40.6,Count40.6to61.0,Count61.0to78.7,CountMore78.7))%>%
  replace(is.na(.), 0)

#Snag Densities by DBH and plot visit
SnagTPH<-Trees %>%
  filter(Status %in% c("D","D1","D2","C","D3"))%>%
  group_by(PV_ID, Veg_Type) %>%
  summarize(CountAllSnags = n(),
     CountLess20.3Snags = sum(DBH<20.3),
     Count20.3to40.6Snags = sum(DBH>=20.3 & DBH<40.6),
     Count40.6to61.0Snags=sum(DBH>=40.6 & DBH<61.0),
     Count61.0to78.7Snags=sum(DBH>=61.0 & DBH<78.7),
     CountMore78.7Snags=sum(DBH>=78.7))%>%
  mutate(SnagsPA = if_else(Veg_Type == "SMC", CountAllSnags*10, CountAllSnags*8.097))%>%
  mutate(SnagsPH = if_else(Veg_Type == "SMC", CountAllSnags*24.71, CountAllSnags*20.01))%>%
  mutate(SnagsPHLess20.3 = if_else(Veg_Type == "SMC", CountLess20.3Snags*24.71, CountLess20.3Snags*20.01))%>%
  mutate(SnagsPH20.3to40.6 = if_else(Veg_Type == "SMC", Count20.3to40.6Snags*24.71, Count20.3to40.6Snags*20.01))%>%
  mutate(SnagsPH40.6to61.0 = if_else(Veg_Type == "SMC", Count40.6to61.0Snags*24.71, Count40.6to61.0Snags*20.01))%>%
  mutate(SnagsPH61.0to78.7 = if_else(Veg_Type == "SMC", Count61.0to78.7Snags*24.71, Count61.0to78.7Snags*20.01))%>%
  mutate(SnagsPHMore78.7 = if_else(Veg_Type == "SMC", CountMore78.7Snags*24.71, CountMore78.7Snags*20.01))%>%
  select(-c(Veg_Type,CountAllSnags,CountLess20.3Snags, Count20.3to40.6Snags,Count40.6to61.0Snags,Count61.0to78.7Snags,CountMore78.7Snags))%>%
  replace(is.na(.), 0)

###TREE SPECIES COMPOSITION BY PLOT (live)
###########################################################
    
 #Did not include species with very few occurrences (SASC, PIFL, etc)
  XTabTrees<-Trees %>%
      filter(Status %in% c("I","L","M"))%>%
      filter(Species %in% c("ABCO", "ABMA","CADE27","JUOC","PICO","PIJE","QUKE","PILA","PIMO","TSME"))%>%
      group_by(Species, PV_ID, Veg_Type) %>%
      tally() %>%
      spread(Species, n)%>%
      replace(is.na(.), 0)%>%
    mutate(ABCO_TPH = if_else(Veg_Type == "SMC", ABCO*24.71, ABCO*20.01))%>%
    mutate(ABMA_TPH = if_else(Veg_Type == "SMC", ABMA*24.71, ABMA*20.01))%>%
    mutate(CADE27_TPH = if_else(Veg_Type == "SMC", CADE27*24.71, CADE27*20.01))%>%
    mutate(JUOC_TPH = if_else(Veg_Type == "SMC", JUOC*24.71, JUOC*20.01))%>%
    mutate(PICO_TPH = if_else(Veg_Type == "SMC", PICO*24.71, PICO*20.01))%>%
    mutate(PIJE_TPH = if_else(Veg_Type == "SMC", PIJE*24.71, PIJE*20.01))%>%
    mutate(QUKE_TPH = if_else(Veg_Type == "SMC", QUKE*24.71, QUKE*20.01))%>%
    mutate(PILA_TPH = if_else(Veg_Type == "SMC", PILA*24.71, PILA*20.01))%>%
    mutate(PIMO_TPH = if_else(Veg_Type == "SMC", PIMO*24.71, PIMO*20.01))%>%
    mutate(TSME_TPH = if_else(Veg_Type == "SMC", TSME*24.71, TSME*20.01))%>%
      select(-c(2:12))
    
###REGEN DENSITIES 
############################################################
  
 #Count Saplings each plot visit. Live versus dead not distinguished in all years, so could not distinguish between these
 #No regen data for RFR plots
 CountSapling<-Regen %>%
       filter(Regen$Seed_or_Sap%in%"Sap")%>%
       group_by(PV_ID) %>%
    summarize(CountSapling = n())
  
#Count Seedlings each plot visit.
  CountSeedling<-Regen %>%
    filter(Regen$Seed_or_Sap%in%"Seed")%>%
    group_by(PV_ID) %>%
    summarize(CountSeedling = sum(Count))
  
  
###TREE HEIGHTS 
############################################################
 #Calculate maximum live tree height per plot visit
  MaxTreeHeight<-Trees %>%
    filter(Status %in% c("I","L","M"))%>%
    mutate_at(c(9),~replace_na(.,0))%>%#replacing NA heights with 0's so that script will calculate max
    group_by(PV_ID)%>%
    summarize(MaxTreeHeight=max(Height))
 
#Calculate maximum live height by crown class (did not include AB or RE coding, no data for RFR plots) 
  XTabCC_Heights<-Trees %>%
    #replace 'OV' and 'SU' with 'OT' so that one consistent code for overtopped trees
    mutate(Crown_Class = replace(Crown_Class, Crown_Class == "OV", "OT")) %>% 
    mutate(Crown_Class = replace(Crown_Class, Crown_Class == "SU", "OT"))  %>% 
    filter(Status %in% c("I","L","M"))%>%
    filter(Crown_Class %in% c("DO", "CO","IN","OT"))%>%
    mutate_at(c(9),~replace_na(.,0))%>%
    group_by(PV_ID, Crown_Class)%>%
    summarise(CC_MaxHeight=max(Height))%>%
    spread(Crown_Class, CC_MaxHeight)%>%
    replace(is.na(.), 0)
  XTabCC_MaxHt<-select(XTabCC_Heights, DO_Max_Ht = DO, CO_Max_Ht=CO, IN_Max_Ht=IN, OT_Max_Ht=OT)

##SPECIES COMPOSITION
#Transform to wide format
  SpCompWide<-SpComp%>%
    group_by(PV_ID, USDACode)%>%
    summarise(TotalCover=sum(PercentCover))%>%
    spread(USDACode,TotalCover)%>%
    replace(is.na(.), 0)
  
  #Remove Unknowns and Species occurring in < 5% of plots (SMC plots only)
  SpCompWideSMC5<-SMCPlotVisits%>%
    left_join(SpCompWide,by=c("PV_ID"))%>%
    select(-c(2:12,381:390))%>%
    purrr::discard(~sum(is.na(.x))/length(.x)* 100 >=95)%>%
    replace(is.na(.), 0)
    
###ATTACH ALL OUTPUTS
 
##Attach these data to Plots as the left side table.  Now this dataset includes ONLY the plot 
##visits designated as either Pre-Caples or Post-Caples visits in avian plots (158 plot visits). Six of
##these are RFR plots with no data for multiple attributes (seedlings, saplings, species comp, crown class)
  Caples_PlotData<-PlotVisits%>% 
    left_join(Veg2,by=c("PV_ID"))%>%
    left_join(LiveTreeTPH,by=c("PV_ID"))%>%
    left_join(SnagTPH, by=c("PV_ID"))%>%
    left_join(CountSapling, by=c("PV_ID"))%>%
    left_join(CountSeedling, by=c("PV_ID"))%>%
    left_join(XTabTrees, by=c("PV_ID"))%>%
    left_join(MaxTreeHeight, by= c("PV_ID"))%>%
    left_join(XTabCC_MaxHt, by  =c("PV_ID"))%>%
    left_join(CWD_cover, by =c("PV_ID"))%>%
    left_join(CWD_vol_ha, by = c("PV_ID"))%>%
    left_join(CWD_density_ha, by =c("PV_ID"))%>%
    left_join(SpCompWideSMC5, by  =c("PV_ID"))
  #Change NA values to zero where data was collected in both SMC and RFR plots (TPA, TPH - size class and species, CWD)
    Caples_PlotData_Zeros<-Caples_PlotData%>%
      mutate_at(c(29:42,45:54,60:74),~replace_na(.,0))
  
  #Change NA values where data not collected for RFR to 0's in SMC plots, but retain NA for RFR (couldn't find more elegant way to do this, although I am sure it exists!)
     Caples_SMC<-Caples_PlotData_Zeros%>%
       filter(Veg_Type=="SMC")%>%
       mutate_at(c(43,44,55:59),~replace_na(.,0))
     Caples_RFR<-Caples_PlotData_Zeros%>%
       filter(Veg_Type=="RFR")
     Caples_PlotData_20220225<-bind_rows(Caples_SMC,Caples_RFR)
     
#Export to Excel and to .RData object
    write_xlsx(Caples_PlotData_20220225,"C:\\Users\\kirstenbovee\\Documents\\FY2022_Tahoe\\Caples\\Caples_2021\\Caples_Tidyverse\\Output\\Caples_PlotData_20220225.xlsx")
    save(Caples_PlotData_20220225, file = "C:\\Users\\kirstenbovee\\Documents\\FY2022_Tahoe\\Caples\\Caples_2021\\Caples_Tidyverse\\Output\\Caples_PlotData_20220225.RData")

#CWD visualizations
library(ggplot2)

#select CWD vars and reshape from wide to long format
CWDvol_long<- Caples_PlotData_20220221%>%
  select(c(1,4,10,61:66))%>%
  gather(CWD_vol_class,CWD_vol_value,c(4:9))%>%
CWDdensity_long<- Caples_PlotData_20220221%>%
  select(c(1,4,10,68:73))%>%
  gather(CWD_density_class,CWD_density_value,c(4:9))

ggplot(CWDvol_long, aes(x = Treatment_Interval, y = CWD_vol_value, fill = CWD_vol_class)) +
  geom_bar(stat = "identity")

