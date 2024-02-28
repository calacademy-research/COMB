library(tidyverse)
library(lubridate)
library(chron)

ml <- read_csv("models/input/file_logit_agg.csv")

ml$Time <- times(ml$Date_Time)

bbwo <- ml %>% filter(species=="BBWO", year(Date_Time)==2020)
class(bbwo$Date_Time)

bbwo$Time <- format(as.POSIXct(bbwo$Date_Time), format = "%H:%M:%S")

head(bbwo)

hits0 <- bbwo %>% filter(max_logit > 0)
hitsn1 <- bbwo %>% filter(max_logit > -1)

n_distinct(hits0$point)
unique(hits0$point)

n_distinct(hitsn1$point)
unique(hitsn1$point)

bytime <- hits0 %>% group_by(Time, point) %>% summarise(ncalls = n())
ggplot(bytime) + 
  geom_col(aes(x=Time, y=ncalls, color=as.factor(point))) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=45, hjust=1), 
        legend.position = "none")

bytime %>% filter(Time >= "06:00:00" & Time < "10:00:00") %>%
  ggplot() + 
  geom_col(aes(x=Time, y=ncalls, color=as.factor(point))) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=45, hjust=1))#, 
  #      legend.position = "none")

bbwo %>% filter(point==800)

dethist <- hits0 %>% filter(Time >= "06:00:00" & Time < "10:00:00", !is.na(point)) %>% group_by(point) %>% summarise(nfiles = n())

max(dethist$nfiles)
dethist %>% filter(nfiles==1)
hits0 %>% filter(point==893)

rowSums(data$y.aru, na.rm=T) # of 




amro <- ml %>% filter(species=="AMRO", year(Date_Time)==2020)
class(amro$Date_Time)

amro$Time <- format(as.POSIXct(amro$Date_Time), format = "%H:%M:%S")

head(amro)
amro0 <- amro %>% filter(max_logit > 0)


amrotime <- amro0 %>% group_by(Time, point) %>% summarise(ncalls = n())
ggplot(amrotime) + 
  geom_col(aes(x=Time, y=ncalls, color=as.factor(point))) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=45, hjust=1), 
        legend.position = "none")

nawa0 <- ml %>% filter(species=="NAWA", year(Date_Time)==2020, max_logit > 0)
class(nawa0$Date_Time)

nawa0$Time <- format(as.POSIXct(nawa0$Date_Time), format = "%H:%M:%S")

nawatime <- nawa0 %>% group_by(Time, point) %>% summarise(ncalls = n())
ggplot(nawatime) + 
  geom_col(aes(x=Time, y=ncalls, color=as.factor(point))) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=45, hjust=1), 
        legend.position = "none")
