##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
#This script calculates indicators of changing winters that are relevant for the ecosystems and people of
#the northern forest region of northeastern North America for the manuscript "Northern forest winters have lost cold, snowy conditions that are important for 
#ecosystems and human communities" published in Ecological Applications (doi:##################)
#Definitions for each indicator are located in "Northern_Forest_Metadata.xlsx"


#There are two main steps

#1. Calculate indicators of changing temperature, snow, and temperature plus snow conditions
#2. For each indicator, summarize frequency of occurrence within each site and year 

#Code was developed by N. Casson, A. Contosta, and S. Nelson

#the dataset required to run this script ("metfin.csv") is located in this repository (with associated metadata)
#and was produced by running the "Northern_Forest_SWE_Modeling" script along with the file
#"daily_station_data.csv."

####################################################################################
#Initial set up
####################################################################################

#call libraries
library(dplyr)
library(stringr)
library(tidyr)
library(varhandle)

#read data
metfin <- read.table("metfin.csv", head = TRUE, sep = ",", na.strings=c("NA", "NAN"))

####################################################################################
#1. Calculate winter climate change indicators related to changes in temperature
#conditions, snow conditions, or both
####################################################################################
#Temperature conditions

#thaw days
metfin$thawday = as.numeric(ifelse(metfin$TMAXfin > 0, 1, 0))

#frost days
metfin$frostday = ifelse(metfin$TMINfin < 0, 1, 0)

#freeze days
metfin$freezeday = ifelse(metfin$TMAXfin < 0, 1, 0)

#extreme cold days / pine beetle kills days
metfin$extremecold = ifelse(metfin$TMINfin < -18, 1, 0)

#hemlock woolly adelgid kills days
metfin$hwaday = ifelse(metfin$TMINfin < -30, 1, 0)

#snowmaking days / mosquito kill days
#snowmaking days determined at -5 and -2 degree thresholds(to account for changes in technology) 
# and at both annual timesteps (AN) and prior to Christmas (CH)
#note that snowmaking typically ends around 2/28
metfin$SM5AN = ifelse(metfin$TMIN < -5 & metfin$doy2 < 120,  1, 0)
metfin$SM5CH = ifelse(metfin$TMIN < -5 & metfin$Month <=12 & metfin$Day <= 25, 1, 0)

metfin$SM2AN = ifelse(metfin$TMIN < -2 & metfin$doy2 < 120, 1, 0)
metfin$SM2CH = ifelse(metfin$TMIN < -2 & metfin$Month <=12 & metfin$Day <= 22, 1, 0)

####################################################################################
#Snow conditions

#snowfall days
metfin$modSFD = ifelse(metfin$modSNOW > 0, 1, 0)

#non-snowfall days (liquid precip)
metfin$modLIQD = ifelse(metfin$modLIQ > 0, 1, 0)

#SCDs (snow covered days)
metfin$modSCD = ifelse(metfin$modSWE > 0, 1, 0)

#rain on snow days
metfin$modROSD = ifelse(metfin$modSCD > 0 & metfin$modLIQD > 0, 1, 0)

#bare ground days
metfin$modBGD = ifelse(metfin$modSWE == 0, 1, 0)

####################################################################################
#Temperature plus snow conditions

#bare ground plus thaw days
metfin$modBGTD = ifelse(metfin$modSWE == 0 & metfin$TMAXfin > 0, 1, 0)

#bare ground plus frost days
metfin$modBGFROD = ifelse(metfin$modSWE == 0 & metfin$TMINfin < 0, 1, 0)

####################################################################################
#2. For each indicator, summarize frequency of occurrence within each site and year 
####################################################################################
#Create blank data.frame to which summary stats will be added

WYear = unique(metfin$WYear)
Site_ID = unique(metfin$Site_ID)

#paste LAT and LON together because there are duplicates of each
LAT_LON = unique(paste(metfin$LAT, metfin$LON, sep = ","))

#split the string to obtain correct number of sites (n = 37) for LAT and LON
spl_LAT_LON <- data.frame(do.call(rbind, str_split(LAT_LON, ",")))
names(spl_LAT_LON) <- c("LAT", "LON")

#unfactor the LAT and LON variables (convert categorical / factor variable back into numeric format)
LAT = unfactor(spl_LAT_LON$LAT)
LON = unfactor(spl_LAT_LON$LON)

#create empty data.frame into which summary statistics will be added
calc.ind = data.frame(matrix(data = NA, nrow =  length(WYear) * length(Site_ID), ncol = 4))
names(calc.ind) = c("WYear", "Site_ID", "LAT", "LON")

#populate data.frame with values for WYear (winter year)
calc.ind$WYear = rep(WYear, nrow(calc.ind) / length(WYear))
#sort data by WYear
calc.ind = calc.ind[order(calc.ind$WYear), ]
#populate data.frame with site ID, LAT, and LON
calc.ind$Site_ID = rep(Site_ID, nrow(calc.ind) / length(Site_ID))
calc.ind$LAT = rep(LAT, nrow(calc.ind) / length(LAT))
calc.ind$LON = rep(LON, nrow(calc.ind) / length(LON))

#create column for unique site x winter year combination
calc.ind$SiteYr = paste(calc.ind$Site_ID, calc.ind$WYear, sep = " ")

#order by winter year and site 
calc.ind = calc.ind[order(calc.ind$WYear, calc.ind$Site_ID), ]

####################################################################################
#Calculate summary sites by SiteYr of condition, timing, or duration of indicators

#######
#THAWS#
#######

thaw = metfin[ , c("SiteYr", "thawday")]
thaw <- na.omit(thaw)

#count over entire "dormant" season, which is November 1 to May 31
count.thaw = thaw %>%
  group_by(SiteYr) %>%
  summarize(count.thaw = sum(thawday))

#add to calc.ind
calc.ind.1 = full_join(calc.ind, count.thaw, by = "SiteYr")

########
#FROSTS#
########

frost = metfin[ , c("SiteYr", "frostday")]
frost <- na.omit(frost)

#count over entire "dormant" season
count.frost = frost %>%
  group_by(SiteYr) %>%
  summarize(count.frost = sum(frostday))

#add to calc.ind
calc.ind.2 = full_join(calc.ind.1, count.frost, by = "SiteYr")

#########
#FREEZES#
#########

freeze = metfin[ , c("SiteYr", "freezeday")]
freeze <- na.omit(freeze)

#count over entire "dormant" season
count.freeze = freeze %>%
  group_by(SiteYr) %>%
  summarize(count.freeze = sum(freezeday))

#add to calc.ind
calc.ind.3 = full_join(calc.ind.2, count.freeze, by = "SiteYr")

############################################
#EXTREME COLD DAYS / PINE BEETLE KILLS DAYS#
############################################

extremecold = metfin[ , c("SiteYr", "extremecold")]
extremecold <- na.omit(extremecold)

#count over entire "dormant" season
count.extremecold = extremecold %>%
  group_by(SiteYr) %>%
  summarize(count.extremecold = sum(extremecold))

#add to calc.ind
calc.ind.4 = full_join(calc.ind.3, count.extremecold, by = "SiteYr")

#######################################
#Hemlock Wooly Adelgid (HWA) KILL DAYS#
#######################################

hwaday = metfin[ , c("SiteYr", "hwaday")]
hwaday <- na.omit(hwaday)

#count over entire "dormant" season
count.hwaday = hwaday %>%
  group_by(SiteYr) %>%
  summarize(count.hwaday = sum(hwaday))

#add to calc.ind
calc.ind.5 = full_join(calc.ind.4, count.hwaday, by = "SiteYr")

######################################
#SNOWMAKING DAYS / MOSQUITO KILL DAYS#
######################################

SMD = metfin[ , c("SiteYr", "SM5AN", "SM5CH", "SM2AN", "SM2CH")]
SMD <- na.omit(SMD)

#count occurrences before Christmas and for entire season with two different temperature thresholds
count.SM5AN = SMD %>%
  group_by(SiteYr) %>%
  summarize(count.SM5AN = sum(SM5AN))

count.SM5CH = SMD %>%
  group_by(SiteYr) %>%
  summarize(count.SM5CH = sum(SM5CH))

count.SM2AN = SMD %>%
  group_by(SiteYr) %>%
  summarize(count.SM2AN = sum(SM2AN))

count.SM2CH = SMD %>%
  group_by(SiteYr) %>%
  summarize(count.SM2CH = sum(SM2CH))

SMD.1 = full_join(count.SM5AN, count.SM5CH, by = "SiteYr")
SMD.2 = full_join(count.SM2AN, count.SM2CH, by = "SiteYr")
SMD.3 = full_join(SMD.1, SMD.2, by = "SiteYr")

#add to calc.ind
calc.ind.6 = full_join(calc.ind.5, SMD.3, by = "SiteYr")

#####################################################
#MODELED SNOW-COVERED DAYS simulated from hydromad in
#Northern_Forest_SWE_modeling.R script
#####################################################

modSCD = metfin[ , c("SiteYr", "modSCD")]
modSCD <- na.omit(modSCD)

#count over entire "dormant" season
count.modSCD = modSCD %>%
  group_by(SiteYr) %>%
  summarize(count.modSCD = sum(modSCD))

#add to calc.ind
calc.ind.7 = full_join(calc.ind.6, count.modSCD, by = "SiteYr")

####################################################
#MODELED BARE GROUND DAYS simulated from hydromad in
#Northern_Forest_SWE_modeling.R script
####################################################

modBGD = metfin[ , c("SiteYr", "modBGD")]
modBGD <- na.omit(modBGD)

#count over entire "dormant" season
count.modBGD = modBGD %>%
  group_by(SiteYr) %>%
  summarize(count.modBGD = sum(modBGD))

#add to calc.ind
calc.ind.8 = full_join(calc.ind.7, count.modBGD, by = "SiteYr")

####################################################
#MODELED RAIN-ON-SNOW DAYSsimulated from hydromad in
#Northern_Forest_SWE_modeling.R script
####################################################

modROSD = metfin[ , c("SiteYr", "modROSD")]
modROSD <- na.omit(modROSD)

#count over entire "dormant" season
count.modROSD = modROSD %>%
  group_by(SiteYr) %>%
  summarize(count.modROSD = sum(modROSD))

#add to calc.ind
calc.ind.9 = full_join(calc.ind.8, count.modROSD, by = "SiteYr")

#########################################################
#MODELED BARE GROUND THAW DAYS simulated from hydromad in
#Northern_Forest_SWE_modeling.R script
#########################################################

modBGTD = metfin[ , c("SiteYr", "modBGTD")]
modBGTD <- na.omit(modBGTD)

#count over entire "dormant" season
count.modBGTD = modBGTD %>%
  group_by(SiteYr) %>%
  summarize(count.modBGTD = sum(modBGTD))

#add to calc.ind
calc.ind.10 = full_join(calc.ind.9, count.modBGTD, by = "SiteYr")

##########################################################
#MODELED BARE GROUND FROST DAYS simulated from hydromad in
#Northern_Forest_SWE_modeling.R script
##########################################################

modBGFROD = metfin[ , c("SiteYr", "modBGFROD")]
modBGFROD <- na.omit(modBGFROD)

#count over entire "dormant" season
count.modBGFROD = modBGFROD %>%
  group_by(SiteYr) %>%
  summarize(count.modBGFROD = sum(modBGFROD))

#add to calc.ind
calc.ind.11 = full_join(calc.ind.10, count.modBGFROD, by = "SiteYr")

####################################################################################
#write calc.ind.11 table
write.table(calc.ind.11, "calcfin.csv", sep = ",", na="NA", append=FALSE, col.names=TRUE, row.names = FALSE)

