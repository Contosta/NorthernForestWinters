##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
#This script calculates uses the hydromad model in the R package http://hydromad.catchment.org/
#to estimate snowater equilavent (SWE) from temperature and precipitation data
#and then to use temperature thresholds that best predicted SWE to partition total precipitation into
#snowfall and rainfall

#There are four main steps: 

#1. Gap fill temperature and precipitation data since the model won't run with NAs
#2. Model set-up
#3. Run model and merge results with metsub data.frame
#4. Validate model and partition precipitation into rainfall and snowfall

#Code was developed by A. Contosta as part of developing indicators of changing winters in the northern forest
#of northeastern North America for the manuscript "Northern forest winters have lost cold, snowy conditions that are important for 
#ecosystems and human communities" published in Ecological Applications (doi:##################)

#the dataset required to run this script("daily_station_data.csv") is located in this repository 
#(with associated metadata) and was produced by screening stations for daily total precipitation (PRCPfin), 
#maximum air temeprature (TMAXfin), and minimum air temperature (TMINfin) for record completeness
#the dataset also contains columns for snowfall (SNOWfin) and snow depth (SNWDfin), but these were not part of 
#record completeness screening

####################################################################################
#Initial set up                                                                   
####################################################################################

#call libraries

library(plyr)
library(hydromad)
library(matrixStats)
library(zoo)

#read data
metsub <- read.table("daily_station_data.csv", head = TRUE, sep = ",", na.strings=c("NA", "NAN"))

####################################################################################
#1. Gap fill temperature and precipitation data since the model won't run with NAs
####################################################################################

#create new columns for nation and SNOWFALL as SWE,
#where Canadian snowfall data are SWE and US data are converted from total precip but constrained by 
#whether or not TMEAN was < 0

metsub$NAT = ifelse(metsub$Site_ID > 1000000, "CA", "US")

#calculate TMEAN for original measurements (if they had been dropped during record completeness screening)
#then interpolate between dates where data are still missing

metsub$TMEAN2 = ifelse(is.na(metsub$TMEAN) == T & is.na(metsub$TMAX) == F & is.na(metsub$TMIN) == F,
                (metsub$TMAX + metsub$TMIN) / 2, metsub$TMEAN)

#interpolate between measurements that still have missing values

metsub$TMEAN3 = ifelse(is.na(metsub$TMEAN2) == T, na.approx(metsub$TMEAN2, na.rm = F), metsub$TMEAN2)

#fill PRCP data from original measurements (if they had been dropped during record completeness screening)

metsub$PRCP2 = ifelse(is.na(metsub$PRCPfin) == T, metsub$PRCP, metsub$PRCPfin)

#estimate precipitation from snowfall records, using raw SWE if NAT = CAN or SNOW / 10 if NAT = US 
#(assumes a snow density of 100 kg / m3)

metsub$PRCP3 = ifelse(is.na(metsub$PRCP) == T & is.na(metsub$SNOW) == F & metsub$NAT == "US", metsub$SNOW / 10,
               ifelse(is.na(metsub$PRCP) == T & is.na(metsub$SNOW) == F & metsub$NAT == "CAN", metsub$SNOW,
               metsub$PRCP))

#fill any additional missing PRCP values as 0
metsub$PRCP4 = ifelse(is.na(metsub$PRCP3) == T, 0, metsub$PRCP3)

#subset to only include complete cases (there are still some instances where TMEAN3 = NA)
metsubb = metsub[complete.cases(metsub[ , "TMEAN3"]), ]

#order by Site, Year, Month, Day
metsubb = metsubb[order(metsubb$Site_ID, metsubb$Year, metsubb$Month, metsubb$Day), ]

####################################################################################
#2. Model set-up
####################################################################################

#omit duplicate SiteYrs from the beginning and the end of the dataframe to determine when to 
#start / end model subsets
start.rows <- !duplicated(metsubb$SiteYr)
end.rows <- !duplicated(metsubb$SiteYr, fromLast = TRUE)

#extract original row indexes associated with the start and end of an SiteYr
sr.ind <- seq_along(metsubb$SiteYr)[start.rows]
er.ind <- seq_along(metsubb$SiteYr)[end.rows]

#check that the number of events is the same for the start and end
print(length(sr.ind))
print(length(er.ind))

#create vectors with all possible values for temperature thresholds for rainfall (Tmax) and snowmelt (Tmin)
#and all possible values for melt coefficient (kd)
Tmax = seq(-2, 2, by = 1)
Tmin = seq(-2, 2, by = 1)
kd = seq(0, 4, by = 1)

#compute all possible combinations of Tmax and Tmin
tpar = expand.grid(Tmax, Tmin)
names(tpar) = c("Tmax", "Tmin")

#select combinations where Tmax is > Tmin (model will not run if Tmax < Tmin
tpar = tpar[tpar$Tmax > tpar$Tmin, ]

#extract Tmax and Tmin for merging with kd
Tmax = tpar$Tmax
Tmin = tpar$Tmin

#compute all possible combinations of Tmax, Tmin, and kd
pars = data.frame(expand.grid(Tmax, Tmin, kd))
names(pars) = c("Tmax", "Tmin", "kd")

#select combinations where Tmax is > Tmin
#this step is critical because the model will not run if Tmax < Tmin!
pars = pars[pars$Tmax > pars$Tmin, ]

#make id column for removing duplicates
pars$id = paste(pars$Tmax, pars$Tmin, pars$kd, sep = " ")

#remove duplicate parameters
pars = subset(pars, !duplicated(pars$id))

#subset to only include Tmax, Tmin, and kd
pars = pars[ , c("Tmax", "Tmin", "kd")]

#transpose so that columns become rows (makes it easier to fill WEI)
pars = t(pars)

#create container to hold results of model. WEI = water equivalent of ice
WEI = data.frame(matrix(nrow = nrow(metsubb), ncol = ncol(pars)))

####################################################################################
#3. Run model and merge results with metsub data.frame
####################################################################################

#begin loop
for (h in seq(1,length(sr.ind))) {

P = ts(metsubb$PRCP4[sr.ind[h]:er.ind[h]])
E = ts(metsubb$TMEAN3[sr.ind[h]:er.ind[h]])
met.1 = ts.union(P,E)

for (i in 1:ncol(pars)) {

Tmax = pars[1, i]
Tmin = pars[2, i]
kd = pars[3, i]

metmod = snow.sim(met.1, Tmax = Tmax, Tmin = Tmin, cr = 1, cs = 1,
               kd = kd, kf = 0, rcap = 0,
               d = 0, f = 0, e = 0, LSWE_0 = 0, ISWE_0 = 0, return_state = T)

WEI[[i]][sr.ind[h]:er.ind[h]] =  metmod[ , 'SWE']

       }
}

#merge WEI with metsubb
metsubb.1 = cbind(metsubb, WEI)

#remove predicted snow values from data table when there had been NAs for TMEAN and PRCPfin
#to remove potentially spurious predictions from interpolated / filled values generated
#so that the model would not crash
#using lapply across columns
metsubb.1[32:81] = lapply(metsubb.1[32:81], function(x) ifelse(is.na(metsubb.1$PRCPfin) == T | is.na(metsubb.1$TMEAN) == T, NA, x))

#add back in variables that were not part of modeling (63 records that could not be gap filled)
metsubb.2 = metsub[!complete.cases(metsub[ , "TMEAN3"]), ]
WEI.1 = data.frame(matrix(nrow = nrow(metsubb.2), ncol = ncol(pars)))
metsubb.2 = cbind(metsubb.2, WEI.1)
names(metsubb.2) = names(metsubb.1)

#merge metsubb.1 and metsubb.2
metsubb.3 = rbind(metsubb.1, metsubb.2)

#change column names of predicted values to something more intuitive and consistent
pred = data.frame("predno" = c(1:50))
pred$prednam = ifelse(pred$predno < 10, paste("pred", 0, pred$predno, sep = ""), paste("pred", pred$predno, sep = ""))
names(metsubb.3)[32:81] <- pred$prednam

####################################################################################
#4. Validate model and partition precipitation into rainfall and snowfall
####################################################################################

#omit all SNOW zero values since many zero values can result from people recording data
snoz.1 = metsubb.3[metsubb.3$SNWD > 0, ]

#omit all NAs
snoz.2 = snoz.1[complete.cases(snoz.1[ , "SNWD"]) & complete.cases(snoz.1[ , "pred01"]), ]

#calculate the r2 between observed SNWD and estimated WEI for each site
func <- function(snoz.2)
{
return(data.frame(unlist(cor(snoz.2[ , "SNWD"], snoz.2[32:81]))^2))
}

WEIcor.s = ddply(snoz.2, .(Site_ID), func)

#identify the max correlation to select best WEI predictions

#Site
WEIcor.1s =  as.matrix(WEIcor.s[ , -1])
WEImod.1s = data.frame(cbind("rmax.s" = as.numeric(rowMaxs(WEIcor.1s)), "rmin.s" = as.numeric(rowMins(WEIcor.1s)),
         "Site_ID" = WEIcor.s$Site_ID))
WEImod.1s$mod_ID.s = colnames(WEIcor.1s)[max.col(WEIcor.1s,ties.method="first")]

#subset to only include rmax.s, Site_ID.s, and mod_ID.s
WEImod.2s = WEImod.1s[ , c("rmax.s", "Site_ID", "mod_ID.s")]

#merge WEImod.2s with metsub.3
metsubb.4 = merge(metsubb.3, WEImod.2s, by.x = "Site_ID", by.y = "Site_ID", all.x = TRUE, all.y = TRUE)

#create WEI with best model parameter sets for each site .
metsubb.4$pred.s = as.numeric(metsubb.4[cbind(seq_len(nrow(metsubb.4)), match(as.character(metsubb.4$mod_ID.s), colnames(metsubb.4)))])

#remove  all the columns not needed
metsubb.5 = metsubb.4[ , -(32:81)]

#create modeled snowfall columns with different model parameter sets for each site / year / month by model combination.

#make pars names the same as pred
pars = t(pars)
pars = data.frame(pars)
pars$pred = pred$prednam

#merge r2 values for with associated model parameters into single dataframe
metsubb.6 = merge(metsubb.5, pars, by.x = "mod_ID.s", by.y = "pred", all.x = T, all.y = F)

#sort by Site_ID, Winter Year, and doy2
metsubb.6 = metsubb.6[order(metsubb.6$Site_ID, metsubb.6$WYear, metsubb.6$doy2), ]

#calculate modeled snowfall
metsubb.6$mLIQ = ifelse(is.na(metsubb.6$PRCPfin) == T | is.na(metsubb.6$TMEAN) == T, NA,
                    ifelse(metsubb.6$TMEAN > metsubb.6$Tmax, metsubb.6$PRCPfin, 0))
metsubb.6$mSNOW = ifelse(is.na(metsubb.6$PRCPfin) == T | is.na(metsubb.6$TMEAN) == T, NA,
                     ifelse(metsubb.6$TMEAN < metsubb.6$Tmin, metsubb.6$PRCPfin, 0))

#compare modeled versus actual by nation (US = solid snow; Canada = liquid SWE)
can = metsubb.6[metsubb.6$NAT == "CA" & metsubb.6$SNOW > 0 & metsubb.6$mSNOW > 0, ]
us = metsubb.6[metsubb.6$NAT == "US" & metsubb.6$SNOW > 0 & metsubb.6$mSNOW > 0, ]

#Site, Year, and Month all perform similarly for estimating snowfall

#determine when the model predicts no snow when snow is actually on the ground (presence / absence) 
#or when the model predicts that
#precip falls as rain when in fact it falls as snow
metsubb.6$PASNWD = ifelse(metsubb.6$pred.s == 0 & metsubb.6$SNWD > 0, 1, 0)
metsubb.6$PASNOW = ifelse(metsubb.6$mSNOW.s == 0 & metsubb.6$SNOW > 0, 1, 0)

#remove all columns not needed in final met file
metfin = metsubb.6[ , c("STA_PROV", "STNAME", "LAT", "LON", "ELEV", "Site_ID", "DATETIME", "Year", "Month",   
                        "Day", "doy", "doy2", "WYear", "SiteYr", 
                        "PRCPfin", "SNWDfin", "SNOWfin", "TMAXfin", "TMINfin", "TMEAN",
                         "NAT",  "pred.s", "mSNOW", "mLIQ")]

names(metfin) = c("STA_PROV", "STNAME", "LAT", "LON", "ELEV",  "Site_ID", "DATETIME", "Year", "Month",   
                  "Day", "doy", "doy2", "WYear", "SiteYr", 
                  "PRCPfin", "SNWDfin", "SNOWfin", "TMAXfin", "TMINfin", "TMEAN",
                  "NAT",  "modSWE", "modSNOW", "modLIQ")

#write metfin table 
write.table(metfin, file = paste("metfin.csv"), sep = ",", na="NA", append=FALSE, col.names=TRUE, row.names = FALSE)

##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################

