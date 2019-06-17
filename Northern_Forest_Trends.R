##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
#This script analyzes change over time in the frequency of winter climate change indicators for the northern
#forest region of northeastern North America for the manuscript "Northern forest winters have lost cold, snowy conditions that are important for 
#ecosystems and human communities" published in Ecological Applications (doi:##################)


#There are three main steps:
#1. Analyze trends for each site
#2. Analyze trends across the entire region
#3. Analyze trends within each subregion

#Code was developed by N. Casson, A. Contosta, and S. Nelson

#the dataset required to run this script ("calcfin.csv") is located in this repository (with associated metadata)
#and was produced by running the "Northern_Forest_Winter_Indicators" script along with the file
#"metfin.csv."

####################################################################################
#Initial set up
####################################################################################

#call libraries
library(Kendall)
library(rkt)
library(trend)

#read data
calc.fin <- read.table("calcfin.csv", head = TRUE, sep = ",", na.strings=c("NA", "NAN"))

#write function for sen slope Sen, P (1968). Estimated of the regression coefficient based on Kendall’s Tau. J Am Stat Assoc 39:1379-1389
sen = function(x, y){
  xx = outer(x, x, "-")
  yy = outer(y, y, "-")
  zz = yy/xx

  s1 = zz[lower.tri(zz)]
  s2 = subset(s1, s1!="NA" & s1!=Inf & s1!=-Inf)
  slope = median(s2)

  i1 = y - slope*x
  i2 = subset(i1, i1!="NA")
  intercept = median(i2)
  gi = data.frame(slope, intercept)
  return(gi)
}

####################################################################################
#1. Analyze trends for each site
####################################################################################

#order by Site_ID and WYear
calc.fin = calc.fin[order(calc.fin$Site_ID, calc.fin$WYear), ]

#omit duplicate SiteYrs from the beginning and the end of the dataframe to determine when 
#to start / end model subsets

start.rows <- !duplicated(calc.fin$Site_ID)
end.rows <- !duplicated(calc.fin$Site_ID, fromLast = TRUE)

#extract original row indexes associated with the start and end of an Site_ID
#add one row onto the start and subtract one row from the end to remove potential artifacts with 
#the first and last year of the record

sr.ind <- seq_along(calc.fin$Site_ID)[start.rows] + 1
er.ind <- seq_along(calc.fin$Site_ID)[end.rows] - 1

#check that the number of events is the same for the start and end
print(length(sr.ind))
print(length(er.ind))

#extract response variables (indicators or inds), predictor (WYear), and Site_ID, where WYear = Winter Year, or 
#November 1 to October 31
inds = calc.fin[ , -(1:5)]

pred = calc.fin[ , "WYear"]

vars = (calc.fin[ , c(2:4)])
vars.1 = vars[!duplicated(vars[1]), ]

#create containers to hold results of analysis
tau = data.frame(matrix(nrow = length(sr.ind), ncol = ncol(inds)))
p.value = data.frame(matrix(nrow = length(sr.ind), ncol = ncol(inds)))
slope = data.frame(matrix(nrow = length(sr.ind), ncol = ncol(inds)))

#begin loop over sites
for (h in seq(1,length(sr.ind))) {

#begin loop over indicators
for (i in 1:ncol(inds)) {

#make a data frame for each site x indicator 
nn = data.frame(inds[sr.ind[h]:er.ind[h],i], pred[sr.ind[h]:er.ind[h]])
#remove NA values (otherwise mk.test and sen won't work)
nn. = nn[complete.cases(nn),]
		if(nrow(nn.)>0 ){
		  #calculate Kendall's tau and associated p, and Sen slope and write to tau, p.value, and slope containers
      tau[h,i] = mk.test(nn.[,1])$estimates[3]
      p.value[h,i] = mk.test(nn.[,1])$p.value
      slope[h,i] = sen(nn.[,2], nn.[,1])$slope
    }
    }
    }

names(tau) = paste("tau", names(inds), "all", sep = "_")
names(p.value) = paste("pval", names(inds), "all", sep = "_")
names(slope) = paste("slope", names(inds), "all", sep = "_")

#compile tau, p, and sen slope into one master summary table
sum.tab.all = cbind(vars.1, tau, p.value, slope)

#write table
write.table(sum.tab.all, file = paste("sumtab_all.csv"), sep = ",", na="NA", append=FALSE, col.names=TRUE, row.names = FALSE)

####################################################################################
#2. Analyze trends across the entire northern forest region 
    #using regional Mann Kendall and Sen slope analyses
    #"Northern_Forest_Metadata.xlsx" includes defintions for each indicator
####################################################################################

#make Site_ID the blocking factor
block = calc.fin[ , "Site_ID"]

#run analyses for each indicator (brute force method; no fancy loops)
rkt.thaw = rkt(calc.fin$WYear, calc.fin$count.thaw, block = calc.fin$Site_ID)
rkt.frost = rkt(calc.fin$WYear, calc.fin$count.frost, block = calc.fin$Site_ID)
rkt.freeze = rkt(calc.fin$WYear, calc.fin$count.freeze, block = calc.fin$Site_ID)
rkt.extremecold = rkt(calc.fin$WYear, calc.fin$count.extremecold, block = calc.fin$Site_ID)
rkt.hwaday = rkt(calc.fin$WYear, calc.fin$count.hwaday, block = calc.fin$Site_ID)
rkt.SM2AN = rkt(calc.fin$WYear, calc.fin$count.SM2AN, block = calc.fin$Site_ID)
rkt.SM2CH = rkt(calc.fin$WYear, calc.fin$count.SM2CH, block = calc.fin$Site_ID)
rkt.SM5AN = rkt(calc.fin$WYear, calc.fin$count.SM5AN, block = calc.fin$Site_ID)
rkt.SM5CH = rkt(calc.fin$WYear, calc.fin$count.SM5CH, block = calc.fin$Site_ID)
rkt.modSCD = rkt(calc.fin$WYear, calc.fin$count.modSCD, block = calc.fin$Site_ID)
rkt.modBGD = rkt(calc.fin$WYear, calc.fin$count.modBGD, block = calc.fin$Site_ID)
rkt.modROSD = rkt(calc.fin$WYear, calc.fin$count.modROSD, block = calc.fin$Site_ID)
rkt.modBGFROD = rkt(calc.fin$WYear, calc.fin$count.modBGFROD.x, block = calc.fin$Site_ID)
rkt.modBGTD = rkt(calc.fin$WYear, calc.fin$count.modBGTD, block = calc.fin$Site_ID)

####################################################################################
#3. Analyze trends across three subregions: west, central, and east
####################################################################################
  
#subset data into subregions defined by longitude
west = calc.fin[calc.fin$LON < -87, ]
central = calc.fin[calc.fin$LON >= -87 & calc.fin$LON < -78, ]
east = calc.fin[calc.fin$LON >= -78, ]

#run analyses for each indicator and each subregion (brute force method; no fancy loops)

####################################################################################
#west
west.thaw = rkt(west$WYear, west$count.thaw, block = west$Site_ID)
west.frost = rkt(west$WYear, west$count.frost, block = west$Site_ID)
west.freeze = rkt(west$WYear, west$count.freeze, block = west$Site_ID)
west.extremecold = rkt(west$WYear, west$count.extremecold, block = west$Site_ID)
west.hwaday = rkt(west$WYear, west$count.hwaday, block = west$Site_ID)
west.SM2AN = rkt(west$WYear, west$count.SM2AN, block = west$Site_ID)
west.SM2CH = rkt(west$WYear, west$count.SM2CH, block = west$Site_ID)
west.SM5AN = rkt(west$WYear, west$count.SM5AN, block = west$Site_ID)
west.SM5CH = rkt(west$WYear, west$count.SM5CH, block = west$Site_ID)
west.modSCD = rkt(west$WYear, west$count.modSCD, block = west$Site_ID)
west.modBGD = rkt(west$WYear, west$count.modBGD, block = west$Site_ID)
west.modROSD = rkt(west$WYear, west$count.modROSD, block = west$Site_ID)
west.modBGFROD = rkt(west$WYear, west$count.modBGFROD.x, block = west$Site_ID)
west.modBGTD = rkt(west$WYear, west$count.modBGTD, block = west$Site_ID)

####################################################################################
#central
central.thaw = rkt(central$WYear, central$count.thaw, block = central$Site_ID)
central.frost = rkt(central$WYear, central$count.frost, block = central$Site_ID)
central.freeze = rkt(central$WYear, central$count.freeze, block = central$Site_ID)
central.extremecold = rkt(central$WYear, central$count.extremecold, block = central$Site_ID)
central.hwaday = rkt(central$WYear, central$count.hwaday, block = central$Site_ID)
central.SM2AN = rkt(central$WYear, central$count.SM2AN, block = central$Site_ID)
central.SM2CH = rkt(central$WYear, central$count.SM2CH, block = central$Site_ID)
central.SM5AN = rkt(central$WYear, central$count.SM5AN, block = central$Site_ID)
central.SM5CH = rkt(central$WYear, central$count.SM5CH, block = central$Site_ID)
central.modSCD = rkt(central$WYear, central$count.modSCD, block = central$Site_ID)
central.modBGD = rkt(central$WYear, central$count.modBGD, block = central$Site_ID)
central.modROSD = rkt(central$WYear, central$count.modROSD, block = central$Site_ID)
central.modBGFROD = rkt(central$WYear, central$count.modBGFROD.x, block = central$Site_ID)
central.modBGTD = rkt(central$WYear, central$count.modBGTD, block = central$Site_ID)

####################################################################################
#east
east.thaw = rkt(east$WYear, east$count.thaw, block = east$Site_ID)
east.frost = rkt(east$WYear, east$count.frost, block = east$Site_ID)
east.freeze = rkt(east$WYear, east$count.freeze, block = east$Site_ID)
east.extremecold = rkt(east$WYear, east$count.extremecold, block = east$Site_ID)
east.hwaday = rkt(east$WYear, east$count.hwaday, block = east$Site_ID)
east.SM2AN = rkt(east$WYear, east$count.SM2AN, block = east$Site_ID)
east.SM2CH = rkt(east$WYear, east$count.SM2CH, block = east$Site_ID)
east.SM5AN = rkt(east$WYear, east$count.SM5AN, block = east$Site_ID)
east.SM5CH = rkt(east$WYear, east$count.SM5CH, block = east$Site_ID)
east.modSCD = rkt(east$WYear, east$count.modSCD, block = east$Site_ID)
east.modBGD = rkt(east$WYear, east$count.modBGD, block = east$Site_ID)
east.modROSD = rkt(east$WYear, east$count.modROSD, block = east$Site_ID)
east.modBGFROD = rkt(east$WYear, east$count.modBGFROD.x, block = east$Site_ID)
east.modBGTD = rkt(east$WYear, east$count.modBGTD, block = east$Site_ID)
