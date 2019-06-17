##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
#This script produces the figures for the manuscript "Northern forest winters have lost cold, snowy conditions that are important for 
#ecosystems and human communities" published in Ecological Applications (doi:##################)

#Code was developed by N. Casson, A. Contosta, and S. Nelson

#the datasets required to run this script ("calcfin.csv" and "sumtab.csv") is located in this repository (with associated metadata)
#and was produced by running the "Northern_Forest_Winter_Indicators_Eco_Apps" and "Northern_Forest_Trends" script along with the file
#"daily_station_data.csv."

#Three spatial data files are also required and are located in this repository. They are:
# "gpr_000b11a_e.shp"; Source:Boundary Files, 2011 Census. Statistics Canada Catalogue no. 92-160-X.
# "cb_2016_us_state_500k.shp"; Source: US Census Bureau, Geography Division
# "NA_CEC_Eco_Level3.shp"; Source: US Environmental Protection Agency, 2010

#Two image files are required to create the map legends and are located in this repository. They are "legend.png" and "legend_inv.png"

####################################################################################
#Initial set up
####################################################################################

#call libraries
library(tidyverse)
library(raster)
library(rgdal)
library(cowplot)
library(png)
library(grid)
library(rgeos)
library(maptools)
#read in data
calcfin<-read.csv("calcfin.csv")
sumtab<-read.csv("sumtab_all.csv")

#classify sites into three zones
calcfin$ZONE<-ifelse(calcfin$LON< -87, "West", 
                      ifelse(calcfin$LON> -78, "East", "Central"))

calcfin$ZONE<-factor(calcfin$ZONE, levels=c("West", "Central", "East"))


#order zones from west to east
col_idx <- grep("ZONE", names(calcfin))
calcfin <- calcfin[, c(col_idx, (1:ncol(calcfin))[-col_idx])]


####################################################################################
#1. Format data for making figures
####################################################################################

##convert data to long format and merge (NOTE - WATCH COLUMN TITLES)
calcfin_long <- gather(calcfin, metric.name, value, count.thaw:count.modBGFROD, factor_key=TRUE)
sumtab_long<-gather(sumtab, metric.name, value, tau_count.thaw_all:slope_count.modBGFROD_all, factor_key=TRUE)
sumtab_long2<-separate(data = sumtab_long, col = metric.name, into=c("type", "metric.name","period"), sep = "_")

#merg calcfin and sumtab to format data for making plots
plotdata<-merge(calcfin_long, sumtab_long2, by=c("metric.name", "Site_ID", "LAT", "LON"))

plotdata_split<-spread(plotdata, type, value.y)


plotdata_split$bp<-ifelse(plotdata_split$pval<0.05&plotdata_split$tau>0, "pos", 
                          ifelse(plotdata_split$pval<0.05&plotdata_split$tau<0, "neg", 
                                 "ns"))
plotdata_split$slope<-ifelse(plotdata_split$bp=="ns", 0, plotdata_split$slope)


plotdata_split$ZONE<-as.factor(plotdata_split$ZONE)

plotdata_split<-transform(plotdata_split,
                           bp=factor(bp,levels=c("neg", "ns","pos")))



###read in and format map data
#read in canadian and us data and crop to study extent
can<-readOGR(dsn=".", layer = "gpr_000b11a_e")
can_crop <- crop(can, extent(-100, -62, 37.5, 53))


us<-readOGR(dsn=".", layer = "cb_2016_us_state_500k")
us_crop <- crop(us, extent(-100, -62, 37.5, 53))

#read in ecoregions, transform to match canadian and us data and crop to study extent
ecoregions<-readOGR(dsn=".", layer="NA_CEC_Eco_Level3")
ecoregions.2 <- spTransform(ecoregions, CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))
ecoregions_crop<-crop(ecoregions.2, extent(-100, -62, 37.5, 53))



#convert to ggplot object
us_df<-fortify(us_crop)
can_df<-fortify(can_crop)
ecoregions_df<-fortify(ecoregions_crop, region = "NA_L2NAME")


#identify ecoregions used in study
ecoregions_df$id2<-ifelse(ecoregions_df$id=="MIXED WOOD SHIELD"| ecoregions_df$id=="ATLANTIC HIGHLANDS"
                          | ecoregions_df$id=="MIXED WOOD PLAINS" |ecoregions_df$id=="CENTRAL USA PLAINS", "Study Region", "")

ecoregions_levels<-c("MIXED WOOD SHIELD","ATLANTIC HIGHLANDS","MIXED WOOD PLAINS","CENTRAL USA PLAINS")

##########################
########FIGURE 1##########
##########################

study_site_map=ggplot()+
  geom_polygon(data=ecoregions_df, aes(x=long, y=lat, group=group, fill=id2))+
  geom_point(data=plotdata_split, aes(x = LON, y = LAT), color="#474747", fill="#474747", pch=21, size=4)+
  scale_fill_manual(values = c("white", "#9ec2a8"))+
  theme_void()+theme(legend.position="top",legend.title = element_blank(), plot.margin=margin(t=0, r=0, l=0, b=0, "cm"))+
  geom_polygon(data=us_df, alpha=0.9,fill=NA, colour="black", size=0.3, aes(x=long, y=lat, group=group))+
  geom_polygon(data=can_df, alpha=0.9,fill=NA, colour="black", size=0.3, aes(x=long, y=lat, group=group))+
  geom_vline(xintercept = -78, color="black", size=1)+geom_vline(xintercept = -87, color="black", size=1)
  

jpeg('fig1.jpg', width = 6.0, height=4.125, units="in", res=600 )
study_site_map
dev.off()


##########################
########FIGURE 2##########
##########################


legend <- readPNG(source = "legend.png")
legend_inv <- readPNG(source = "legend_inv.png")
legend_g <- rasterGrob(legend, interpolate=TRUE)
legend_inv_g <- rasterGrob(legend_inv, interpolate=TRUE)

count_frost<-plotdata_split[which(plotdata_split$metric.name=="count.frost"),]
count_frost_plot<-ggplot(data = count_frost) + aes(x = WYear, y = value.x, group=Site_ID, colour=bp) + theme_bw()+ facet_wrap(~ZONE)+
  geom_line(aes(WYear, value.x, group=Site_ID), alpha=0.1, size=0.5)+
  scale_color_manual(values = alpha(c("#d65555", "grey","#595bcc"), 0.8), guide=F)+
  ylab("# Frost Days")+
  xlab("")+
  stat_smooth(data=count_frost,method="lm", se=F, alpha=0.4, size=0.75)+
  theme(text = element_text(size=6),axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_x_continuous(breaks=c(1920, 1960, 2000))+ylim(0,212)

count_freeze<-plotdata_split[which(plotdata_split$metric.name=="count.freeze"),]
count_freeze_plot<-ggplot(data = count_freeze) + aes(x = WYear, y = value.x, group=Site_ID, colour=bp) + theme_bw()+ facet_wrap(~ZONE)+
  geom_line(aes(WYear, value.x, group=Site_ID), alpha=0.1, size=0.5)+
  scale_color_manual(values = alpha(c("#d65555", "grey","#595bcc"), 0.8), guide=F)+
  ylab("# Ice Days")+
  xlab("")+
  stat_smooth(data=count_freeze,method="lm", se=F, alpha=0.4, size=0.75)+
  theme(text = element_text(size=6),axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_x_continuous(breaks=c(1920, 1960, 2000))+ylim(0,212)

count_extremecold<-plotdata_split[which(plotdata_split$metric.name=="count.extremecold"),]
count_extremecold_plot<-ggplot(data = count_extremecold) + aes(x = WYear, y = value.x, group=Site_ID, colour=bp) + theme_bw()+ facet_wrap(~ZONE)+
  geom_line(aes(WYear, value.x, group=Site_ID), alpha=0.1, size=0.5)+
  scale_color_manual(values = alpha(c("#d65555", "grey","#595bcc"), 0.8), guide=F)+
  ylab("# Extreme Cold Days")+
  xlab("Year")+stat_smooth(data=count_extremecold,method="lm", se=F, alpha=0.4, size=0.75)+
  theme(text = element_text(size=6))+
  scale_x_continuous(breaks=c(1920, 1960, 2000))+ylim(0,212)


count_frost_map<-ggplot()+
  geom_polygon(data=us_df, alpha=0.5,fill="#d8d0d0", colour="black", size=0.3, aes(x=long, y=lat, group=group))+
  geom_polygon(data=can_df, alpha=0.5,fill="#d8d0d0", colour="black", size=0.3, aes(x=long, y=lat, group=group))+
  geom_point(data=count_frost, aes(x = LON, y = LAT, fill=bp, size=abs(slope*4.5)), color="black", pch=21,alpha=0.5)+
  scale_fill_manual(values = alpha(c( "#d65555", "dark grey","#595bcc"),0.5), guide=F)+
  scale_size(breaks=c(0, 0.45, 0.9,1.35), range = c(0.9, 3.6))+
  theme_void()+theme(legend.position ="none")+
  geom_vline(xintercept = -78, color="black", size=1)+geom_vline(xintercept = -87, color="black", size=1)+
  annotation_custom(legend_g, xmin=-76, xmax=-64.47, ymin=49, ymax=53)

count_freeze_map<-ggplot()+
  geom_polygon(data=us_df, alpha=0.5,fill="#d8d0d0", colour="black", size=0.3, aes(x=long, y=lat, group=group))+
  geom_polygon(data=can_df, alpha=0.5,fill="#d8d0d0", colour="black", size=0.3, aes(x=long, y=lat, group=group))+
  geom_point(data=count_freeze, aes(x = LON, y = LAT, fill=bp, size=abs(slope*4.5)), color="black", pch=21,alpha=0.5)+
  scale_fill_manual(values = alpha(c( "#d65555", "dark grey","#595bcc"),0.5), guide=F)+
  scale_size(breaks=c(0, 0.45, 0.9,1.35), range = c(0.9, 3.6))+
  theme_void()+theme(legend.position ="none")+
  geom_vline(xintercept = -78, color="black", size=1)+geom_vline(xintercept = -87, color="black", size=1)+
  annotation_custom(legend_g, xmin=-76, xmax=-64.47, ymin=49, ymax=53)


count_extremecold_map<-ggplot()+
  geom_polygon(data=us_df, alpha=0.5,fill="#d8d0d0", colour="black", size=0.3, aes(x=long, y=lat, group=group))+
  geom_polygon(data=can_df, alpha=0.5,fill="#d8d0d0", colour="black", size=0.3, aes(x=long, y=lat, group=group))+
  geom_point(data=count_extremecold, aes(x = LON, y = LAT, fill=bp, size=abs(slope*4.5)), color="black", pch=21,alpha=0.5)+
  scale_fill_manual(values = alpha(c( "#d65555", "dark grey","#595bcc"),0.5), guide=F)+
  scale_size(breaks=c(0, 0.45, 0.9,1.35), range = c(0.9, 3.6))+
  theme_void()+theme(legend.position ="none")+
  geom_vline(xintercept = -78, color="black", size=1)+geom_vline(xintercept = -87, color="black", size=1)+
  annotation_custom(legend_g, xmin=-76, xmax=-64.47, ymin=49, ymax=53)

fig2<-plot_grid(count_frost_plot, count_frost_map, count_freeze_plot, count_freeze_map, count_extremecold_plot,
                count_extremecold_map, ncol=2, 
                align = 'h', axis = 'b',labels = c("A", "D", "B", "E", "C", "F"), label_size = 8)

jpeg('fig2.jpg', width=6, height=6, units="in", res=600)
fig2
dev.off()


##########################
########FIGURE 3##########
##########################



count_modSCD<-plotdata_split[which(plotdata_split$metric.name=="count.modSCD"),]
count_modSCD_plot<-ggplot(data = count_modSCD) + aes(x = WYear, y = value.x, group=Site_ID, colour=bp) + theme_bw()+ facet_wrap(~ZONE)+
  geom_line(aes(WYear, value.x, group=Site_ID), alpha=0.1, size=0.5)+
  scale_color_manual(values = alpha(c("#d65555", "grey","#595bcc"), 0.8), guide=F)+
  ylab("# Snow Covered Days")+
  xlab("")+stat_smooth(data=count_modSCD,method="lm", se=F, alpha=0.4, size=0.75)+
  theme(text = element_text(size=6),axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_x_continuous(breaks=c(1920, 1960, 2000))+ylim(0,212)

count_modBGD<-plotdata_split[which(plotdata_split$metric.name=="count.modBGD"),]
count_modBGD_plot<-ggplot(data = count_modBGD) + aes(x = WYear, y = value.x, group=Site_ID, colour=bp) + theme_bw()+ facet_wrap(~ZONE)+
  geom_line(aes(WYear, value.x, group=Site_ID), alpha=0.1, size=0.5)+
  scale_color_manual(values = alpha(c("#595bcc", "grey","#d65555"), 0.8), guide=F)+
  ylab("# Bare Ground Days")+
  xlab("Year")+stat_smooth(data=count_modBGD,method="lm", se=F, alpha=0.4, size=0.75)+
  theme(text = element_text(size=6))+
  scale_x_continuous(breaks=c(1920, 1960, 2000))+ylim(0,212)

count_modSCD_map<-ggplot()+
  geom_polygon(data=us_df, alpha=0.5,fill="#d8d0d0", colour="black", size=0.3, aes(x=long, y=lat, group=group))+
  geom_polygon(data=can_df, alpha=0.5,fill="#d8d0d0", colour="black", size=0.3, aes(x=long, y=lat, group=group))+
  geom_point(data=count_modSCD, aes(x = LON, y = LAT, fill=bp, size=abs(slope*4.5)), color="black", pch=21,alpha=0.5)+
  scale_fill_manual(values = alpha(c( "#d65555", "dark grey","#595bcc"),0.5), guide=F)+
  scale_size(breaks=c(0, 0.45, 0.9,1.35), range = c(0.9, 3.6))+
  theme_void()+theme(legend.position ="none")+
  geom_vline(xintercept = -78, color="black", size=1)+geom_vline(xintercept = -87, color="black", size=1)+
  annotation_custom(legend_g, xmin=-76, xmax=-64.47, ymin=49, ymax=53)

count_modBGD_map<-ggplot()+
  geom_polygon(data=us_df, alpha=0.5,fill="#d8d0d0", colour="black", size=0.3, aes(x=long, y=lat, group=group))+
  geom_polygon(data=can_df, alpha=0.5,fill="#d8d0d0", colour="black", size=0.3, aes(x=long, y=lat, group=group))+
  geom_point(data=count_modBGD, aes(x = LON, y = LAT, fill=bp, size=abs(slope*4.5)), color="black", pch=21,alpha=0.5)+
  scale_fill_manual(values = alpha(c( "#595bcc", "dark grey","#d65555"),0.5), guide=F)+
  scale_size(breaks=c(0, 0.45, 0.9,1.35), range = c(0.9, 3.6))+
  theme_void()+theme(legend.position ="none")+
  geom_vline(xintercept = -78, color="black", size=1)+geom_vline(xintercept = -87, color="black", size=1)+
  annotation_custom(legend_inv_g, xmin=-76, xmax=-64.47, ymin=49, ymax=53)





fig3<-plot_grid(count_modSCD_plot,count_modSCD_map, count_modBGD_plot, count_modBGD_map, 
                ncol=2, align = 'h', axis = 'b', labels = c("A", "C", "B", "D"), label_size = 8 )

jpeg('fig3.jpg', width=6, height=4, units="in", res=600)
fig3
dev.off()


##########################
########FIGURE 4##########
##########################


count.SM2AN<-plotdata_split[which(plotdata_split$metric.name=="count.SM2AN"),]
count.SM2AN_plot<-ggplot(data = count.SM2AN) + aes(x = WYear, y = value.x, group=Site_ID, colour=bp) + theme_bw()+ facet_wrap(~ZONE)+
  geom_line(aes(WYear, value.x, group=Site_ID), alpha=0.1, size=0.5)+
  scale_color_manual(values = alpha(c("#d65555", "grey","#595bcc"), 0.8), guide=F)+
  ylab("# Snow Making Days/\nMosquito Kill Days")+
  xlab("Year")+stat_smooth(data=count.SM2AN,method="lm", se=F, alpha=0.4, size=0.75)+
  theme(text = element_text(size=6))+
  scale_x_continuous(breaks=c(1920, 1960, 2000))+ylim(0,212)

count.SM2AN_map<-ggplot()+
  geom_polygon(data=us_df, alpha=0.5,fill="#d8d0d0", colour="black", size=0.3, aes(x=long, y=lat, group=group))+
  geom_polygon(data=can_df, alpha=0.5,fill="#d8d0d0", colour="black", size=0.3, aes(x=long, y=lat, group=group))+
  geom_point(data=count.SM2AN, aes(x = LON, y = LAT, fill=bp, size=abs(slope*4.5)), color="black", pch=21,alpha=0.5)+
  scale_fill_manual(values = alpha(c( "#d65555", "dark grey","#595bcc"),0.5), guide=F)+
  scale_size(breaks=c(0, 0.45, 0.9,1.35), range = c(0.9, 3.6))+
  theme_void()+theme(legend.position ="none")+
  geom_vline(xintercept = -78, color="black", size=1)+geom_vline(xintercept = -87, color="black", size=1)+
  annotation_custom(legend_g, xmin=-76, xmax=-64.47, ymin=49, ymax=53)

fig4<-plot_grid(count.SM2AN_plot,count.SM2AN_map, 
                ncol=2, align = 'h', axis = 'b', labels = "AUTO", label_size = 8 )

jpeg('fig4.jpg', width=6, height=2, units="in", res=600)
fig4
dev.off()

##########################
########FIGURE 5##########
##########################


count.modBGTD<-plotdata_split[which(plotdata_split$metric.name=="count.modBGTD"),]
count.modBGTD_plot<-ggplot(data = count.modBGTD) + aes(x = WYear, y = value.x, group=Site_ID, colour=bp) + theme_bw()+ facet_wrap(~ZONE)+
  geom_line(aes(WYear, value.x, group=Site_ID), alpha=0.1, size=0.5)+
  scale_color_manual(values = alpha(c("#595bcc", "dark grey","#d65555"), 0.8), guide=F)+
  ylab("# Mud Days")+
  xlab("Year")+stat_smooth(data=count.modBGTD,method="lm", se=F, alpha=0.4, size=0.75)+
  theme(text = element_text(size=6))+
  scale_x_continuous(breaks=c(1920, 1960, 2000))+ylim(0,212)

count.modBGTD_map<-ggplot()+
  geom_polygon(data=us_df, alpha=0.5,fill="#d8d0d0", colour="black", size=0.3, aes(x=long, y=lat, group=group))+
  geom_polygon(data=can_df, alpha=0.5,fill="#d8d0d0", colour="black", size=0.3, aes(x=long, y=lat, group=group))+
  geom_point(data=count.modBGTD, aes(x = LON, y = LAT, fill=bp, size=abs(slope*4.5)), color="black", pch=21,alpha=0.5)+
  scale_fill_manual(values = alpha(c("#595bcc", "dark grey","#d65555"),0.5), guide=F)+
  scale_size(breaks=c(0, 0.45, 0.9,1.35), range = c(0.9, 3.6))+
  theme_void()+theme(legend.position ="none")+
  geom_vline(xintercept = -78, color="black", size=1)+geom_vline(xintercept = -87, color="black", size=1)+
  annotation_custom(legend_inv_g, xmin=-76, xmax=-64.47, ymin=49, ymax=53)

fig5<-plot_grid(count.modBGTD_plot,count.modBGTD_map, 
                ncol=2, align = 'h', axis = 'b', labels = "AUTO", label_size = 8 )

jpeg('fig5.jpg', width=6, height=2, units="in", res=600)
fig5
dev.off()
