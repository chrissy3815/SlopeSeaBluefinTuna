### Calculations of larval abundance with different temporal and spatial parameters
# for manuscript revisions
# Chrissy Hernandez
# May 5 2021

# packages needed:
library(here)
library(sp)
library(rgdal)

# load data:
source(here("code", "SlopeSea_GoMex_tunamaps.R"))
# from this other script, we want the tunacounts_SS object, which contains 
# (almost all) the information we need. Just need to add day and month.
tunacounts_SS<- merge(tunacounts_SS, 
                      SS2016_eventdata[,c("Cruise","Station","Operation","Date",'Day','Month')], 
                      by = c("Cruise", "Station", "Operation"), all.x=T)

## Check a few stats on overall bongo sampling:
all_bongo_stns<- SS2016_eventdata[SS2016_eventdata$Operation=="BON/CTD",]
# breakdown of Gunther vs. Bigelow cruises:
sum(all_bongo_stns$Cruise=="GU1608")
sum(all_bongo_stns$Cruise=="HB1603")
#1. How many stations are shallower than 1000 m for each cruise?
# For the Gunther cruise, if we only include the sampling from June 17-20, then 
# there are 24 stations from that cruise, all deeper than 1000 m.
I<- which(all_bongo_stns$BottomDepth<1000 & all_bongo_stns$Cruise=="HB1603")
length(I) #44 stations shallower than 1000 m on Bigelow cruise
I<- which(all_bongo_stns$BottomDepth<1000 & all_bongo_stns$Cruise=="HB1603" 
          & all_bongo_stns$Month=="AUG" & all_bongo_stns$Day<=15)
length(I) # We exclude 12 stations from August with the 1000 m requirement
I<- which(all_bongo_stns$BottomDepth<1000 & all_bongo_stns$Cruise=="HB1603" 
          & all_bongo_stns$Month!="AUG")
length(I) # We exclude 23 stations from June and July with the 1000 m requirement

# so let's make an object that has all the station information, with and without
# positive bluefin catches, for June 17 through end of summer. 
I<- which(all_bongo_stns$Month=="MAY")
all_bongo_stns<- all_bongo_stns[-I,] # drop the May sampling from ECOMON
I<- which(all_bongo_stns$Month=="JUN" & all_bongo_stns$Day<17)
all_bongo_stns<- all_bongo_stns[-I,] # drop the early June sampling from ECOMON
all_bongo_stns<- merge(all_bongo_stns, tunacounts_SS[,c("Cruise", "Station", "Operation", "Nbluefin","Vol_filtered", "Abundance", "Density")], 
                       by=c("Cruise", "Station","Operation"), all.x = T)
all_bongo_stns$Abundance[is.na(all_bongo_stns$Abundance)]<- 0
# some of the day column is a character and it seems to be messing things up.
all_bongo_stns$Day<- as.numeric(all_bongo_stns$Day)

## configurations:

#1. Gunther and Bigelow cruises, June 17-Aug 15, 1000 m and deeper
I<- which(all_bongo_stns$Month=="AUG" & all_bongo_stns$Day>15) 
subdata<- all_bongo_stns[-I,] # drop the late august stations
I<- which(subdata$BottomDepth<1000)
subdata<- subdata[-I,] # drop the shallow stations
meanAbund<- mean(subdata$Abundance)
meanAbund
meanPosStn<- mean(subdata$Abundance[subdata$Abundance>0])
meanPosStn

#2. Bigelow cruise only, June 28-Aug 15, 1000 m and deeper
I<- which(all_bongo_stns$Cruise=="HB1603")
subdata<- all_bongo_stns[I,] # pull out only the Bigelow cruise samples
I<- which(subdata$Month=="AUG" & subdata$Day>15)
subdata<- subdata[-I,] # drop the late august stations
I<- which(subdata$BottomDepth<1000)
subdata<- subdata[-I,] # drop the shallow samples
meanAbund<- mean(subdata$Abundance)
meanAbund
meanPosStn<- mean(subdata$Abundance[subdata$Abundance>0])
meanPosStn

#3. Bigelow cruise only, 42-day window (June 28 - Aug 8)
I<- which(all_bongo_stns$Cruise=="HB1603")
subdata<- all_bongo_stns[I,] # pull out only the Bigelow cruise samples
I<- which(subdata$Month=="AUG" & subdata$Day>8)
subdata<- subdata[-I,] # drop the later august stations
meanAbund<- mean(subdata$Abundance)
meanAbund
meanPosStn<- mean(subdata$Abundance[subdata$Abundance>0])
meanPosStn

#4. Bigelow cruise only, 42-day window (June 28 - Aug 8), 1000 m or deeper
I<- which(all_bongo_stns$Cruise=="HB1603")
subdata<- all_bongo_stns[I,] # pull out only the Bigelow cruise samples
I<- which(subdata$Month=="AUG" & subdata$Day>8)
subdata<- subdata[-I,] # drop the later august stations
I<- which(subdata$BottomDepth<1000)
subdata<- subdata[-I,] # drop the shallow stations
meanAbund<- mean(subdata$Abundance)
meanAbund
meanPosStn<- mean(subdata$Abundance[subdata$Abundance>0])
meanPosStn

#5. Bigelow cruise only, include all stations, stratified mean
# load strata:
offshore_polygon<- read.csv(here("data","AMAPPS_Strata","OffshelfStrataByMammal.csv"),
                            header=F)
names(offshore_polygon)<- c("Latitude", "Longitude")
shelfbreak_polygon<- read.csv(here("data","AMAPPS_Strata","ShelfbreakStrataByMammal.csv"),
                            header=F)
names(shelfbreak_polygon)<- c("Latitude", "Longitude")
# need to swap the order so that it's (x,y) 
shelfbreak_polygon<- cbind(shelfbreak_polygon$Longitude, shelfbreak_polygon$Latitude)
offshore_polygon<- cbind(offshore_polygon$Longitude, offshore_polygon$Latitude)
# convert offshore polygon to spatial polygon object
p<- Polygon(offshore_polygon)
ps<- Polygons(list(p),1)
offshore_sps<- SpatialPolygons(list(ps))
plot(offshore_sps)
proj4string(offshore_sps)<- CRS("+proj=longlat")
# now repeat for shelfbreak polygon
p<- Polygon(shelfbreak_polygon)
ps<- Polygons(list(p),1)
shelfbreak_sps<- SpatialPolygons(list(ps))
plot(shelfbreak_sps, add=T)
proj4string(shelfbreak_sps)<- CRS("+proj=longlat")
# make a spatial points object of all the stations
pointdata<- all_bongo_stns[all_bongo_stns$Cruise=="HB1603",]
points_sps<- SpatialPointsDataFrame(pointdata[,c("Longitude","Latitude")], pointdata)
proj4string(points_sps)<- CRS("+proj=longlat")
# extract offshore points and calculate mean:
offshorepoints<- points_sps[!is.na(over(points_sps, offshore_sps)),]
offshoremean<- mean(offshorepoints$Abundance)
offshoremean
offshoremeanPos<- mean(offshorepoints$Abundance[offshorepoints$Abundance>0])
offshoremeanPos
# extract shelfbreak points and calculate mean:
shelfbreakpoints<- points_sps[!is.na(over(points_sps, shelfbreak_sps)),]
shelfbreakmean<- mean(shelfbreakpoints$Abundance)
shelfbreakmean
shelfbreakmeanPos<- mean(shelfbreakpoints$Abundance[shelfbreakpoints$Abundance>0])
shelfbreakmeanPos
# calculate the polygon areas:
shelfbreak_sps<- spTransform(shelfbreak_sps, CRS=("+proj=utm +zone=18 +datum=WGS84 +units=km"))
offshore_sps<- spTransform(offshore_sps, CRS=("+proj=utm +zone=18 +datum=WGS84 +units=km"))
shelfbreak_area<- shelfbreak_sps@polygons[[1]]@area
offshore_area<- offshore_sps@polygons[[1]]@area
# calculate the stratified means:
stratMean<- 1/(shelfbreak_area+offshore_area)*(shelfbreakmean*shelfbreak_area+offshoremean*offshore_area)
stratMean
stratMeanPos<- 1/(shelfbreak_area+offshore_area)*(shelfbreakmeanPos*shelfbreak_area+offshoremeanPos*offshore_area)
stratMeanPos

#6. Bigelow cruise only, 42 day window (June 28-Aug 8) include all stations, stratified mean
# make a spatial points object of the subset of the stations
pointdata<- all_bongo_stns[all_bongo_stns$Cruise=="HB1603",]
I<- which(pointdata$Month=="AUG" & pointdata$Day>8)
pointdata<- pointdata[-I,] # drop the stations after August 8
points_sps<- SpatialPointsDataFrame(pointdata[,c("Longitude","Latitude")], pointdata)
proj4string(points_sps)<- CRS("+proj=longlat")
points_sps<- spTransform(points_sps, CRS=("+proj=utm +zone=18 +datum=WGS84 +units=km"))
# extract offshore points and calculate mean:
offshorepoints<- points_sps[!is.na(over(points_sps, offshore_sps)),]
offshoremean<- mean(offshorepoints$Abundance)
offshoremean
offshoremeanPos<- mean(offshorepoints$Abundance[offshorepoints$Abundance>0])
offshoremeanPos
# extract shelfbreak points and calculate mean:
shelfbreakpoints<- points_sps[!is.na(over(points_sps, shelfbreak_sps)),]
shelfbreakmean<- mean(shelfbreakpoints$Abundance)
shelfbreakmean
shelfbreakmeanPos<- mean(shelfbreakpoints$Abundance[shelfbreakpoints$Abundance>0])
shelfbreakmeanPos
# calculate the stratified means:
stratMean<- 1/(shelfbreak_area+offshore_area)*(shelfbreakmean*shelfbreak_area+offshoremean*offshore_area)
stratMean
stratMeanPos<- 1/(shelfbreak_area+offshore_area)*(shelfbreakmeanPos*shelfbreak_area+offshoremeanPos*offshore_area)
stratMeanPos

#7. Bigelow cruise only, 31 day window (June 28-July28) include all stations, stratified mean
# make a spatial points object of the subset of the stations
pointdata<- all_bongo_stns[all_bongo_stns$Cruise=="HB1603",]
I<- which(pointdata$Month=="AUG")
pointdata<- pointdata[-I,] # drop the stations in August
I<- which(pointdata$Month=="JUL" & pointdata$Day>28)
pointdata<- pointdata[-I,] # drop the stations after July 28
points_sps<- SpatialPointsDataFrame(pointdata[,c("Longitude","Latitude")], pointdata)
proj4string(points_sps)<- CRS("+proj=longlat")
points_sps<- spTransform(points_sps, CRS=("+proj=utm +zone=18 +datum=WGS84 +units=km"))
# extract offshore points and calculate mean:
offshorepoints<- points_sps[!is.na(over(points_sps, offshore_sps)),]
offshoremean<- mean(offshorepoints$Abundance)
offshoremean
offshoremeanPos<- mean(offshorepoints$Abundance[offshorepoints$Abundance>0])
offshoremeanPos
# extract shelfbreak points and calculate mean:
shelfbreakpoints<- points_sps[!is.na(over(points_sps, shelfbreak_sps)),]
shelfbreakmean<- mean(shelfbreakpoints$Abundance)
shelfbreakmean
shelfbreakmeanPos<- mean(shelfbreakpoints$Abundance[shelfbreakpoints$Abundance>0])
shelfbreakmeanPos
# calculate the stratified means:
stratMean<- 1/(shelfbreak_area+offshore_area)*(shelfbreakmean*shelfbreak_area+offshoremean*offshore_area)
stratMean
stratMeanPos<- 1/(shelfbreak_area+offshore_area)*(shelfbreakmeanPos*shelfbreak_area+offshoremeanPos*offshore_area)
stratMeanPos

#8. SEAMAP data:
# seamap_bongos object is already loaded in from running the tuna maps script.
# this dataset that I have covers April 30 to May 30, 31 days.
seamapmean<- mean(seamap_bongos$Abundance)
seamapmean
seamapmeanPos<- mean(seamap_bongos$Abundance[seamap_bongos$Abundance>0])
seamapmeanPos
