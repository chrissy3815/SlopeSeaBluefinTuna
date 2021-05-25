### Calculations of larval abundance with different temporal and spatial parameters
# for manuscript revisions
# Chrissy Hernandez
# May 5 2021

# packages needed:
library(here)
library(sp)
library(rgdal)

# load data:
source(here("code", "SlopeSeaAbundProcess.R"))

# load SEAMAP data from 2016:
source(here("code", "GoMexAbundProcess.R"))

# initialize the table:
larvAbundSensTbl<- data.frame(Configuration = vector(),
                              NDays = vector(),
                              MeanAbund = vector(),
                              MeanAbundPosStn = vector())

## configurations:

#1. Gunther and Bigelow cruises, June 17-Aug 15, 1000 m and deeper
I<- which(all_bongo_stns$Month=="AUG" & all_bongo_stns$Day>15) 
subdata<- all_bongo_stns[-I,] # drop the late august stations
I<- which(subdata$BottomDepth<1000)
subdata<- subdata[-I,] # drop the shallow stations
meanAbund<- mean(subdata$Abundance)
meanPosStn<- mean(subdata$Abundance[subdata$Abundance>0])
# add row to table:
newrow<- data.frame(Configuration = "GU1608+HB1603, June 17-Aug 15, 1000m and deeper",
                    NDays = (30-16)+31+15,
                    MeanAbund = meanAbund,
                    MeanAbundPosStn = meanPosStn)
larvAbundSensTbl<- rbind(larvAbundSensTbl, newrow)

#2. Bigelow cruise only, June 28-Aug 15, 1000 m and deeper
I<- which(all_bongo_stns$Cruise=="HB1603")
subdata<- all_bongo_stns[I,] # pull out only the Bigelow cruise samples
I<- which(subdata$Month=="AUG" & subdata$Day>15)
subdata<- subdata[-I,] # drop the late august stations
I<- which(subdata$BottomDepth<1000)
subdata<- subdata[-I,] # drop the shallow samples
meanAbund<- mean(subdata$Abundance)
meanPosStn<- mean(subdata$Abundance[subdata$Abundance>0])
# add row to table:
newrow<- data.frame(Configuration = "HB1603, June 28-Aug 15, 1000m and deeper",
                    NDays = (30-27)+31+15,
                    MeanAbund = meanAbund,
                    MeanAbundPosStn = meanPosStn)
larvAbundSensTbl<- rbind(larvAbundSensTbl, newrow)

#3. Bigelow cruise only, 42-day window (June 28 - Aug 8)
I<- which(all_bongo_stns$Cruise=="HB1603")
subdata<- all_bongo_stns[I,] # pull out only the Bigelow cruise samples
I<- which(subdata$Month=="AUG" & subdata$Day>8)
subdata<- subdata[-I,] # drop the later august stations
meanAbund<- mean(subdata$Abundance)
meanPosStn<- mean(subdata$Abundance[subdata$Abundance>0])
# add row to table:
newrow<- data.frame(Configuration = "HB1603, June 28-Aug 8",
                    NDays = (30-27)+31+8,
                    MeanAbund = meanAbund,
                    MeanAbundPosStn = meanPosStn)
larvAbundSensTbl<- rbind(larvAbundSensTbl, newrow)

#4. Bigelow cruise only, 42-day window (June 28 - Aug 8), 1000 m or deeper
I<- which(all_bongo_stns$Cruise=="HB1603")
subdata<- all_bongo_stns[I,] # pull out only the Bigelow cruise samples
I<- which(subdata$Month=="AUG" & subdata$Day>8)
subdata<- subdata[-I,] # drop the later august stations
I<- which(subdata$BottomDepth<1000)
subdata<- subdata[-I,] # drop the shallow stations
meanAbund<- mean(subdata$Abundance)
meanPosStn<- mean(subdata$Abundance[subdata$Abundance>0])
# add row to table:
newrow<- data.frame(Configuration = "HB1603, June 28-Aug 8, 1000m and deeper",
                    NDays = (30-27)+31+8,
                    MeanAbund = meanAbund,
                    MeanAbundPosStn = meanPosStn)
larvAbundSensTbl<- rbind(larvAbundSensTbl, newrow)

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
# add row to table:
newrow<- data.frame(Configuration = "HB1603, June 28-Aug 24, stratified mean",
                    NDays = (30-27)+31+24,
                    MeanAbund = stratMean,
                    MeanAbundPosStn = stratMeanPos)
larvAbundSensTbl<- rbind(larvAbundSensTbl, newrow)


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
stratMeanPos<- 1/(shelfbreak_area+offshore_area)*(shelfbreakmeanPos*shelfbreak_area+offshoremeanPos*offshore_area)
# add row to table:
newrow<- data.frame(Configuration = "HB1603, June 28-Aug 8, stratified mean",
                    NDays = (30-27)+31+8,
                    MeanAbund = stratMean,
                    MeanAbundPosStn = stratMeanPos)
larvAbundSensTbl<- rbind(larvAbundSensTbl, newrow)

#7. Bigelow cruise only, 31 day window (June 28-July 28) include all stations, stratified mean
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
# add row to table:
newrow<- data.frame(Configuration = "HB1603, June 28-July 28, stratified mean",
                    NDays = (30-27)+28,
                    MeanAbund = stratMean,
                    MeanAbundPosStn = stratMeanPos)
larvAbundSensTbl<- rbind(larvAbundSensTbl, newrow)

#8. SEAMAP data:
# this dataset that I have covers April 30 to May 30, 31 days.
seamapmean<- mean(seamap_bluefin$Abundance)
seamapmean
seamapmeanPos<- mean(seamap_bluefin$Abundance[seamap_bluefin$Abundance>0])
seamapmeanPos
# add row to table:
newrow<- data.frame(Configuration = "SEAMAP 2016, April 30-May 30",
                    NDays = 31,
                    MeanAbund = seamapmean,
                    MeanAbundPosStn = seamapmeanPos)
larvAbundSensTbl<- rbind(larvAbundSensTbl, newrow)

#9. 2013 Slope Sea data: June 21 to Aug 18, 1000 m or deeper:
# Note: the Gunther cruise didn't get into deep water until June 21
I<- which((-1*all_bongos_2013$MODEL_DEPTH)>=1000)
subdata<- all_bongos_2013[I,]
meanAbund<- mean(subdata$Abundance)
meanPosStn<- mean(subdata$Abundance[subdata$Abundance>0])
# add row to table:
newrow<- data.frame(Configuration = "GU1308+HB1303, June 21-Aug 18, 1000 m or deeper",
                    NDays = (30-20)+31+18,
                    MeanAbund = meanAbund,
                    MeanAbundPosStn = meanPosStn)
larvAbundSensTbl<- rbind(larvAbundSensTbl, newrow)

# 10. 2013 Bigelow cruise all stations, July 2-Aug 18
I<- which(all_bongos_2013$CRUISE_NAME=="HB1303")
subdata<- all_bongos_2013[I,] # keep only the bigelow cruise
meanAbund<- mean(subdata$Abundance)
meanPosStn<- mean(subdata$Abundance[subdata$Abundance>0])
# add row to table:
newrow<- data.frame(Configuration = "HB1303, July 2-Aug 18, all stations",
                    NDays = (31-1)+18,
                    MeanAbund = meanAbund,
                    MeanAbundPosStn = meanPosStn)
larvAbundSensTbl<- rbind(larvAbundSensTbl, newrow)

# 11. 2013 Bigelow cruise, July 2-Aug 12, 1000 m or deeper
I<- which(all_bongos_2013$CRUISE_NAME=="HB1303")
subdata<- all_bongos_2013[I,] # keep only the bigelow cruise
I<- which((-1*subdata$MODEL_DEPTH)>=1000)
subdata<- subdata[I,] # keep only deep stations
I<- which(subdata$Month=="Aug" & subdata$Day>12)
subdata<- subdata[-I,] #drop the late August stations
meanAbund<- mean(subdata$Abundance)
meanPosStn<- mean(subdata$Abundance[subdata$Abundance>0])
# add row to table:
newrow<- data.frame(Configuration = "HB1303, July2-Aug 12, 1000 m or deeper",
                    NDays = (31-1)+12,
                    MeanAbund = meanAbund,
                    MeanAbundPosStn = meanPosStn)
larvAbundSensTbl<- rbind(larvAbundSensTbl, newrow)

# 12. 2013 Bigelow cruise, all stations (July 2-Aug 18), stratified mean 
# make a spatial points object of the subset of the stations
pointdata<- all_bongos_2013[all_bongos_2013$CRUISE_NAME=="HB1303",]
points_sps<- SpatialPointsDataFrame(pointdata[,c("LONGITUDE","LATITUDE")], pointdata)
proj4string(points_sps)<- CRS("+proj=longlat")
points_sps<- spTransform(points_sps, CRS=("+proj=utm +zone=18 +datum=WGS84 +units=km"))
# extract offshore points and calculate mean:
offshorepoints<- points_sps[!is.na(over(points_sps, offshore_sps)),]
offshoremean<- mean(offshorepoints$Abundance)
offshoremeanPos<- mean(offshorepoints$Abundance[offshorepoints$Abundance>0])
# extract shelfbreak points and calculate mean:
shelfbreakpoints<- points_sps[!is.na(over(points_sps, shelfbreak_sps)),]
shelfbreakmean<- mean(shelfbreakpoints$Abundance)
shelfbreakmeanPos<- mean(shelfbreakpoints$Abundance[shelfbreakpoints$Abundance>0])
# calculate the stratified means:
stratMean<- 1/(shelfbreak_area+offshore_area)*(shelfbreakmean*shelfbreak_area+offshoremean*offshore_area)
stratMeanPos<- 1/(shelfbreak_area+offshore_area)*(shelfbreakmeanPos*shelfbreak_area+offshoremeanPos*offshore_area)
# add row to table:
newrow<- data.frame(Configuration = "HB1303, July 2-Aug 18, stratified mean",
                    NDays = (31-1)+18,
                    MeanAbund = stratMean,
                    MeanAbundPosStn = stratMeanPos)
larvAbundSensTbl<- rbind(larvAbundSensTbl, newrow)

# 12. 2013 Bigelow cruise, all stations (July 2-Aug 1), stratified mean 
# make a spatial points object of the subset of the stations
pointdata<- all_bongos_2013[all_bongos_2013$CRUISE_NAME=="HB1303",]
I<- which(pointdata$Month=="Aug" & pointdata$Day>1)
pointdata<- pointdata[-I,] # drop most of the August data
points_sps<- SpatialPointsDataFrame(pointdata[,c("LONGITUDE","LATITUDE")], pointdata)
proj4string(points_sps)<- CRS("+proj=longlat")
points_sps<- spTransform(points_sps, CRS=("+proj=utm +zone=18 +datum=WGS84 +units=km"))
# extract offshore points and calculate mean:
offshorepoints<- points_sps[!is.na(over(points_sps, offshore_sps)),]
offshoremean<- mean(offshorepoints$Abundance)
offshoremeanPos<- mean(offshorepoints$Abundance[offshorepoints$Abundance>0])
# extract shelfbreak points and calculate mean:
shelfbreakpoints<- points_sps[!is.na(over(points_sps, shelfbreak_sps)),]
shelfbreakmean<- mean(shelfbreakpoints$Abundance)
shelfbreakmeanPos<- mean(shelfbreakpoints$Abundance[shelfbreakpoints$Abundance>0])
if (is.na(shelfbreakmeanPos)){shelfbreakmeanPos<- 0}
# calculate the stratified means:
stratMean<- 1/(shelfbreak_area+offshore_area)*(shelfbreakmean*shelfbreak_area+offshoremean*offshore_area)
stratMeanPos<- 1/(shelfbreak_area+offshore_area)*(shelfbreakmeanPos*shelfbreak_area+offshoremeanPos*offshore_area)
# add row to table:
newrow<- data.frame(Configuration = "HB1303, July 2-Aug 1, stratified mean",
                    NDays = (31-1)+1,
                    MeanAbund = stratMean,
                    MeanAbundPosStn = stratMeanPos)
larvAbundSensTbl<- rbind(larvAbundSensTbl, newrow)


### Write out the table as a csv, and then can edit the text in Excel:
write.csv(larvAbundSensTbl, 
          file = here("results", "LarvalAbundanceSensitivity.csv"), 
          row.names=F)

