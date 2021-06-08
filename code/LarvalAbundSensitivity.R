### Calculations of larval abundance with different temporal and spatial parameters
# for manuscript revisions
# Chrissy Hernandez
# May 5 2021

# packages needed:
library(here)
library(sp)
library(rgdal)
library(oce)
library(ocedata)
library(ncdf4)

# A function for calculating the convex hull area in km2:
sample_area<- function(subdata, zone, buffer=F){
  points<- subdata[,c("Longitude","Latitude")]
  this_hull<- chull(points)
  this_polygon<- subdata[this_hull,c("Longitude","Latitude")]
  # convert polygon to spatial polygon object
  p<- Polygon(this_polygon)
  ps<- Polygons(list(p),1)
  this_sps<- SpatialPolygons(list(ps))
  proj4string(this_sps)<- CRS("+proj=longlat")
  this_CRS<- paste("+proj=utm +zone=",zone," +datum=WGS84 +units=km", sep='')
  this_sps<- spTransform(this_sps, CRS=this_CRS)
  this_area<- this_sps@polygons[[1]]@area
  if (buffer==T){
    points_to_buffer<- this_sps@polygons[[1]]@Polygons[[1]]@coords
    centroid<- this_sps@polygons[[1]]@labpt
    new_coords<- matrix(data=NA, nrow=dim(points_to_buffer)[1], ncol=2)
    for (i in 1:dim(points_to_buffer)[1]){
      sub_vec<- points_to_buffer[i,]-centroid
      lengthi<- sqrt(sum(sub_vec^2))
      sub_vec<- sub_vec/lengthi
      new_point<- centroid+sub_vec*(lengthi+28)
      new_coords[i,]<- new_point
    }
    p<- Polygon(new_coords)
    ps<- Polygons(list(p),1)
    this_sps<- SpatialPolygons(list(ps))
    this_CRS<- paste("+proj=utm +zone=",zone," +datum=WGS84 +units=km", sep='')
    proj4string(this_sps)<- this_CRS
    this_area<- this_sps@polygons[[1]]@area
  }
  return(this_area)
}

# load data:
source(here("code", "SlopeSeaAbundProcess.R"))

# load SEAMAP data from 2016:
source(here("code", "GoMexAbundProcess.R"))

# initialize the table:
larvAbundSensTbl<- data.frame(Configuration = vector(),
                              NDays = vector(),
                              Area = vector(),
                              NPerTow = vector(),
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
catchrate<- sum(subdata$Nbluefin, na.rm=T)/length(subdata$Station)
sum(subdata$DayNight=="Day")/length(subdata$DayNight)
# add row to table:
newrow<- data.frame(Configuration = "GU1608+HB1603, June 17-Aug 15, 1000m and deeper",
                    NDays = (30-16)+31+15,
                    Area = sample_area(subdata, 18, buffer=T),
                    NPerTow = catchrate,
                    MeanAbund = meanAbund,
                    MeanAbundPosStn = meanPosStn)
larvAbundSensTbl<- rbind(larvAbundSensTbl, newrow)


#1a. Gunther and Bigelow cruises, June 17-Aug 15, all stations
I<- which(all_bongo_stns$Month=="AUG" & all_bongo_stns$Day>15) 
subdata<- all_bongo_stns[-I,] # drop the late august stations
meanAbund<- mean(subdata$Abundance)
meanPosStn<- mean(subdata$Abundance[subdata$Abundance>0])
catchrate<- sum(subdata$Nbluefin, na.rm=T)/length(subdata$Station)
# add row to table:
newrow<- data.frame(Configuration = "GU1608+HB1603, June 17-Aug 15, all stations",
                    NDays = (30-16)+31+15,
                    Area = sample_area(subdata, 18, buffer=T),
                    NPerTow = catchrate,
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
catchrate<- sum(subdata$Nbluefin, na.rm=T)/length(subdata$Station)
# add row to table:
newrow<- data.frame(Configuration = "HB1603, June 28-Aug 15, 1000m and deeper",
                    NDays = (30-27)+31+15,
                    Area = sample_area(subdata, 18, buffer=T),
                    NPerTow = catchrate,
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
catchrate<- sum(subdata$Nbluefin, na.rm=T)/length(subdata$Station)
# add row to table:
newrow<- data.frame(Configuration = "HB1603, June 28-Aug 8",
                    NDays = (30-27)+31+8,
                    Area = sample_area(subdata, 18, buffer=T),
                    NPerTow = catchrate,
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
catchrate<- sum(subdata$Nbluefin, na.rm=T)/length(subdata$Station)
# add row to table:
newrow<- data.frame(Configuration = "HB1603, June 28-Aug 8, 1000m and deeper",
                    NDays = (30-27)+31+8,
                    Area = sample_area(subdata, 18, buffer=T),
                    NPerTow = catchrate,
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
offshoremeanPos<- mean(offshorepoints$Abundance[offshorepoints$Abundance>0])
offshorecatchrate<- sum(offshorepoints$Nbluefin, na.rm=T)/length(offshorepoints$Station)
# extract shelfbreak points and calculate mean:
shelfbreakpoints<- points_sps[!is.na(over(points_sps, shelfbreak_sps)),]
shelfbreakmean<- mean(shelfbreakpoints$Abundance)
shelfbreakmeanPos<- mean(shelfbreakpoints$Abundance[shelfbreakpoints$Abundance>0])
shelfbreakcatchrate<- sum(shelfbreakpoints$Nbluefin, na.rm=T)/length(shelfbreakpoints$Station)
# calculate the polygon areas:
shelfbreak_sps<- spTransform(shelfbreak_sps, CRS=("+proj=utm +zone=18 +datum=WGS84 +units=km"))
offshore_sps<- spTransform(offshore_sps, CRS=("+proj=utm +zone=18 +datum=WGS84 +units=km"))
shelfbreak_area<- shelfbreak_sps@polygons[[1]]@area
offshore_area<- offshore_sps@polygons[[1]]@area
# calculate the stratified means:
stratMean<- 1/(shelfbreak_area+offshore_area)*(shelfbreakmean*shelfbreak_area+offshoremean*offshore_area)
stratMeanPos<- 1/(shelfbreak_area+offshore_area)*(shelfbreakmeanPos*shelfbreak_area+offshoremeanPos*offshore_area)
stratcatchrate<- 1/(shelfbreak_area+offshore_area)*(shelfbreakcatchrate*shelfbreak_area+offshorecatchrate*offshore_area)
# add row to table:
newrow<- data.frame(Configuration = "HB1603, June 28-Aug 24, stratified mean",
                    NDays = (30-27)+31+24,
                    Area = shelfbreak_area+offshore_area,
                    NPerTow = stratcatchrate,
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
offshoremeanPos<- mean(offshorepoints$Abundance[offshorepoints$Abundance>0])
offshorecatchrate<- sum(offshorepoints$Nbluefin, na.rm=T)/length(offshorepoints$Station)
# extract shelfbreak points and calculate mean:
shelfbreakpoints<- points_sps[!is.na(over(points_sps, shelfbreak_sps)),]
shelfbreakmean<- mean(shelfbreakpoints$Abundance)
shelfbreakmeanPos<- mean(shelfbreakpoints$Abundance[shelfbreakpoints$Abundance>0])
shelfbreakcatchrate<- sum(shelfbreakpoints$Nbluefin, na.rm=T)/length(shelfbreakpoints$Station)
# calculate the stratified means:
stratMean<- 1/(shelfbreak_area+offshore_area)*(shelfbreakmean*shelfbreak_area+offshoremean*offshore_area)
stratMeanPos<- 1/(shelfbreak_area+offshore_area)*(shelfbreakmeanPos*shelfbreak_area+offshoremeanPos*offshore_area)
stratcatchrate<- 1/(shelfbreak_area+offshore_area)*(shelfbreakcatchrate*shelfbreak_area+offshorecatchrate*offshore_area)
# add row to table:
newrow<- data.frame(Configuration = "HB1603, June 28-Aug 8, stratified mean",
                    NDays = (30-27)+31+8,
                    Area = offshore_area+shelfbreak_area,
                    NPerTow = stratcatchrate,
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
offshoremeanPos<- mean(offshorepoints$Abundance[offshorepoints$Abundance>0])
offshorecatchrate<- sum(offshorepoints$Nbluefin, na.rm=T)/length(offshorepoints$Station)
# extract shelfbreak points and calculate mean:
shelfbreakpoints<- points_sps[!is.na(over(points_sps, shelfbreak_sps)),]
shelfbreakmean<- mean(shelfbreakpoints$Abundance)
shelfbreakmeanPos<- mean(shelfbreakpoints$Abundance[shelfbreakpoints$Abundance>0])
shelfbreakcatchrate<- sum(shelfbreakpoints$Nbluefin, na.rm=T)/length(shelfbreakpoints$Station)
# calculate the stratified means:
stratMean<- 1/(shelfbreak_area+offshore_area)*(shelfbreakmean*shelfbreak_area+offshoremean*offshore_area)
stratMeanPos<- 1/(shelfbreak_area+offshore_area)*(shelfbreakmeanPos*shelfbreak_area+offshoremeanPos*offshore_area)
stratcatchrate<- 1/(shelfbreak_area+offshore_area)*(shelfbreakcatchrate*shelfbreak_area+offshorecatchrate*offshore_area)
# add row to table:
newrow<- data.frame(Configuration = "HB1603, June 28-July 28, stratified mean",
                    NDays = (30-27)+28,
                    Area = shelfbreak_area+offshore_area,
                    NPerTow = stratcatchrate,
                    MeanAbund = stratMean,
                    MeanAbundPosStn = stratMeanPos)
larvAbundSensTbl<- rbind(larvAbundSensTbl, newrow)

#8. SEAMAP data:
# this dataset that I have covers April 30 to May 30, 31 days.
seamapmean<- mean(seamap_bluefin$Abundance)
seamapmeanPos<- mean(seamap_bluefin$Abundance[seamap_bluefin$Abundance>0])
catchrate<- sum(seamap_bluefin$TOT_LARVAE, na.rm=T)/length(seamap_bluefin$P_STA_NO)
# lat/lon for area:
seamap_stations<- seamap_bluefin[,c("STA_LAT","STA_LON")]
names(seamap_stations)<- c("Latitude","Longitude")
seamap_xx<- c(28.25,   28.25,  29.75,  30.25,  28.5,  26.25, 24.75, 23.75,  23.75,  24.25,  25.75,  25.75,  28.25)
seamap_yy<- c(-96.25, -88.25, -88.25, -87.00, -84.75, -83.75, -83.25, -83.25, -84.25, -85.25, -87.75,  -96.25, -96.25)
this_polygon<- data.frame(Latitude=seamap_yy, Longitude=seamap_xx)
# convert polygon to spatial polygon object
p<- Polygon(this_polygon)
ps<- Polygons(list(p),1)
this_sps<- SpatialPolygons(list(ps))
proj4string(this_sps)<- CRS("+proj=longlat")
this_sps<- spTransform(this_sps, CRS="+proj=utm +zone=16 +datum=WGS84 +units=km")
this_area<- this_sps@polygons[[1]]@area
# add row to table:
newrow<- data.frame(Configuration = "SEAMAP 2016, April 30-May 30",
                    NDays = 31,
                    Area = this_area,
                    NPerTow = catchrate,
                    MeanAbund = seamapmean,
                    MeanAbundPosStn = seamapmeanPos)
larvAbundSensTbl<- rbind(larvAbundSensTbl, newrow)

#9. 2013 Slope Sea data: June 21 to Aug 18, 1000 m or deeper:
# Note: the Gunter cruise didn't get into deep water until June 21
I<- which((-1*all_bongos_2013$MODEL_DEPTH)>=1000)
subdata<- all_bongos_2013[I,]
meanAbund<- mean(subdata$Abundance)
meanPosStn<- mean(subdata$Abundance[subdata$Abundance>0])
catchrate<- sum(subdata$BLUEFIN_1, subdata$BLUEFIN_2, na.rm=T)/length(subdata$STATION)
ss_stations<- subdata[,c("LONGITUDE","LATITUDE")]
names(ss_stations)<- c("Longitude","Latitude")
# add row to table:
newrow<- data.frame(Configuration = "GU1302+HB1303, June 21-Aug 18, 1000 m or deeper",
                    NDays = (30-20)+31+18,
                    Area = sample_area(ss_stations, 18, buffer=T),
                    NPerTow = catchrate,
                    MeanAbund = meanAbund,
                    MeanAbundPosStn = meanPosStn)
larvAbundSensTbl<- rbind(larvAbundSensTbl, newrow)

# 10. 2013 Bigelow cruise all stations, July 2-Aug 18
I<- which(all_bongos_2013$CRUISE_NAME=="HB1303")
subdata<- all_bongos_2013[I,] # keep only the bigelow cruise
meanAbund<- mean(subdata$Abundance)
meanPosStn<- mean(subdata$Abundance[subdata$Abundance>0])
catchrate<- sum(subdata$BLUEFIN_1, subdata$BLUEFIN_2, na.rm=T)/length(subdata$STATION)
ss_stations<- subdata[,c("LONGITUDE","LATITUDE")]
names(ss_stations)<- c("Longitude","Latitude")
# add row to table:
newrow<- data.frame(Configuration = "HB1303, July 2-Aug 18, all stations",
                    NDays = (31-1)+18,
                    Area = sample_area(ss_stations, 18, buffer=T),
                    NPerTow = catchrate,
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
catchrate<- sum(subdata$BLUEFIN_1, subdata$BLUEFIN_2, na.rm=T)/length(subdata$STATION)
ss_stations<- subdata[,c("LONGITUDE","LATITUDE")]
names(ss_stations)<- c("Longitude","Latitude")
# add row to table:
newrow<- data.frame(Configuration = "HB1303, July2-Aug 12, 1000 m or deeper",
                    NDays = (31-1)+12,
                    Area = sample_area(ss_stations, 18, buffer=T),
                    NPerTow = catchrate,
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
offshorecatchrate<- sum(offshorepoints$BLUEFIN_1, offshorepoints$BLUEFIN_2, na.rm=T)/length(offshorepoints$STATION)
# extract shelfbreak points and calculate mean:
shelfbreakpoints<- points_sps[!is.na(over(points_sps, shelfbreak_sps)),]
shelfbreakmean<- mean(shelfbreakpoints$Abundance)
shelfbreakmeanPos<- mean(shelfbreakpoints$Abundance[shelfbreakpoints$Abundance>0])
shelfbreakcatchrate<- sum(shelfbreakpoints$BLUEFIN_1, shelfbreakpoints$BLUEFIN_2, na.rm=T)/length(shelfbreakpoints$STATION)
# calculate the stratified means:
stratMean<- 1/(shelfbreak_area+offshore_area)*(shelfbreakmean*shelfbreak_area+offshoremean*offshore_area)
stratMeanPos<- 1/(shelfbreak_area+offshore_area)*(shelfbreakmeanPos*shelfbreak_area+offshoremeanPos*offshore_area)
stratcatchrate<- 1/(shelfbreak_area+offshore_area)*(shelfbreakcatchrate*shelfbreak_area+offshorecatchrate*offshore_area)
# add row to table:
newrow<- data.frame(Configuration = "HB1303, July 2-Aug 18, stratified mean",
                    NDays = (31-1)+18,
                    Area = shelfbreak_area+offshore_area,
                    NPerTow = stratcatchrate,
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
offshorecatchrate<- sum(offshorepoints$BLUEFIN_1, offshorepoints$BLUEFIN_2, na.rm=T)/length(offshorepoints$STATION)
# extract shelfbreak points and calculate mean:
shelfbreakpoints<- points_sps[!is.na(over(points_sps, shelfbreak_sps)),]
shelfbreakmean<- mean(shelfbreakpoints$Abundance)
shelfbreakmeanPos<- mean(shelfbreakpoints$Abundance[shelfbreakpoints$Abundance>0])
shelfbreakcatchrate<- sum(shelfbreakpoints$BLUEFIN_1, shelfbreakpoints$BLUEFIN_2, na.rm=T)/length(shelfbreakpoints$STATION)
if (is.na(shelfbreakmeanPos)){shelfbreakmeanPos<- 0}
if (is.na(shelfbreakcatchrate)){shelfbreakcatchrate<- 0}
# calculate the stratified means:
stratMean<- 1/(shelfbreak_area+offshore_area)*(shelfbreakmean*shelfbreak_area+offshoremean*offshore_area)
stratMeanPos<- 1/(shelfbreak_area+offshore_area)*(shelfbreakmeanPos*shelfbreak_area+offshoremeanPos*offshore_area)
stratcatchrate<- 1/(shelfbreak_area+offshore_area)*(shelfbreakcatchrate*shelfbreak_area+offshorecatchrate*offshore_area)
# add row to table:
newrow<- data.frame(Configuration = "HB1303, July 2-Aug 1, stratified mean",
                    NDays = (31-1)+1,
                    Area = shelfbreak_area+offshore_area,
                    NPerTow = stratcatchrate,
                    MeanAbund = stratMean,
                    MeanAbundPosStn = stratMeanPos)
larvAbundSensTbl<- rbind(larvAbundSensTbl, newrow)


### Write out the table as a csv, and then can edit the text in Excel:
write.csv(larvAbundSensTbl, 
          file = here("results", "LarvalAbundanceSensitivity.csv"), 
          row.names=F)
