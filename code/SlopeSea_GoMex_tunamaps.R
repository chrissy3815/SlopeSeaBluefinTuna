
# required packages:
library(oce)
library(ocedata)
library(ncdf4)
library(here)

# run the Slope Sea processing scripts:
source(here("code", "SlopeSeaOtoProcess.R"))
source(here("code","SlopeSeaAbundProcess.R"))
# run the Gulf of Mexico processing script:
source(here('code','GoMexOtoProcess.R'))

## make a map of GoMex aged larvae:
# Load in the coastlines
data(coastlineWorldFine, package="ocedata")
# read in the kickass bathymetry
ncid<- nc_open('data/GoMexBathymetry/GEBCO_2014_2D_-101.9175_17.1845_-76.1408_31.6019.nc')
lat<- ncvar_get(ncid, varid='lat')
lon<- ncvar_get(ncid, varid='lon')
elev<- ncvar_get(ncid, varid='elevation')
nc_close(ncid)

# tuna sampled for otoliths:
OR317<- GOM_oto_data[,c("Fish", "Slat", "Slon")]
head(OR317)
toplot<- aggregate(Fish~Slat+Slon, data=OR317, FUN = length)

setEPS()
postscript('results/GOM2016_aged_map.eps', height=5.5, width=7)
plot(coastlineWorldFine, longitudelim=c(-78, -99), latitudelim=c(25, 30), 
     xlab='Longitude', ylab='Latitude')
contour(lon, lat, elev, levels=c(-100, -200, -1000, -2000), add=TRUE, 
        drawlabels=FALSE, lwd=0.75, col='dark grey')
onefour<- toplot[toplot$Fish>0 & toplot$Fish<=4,]
fivenine<- toplot[toplot$Fish>4 & toplot$Fish<10,]
tenplus<- toplot[toplot$Fish>9,]
points(onefour$Slon, onefour$Slat, cex=1.5, col='blue', lwd=2)
points(fivenine$Slon, fivenine$Slat, cex=3, col='blue', lwd=2)
points(tenplus$Slon, tenplus$Slat, cex=4.5, col='blue', lwd=2)
# add a legend
polygon(c(-82.5, -82.5, -86.5, -86.5), c(22.5, 25.5, 25.5, 22.5), col='white')
legend(-86, 25.5, legend=c('1-4', '5-9', '10+'), 
       col=c('blue', 'blue', 'blue'), pch=c(1, 1, 1), 
       pt.cex=rep(NA, 3), bty='n', 
       title='N for aging')
legend(-86.1, 25, legend=rep(NA, 3), 
       col=c('blue', 'blue', 'blue'), pch=c(1, 1, 1), 
       pt.cex=c(1.5, 3, 4.5), pt.lwd=c(2,2,2), bty='n')
dev.off()

# Make a map of Slope Sea aged larvae:

# aggregate aged larvae by station:
SS_aged_larvae<- aggregate(Fish~Cruise+Station+Gear, data=SS_oto_data, FUN=length)
# add lat and lon:
SS2016_eventdata<- read.csv('data/GU1608HB1603Event.csv')
# pull out the relevant columns:
SS2016_eventdata<- SS2016_eventdata[,c("CRUISE_NAME", "STATION", "OPERATION", "LATITUDE", "LONGITUDE")]
names(SS2016_eventdata)<- c("Cruise", "Station", "Operation", "Latitude", "Longitude")
latitude<- aggregate(Latitude~Cruise+Station, data=SS2016_eventdata, FUN=mean)
longitude<- aggregate(Longitude~Cruise+Station, data=SS2016_eventdata, FUN=mean)
# merge the lat and lon back into tuna data:
SS_aged_larvae<- merge(SS_aged_larvae, latitude, all.x=T, all.y=F)
SS_aged_larvae<- merge(SS_aged_larvae, longitude, all.x=T, all.y=F)
head(SS_aged_larvae)

# Load in the coastlines
data(coastlineWorldFine, package="ocedata")
# read in the kickass bathymetry
ncid<- nc_open('data/SlopeSeaBathymetry/GEBCO_2014_2D_-85.1699_25.7039_-55.3641_44.7816.nc')
lat<- ncvar_get(ncid, varid='lat')
lon<- ncvar_get(ncid, varid='lon')
elev<- ncvar_get(ncid, varid='elevation')
nc_close(ncid)

setEPS()
postscript('results/SS2016_aged_map_ExcludeLapillae.eps', height=5.5, width=7)
plot(coastlineWorldFine, longitudelim=c(-64, -76), latitudelim=c(34, 42))
contour(lon, lat, elev, levels=c(-100, -200, -1000, -2000), add=TRUE, 
        drawlabels=FALSE, lwd=0.75, col='dark grey')
onefour<- SS_aged_larvae[SS_aged_larvae$Fish>0 & SS_aged_larvae$Fish<=4,]
fivenine<- SS_aged_larvae[SS_aged_larvae$Fish>4 & SS_aged_larvae$Fish<10,]
tenplus<- SS_aged_larvae[SS_aged_larvae$Fish>9,]
points(onefour$Longitude, onefour$Latitude, cex=1.5, col='blue', lwd=2)
points(fivenine$Longitude, fivenine$Latitude, cex=3, col='blue', lwd=2)
points(tenplus$Longitude, tenplus$Latitude, cex=4.5, col='blue', lwd=2)
# add a legend
polygon(c(-64, -64, -67.1, -67.1), c(34.75, 37.1, 37.1, 34.75), col='white')
legend(-66.5, 37, legend=c('1-4', '5-9', '10+'), 
       col=c('blue', 'blue', 'blue'), pch=c(1, 1, 1), 
       pt.cex=rep(NA, 3), bty='n', 
       title='N for aging')
legend(-66.8, 36.45, legend=rep(NA, 3), 
       col=c('blue', 'blue', 'blue'), pch=c(1, 1, 1), 
       pt.cex=c(1.5, 3, 4.5), pt.lwd=c(2,2,2), bty='n')
dev.off()


#### Okay, and the main event, a map of Slope Sea larval distributions:

# restrict to bongo stations through Aug 15:
I<- which(all_bongo_stns$Month=="AUG" & all_bongo_stns$Day>15)
all_bongo_stns<- all_bongo_stns[-I,]

# pull out the zero stations:
zerostations<- all_bongo_stns[all_bongo_stns$Abundance==0 & all_bongo_stns$BottomDepth>=1000,]
# pull out the positive stations:
catchstations<- all_bongo_stns[all_bongo_stns$Abundance>0,]

# Load in the coastlines
data(coastlineWorldFine, package="ocedata")
# read in the kickass bathymetry
ncid<- nc_open('data/SlopeSeaBathymetry/GEBCO_2014_2D_-85.1699_25.7039_-55.3641_44.7816.nc')
lat<- ncvar_get(ncid, varid='lat')
lon<- ncvar_get(ncid, varid='lon')
elev<- ncvar_get(ncid, varid='elevation')
nc_close(ncid)

# plot as eps:
setEPS()
postscript('results/SS2016_abund_map_through_Aug15_1000m.eps', height=5.5, width=7)
plot(coastlineWorldFine, longitudelim=c(-64, -76), latitudelim=c(34, 42))
contour(lon, lat, elev, levels=c(-100, -200, -1000, -2000), add=TRUE, 
        drawlabels=FALSE, lwd=0.75, col='dark grey')
points(catchstations$Longitude, catchstations$Latitude,
       cex=sqrt(catchstations$Abundance)/1.2, lwd=2)
points(zerostations$Longitude, zerostations$Latitude, pch=3, lwd=1.5)
# add a legend
legend(-65, 37, legend=c('0', '2', '5', '10', '20'), 
       col=c('black', 'black', 'blue', 'blue', 'blue'), pch=rep(1,5), 
       pt.cex=rep(NA, 4), bty='n', 
       title='N per 10 m2')
legend(-65.3, 36.45, legend=rep(NA, 5), pch=c(3, 1, 1, 1, 1), 
       pt.cex=c(1, sqrt(2)/1.2, sqrt(5)/1.2, sqrt(10)/1.2, sqrt(20)/1.2), pt.lwd=c(1.5,2,2,2,2), bty='n')

# Add an inset map:
par(new=TRUE)
par(mar=c(1,1,1,1))
plotInset("bottomleft", 
          expr={plot(coastlineWorldFine, longitudelim=c(-60, -80), latitudelim=c(20, 49), 
                     inset=TRUE, bg='white', axes=F, lwd=0.5)
                polygon(x=c(-63, -63, -77, -77), y=c(34, 42, 42, 34), 
                        density=NULL, border='blue')
          })
dev.off()


## SEAMAP abundance in 2016:
# run the processing script:
source(here("code","GoMexAbundProcess.R"))

# Plot a map:
# Load in the coastlines
data(coastlineWorldFine, package="ocedata")
# read in the kickass bathymetry
ncid<- nc_open('data/GoMexBathymetry/GEBCO_2014_2D_-101.9175_17.1845_-76.1408_31.6019.nc')
lat<- ncvar_get(ncid, varid='lat')
lon<- ncvar_get(ncid, varid='lon')
elev<- ncvar_get(ncid, varid='elevation')
nc_close(ncid)

setEPS()
postscript('results/GOM2016_abund_map.eps', height=5.5, width=6.3)
plot(coastlineWorldFine, longitudelim=c(-78.5, -99), latitudelim=c(25, 30), 
     xlab='Longitude', ylab='Latitude')
contour(lon, lat, elev, levels=c(-100, -200, -1000, -2000), add=TRUE, 
        drawlabels=FALSE, lwd=0.75, col='dark grey')
I<- which(seamap_bluefin$Abundance>0 & seamap_bluefin$Abundance<75)
points(seamap_bluefin$STA_LON[I], seamap_bluefin$STA_LAT[I],
       cex=sqrt(seamap_bluefin$Abundance[I])/2, lwd=2)
I<- which(seamap_bluefin$Abundance>75)
points(seamap_bluefin$STA_LON[I], seamap_bluefin$STA_LAT[I],
       cex=4.5, lwd=2)
I<- which(seamap_bluefin$Abundance==0)
points(seamap_bluefin$STA_LON[I], seamap_bluefin$STA_LAT[I], pch=3, lwd=1.5, cex=0.75)
# add a legend
legend('topleft', legend=c('0', '5', '10', '25', '>75'), 
       pch=c(3, 1, 1, 1, 1), pt.cex=c(0.75, sqrt(5)/2, sqrt(10)/2, 2.5, 4.5),
       title="N per 10 m2")
dev.off()


## Map of AMAPPS Survey strata:

# Color palette that is colorblind-accessible:
# Grey, orange, light blue, green, yellow, dark blue, red, pink, black
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
               "#D55E00", "#CC79A7", "#000000")

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
proj4string(offshore_sps)<- CRS("+proj=longlat")
# now repeat for shelfbreak polygon
p<- Polygon(shelfbreak_polygon)
ps<- Polygons(list(p),1)
shelfbreak_sps<- SpatialPolygons(list(ps))
proj4string(shelfbreak_sps)<- CRS("+proj=longlat")

# Load in the survey lines:
offshore_lines<- read.csv("data/AMAPPS_Strata/Offshelf_Lines.csv")
offshore_lines<- offshore_lines[,c("Lon","Lat")]
l<- Line(offshore_lines)
ls<- Lines(list(l),1)
offshore_lines_sls<- SpatialLines(list(ls))
proj4string(offshore_lines_sls)<- CRS("+proj=longlat")
shelfbreak_lines<- read.csv("data/AMAPPS_Strata/Shelfbreak_Lines.csv")
shelfbreak_lines<- shelfbreak_lines[,c("Lon","Lat")]
l<- Line(shelfbreak_lines)
ls<- Lines(list(l),1)
shelfbreak_lines_sls<- SpatialLines(list(ls))
proj4string(shelfbreak_lines_sls)<- CRS("+proj=longlat")

# Load in the coastlines
data(coastlineWorldFine, package="ocedata")
# read in the kickass bathymetry
ncid<- nc_open(here('data',"SlopeSeaBathymetry",'GEBCO_2014_2D_-85.1699_25.7039_-55.3641_44.7816.nc'))
lat<- ncvar_get(ncid, varid='lat')
lon<- ncvar_get(ncid, varid='lon')
elev<- ncvar_get(ncid, varid='elevation')
nc_close(ncid)

# plot as eps:
setEPS()
postscript('results/AMAPPS_survey_strata.eps', height=5.5, width=7)
plot(coastlineWorldFine, longitudelim=c(-64, -76), latitudelim=c(34, 42.5))
contour(lon, lat, elev, levels=c(-100, -200, -1000, -2000), add=TRUE, 
        drawlabels=FALSE, lwd=0.75, col='dark grey')
plot(offshore_sps, border=cbPalette[7], lwd=1.5, add=T)
plot(shelfbreak_sps, border=cbPalette[7], lwd=1.5, add=T)
plot(offshore_lines_sls, col=cbPalette[3], lwd=1.25, add=T)
plot(shelfbreak_lines_sls, col=cbPalette[3], lwd=1.25, add=T)
dev.off()

