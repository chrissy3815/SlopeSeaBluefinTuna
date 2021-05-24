
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

seamap_bongos<- read.csv('data/SEAMAP_SPRING2016_TUNA_LARVAE_bongos.csv')
seamap_shallow_bongos<- read.csv('data/SEAMAP_SPRING2016_TUNA_LARVAE_shallow_bongos.csv')
seamap_zeros<- seamap_bongos[seamap_bongos$TAXON=="NO TUNA LARVAE CAUGH",]
seamap_bluefin<- seamap_bongos[seamap_bongos$TAXON=="Thunnus thynnus",]
seamap_bongos<- rbind(seamap_zeros, seamap_bluefin)
seamap_bongos$NetDepth<- 200
I<- which(seamap_bongos$STA_DPTH<200)
seamap_bongos$NetDepth[I]<- seamap_bongos$STA_DPTH[I]-5
seamap_bongos$Abundance<- 10*seamap_bongos$TOT_LARVAE/seamap_bongos$VOL_FILT*seamap_bongos$NetDepth
# remove the repeated rows:
I<- which(seamap_bongos$P_STA_NO==63)
seamap_bongos<- seamap_bongos[-I[1],]
I<- which(seamap_bongos$P_STA_NO==66)
seamap_bongos<- seamap_bongos[-I[1],]

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
I<- which(seamap_bongos$Abundance>0 & seamap_bongos$Abundance<75)
points(seamap_bongos$STA_LON[I], seamap_bongos$STA_LAT[I],
       cex=sqrt(seamap_bongos$Abundance[I])/2, lwd=2)
I<- which(seamap_bongos$Abundance>75)
points(seamap_bongos$STA_LON[I], seamap_bongos$STA_LAT[I],
       cex=4.5, lwd=2)
I<- which(seamap_bongos$Abundance==0)
points(seamap_bongos$STA_LON[I], seamap_bongos$STA_LAT[I], pch=3, lwd=1.5, cex=0.75)
# add a legend
legend('topleft', legend=c('0', '5', '10', '25', '>75'), 
       pch=c(3, 1, 1, 1, 1), pt.cex=c(0.75, sqrt(5)/2, sqrt(10)/2, 2.5, 4.5),
       title="N per 10 m2")
dev.off()

