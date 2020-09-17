
# required packages:
library(oce)
library(ocedata)
library(ncdf4)

# run the Slope Sea processing script:
source('SlopeSeaOtoProcess.R')
# run the Gulf of Mexico processing script:
source('GoMexOtoProcess.R')

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

png(filename='results/GOM2016_map.png', height=5.5, width=7, units='in', res=300)
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

png(filename='results/SS2016_aged_map_ExcludeLapillae.png', height=5.5, width=7, units='in', res=300)
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

# I think the full set of tuna larvae ID'ed from Slope Sea in 2016 is comprised of 
# the larvae that Dave and I ID'ed (compiled in the usalengths)
# and the larvae that were ID'ed in Poland (but in this case I only care about counts)

# read in the poland data again:
polanddata<- read.csv('data/HB1603_GU1608_IchData_7Nov2019.csv')
# but we only want the bluefin tuna:
polanddata<- polanddata[polanddata$TAXA_NAME=="Thunnus thynnus",]
# keep only the relevant columns:
polanddata<- polanddata[,c("CRUISE_NAME", "STATION", "GEAR", 
                           "LATITUDE", "LONGITUDE", "TOTAL_COUNT")]
head(polanddata)
names(polanddata)<- c("Cruise", "Station", "Gear" , "LatDec", "LonDec", "Nbluefin")
# switch longitudes to degW to match the other dataset:
polanddata$LonDec<- abs(polanddata$LonDec)
# eliminate repeated rows:
polanddata<- unique(polanddata)

# aggregate the usalengths data to get N per net:
usadata<- aggregate(Fish~Cruise+Station+Gear, data=usalengths, FUN=length)
# rename the "Fish" column:
names(usadata)[length(names(usadata))]<- "Nbluefin"
# do a merge to get the lat and lon back:
usadata<- merge(usadata, 
                unique(usalengths[,c("Cruise", "Station", "Gear", "LatDec", "LonDec")]))
# reorder the columns to match polanddata:
usadata<- usadata[,c("Cruise", "Station", "Gear" , "LatDec", "LonDec", "Nbluefin")]

# merge the usadata and polanddata:
tunacounts_SS<- merge(polanddata, usadata, all=TRUE)
# drop the 2B1 samples because we won't use those for abundance calculations:
I<- which(tunacounts_SS$Gear=="2B1")
tunacounts_SS<- tunacounts_SS[-I,]
# add the "Operation" column:
tunacounts_SS$Operation<- NA
I<- which(tunacounts_SS$Gear=="6B3I" | tunacounts_SS$Gear=="6B3" | tunacounts_SS$Gear=="6B3Z")
tunacounts_SS$Operation[I]<- "BON/CTD"
I<- which(tunacounts_SS$Gear=="2N3")
tunacounts_SS$Operation[I]<- "CTD/IKMT"
# combine the bongos:
tunacounts_SS<- aggregate(Nbluefin~Cruise+Station+Operation, data=tunacounts_SS,
                          FUN=sum)

# need to add the volume filtered:
SS2016_netdata<- read.csv('data/GU1608HB1603Net.csv')
# pull out the relevant columns:
SS2016_netdata<- SS2016_netdata[,c("CRUISE_NAME", "STATION", "GEAR", "GEAR_VOLUME_FILTERED")]
names(SS2016_netdata)<- c("Cruise", "Station", "Gear", "Vol_filtered")
# add the Operation column:
SS2016_netdata$Operation<- NA
I<- which(SS2016_netdata$Gear=="6B3I" | SS2016_netdata$Gear=="6B3" | SS2016_netdata$Gear=="6B3Z")
SS2016_netdata$Operation[I]<- "BON/CTD"
I<- which(SS2016_netdata$Gear=="2N3")
SS2016_netdata$Operation[I]<- "CTD/IKMT"
# we're going to combine the Bongo samples:
SS2016_netdata<- aggregate(Vol_filtered~Cruise+Station+Operation, 
                            data=SS2016_netdata, FUN=sum, na.rm=T)
# merge with the tunacounts_SS:
tunacounts_SS<- merge(tunacounts_SS, SS2016_netdata, all.x=T, all.y=F)

# need to add the sampling depth and lat/lon:
SS2016_eventdata<- read.csv('data/GU1608HB1603Event.csv')
# pull out the relevant columns:
SS2016_eventdata<- SS2016_eventdata[,c("CRUISE_NAME", "STATION", "OPERATION", 
                                       "TOW_MAXIMUM_DEPTH", "BOTTOM_DEPTH_MAX_WIRE_OUT",
                                       "LATITUDE", "LONGITUDE")]
names(SS2016_eventdata)<- c("Cruise", "Station", "Operation", "SamplingDepth",
                            "BottomDepth", "Latitude", "Longitude")
# merge with the tuna data:
tunacounts_SS<- merge(tunacounts_SS, SS2016_eventdata, all.x=T, all.y=F)

# Calculate abundance as N per 10m2:
tunacounts_SS$Abundance<- 10*tunacounts_SS$Nbluefin/tunacounts_SS$Vol_filtered*tunacounts_SS$SamplingDepth
# Calculate density as N per 100 m3:
tunacounts_SS$Density<- 100*tunacounts_SS$Nbluefin/tunacounts_SS$Vol_filtered


# make a map for Slope Sea abundance!        

# find the zero stations:
SS2016_eventdata<- read.csv('data/GU1608HB1603Event.csv')
# pull out the relevant columns:
SS2016_eventdata<- SS2016_eventdata[,c("CRUISE_NAME", "STATION", "OPERATION", 
                                       "TOW_MAXIMUM_DEPTH", "BOTTOM_DEPTH_MAX_WIRE_OUT",
                                       "LATITUDE", "LONGITUDE", "EVENT_DATE")]
names(SS2016_eventdata)<- c("Cruise", "Station", "Operation", "SamplingDepth",
                            "BottomDepth", "Latitude", "Longitude", "Date")
# get day and month out:
SSmoday<- strsplit(SS2016_eventdata$Date, "-")
SS2016_eventdata$Day<- sapply(SSmoday, '[', 1)
SS2016_eventdata$Month<- sapply(SSmoday, '[', 2)
# restrict to the bongo stations:
zerostations<- SS2016_eventdata[SS2016_eventdata$Operation=="BON/CTD",]
I<- which(zerostations$BottomDepth<1000)
zerostations<- zerostations[-I,]
I<- which(zerostations$Month=='AUG' & zerostations$Day>15)
zerostations<- zerostations[-I,]
I<- which(zerostations$Month=="MAY")
zerostations<- zerostations[-I,]
temp<- tunacounts_SS[tunacounts_SS$Operation=="BON/CTD",]
I<- which(zerostations$Station %in% temp$Station)
zerostations<- zerostations[-I,]

# Load in the coastlines
data(coastlineWorldFine, package="ocedata")
# read in the kickass bathymetry
ncid<- nc_open('data/SlopeSeaBathymetry/GEBCO_2014_2D_-85.1699_25.7039_-55.3641_44.7816.nc')
lat<- ncvar_get(ncid, varid='lat')
lon<- ncvar_get(ncid, varid='lon')
elev<- ncvar_get(ncid, varid='elevation')
nc_close(ncid)

# plot as png:
png(filename='results/SS2016_abund_map.png', height=5.5, width=7, units='in', res=300)
plot(coastlineWorldFine, longitudelim=c(-64, -76), latitudelim=c(34, 42))
contour(lon, lat, elev, levels=c(-100, -200, -1000, -2000), add=TRUE, 
        drawlabels=FALSE, lwd=0.75, col='dark grey')
points(tunacounts_SS$Longitude, tunacounts_SS$Latitude,
       cex=tunacounts_SS$Abundance/5, lwd=2)
points(zerostations$Longitude, zerostations$Latitude, pch=3, lwd=1.5)
# add a legend
legend(-65, 37, legend=c('0', '2', '5', '10', '20'), 
       col=c('black', 'black', 'blue', 'blue', 'blue'), pch=rep(1,5), 
       pt.cex=rep(NA, 4), bty='n', 
       title='N per 10 m2')
legend(-65.3, 36.45, legend=rep(NA, 5), pch=c(3, 1, 1, 1, 1), 
       pt.cex=c(1, 2/5, 1, 2, 4), pt.lwd=c(1.5,2,2,2,2), bty='n')

# Add an inset map:
par(new=TRUE)
par(mar=c(1,1,1,1))
plotInset("bottomleft", 
          expr={plot(coastlineWorldFine, longitudelim=c(-60, -80), latitudelim=c(20, 50), 
                     inset=TRUE, bg='white', axes=F, lwd=0.5)
                  polygon(x=c(-63, -63, -77, -77), y=c(34, 42, 42, 34), 
                          density=NULL, border='blue')
          })
dev.off()

# plot as eps:
setEPS()
postscript('results/SS2016_abund_map.eps', height=5.5, width=7)
plot(coastlineWorldFine, longitudelim=c(-64, -76), latitudelim=c(34, 42))
contour(lon, lat, elev, levels=c(-100, -200, -1000, -2000), add=TRUE, 
        drawlabels=FALSE, lwd=0.75, col='dark grey')
points(tunacounts_SS$Longitude, tunacounts_SS$Latitude,
       cex=sqrt(tunacounts_SS$Abundance)/1.2, lwd=2)
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

png(filename='results/GOM2016_abund_map.png', height=5.5, width=6.3, units='in', res=300)
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

