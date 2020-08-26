

# run the Slope Sea processing script:
source('SlopeSeaOtoProcess.R')
# run the Gulf of Mexico processing script:
source('GoMexOtoProcess.R')

## make a map of GoMex aged larvae:
library(oce)
library(ocedata)

# Load in the coastlines
data(coastlineWorldFine, package="ocedata")
# read in the kickass bathymetry
library(ncdf4)
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
SS_aged_larvae<- aggregate(Fish~Cruise+Station, data=SS_oto_data, FUN=length)
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
library(ncdf4)
ncid<- nc_open('data/SlopeSeaBathymetry/GEBCO_2014_2D_-85.1699_25.7039_-55.3641_44.7816.nc')
lat<- ncvar_get(ncid, varid='lat')
lon<- ncvar_get(ncid, varid='lon')
elev<- ncvar_get(ncid, varid='elevation')
nc_close(ncid)

png(filename='results/SS2016_aged_map.png', height=5.5, width=7, units='in', res=300)
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

# make a map for Slope Sea abundance!        

# Load in the coastlines
data(coastlineWorldFine, package="ocedata")
# read in the kickass bathymetry
library(ncdf4)
ncid<- nc_open('data/SlopeSeaBathymetry/GEBCO_2014_2D_-85.1699_25.7039_-55.3641_44.7816.nc')
lat<- ncvar_get(ncid, varid='lat')
lon<- ncvar_get(ncid, varid='lon')
elev<- ncvar_get(ncid, varid='elevation')
nc_close(ncid)

png(filename='results/SS2016_abund_map.png', height=5.5, width=7, units='in', res=300)
plot(coastlineWorldFine, longitudelim=c(-64, -76), latitudelim=c(34, 42))
contour(lon, lat, elev, levels=c(-100, -200, -1000, -2000), add=TRUE, 
        drawlabels=FALSE, lwd=0.75, col='dark grey')
points(tunacounts_SS$Longitude, tunacounts_SS$Latitude,
       cex=tunacounts_SS$Abundance/5, lwd=2)
# add a legend
legend(-65, 37, legend=c('2', '5', '10', '20'), 
       col=c('black', 'blue', 'blue', 'blue'), pch=c(1, 1, 1, 1), 
       pt.cex=rep(NA, 4), bty='n', 
       title='N per 10 m2')
legend(-65.3, 36.45, legend=rep(NA, 4), pch=c(1, 1, 1, 1), 
       pt.cex=c(2/5, 1, 2, 4), pt.lwd=c(2,2,2,2), bty='n')
dev.off()
