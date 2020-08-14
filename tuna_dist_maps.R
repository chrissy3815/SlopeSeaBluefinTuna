## Plot 2013 and 2016 Slope Sea data

setwd('/Users/chrissy/JointProgram/SlopeSea/')
library(oce)
library(ocedata)

# Load in the coastlines
data(coastlineWorldFine, package="ocedata")
# read in the kickass bathymetry
library(ncdf4)
ncid<- nc_open('Bathymetry_GEBCO_2014_2D_-85.1699_25.7039_-55.3641_44.7816.nc')
print(ncid)
lat<- ncvar_get(ncid, varid='lat')
lon<- ncvar_get(ncid, varid='lon')
elev<- ncvar_get(ncid, varid='elevation')
nc_close(ncid)

# 2016 data:
tuna2016<- read.csv('2016_data/2016MeasuredFish.csv', header=TRUE)
head(tuna2016)
tuna2016<- aggregate(Fish~Cruise+Station, data=tuna2016, FUN=length)
# Read in Station Locations
stations2016<- read.csv('2016_data/2016StationLocations.csv', header=TRUE)
# fix the latitude and longitude:
LatitudeDeg<- floor(stations2016$lat/100)
LatitudeMin<- stations2016$lat-LatitudeDeg*100
stations2016$LatDec<- LatitudeDeg+LatitudeMin/60
LongitudeDeg<- floor(stations2016$lon/100)
LongitudeMin<- stations2016$lon-LongitudeDeg*100
stations2016$LonDec<- -(LongitudeDeg+LongitudeMin/60)
names(stations2016)<- c('Cruise', 'Time', 'lat', 'lon', 'Station', 'LatDec', 'LonDec')
# merge the lat/lons with the fish data
toplot<- merge(tuna2016, stations2016[,c('Cruise', 'Station', 'LatDec', 'LonDec')])
# manually remove the duplicate rows:
toplot<- toplot[-c(8, 10,12),]
toplot
# get the other station locations:
uniquestations<- unique(stations2016[,c('Cruise', 'Station')])
I<- which(uniquestations$Station %in% toplot$Station)
J<- as.numeric(row.names(uniquestations))
J<- J[-I]
zeros<- stations2016[J,c('LatDec', 'LonDec')]

plot(coastlineWorldFine, longitudelim=c(-64, -76), latitudelim=c(34, 42))
contour(lon, lat, elev, levels=c(-100, -200, -1000, -2000), add=TRUE, 
        drawlabels=FALSE, lwd=0.75, col='dark grey')
points(zeros$LonDec, zeros$LatDec, cex=0.8, lwd=2)
onefour<- toplot[toplot$Fish>0 & toplot$Fish<=4,]
fivenine<- toplot[toplot$Fish>4 & toplot$Fish<10,]
tenplus<- toplot[toplot$Fish>9,]
points(onefour$LonDec, onefour$LatDec, cex=2, col='blue', lwd=2)
points(fivenine$LonDec, fivenine$LatDec, cex=4, col='blue', lwd=2)
points(tenplus$LonDec, tenplus$LatDec, cex=6, col='blue', lwd=2)
# add a legend
legend(-65, 37, legend=c('unknown', '1-4', '5-9', '10+'), 
       col=c('black', 'blue', 'blue', 'blue'), pch=c(1, 1, 1, 1), 
       pt.cex=rep(NA, 4), bty='n', 
       title='N per tow')
legend(-65.3, 36.45, legend=rep(NA, 4), 
       col=c('black', 'blue', 'blue', 'blue'), pch=c(1, 1, 1, 1), 
       pt.cex=c(0.8, 2, 4, 6), pt.lwd=c(2,2,2,2), bty='n')



# 2013 data:
tuna2013<- read.csv('2013_data/RichardsonPaper_TableS1.csv', header=TRUE)
head(tuna2013)

plot(coastlineWorldFine, longitudelim=c(-64, -76), latitudelim=c(34, 42))
contour(lon, lat, elev, levels=c(-100, -200, -1000, -2000), add=TRUE, 
        drawlabels=FALSE, lwd=0.75, col='dark grey')
zeros<- tuna2013[tuna2013$nBluefin==0,]
points(zeros$Longitude, zeros$Latitude, cex=0.8, lwd=2)
onefour<- tuna2013[tuna2013$nBluefin>0 & tuna2013$nBluefin<=4,]
fivenine<- tuna2013[tuna2013$nBluefin>4 & tuna2013$nBluefin<10,]
tenplus<- tuna2013[tuna2013$nBluefin>9,]
points(onefour$Longitude, onefour$Latitude, cex=2, col='blue', lwd=2)
points(fivenine$Longitude, fivenine$Latitude, cex=4, col='blue', lwd=2)
points(tenplus$Longitude, tenplus$Latitude, cex=6, col='blue', lwd=2)
# add a legend
legend(-65, 37, legend=c('0', '1-4', '5-9', '10+'), 
       col=c('black', 'blue', 'blue', 'blue'), pch=c(1, 1, 1, 1), 
       pt.cex=rep(NA, 4), bty='n', 
       title='N per tow')
legend(-65.3, 36.45, legend=rep(NA, 4), 
       col=c('black', 'blue', 'blue', 'blue'), pch=c(1, 1, 1, 1), 
       pt.cex=c(0.8, 2, 4, 6), pt.lwd=c(2,2,2,2), bty='n')
