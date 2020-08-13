

## make a map of GoMex:
library(oce)
library(ocedata)

# Load in the coastlines
data(coastlineWorldFine, package="ocedata")
# read in the kickass bathymetry
library(ncdf4)
ncid<- nc_open('GoMexBathymetry/GEBCO_2014_2D_-101.9175_17.1845_-76.1408_31.6019.nc')
print(ncid)
lat<- ncvar_get(ncid, varid='lat')
lon<- ncvar_get(ncid, varid='lon')
elev<- ncvar_get(ncid, varid='elevation')
nc_close(ncid)

# tuna sampled for otoliths:
OR317<- toprocess[,c("Fish", "Slat", "Slon")]
head(OR317)
toplot<- aggregate(Fish~Slat+Slon, data=OR317, FUN = length)

plot(coastlineWorldFine, longitudelim=c(-81, -97), latitudelim=c(25, 30))
contour(lon, lat, elev, levels=c(-100, -200, -1000, -2000), add=TRUE, 
        drawlabels=FALSE, lwd=0.75, col='dark grey')
onefour<- toplot[toplot$Fish>0 & toplot$Fish<=4,]
fivenine<- toplot[toplot$Fish>4 & toplot$Fish<10,]
tenplus<- toplot[toplot$Fish>9,]
points(onefour$Slon, onefour$Slat, cex=1.5, col='blue', lwd=2)
points(fivenine$Slon, fivenine$Slat, cex=3, col='blue', lwd=2)
points(tenplus$Slon, tenplus$Slat, cex=4.5, col='blue', lwd=2)
# add a legend
polygon(c(-82.5, -82.5, -86.5, -86.5), c(20, 25.5, 25.5, 20), col='white')
legend(-86, 25.5, legend=c('1-4', '5-9', '10+'), 
       col=c('blue', 'blue', 'blue'), pch=c(1, 1, 1), 
       pt.cex=rep(NA, 3), bty='n', 
       title='N for aging')
legend(-86.1, 25, legend=rep(NA, 3), 
       col=c('blue', 'blue', 'blue'), pch=c(1, 1, 1), 
       pt.cex=c(1.5, 3, 4.5), pt.lwd=c(2,2,2), bty='n')