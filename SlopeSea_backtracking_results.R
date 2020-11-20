## Plot the results from Irina's backtracking work:

# required packages:
library(rmatio)
library(oce)
library(ocedata)
library(ncdf4)
library(rgdal)

# Color palette that is colorblind-accessible:
# Grey, orange, light blue, green, yellow, dark blue, red, pink, black
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
               "#D55E00", "#CC79A7", "#000000")

# function for getting matlab dates out:
Matlab2Rdate <- function(val) as.Date(val - 1, origin = '0000-01-01') 

# read in one of the trajectory files:
test<- read.mat('data/MABGOM2_2016_tr/MABGOM2_2016_tr1.mat')
dates<- Matlab2Rdate(test$tmesh_b)
test2<- read.mat('data/MABGOM2_2016_tr/MABGOM2_2016_tr10.mat')

# read in the backtracking input:
sim_input_aged<- read.csv('results/SlopeSea2016_agedlarvae_forbacktracking_SENT.csv')
sim_input_est<- read.csv('results/SlopeSea2016_estimatedlarvae_forbacktracking_SENT.csv')

# Load in the coastlines
data(coastlineWorldFine, package="ocedata")
# read in the kickass bathymetry
ncid<- nc_open('data/SlopeSeaBathymetry/GEBCO_2014_2D_-85.1699_25.7039_-55.3641_44.7816.nc')
lat<- ncvar_get(ncid, varid='lat')
lon<- ncvar_get(ncid, varid='lon')
elev<- ncvar_get(ncid, varid='elevation')
nc_close(ncid)
# load in the Slope Sea boundary:
SS_polygon<- read.mat('data/LonLatSlopeSea.mat') # Irina's SlopeSea boundary
SS_shp<- readOGR('data/gazetteer_polygon/gazetteer_polygon.shp') #Dave's SlopeSea boundary

# make a plot of the backtracked locations:
setEPS()
postscript('results/SS2016_trajectory_map_backwards.eps', height=5.5, width=7)
plot(coastlineWorldFine, longitudelim=c(-60, -77), latitudelim=c(33, 43), 
     xlab='Longitude', ylab='Latitude')
contour(lon, lat, elev, levels=c(-100, -200, -1000, -2000), add=TRUE, 
        drawlabels=FALSE, lwd=0.75, col=cbPalette[1])
# polygon(SS_polygon$lonss, SS_polygon$latss, border= cbPalette[5], lwd=1.5) # for adding Irina's boundary
plot(SS_shp, border=cbPalette[2], lwd=1.5, add=TRUE)

for (i in 1:68){
  thisfile<- paste('data/MABGOM2_2016_tr/MABGOM2_2016_tr',i,'.mat', sep='')
  traj<- read.mat(thisfile)
  points(traj$lontr_b[1], traj$lattr_b[1], pch=19)
  lines(traj$lontr_b, traj$lattr_b, col=cbPalette[6])
  n<- length(traj$lontr_b)
  points(traj$lontr_b[n], traj$lattr_b[n], pch=24, col=cbPalette[6], bg=cbPalette[6])
}
legend("bottomright", legend=c('collection', 'backtracked origin'),
       col=c('black', cbPalette[6]), pch=c(19,24), 
       pt.bg=c('black', cbPalette[6]))
dev.off()

# make a plot of the forward-tracked locations:
setEPS()
postscript('results/SS2016_trajectory_map_forwards.eps', height=5.5, width=7)
plot(coastlineWorldFine, longitudelim=c(-60, -77), latitudelim=c(33, 43), 
     xlab='Longitude', ylab='Latitude')
contour(lon, lat, elev, levels=c(-100, -200, -1000, -2000), add=TRUE, 
        drawlabels=FALSE, lwd=0.75, col=cbPalette[1])
# polygon(SS_polygon$lonss, SS_polygon$latss, border= cbPalette[5], lwd=1.5) # for adding Irina's boundary
plot(SS_shp, border=cbPalette[2], lwd=1.5, add=TRUE)

for (i in 1:68){
  thisfile<- paste('data/MABGOM2_2016_tr/MABGOM2_2016_tr',i,'.mat', sep='')
  traj<- read.mat(thisfile)
  points(traj$lontr_b[1], traj$lattr_b[1], pch=19)
  lines(traj$lontr_f, traj$lattr_f, col=cbPalette[3], lwd=0.5)
  n<- length(traj$lontr_f)
  if (is.na(traj$lattr_f[n])){
    I<- which(!is.na(traj$lattr_f))
    n<- I[length(I)]
  }
  points(traj$lontr_f[n], traj$lattr_f[n], pch=25, col=cbPalette[3], bg=cbPalette[3])
}
legend("bottomright", legend=c('collection', '25 days post spawn'),
       col=c('black', cbPalette[3]), pch=c(19,25), 
       pt.bg=c('black', cbPalette[3]))
dev.off()


# the ones that left the domain:
leftdomain<- c(25, 44, 45, 46, 47, 48, 49)
setEPS()
postscript('results/SS2016_traj_domainEdge.eps', height=5.5, width=7)
plot(coastlineWorldFine, longitudelim=c(-60, -77), latitudelim=c(33, 43),
     xlab='Longitude', ylab='Latitude')
contour(lon, lat, elev, levels=c(-100, -200, -1000, -2000), add=TRUE, 
        drawlabels=FALSE, lwd=0.75, col=cbPalette[1])
polygon(SS_polygon$lonss, SS_polygon$latss, border= cbPalette[5], lwd=1.5)

for (i in leftdomain){
  thisfile<- paste('data/MABGOM2_2016_tr/MABGOM2_2016_tr',i,'.mat', sep='')
  traj<- read.mat(thisfile)
  points(traj$lontr_b[1], traj$lattr_b[1], pch=19)
  lines(traj$lontr_b, traj$lattr_b, col=cbPalette[6])
  n<- length(traj$lontr_b)
  points(traj$lontr_b[n], traj$lattr_b[n], pch=24, col=cbPalette[6], bg=cbPalette[6])
  lines(traj$lontr_f, traj$lattr_f, col=cbPalette[3], lwd=0.5)
  n<- length(traj$lontr_f)
  if (is.na(traj$lattr_f[n])){
    I<- which(!is.na(traj$lattr_f))
    n<- I[length(I)]
  }
  points(traj$lontr_f[n], traj$lattr_f[n], pch=25, col=cbPalette[3], bg=cbPalette[3])
}
legend("bottomright", legend=c('collection', 'backtracked origin', '25 days post spawn'),
       col=c('black', cbPalette[6], cbPalette[3]), pch=c(19,24,25), 
       pt.bg=c('black', cbPalette[6], cbPalette[3]))
dev.off()

