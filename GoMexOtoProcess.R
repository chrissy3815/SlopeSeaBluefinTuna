## GoMex otoliths from Katie:

setwd('/Users/chrissy/JointProgram/SlopeSea/GoMex_oto_fromKatie/')

# load the data from the first read
read1<- read.csv("GoMexOtolithReads1.csv")

# load the data from the second read
read2<- read.csv("GoMexOtolithReads2.csv")

# Check which fish have Read1 and Read2 agreeing within 1d
agecheck<- merge(read1[,c("Fish","Age")], read2[,c("Fish", "Age")], by="Fish")
agecheck$d <- abs(agecheck$Age.y-agecheck$Age.x)
agecheck$within1 <- 0
agecheck$within1[agecheck$d<2] <- 1

# If the first and second read are within 1, keep second:
# (We are assuming that skill increased with time spent reading)
keepreads <- agecheck$Fish[agecheck$within1==1]
toprocess <- read2[read2$Fish %in% keepreads, ]
# correct the Age to be Increments:
toprocess$Increments<- toprocess$Age-1

# read in size data
require(xlsx)
seamapdata<- read.xlsx('BFT2016_SEAMAP_Metadata_forCHernandez_27Mar2018.xlsx', 1)
seamapdata<- seamapdata[,c('ELH_ID', 'SL_mm_EtOH', 'Slat', 'Slon')]
head(seamapdata)

# make matching columns in toprocess and seamapdata that contain the 4-digit fish ID:
toprocess$Fish<- as.character(toprocess$Fish)
sampleid<- strsplit(toprocess$Fish, '_')
Cruise<- sapply(sampleid, '[', 1)
Fish<- sapply(sampleid, '[', 2)
toprocess$Fish<- sapply(Fish, substr, start=1, stop=4)

seamapdata$Fish<- sapply(seamapdata$ELH_ID, substr, start=6, stop=9)

# append size and lat/lon to the otolith reads:
toprocess<- merge(toprocess, seamapdata, by='Fish')
toprocess<- toprocess[,-2] # get rid of the 'n' column
head(toprocess) # check that all the otolith data is still there.

# plot of size at radius (to R graphics device)
plot(toprocess$SL_mm_EtOH, toprocess$Radius, pch=19)

# plot of size at age (to file)
png(filename='GoMex_oto_fromKatie/GoM2016SizeAtAge_20180517_dashedonly.png', height=6.5, width=7.5, units= 'in', res=300)
toplot<- toprocess[toprocess$Increments>0,]
plot(toplot$Increments, toplot$SL_mm_EtOH, xlab='Daily Increments', 
     ylab='Standard Length (mm)', pch=19, cex=1.25, cex.lab=1.5, cex.axis=1.5)
model<- lm(SL_mm_EtOH~Increments, data=toplot)
summary(model)
exes<- 1:14
whys<- summary(model)$coefficients[1,1]+summary(model)$coefficients[2,1]*exes
lines(exes, whys, lty=2)
text(1, 6.75, "SL=2.5+0.36*DI", pos=4)
# line for fish up to 8 increments, to match with the SS data:
model2<- lm(SL_mm_EtOH~Increments, data=toplot[toplot$Increments<9,])
summary(model2)
exes2<- 1:8
whys2<- summary(model2)$coefficients[1,1]+summary(model2)$coefficients[2,1]*exes2
lines(exes2, whys2, lwd=2)
text(1, 7.5, "SL=2.47+0.46*DI", pos=4, cex=1.5)
dev.off()

# radius at age (to get rid of shrinkage differences)
png(filename='GoMex_oto_fromKatie/GoM2016RadiusAtAge_20180515.png', height=6.5, width=7.5, units= 'in', res=300)
toplot<- toprocess[toprocess$Increments>0,]
plot(toplot$Increments, toplot$Radius, xlab='Daily Increments', 
     ylab='Otolith Radius (um)', pch=19)
model<- lm(Radius~Increments, data=toplot)
summary(model)
exes<- 0:14
whys<- summary(model)$coefficients[1,1]+summary(model)$coefficients[2,1]*exes
lines(exes, whys)
text(1, 55, "Radius=8.13+3.35*DI", pos=4, cex=1.5)
dev.off()

# line for fish up to 8 increments, to match with the SS data:
model2<- lm(Radius~Increments, data=toplot[toplot$Increments<9,])
summary(model2)
exes2<- 1:8
whys2<- summary(model2)$coefficients[1,1]+summary(model2)$coefficients[2,1]*exes2
lines(exes2, whys2, lty=2)
text(1, 50, "Radius=7.414+3.50*DI", pos=4)

# calculate increment widths:
I<- which(names(toprocess)=='ToRing14') #bc the max age is 13
J<- which(names(toprocess)=='ToRing2') #bc I mark the edge of the core
# radii to outer and inner edges of each increment
outer<- toprocess[,J:I]
inner<- toprocess[,(J-1):(I-1)]
# increment width
incwidthGOM<- outer-inner
incwidthGOM<- cbind(toprocess$Fish, incwidthGOM)
names(incwidthGOM)<- c("Fish", "Inc1", "Inc2", "Inc3", "Inc4", "Inc5", "Inc6", "Inc7",
                      "Inc8", "Inc9", "Inc10", "Inc11", "Inc12", "Inc13")

# Make a second increment width dataframe with only fish that have 8 increments or fewer
subdata<- toprocess[toprocess$Increments<9,]
# calculate increment widths:
I<- which(names(subdata)=='ToRing9') #bc the max age is 8
J<- which(names(subdata)=='ToRing2') #bc I mark the edge of the core
# radii to outer and inner edges of each increment
outer<- subdata[,J:I]
inner<- subdata[,(J-1):(I-1)]
# increment width
incwidthGOM2<- outer-inner
incwidthGOM2<- cbind(subdata$Fish, incwidthGOM2)
names(incwidthGOM2)<- c("Fish", "Inc1", "Inc2", "Inc3", "Inc4", "Inc5", "Inc6", "Inc7",
                       "Inc8")

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
