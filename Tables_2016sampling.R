## This is for making tables for my chapter/manuscript on bluefin tuna:
library(readxl)

# Larvae that were ID'ed in the US:
davelengths<- read.csv('data/2016MeasuredFish_CMH200811.csv')
davecounts<- aggregate(Fish~Cruise+Gear+Station, data=davelengths, FUN=length)
davecounts$Gear[davecounts$Gear=='2n3']<- '2N3'

chrissylengths<- read_excel('data/HB1603_6B3I_BFTlengths_20200811.xlsx', sheet = 1)
chrissycounts<- aggregate(Fish~Cruise+Station+Gear, data=chrissylengths, FUN=length)

usacounts<- rbind(davecounts, chrissycounts)
usacounts<- aggregate(Fish~Cruise+Gear+Station, data=usacounts, FUN=sum)
# rename the "Fish" column:
names(usacounts)[names(usacounts)=="Fish"]<- "Nbluefin"

# read in the poland data:
polandcounts<- read.csv('data/HB1603_GU1608_IchData_7Nov2019.csv')
# but we only want the bluefin tuna:
polandcounts<- polandcounts[polandcounts$TAXA_NAME=="Thunnus thynnus",]
# keep only the relevant columns:
polandcounts<- polandcounts[,c("CRUISE_NAME", "GEAR", "STATION", "TOTAL_COUNT")]
head(polandcounts)
names(polandcounts)<- c("Cruise", "Gear", "Station", "Nbluefin")
polandcounts<- unique(polandcounts)

# merge the usacounts and polandcounts:
tunacounts_SS<- merge(polandcounts, usacounts, all=TRUE)
# add the "Operation" column:
tunacounts_SS$Operation<- NA
I<- which(tunacounts_SS$Gear %in% c("6B3I", "6B3", "6B3Z"))
tunacounts_SS$Operation[I]<- "BON/CTD"
I<- which(tunacounts_SS$Gear=="2N3")
tunacounts_SS$Operation[I]<- "CTD/IKMT"
I<- which(tunacounts_SS$Gear=='2B1')
tunacounts_SS$Operation[I]<- "BB/BON"

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

# add lat, lon, sampling depth, and SST:
SS2016_eventdata<- read.csv('data/GU1608HB1603Event.csv')
# pull out the relevant columns:
SS2016_eventdata<- SS2016_eventdata[,c("CRUISE_NAME", "STATION", "OPERATION", 
                                       "TOW_MAXIMUM_DEPTH", "BOTTOM_DEPTH_MAX_WIRE_OUT",
                                       "LATITUDE", "LONGITUDE", "EVENT_DATE", 
                                       "SURFACE_TEMPERATURE")]
names(SS2016_eventdata)<- c("Cruise", "Station", "Operation", "SamplingDepth",
                            "BottomDepth", "Latitude", "Longitude", "Date", "SST")
# get day and month out:
SSmoday<- strsplit(SS2016_eventdata$Date, "-")
SS2016_eventdata$Day<- sapply(SSmoday, '[', 1)
SS2016_eventdata$Month<- sapply(SSmoday, '[', 2)
# merge with the tuna data:
tunacounts_SS<- merge(tunacounts_SS, SS2016_eventdata, all.x=T, all.y=F)

# Calculate abundance as N per 10m2:
tunacounts_SS$Abundance<- 10*tunacounts_SS$Nbluefin/tunacounts_SS$Vol_filtered*tunacounts_SS$SamplingDepth

# find the zero bongo stations:
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
# get the columns we want 
zerostations<- zerostations[, c("Month", "Day", "Cruise", "Station", "Operation", 
                                "Latitude", "Longitude", "BottomDepth", "SST")]
zerostations$Nbluefin<- 0
zerostations$Abundance<- 0

# table to write to a csv:
table_SS<- tunacounts_SS[,c("Date", "Cruise", "Station", "Operation", 
                            "Latitude", "Longitude", "BottomDepth", "SST", 
                            "Nbluefin", "Abundance")]
write.csv(table_SS, file='results/SlopeSeaSampling_Table.csv')

