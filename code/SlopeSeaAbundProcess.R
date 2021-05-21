## Code to calculate the Slope Sea station abundances for all bongo samples:

library(here)

# run the otolith processing script, which generates the list of all lengths:
source(here("code","SlopeSeaOtoProcess.R"))

# I think the full set of tuna larvae ID'ed from Slope Sea in 2016 is comprised of 
# the larvae that Dave and I ID'ed (compiled in the usalengths)
# and the larvae that were ID'ed in Poland (but in this case I only care about counts)

# December 7 2020: I updated the SlopeSeaOtoProcess script to include ALL larvae 
# in the all_lengths_SS object, including 3 larvae that we do not have length 
# data for AND the newly identified bluefin that were mislabeled as Scomber scombrus. 
# Rather than read in all the data a second time in this script, I updated this 
# section to use all_lengths_SS.

# Count how many rows we have for each net in the all_lengths_SS data to get Nbluefin per net:
# *** note that I'm using SST here so that "aggregate" counts the rows with NA for the 
# fish number and the rows with NA for the length...SST is just a column that's 
# guaranteed to have an entry so that I can get row counts.
tunacounts_SS<- aggregate(SST~Cruise+Station+Gear, data=all_lengths_SS, FUN=length)
# rename the "SST" column:
names(tunacounts_SS)[length(names(tunacounts_SS))]<- "Nbluefin"
# do a merge to get the lat and lon back:
tunacounts_SS<- merge(tunacounts_SS, 
                      unique(all_lengths_SS[,c("Cruise", "Station", "Gear", "LatDec", "LonDec")]))
# reorder the columns:
tunacounts_SS<- tunacounts_SS[,c("Cruise", "Station", "Gear" , "LatDec", "LonDec", "Nbluefin")]

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
# exclude the volume filtered for a few samples that weren't sorted:
# GU1608-st240; HB1603-sts 16,36,121,125
not_sorted<- data.frame(Cruise=c("GU1608", rep("HB1603",4)), 
                        Station=c(240, 16, 36, 121, 125), 
                        Gear = rep("6B3I", 5))
to_exclude<- vector()
for (i in 1:length(not_sorted$Cruise)){
  J<- which(SS2016_netdata$Cruise==not_sorted$Cruise[i] &
              SS2016_netdata$Station==not_sorted$Station[i] &
              SS2016_netdata$Gear==not_sorted$Gear[i])
  to_exclude<- c(to_exclude,J)
}
SS2016_netdata<- SS2016_netdata[-to_exclude,]
# now we're going to combine the volume filtered for paired Bongo samples:
SS2016_netdata<- aggregate(Vol_filtered~Cruise+Station+Operation, 
                           data=SS2016_netdata, FUN=sum, na.rm=T)
# merge with the tunacounts_SS:
tunacounts_SS<- merge(tunacounts_SS, SS2016_netdata, all.x=T, all.y=F)

# need to add the sampling depth, lat/lon, and date:
SS2016_eventdata<- read.csv(here('data','GU1608HB1603Event.csv'))
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

# merge with the tuna data:
tunacounts_SS<- merge(tunacounts_SS, SS2016_eventdata, all.x=T, all.y=F)

# Calculate abundance as N per 10m2:
tunacounts_SS$Abundance<- 10*tunacounts_SS$Nbluefin/tunacounts_SS$Vol_filtered*tunacounts_SS$SamplingDepth
# Calculate density as N per 100 m3:
tunacounts_SS$Density<- 100*tunacounts_SS$Nbluefin/tunacounts_SS$Vol_filtered

## Check a few stats on overall bongo sampling:
all_bongo_stns<- SS2016_eventdata[SS2016_eventdata$Operation=="BON/CTD",]
# breakdown of Gunther vs. Bigelow cruises:
sum(all_bongo_stns$Cruise=="GU1608")
sum(all_bongo_stns$Cruise=="HB1603")
#1. How many stations are shallower than 1000 m for each cruise?
# For the Gunther cruise, if we only include the sampling from June 17-20, then 
# there are 24 stations from that cruise, all deeper than 1000 m.
I<- which(all_bongo_stns$BottomDepth<1000 & all_bongo_stns$Cruise=="HB1603")
length(I) #44 stations shallower than 1000 m on Bigelow cruise
I<- which(all_bongo_stns$BottomDepth<1000 & all_bongo_stns$Cruise=="HB1603" 
          & all_bongo_stns$Month=="AUG" & all_bongo_stns$Day<=15)
length(I) # We exclude 12 stations from August with the 1000 m requirement
I<- which(all_bongo_stns$BottomDepth<1000 & all_bongo_stns$Cruise=="HB1603" 
          & all_bongo_stns$Month!="AUG")
length(I) # We exclude 23 stations from June and July with the 1000 m requirement

# so let's make an object that has all the station information, with and without
# positive bluefin catches, for June 17 through end of summer. 
I<- which(all_bongo_stns$Month=="MAY")
all_bongo_stns<- all_bongo_stns[-I,] # drop the May sampling from ECOMON
I<- which(all_bongo_stns$Month=="JUN" & all_bongo_stns$Day<17)
all_bongo_stns<- all_bongo_stns[-I,] # drop the early June sampling from ECOMON
all_bongo_stns<- merge(all_bongo_stns, tunacounts_SS[,c("Cruise", "Station", "Operation", "Nbluefin","Vol_filtered", "Abundance", "Density")], 
                       by=c("Cruise", "Station","Operation"), all.x = T)
all_bongo_stns$Abundance[is.na(all_bongo_stns$Abundance)]<- 0
# some of the day column is a character and it seems to be messing things up.
all_bongo_stns$Day<- as.numeric(all_bongo_stns$Day)

