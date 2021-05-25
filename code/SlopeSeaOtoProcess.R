# Slope Sea otolith processing:

# required packages:
library(readxl)
library(here)

## First, we need to shuffle the read order for read 2:
read1<- read_excel(here('data', 'SlopeSeaOto2016_FullRead_1.xlsx'))
filenames<- read1$Fish
shuffledorder<- sample(filenames)
write.csv(shuffledorder, file='results/SlopeSeaOto_2020Reads_read2order.csv')

## Next, some QC checks on reads 1 and 2:
# bring in the data from otoliths read 1:
read1<- read_excel(here('data','SlopeSeaOto2016_FullRead_1.xlsx'), sheet = 1)
# bring in the data from otoliths read 2:
read2<- read_excel(here('data','SlopeSeaOto2016_FullRead_2_retakes.xlsx'), sheet = 1)
names(read2)[length(names(read2))]<- "Age2"

# use merge to combine both reads:
bothreads<- merge(read1[,c('Fish', 'Age')], read2[,c('Fish', 'Age2')], by='Fish', all.y=T)
# difference between the two reads:
bothreads$diff = abs(bothreads$Age-bothreads$Age2)
summary(bothreads$diff)
# print the ones that differ by more than 1 day:
bothreads[bothreads$diff>1,c("Fish", "diff")]

## Shuffle the order for read 3:
shuffledorder2<- sample(bothreads$Fish[bothreads$diff>1])
write.csv(shuffledorder2, file='results/SlopeSeaOto_2020Reads_read3order.csv')

##Check read 3:
# bring in the read3 data:
read3<- read_excel(here('data','SlopeSeaOto2016_FullRead_3.xlsx'), sheet = 1)
names(read3)[length(names(read3))]<- "Age3"
all3reads<- merge(bothreads, read3[,c("Fish", "Age3")])
all3reads$diff1<- abs(all3reads$Age-all3reads$Age3)
all3reads$diff2<- abs(all3reads$Age2-all3reads$Age3)
all3reads
# We'll lose 2 of the 56 fish from the analysis because their otoliths didn't pass the test on read 3. 

## Okay, now we can treat this data as final and clean it up for use in figures, etc:
# If Read 1 and Read 2 agree within 1 day, use Read 2
read2_tokeep<- bothreads[bothreads$diff<=1,"Fish"]
SS_oto_data<- read2[read2$Fish %in% read2_tokeep,]
names(SS_oto_data)[length(names(SS_oto_data))]<- "Age"
# If Read 3 agrees with either Read 1 or Read 2 to within 1 day, use Read 3:
read3_tokeep<- all3reads[all3reads$diff1<=1 | all3reads$diff2<=1,"Fish"]
names(read3)[length(names(read3))]<- "Age"
SS_oto_data<- rbind(SS_oto_data, read3[read3$Fish %in% read3_tokeep,])

# correct the ages (I mark the edge of the core, which is not a daily increment)
SS_oto_data$Increments<- SS_oto_data$Age-1
# there are a couple that I didn't even mark the edge of the core- these should still be recorded as 0, not -1
SS_oto_data$Increments[SS_oto_data$Increments<0]<- 0

# Let's also get rid of the empty columns:
NAcolumns<- apply(SS_oto_data, MARGIN=2, function(x){sum(is.na(x))})
I<- which(NAcolumns==54)
SS_oto_data<- SS_oto_data[,-I]
# rename the "Fish" column as "Image"
I<- which(names(SS_oto_data)=="Fish")
names(SS_oto_data)[I]<- "Image"
# Also, drop the "n" column which is a carryover from ImageJ and doesn't mean anything
I<- which(names(SS_oto_data)=="n")
if(length(I)>0){
  SS_oto_data<- SS_oto_data[,-I]
}

# break up the "Image" column to get cruise, station, gear, fish:
fishID<- strsplit(SS_oto_data$Image, "-")
# get out the cruise IDs:
cruise<- sapply(fishID, '[', 1)
head(cruise)
SS_oto_data$Cruise<- sapply(strsplit(cruise, "Stn"), '[', 1)
# get out the station numbers:
station<- sapply(strsplit(SS_oto_data$Image, "Stn"), '[', 2)
SS_oto_data$Station<- as.numeric(sapply(strsplit(station, "-"), '[', 1))
# get the gear IDs:
SS_oto_data$Gear<- sapply(strsplit(station, "-"), '[', 2)
# get the fish numbers (actually, need these as numeric to match length data):
fishNum<- sapply(strsplit(station, "-"), '[', 3) # extracts the FX portion
fishNum<- sapply(strsplit(fishNum, 'F'), '[', 2) # cleaves off the "F"
fishNum<- sapply(strsplit(fishNum, '.tif'), '[', 1) # cleaves off ".tif" on any that have that
SS_oto_data$Fish<- as.numeric(fishNum)
head(SS_oto_data)

## Other metadata from the tows:
slopeseaoperations<- read_excel(here('data','SlopeSeaOperations.xlsx'), sheet = 2)
metadata<- slopeseaoperations[,c("cruiseid", "siteid", "event time", 
                                 "deployment", "lat", "lon", "max ctd depth", 
                                 "gm_1", "tot_1", "gm_2", "tot_2", "sfc t")]
names(metadata)<- c("Cruise", "Station", "DateTime", "GearType", "lat", "lon", 
                    "MaxDepth", "Bongo1", "VolumeFiltered_B1", "Bongo2", 
                    "VolumeFiltered_B2", "SST")
metadata$LatDec<- floor(metadata$lat/100)+(metadata$lat-floor(metadata$lat/100)*100)/60
metadata$LonDec<- floor(metadata$lon/100)+(metadata$lon-floor(metadata$lon/100)*100)/60
head(metadata)

## Wrangling the length data:
davelengths<- read.csv(here('data','2016MeasuredFish_CMH200811.csv'))
chrissylengths<- read_excel(here('data','HB1603_6B3I_BFTlengths_20200811.xlsx'), sheet = 1)
names(chrissylengths)[length(names(chrissylengths))]<- "Length"
polanddata<- read.csv(here('data','HB1603_GU1608_IchData_7Nov2019.csv'))

# pull out the correct columns to be able to rbind davelengths and chrissylengths
davelengths<- davelengths[,c('Cruise', 'Station', 'Gear', 'Fish', 'Length')]
chrissylengths<- chrissylengths[,c('Cruise', 'Station', 'Gear', 'Fish', 'Length')]
usalengths<- rbind(davelengths, chrissylengths)
rm(davelengths, chrissylengths)
# fix the typo in geartype:
usalengths$Gear[usalengths$Gear=='2n3']<- "2N3"

# add lat, lon, and date to usalengths
I<- which(usalengths$Gear=="2N3")
query_stations<- unique(usalengths[I,c("Cruise", "Station")])
J<- which(metadata$GearType=="CTD/IKMT Oblique")
framenet_samples<- merge(query_stations, metadata[J, c("Cruise", "Station", "LatDec", "LonDec", "DateTime", "SST")])
framenet_samples$Gear<- "2N3"
# Next, for the baby Bongo samples:
I<- which(usalengths$Gear=="2B1")
query_stations<- unique(usalengths[I,c("Cruise", "Station")])
J<- which(metadata$GearType=="CTD/Bongo Oblique")
babybongo_samples<- merge(query_stations, metadata[J, c("Cruise", "Station", "LatDec", "LonDec", "DateTime", "SST")])
babybongo_samples$Gear<- "2B1"
# Next, for the 6B3I samples:
I<- which(usalengths$Gear=="6B3I")
query_stations<- unique(usalengths[I,c("Cruise", "Station")])
J<- which(metadata$GearType=="CTD/Bongo Oblique")
bongoI_samples<- merge(query_stations, metadata[J, c("Cruise", "Station", "LatDec", "LonDec", "DateTime", "SST")])
bongoI_samples$Gear<- "6B3I"
# Next, for the 6B3Z samples:
I<- which(usalengths$Gear=="6B3Z")
query_stations<- unique(usalengths[I,c("Cruise", "Station")])
J<- which(metadata$GearType=="CTD/Bongo Oblique")
bongoZ_samples<- merge(query_stations, metadata[J, c("Cruise", "Station", "LatDec", "LonDec", "DateTime", "SST")])
bongoZ_samples$Gear<- "6B3Z"
# paste all these rows together
to_merge<- rbind(framenet_samples, babybongo_samples, bongoI_samples, bongoZ_samples)
# merge metadata to usalengths:
usalengths<- merge(usalengths, to_merge)
head(usalengths)

# Organize the Poland length data to match:
polandlengths<- polanddata[polanddata$TAXA_NAME=="Thunnus thynnus",]
polandlengths<- polandlengths[,c("CRUISE_NAME", "STATION", "GEAR", "COUNT_AT_LENGTH", "LENGTH")]
polandlengths_long <- as.data.frame(lapply(polandlengths, rep, polandlengths$COUNT_AT_LENGTH))
I<- which(names(polandlengths_long)=="COUNT_AT_LENGTH")
polandlengths_long<- polandlengths_long[,-I]
names(polandlengths_long)<- c("Cruise", "Station", "Gear", "Length")
polandlengths_long$Fish<- NA
# Find the 3 fish with missing length data:
polanddata_agg<- polanddata[polanddata$TAXA_NAME=="Thunnus thynnus",]
polanddata_agg<- polanddata_agg[,c("CRUISE_NAME", "STATION", "GEAR", "TOTAL_COUNT",
                                 "COUNT_AT_LENGTH", "LENGTH")]
polanddata_agg<- aggregate(COUNT_AT_LENGTH~CRUISE_NAME+STATION+TOTAL_COUNT, 
                           data=polanddata_agg, FUN = sum)
polanddata_agg$REP<- polanddata_agg$TOTAL_COUNT-polanddata_agg$COUNT_AT_LENGTH
I<- which(polanddata_agg$REP>0)
for (i in I){
  J<- which(polandlengths_long$Cruise==polanddata_agg$CRUISE_NAME[i] & 
              polandlengths_long$Station==polanddata_agg$STATION[i])[1]
  newrow<- as.data.frame(lapply(polandlengths_long[J,], rep, polanddata_agg$REP[i]))
  newrow$Length<- NA
  polandlengths_long<- rbind(polandlengths_long, newrow)
}

# Bring in the lat, lon, and datetime fields from the metadata
query_stations<- unique(polandlengths_long[,c("Cruise", "Station")])
J<- which(metadata$GearType=="CTD/Bongo Oblique")
bongoI_samples<- merge(query_stations, metadata[J, c("Cruise", "Station", "LatDec", "LonDec", "DateTime", "SST")])
# merge back together:
polandlengths_long<- merge(polandlengths_long,bongoI_samples)
# reorder the columns to match usalengths so that i can rbind:
polandlengths_long<- polandlengths_long[,names(usalengths)]

# need to read in and add the scomber scombrus corrections from Dave:
scomber<- read_excel(here('data','FishData_HB1603_ScomberCheck.xlsx'), sheet=2)
scomber<- scomber[scomber$`Richardson ID`=="bluefin",]
scomber<- scomber[,c("CRUISE_NAME", "STATION", "GEAR", "COUNT_AT_LENGTH", "LENGTH")]
scomber_long <- as.data.frame(lapply(scomber, rep, scomber$COUNT_AT_LENGTH))
I<- which(names(scomber_long)=="COUNT_AT_LENGTH")
scomber_long<- scomber_long[,-I]
names(scomber_long)<- c("Cruise", "Station", "Gear", "Length")
scomber_long$Fish<- NA
# Bring in the lat, lon, and datetime fields from the metadata
query_stations<- unique(scomber_long[,c("Cruise", "Station")])
J<- which(metadata$GearType=="CTD/Bongo Oblique")
bongoI_samples<- merge(query_stations, metadata[J, c("Cruise", "Station", "LatDec", "LonDec", "DateTime", "SST")])
# merge back together:
scomber_long<- merge(scomber_long,bongoI_samples)
# reorder the columns to match usalengths so that i can rbind:
scomber_long<- scomber_long[,names(usalengths)]

all_lengths_SS<- rbind(usalengths, polandlengths_long, scomber_long)

# Join the otolith data and length data:
SS_oto_data<- merge(SS_oto_data, all_lengths_SS)
dim(SS_oto_data)

# calculate increment widths for Slope Sea:
I<- which(names(SS_oto_data)=='ToRing9') #bc the max age is 8
J<- which(names(SS_oto_data)=='ToRing2') #bc I mark the edge of the core
# radii to outer and inner edges of each increment
outer<- SS_oto_data[,(J+1):I]
inner<- SS_oto_data[,(J):(I-1)]
# increment width
incwidthSS<- outer-inner
incwidthSS<- cbind(SS_oto_data[,c("Cruise", "Station", "Gear","Fish")], incwidthSS)
names(incwidthSS)<- c("Cruise", "Station", "Gear", "Fish", "Inc1", "Inc2", "Inc3", "Inc4", "Inc5", "Inc6", "Inc7")

### Check which of the possible lapillae need to be excluded from the data:
to_exclude<- read_xlsx(here('data','SlopeSea_PossibleLapillae.xlsx'))
plot(SS_oto_data$Increments, SS_oto_data$Radius, pch=19, xlab='Increments', ylab='Radius')
for (i in 1:length(to_exclude$Cruise)){
  I<- which(SS_oto_data$Station==to_exclude$Station[i] &
              SS_oto_data$Fish==to_exclude$Fish[i])
  points(SS_oto_data$Increments[I], SS_oto_data$Radius[I], pch=19, col='red')
}
## From this plot, we determined that the only ones that are likely to be lapillae
## and that, therefore, should be excluded, are the ones with 2 increments because
## they fall at the bottom of the distribution of otolith radius for otoliths with 
## 2 increments.

# So, find the ones that have 2 increments:
to_exclude<- merge(to_exclude, SS_oto_data[,c("Cruise", "Station", "Gear", "Fish", "Increments")])
to_exclude<- to_exclude[to_exclude$Increments==2,]
# And now exclude them from the data:
for (i in 1:length(to_exclude$Cruise)){
  I<- which(SS_oto_data$Station==to_exclude$Station[i] &
              SS_oto_data$Fish==to_exclude$Fish[i])
  SS_oto_data<- SS_oto_data[-I,]
  I<- which(incwidthSS$Station==to_exclude$Station[i] &
              incwidthSS$Fish==to_exclude$Fish[i])
  incwidthSS<- incwidthSS[-I,]
}

# Clean up the workspace:
rm(inner, outer, to_exclude)
rm(babybongo_samples, query_stations, bongoZ_samples, bongoI_samples, framenet_samples, to_merge)
rm(polandlengths, polanddata, usalengths, polanddata_agg)
rm(all3reads, bothreads, fishID, read1, read2, read3)
rm(read2_tokeep, read3_tokeep, cruise, filenames, fishNum, I, J, NAcolumns, shuffledorder, shuffledorder2, station)
rm(slopeseaoperations, polandlengths_long, scomber, scomber_long, newrow, metadata)

