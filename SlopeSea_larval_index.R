
# run the Slope Sea processing script:
source('SlopeSeaOtoProcess.R')

# we'll use only the bongo samples:
# drop the 2B1 samples:
I<- which(all_lengths_SS$Gear=="2B1")
all_lengths_SS<- all_lengths_SS[-I,]
# drop the 2N3 samples as well:
I<- which(all_lengths_SS$Gear=="2N3")
all_lengths_SS<- all_lengths_SS[-I,]
# add the "Operation" column:
all_lengths_SS$Operation<- NA
I<- which(all_lengths_SS$Gear=="6B3I" | all_lengths_SS$Gear=="6B3" | all_lengths_SS$Gear=="6B3Z")
all_lengths_SS$Operation[I]<- "BON/CTD"

# need to add the volume filtered:
SS2016_netdata<- read.csv('data/GU1608HB1603Net.csv')
# pull out the relevant columns:
SS2016_netdata<- SS2016_netdata[,c("CRUISE_NAME", "STATION", "GEAR", "GEAR_VOLUME_FILTERED")]
names(SS2016_netdata)<- c("Cruise", "Station", "Gear", "Vol_filtered")
# add the Operation column:
SS2016_netdata$Operation<- NA
I<- which(SS2016_netdata$Gear=="6B3I" | SS2016_netdata$Gear=="6B3" | SS2016_netdata$Gear=="6B3Z")
SS2016_netdata$Operation[I]<- "BON/CTD"
# exclude the volume filtered for a few samples that weren't sorted:
not_sorted<- data.frame(Cruise=c(rep("HB1603",2),"GU1608"), 
                        Station=c(36,125,240), 
                        Gear = c(rep("6B3I", 2), "6B3Z"))
to_exclude<- vector()
for (i in 1:length(not_sorted$Cruise)){
  J<- which(SS2016_netdata$Cruise==not_sorted$Cruise[i] &
              SS2016_netdata$Station==not_sorted$Station[i] &
              SS2016_netdata$Gear==not_sorted$Gear[i])
  to_exclude<- c(to_exclude,J)
}
SS2016_netdata<- SS2016_netdata[-to_exclude,]
# we're going to combine the Bongo samples:
SS2016_netdata<- aggregate(Vol_filtered~Cruise+Station+Operation, 
                           data=SS2016_netdata, FUN=sum, na.rm=T)
# merge with the all_lengths_SS:
all_lengths_SS<- merge(all_lengths_SS, SS2016_netdata, all.x=T, all.y=F)

# need to add the sampling depth and lat/lon:
SS2016_eventdata<- read.csv('data/GU1608HB1603Event.csv')
# pull out the relevant columns:
SS2016_eventdata<- SS2016_eventdata[,c("CRUISE_NAME", "STATION", "OPERATION", 
                                       "TOW_MAXIMUM_DEPTH", "BOTTOM_DEPTH_MAX_WIRE_OUT",
                                       "LATITUDE", "LONGITUDE", "EVENT_DATE")]
names(SS2016_eventdata)<- c("Cruise", "Station", "Operation", "SamplingDepth",
                            "BottomDepth", "Latitude", "Longitude", "Date")
# we want only the bongo events:
SS2016_eventdata<- SS2016_eventdata[SS2016_eventdata$Operation=="BON/CTD",]
# get day and month out:
SSmoday<- strsplit(SS2016_eventdata$Date, "-")
SS2016_eventdata$Day<- sapply(SSmoday, '[', 1)
SS2016_eventdata$Month<- sapply(SSmoday, '[', 2)
# there are two stations that weren't processed in this date range, 
# so they should be dropped from the analysis:
not_sorted<- data.frame(Cruise=c("HB1603", "GU1608"), Station=c(21, 232))
for (i in 1:length(not_sorted$Cruise)){
  J<- which(SS2016_eventdata$Cruise==not_sorted$Cruise[i] &
            SS2016_eventdata$Station==not_sorted$Station[i])
  if (length(J)>0){
    SS2016_eventdata<- SS2016_eventdata[-J,]
  }
}

#count how many events fit our criteria of June 15-August 15, in water 1000 m or deeper
june<- which(SS2016_eventdata$Month=="JUN" & 
             SS2016_eventdata$BottomDepth>=1000 &
             SS2016_eventdata$Day>=15)
july<- which(SS2016_eventdata$Month=="JUL" & 
               SS2016_eventdata$BottomDepth>=1000)
aug<- which(SS2016_eventdata$Month=="AUG" & 
              SS2016_eventdata$BottomDepth>=1000 &
              SS2016_eventdata$Day<=15)
# ny is the number of total stations sampled:
n2016<- length(june) + length(july) + length(aug)

# merge with the tuna data:
all_lengths_SS<- merge(all_lengths_SS, SS2016_eventdata, all.x=T, all.y=F)

# drop station 68 because it's in shallow water:
I<- which(all_lengths_SS$BottomDepth<1000)
all_lengths_SS<- all_lengths_SS[-I,]

## I'll just use the deterministic age-length relationship and round to nearest increment
SS_agelength_inverse<- lm(Increments~Length, data=SS_oto_data)
all_lengths_SS$DI<- SS_agelength_inverse$coefficients[1] +
                    all_lengths_SS$Length*SS_agelength_inverse$coefficients[2]
all_lengths_SS$DI<- round(all_lengths_SS$DI)
# I<- which(all_lengths_SS$DI<=0)
#all_lengths_SS<- all_lengths_SS[-I,]
#all_lengths_SS$DI[all_lengths_SS$DI<0]<- 0
all_lengths_SS$DI[all_lengths_SS$DI<0]<- 1

unique_stations<- unique(all_lengths_SS[,c("Cruise","Station")])

# my is the number of positive stations:
m2016<- length(unique_stations$Cruise)

Isy<- rep(0,m2016)

for (i in 1:m2016){
  station_i<- unique_stations$Station[i]
  subdata<- all_lengths_SS[all_lengths_SS$Station==station_i,]
  subdata$num<-exp(-0.2*(subdata$DI-1))
  Asy<- subdata$Vol_filtered[1]/subdata$SamplingDepth[1]
  Isy[i]<- sum(subdata$num/Asy)
}

T2016<- mean(Isy, na.rm=T)
s20162<- var(Isy, na.rm = T)

check<- 1
Gterm<- 1+(m2016-1)/m2016*s20162/2
denom<- m2016+1
j<- 2
while (check>0.001){
  if (j>2){
    denom<- denom*(m2016+2*j-3)
  }
  Gterm_prev<- Gterm
  Gterm<- Gterm_prev + (m2016-1)^(2*j-1)/(m2016^j*denom) * (s20162/2)^j/factorial(j)
  j<- j+1
  check<- abs(Gterm-Gterm_prev)
}
Gterm

Iy<- m2016/n2016*exp(T2016)*Gterm

## a couple little analyses of the GoMex index:
gom_index<- read.csv('data/Ingram2018_GoMindexTS.csv')
gom_index_1990_2015<- gom_index[gom_index$Survey.Year>=1990 & gom_index$Survey.Year<2015,]
min(gom_index_1990_2015$ZIDL.Index, na.rm=T)
max(gom_index_1990_2015$ZIDL.Index, na.rm=T)
mean(gom_index_1990_2015$ZIDL.Index, na.rm=T)

gom_index_1981_2015<- gom_index[gom_index$Survey.Year>=1981 & gom_index$Survey.Year<2015,]
min(gom_index_1981_2015$ZIDL.Index, na.rm=T)
max(gom_index_1981_2015$ZIDL.Index, na.rm=T)
mean(gom_index_1981_2015$ZIDL.Index, na.rm=T)
