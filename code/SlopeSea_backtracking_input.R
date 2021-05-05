##### Data to give to Irina to run backtracking:

# first, need to run the Slope Sea otolith processing:
source('SlopeSeaOtoProcess.R')

###  1. For aged larvae, we add 2 days post hatch before first ring (Yufera et al 2014), and 2 days of egg duration 
#      (Reglero et al 2018) to get the "days post spawning."

# first, take the ones that have ages:
for_backtracking<- SS_oto_data[,c("Cruise", "Station", "Gear", "Fish", "Increments")]
head(for_backtracking)

## Add latitude and longitude of collection location, and date/time of collection
# First for the 2N3 samples:
I<- which(SS_oto_data$Gear=="2N3")
query_stations<- unique(SS_oto_data[I,c("Cruise", "Station")])
J<- which(metadata$GearType=="CTD/IKMT Oblique")
framenet_samples<- merge(query_stations, metadata[J, c("Cruise", "Station", "LatDec", "LonDec", "DateTime")])
framenet_samples$Gear<- "2N3"
# Next, for the baby Bongo samples:
I<- which(for_backtracking$Gear=="2B1")
query_stations<- unique(SS_oto_data[I,c("Cruise", "Station")])
J<- which(metadata$GearType=="CTD/Bongo Oblique")
babybongo_samples<- merge(query_stations, metadata[J, c("Cruise", "Station", "LatDec", "LonDec", "DateTime")])
babybongo_samples$Gear<- "2B1"
# Next, for the 6B3I samples:
I<- which(for_backtracking$Gear=="6B3I")
query_stations<- unique(SS_oto_data[I,c("Cruise", "Station")])
J<- which(metadata$GearType=="CTD/Bongo Oblique")
bongoI_samples<- merge(query_stations, metadata[J, c("Cruise", "Station", "LatDec", "LonDec", "DateTime")])
bongoI_samples$Gear<- "6B3I"
# Next, for the 6B3Z samples:
I<- which(for_backtracking$Gear=="6B3Z")
query_stations<- unique(SS_oto_data[I,c("Cruise", "Station")])
J<- which(metadata$GearType=="CTD/Bongo Oblique")
bongoZ_samples<- merge(query_stations, metadata[J, c("Cruise", "Station", "LatDec", "LonDec", "DateTime")])
bongoZ_samples$Gear<- "6B3Z"
# paste all these rows together
to_merge<- rbind(framenet_samples, babybongo_samples, bongoI_samples, bongoZ_samples)

# merge metadata to ages:
for_backtracking<- merge(for_backtracking, to_merge)

# Correct otolith "ages" to days post spawning: add 2 for days post hatch (Malca et al. 2017 references Yufera et al 2014) and add 2 days for egg duration to convert Increments to Days Post Spawning
for_backtracking$DaysPostSpawn<- for_backtracking$Increments+4

# There's no use in running repeats, so take unique rows of the important columns:
for_backtracking_unique<- unique(for_backtracking[,c("Cruise", "Station", "DateTime", "LatDec", "LonDec", "DaysPostSpawn")])
# but I want to keep track of how many larvae each row corresponds to:
for (i in 1:length(for_backtracking_unique$Cruise)){
  x<- length(which(for_backtracking$Cruise==for_backtracking_unique$Cruise[i] &
                   for_backtracking$Station==for_backtracking_unique$Station[i] &
                   for_backtracking$DateTime==for_backtracking_unique$DateTime[i] &
                   for_backtracking$DaysPostSpawn==for_backtracking_unique$DaysPostSpawn[i]))
  for_backtracking_unique$Nlarvae[i]<- x
}

# write the csv out:
write.csv(for_backtracking_unique, 
          file='results/SlopeSea2016_agedlarvae_forbacktracking_121520.csv', 
          row.names=F)

### 2. For larvae that we only have length, need to follow the same steps as from the PIPA project.
#     a. Construct age-length relationship from my Slope Sea reads
#     b. Use the residuals to define a normal distribution spreading around the best-fit line for age-length.
#     c. Sample from this distribution to define the age for each larva that we have a measured length.
#     d. Add 2 days post hatch before first ring (Yufera et al 2014), and 2 days of egg duration 
#        (Reglero et al 2018) to get the "days post spawning."

# There are also 3 larvae without length data. For these 3, we will use the 
# average length of the other larvae collected at that station. 

## Estimating ages for the rest of the larvae that have lengths:
# collect all the lengths that don't have corresponding ages:
no_age<- merge(all_lengths_SS, SS_oto_data[,c("Cruise", "Station", "Gear", "Fish", "Age")], 
               by=c("Cruise", "Station", "Gear", "Fish"), all.x=T)
I<- which(!is.na(no_age$Age))
no_age<- no_age[-I,c("Cruise","Station", "Gear", "Fish", "Length")]
head(no_age)

# Assign the average length of larvae from the station to the larvae without lengths.
avg_length_station<- aggregate(Length~Cruise+Station, data=all_lengths_SS, 
                               FUN=mean, na.rm=T)
I<- which(is.na(no_age$Length))
for (i in I){
  J<- which(avg_length_station$Station==no_age$Station[i])
  no_age$Length[i]<- avg_length_station$Length[J]
}

# Build the relationship with the residuals of the inverted age-length relationship:
model2<- lm(Increments~Length, data=SS_oto_data)
summary(model2)
# check visually if variance changes with time:
plot(SS_oto_data$Length, model2$residuals, pch=19)
lines(c(-1, 13), c(0,0))
# check visually if residuals are normally distributed:
hist(model2$residuals, breaks = c(-3.5, -2.5, -1.5, -.5, .5, 1.5, 2.5, 3.5, 4.5))
# mean of residuals:
mean(model2$residuals) ## should always be very close to 0!


# For all fish, need to estimate an age based on the length
no_age$EstDaysPostSpawn<-NA
# pull out the slope and intercept for the model that predicts age from length
slope_agelength<- model2$coefficients[2]
int_agelength<- model2$coefficients[1]
# standard deviation of residuals of length:
sd_agelength<- sd(model2$residuals)
# go through the list of lengths:
for (i in 1:length(no_age$Cruise)){
  ilength<- no_age$Length[i]
  # estimate age on regression line
  # also use the SD for the genus to pull 1 random number for a normal distribution
  agest<- ilength*slope_agelength+int_agelength
  agerand<- rnorm(1, agest, sd_agelength)
  
  # Round the randomly generated age to the nearest 1 day, and add 4 days to make it post spawn
  if (agerand<0){
    agerand<-0
  }
  no_age$EstDaysPostSpawn[i]<- round(agerand)+4
}

## Add latitude and longitude of collection location, and date/time of collection
# First for the 2N3 samples:
I<- which(no_age$Gear=="2N3")
query_stations<- unique(no_age[I,c("Cruise", "Station")])
J<- which(metadata$GearType=="CTD/IKMT Oblique")
framenet_samples<- merge(query_stations, metadata[J, c("Cruise", "Station", "LatDec", "LonDec", "DateTime")])
framenet_samples$Gear<- "2N3"
# Next, for the baby Bongo samples:
I<- which(no_age$Gear=="2B1")
query_stations<- unique(no_age[I,c("Cruise", "Station")])
J<- which(metadata$GearType=="CTD/Bongo Oblique")
babybongo_samples<- merge(query_stations, metadata[J, c("Cruise", "Station", "LatDec", "LonDec", "DateTime")])
babybongo_samples$Gear<- "2B1"
# Next, for the 6B3I samples:
I<- which(no_age$Gear=="6B3I")
query_stations<- unique(no_age[I,c("Cruise", "Station")])
J<- which(metadata$GearType=="CTD/Bongo Oblique")
bongoI_samples<- merge(query_stations, metadata[J, c("Cruise", "Station", "LatDec", "LonDec", "DateTime")])
bongoI_samples$Gear<- "6B3I"
# Next, for the 6B3Z samples:
I<- which(no_age$Gear=="6B3Z")
query_stations<- unique(no_age[I,c("Cruise", "Station")])
J<- which(metadata$GearType=="CTD/Bongo Oblique")
bongoZ_samples<- merge(query_stations, metadata[J, c("Cruise", "Station", "LatDec", "LonDec", "DateTime")])
bongoZ_samples$Gear<- "6B3Z"
# Next, for the 6B3 samples in the poland dataset:
I<- which(no_age$Gear=="6B3")
query_stations<- unique(no_age[I,c("Cruise", "Station")])
J<- which(metadata$GearType=="CTD/Bongo Oblique")
bongoY_samples<- merge(query_stations, metadata[J, c("Cruise", "Station", "LatDec", "LonDec", "DateTime")])
bongoY_samples$Gear<- "6B3"
# paste all these rows together
to_merge<- rbind(framenet_samples, babybongo_samples, bongoI_samples, 
                 bongoZ_samples, bongoY_samples)

no_age<- merge(no_age, to_merge)

# There's no use in running repeats, so take unique rows of the important columns:
no_age_unique<- unique(no_age[,c("Cruise", "Station", "DateTime", "LatDec", "LonDec", "EstDaysPostSpawn")])
# but I want to keep track of how many larvae each row corresponds to:
for (i in 1:length(no_age_unique$Cruise)){
  x<- length(which(no_age$Cruise==no_age_unique$Cruise[i] &
                   no_age$Station==no_age_unique$Station[i] &
                   no_age$DateTime==no_age_unique$DateTime[i] &
                   no_age$EstDaysPostSpawn==no_age_unique$EstDaysPostSpawn[i]))
  no_age_unique$Nlarvae[i]<- x
}

# write the csv out:
write.csv(no_age_unique, 
          file='results/SlopeSea2016_estimatedlarvae_forbacktracking_121520.csv',
          row.names=F)

### As of December 7, 2020, I think it would make the most sense to combine these
# two dataframes and only send a single one to Irina. 

# rename the column in for_backtracking_unique to "EstDaysPostSpawn"
J<- which(names(for_backtracking_unique)=="DaysPostSpawn")
names(for_backtracking_unique)[J]<- "EstDaysPostSpawn"
# rbind the two dataframes together:
all_backtracking<- rbind(no_age_unique, for_backtracking_unique)
# take the unique again:
all_backtracking_unique<- unique(all_backtracking[,c("Cruise", "Station", "DateTime", "LatDec", "LonDec", "EstDaysPostSpawn")])
# but I want to keep track of how many larvae each row corresponds to:
for (i in 1:length(all_backtracking_unique$Cruise)){
  x<- which(all_backtracking$Cruise==all_backtracking_unique$Cruise[i] &
            all_backtracking$Station==all_backtracking_unique$Station[i] &
            all_backtracking$DateTime==all_backtracking_unique$DateTime[i] &
            all_backtracking$EstDaysPostSpawn==all_backtracking_unique$EstDaysPostSpawn[i])
  all_backtracking_unique$Nlarvae[i]<- sum(all_backtracking$Nlarvae[x])
}

# write the csv!
write.csv(all_backtracking_unique, 
          file='results/SlopeSea2016_allLarvae_forbacktracking_121520.csv',
          row.names=F)