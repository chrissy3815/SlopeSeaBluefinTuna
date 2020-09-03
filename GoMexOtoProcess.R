## GoMex otoliths from Katie:

library(readxl)

# load the data from the first read
read1<- read.csv("data/GoMexOtolithReads1.csv")

# load the data from the second read
read2<- read.csv("data/GoMexOtolithReads2.csv")

# Check which fish have Read1 and Read2 agreeing within 1d
agecheck<- merge(read1[,c("Fish","Age")], read2[,c("Fish", "Age")], by="Fish")
agecheck$d <- abs(agecheck$Age.y-agecheck$Age.x)
agecheck$within1 <- 0
agecheck$within1[agecheck$d<2] <- 1
# print the ones that differ by more than 1 day:
agecheck[agecheck$within1==0, "Fish"]
# shuffle the order and print to a csv file:
shuffledorder<- sample(agecheck[agecheck$within1==0, "Fish"])
write.csv(shuffledorder, file='results/GoMexOto2016_read3order.csv')

# bring in the third read and check those:
# bring in the read3 data:
read3<- read_excel('data/GoMexOtolithReads3.xlsx', sheet = 1)
names(read3)[length(names(read3))]<- "Age3"
all3reads<- merge(agecheck, read3[,c("Fish", "Age3")])
all3reads$diff1<- abs(all3reads$Age.x-all3reads$Age3)
all3reads$diff2<- abs(all3reads$Age.y-all3reads$Age3)
all3reads
# This recovers all but the largest otolith into the dataset
#change the column name back:
names(read3)[length(names(read3))]<- "Age"

# If the first and second read are within 1, keep second:
# (We are assuming that skill increased with time spent reading)
keepreads<- agecheck$Fish[agecheck$within1==1]
GOM_oto_data<- read2[read2$Fish %in% keepreads, ]
# If the third read matches within 1 day to either the 
# first or the second read, keep the third:
keepreads<- which(all3reads$diff1<=1 | all3reads$diff2<=1)
GOM_oto_data<- rbind(GOM_oto_data, read3[keepreads,])

# correct the Age to be Increments:
GOM_oto_data$Increments<- GOM_oto_data$Age-1
# if increments were already 0, should not be -1:
GOM_oto_data$Increments[GOM_oto_data$Increments<0]<- 0

# read in size data
seamapdata<- read_excel('data/BFT2016_SEAMAP_Metadata_forCHernandez_27Mar2018.xlsx', 1)
seamapdata<- seamapdata[,c('ELH_ID', 'SL_mm_EtOH', 'Slat', 'Slon')]
head(seamapdata)

# make matching columns in GOM_oto_data and seamapdata that contain the 4-digit fish ID:
GOM_oto_data$Fish<- as.character(GOM_oto_data$Fish)
sampleid<- strsplit(GOM_oto_data$Fish, '_')
Cruise<- sapply(sampleid, '[', 1)
Fish<- sapply(sampleid, '[', 2)
GOM_oto_data$Fish<- sapply(Fish, substr, start=1, stop=4)

seamapdata$Fish<- sapply(seamapdata$ELH_ID, substr, start=6, stop=9)

# append size and lat/lon to the otolith reads:
GOM_oto_data<- merge(GOM_oto_data, seamapdata, by='Fish')
GOM_oto_data<- GOM_oto_data[,-2] # get rid of the 'n' column
head(GOM_oto_data) # check that all the otolith data is still there.

# plot of size at radius
plot(GOM_oto_data$SL_mm_EtOH, GOM_oto_data$Radius, pch=19)

# calculate increment widths:
I<- which(names(GOM_oto_data)=='ToRing14') #bc the max age is 13
J<- which(names(GOM_oto_data)=='ToRing2') #bc I mark the edge of the core
# radii to outer and inner edges of each increment
outer<- GOM_oto_data[,(J+1):I]
inner<- GOM_oto_data[,(J):(I-1)]
# increment width
incwidthGOM<- outer-inner
incwidthGOM<- cbind(GOM_oto_data$Fish, incwidthGOM)
names(incwidthGOM)<- c("Fish", "Inc1", "Inc2", "Inc3", "Inc4", "Inc5", "Inc6", "Inc7",
                      "Inc8", "Inc9", "Inc10", "Inc11", "Inc12")

# Make a second increment width dataframe with only fish that have 8 increments or fewer
subdata<- GOM_oto_data[GOM_oto_data$Increments<9,]
# calculate increment widths:
I<- which(names(subdata)=='ToRing9') #bc the max age is 8
J<- which(names(subdata)=='ToRing2') #bc I mark the edge of the core
# radii to outer and inner edges of each increment
outer<- subdata[,(J+1):I]
inner<- subdata[,(J):(I-1)]
# increment width
incwidthGOM_sub8inc<- outer-inner
incwidthGOM_sub8inc<- cbind(subdata$Fish, incwidthGOM_sub8inc)
names(incwidthGOM_sub8inc)<- c("Fish", "Inc1", "Inc2", "Inc3", "Inc4", "Inc5", "Inc6", "Inc7")

rm(inner, outer, agecheck, all3reads, read1, read2, read3, sampleid, 
   subdata, Cruise, Fish, I, J, keepreads, shuffledorder)
