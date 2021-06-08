## SEAMAP data processing:

library(here)

# read in the data that Glenn Zapfe sent me for 2016:
seamap_bongos<- read.csv(here('data','SEAMAP_SPRING2016_TUNA_LARVAE_bongos_ENVIRO.csv'))
# seamap_shallow_bongos<- read.csv(here('data','SEAMAP_SPRING2016_TUNA_LARVAE_shallow_bongos_ENVIRO.csv'))
seamap_neuston<- read.csv(here('data','SEAMAP_SPRING2016_TUNA_LARVAE_neuston_ENVIRO.csv'))
# I'm not going to use the "Shallow Bongos" dataset because they seem to go to similar depths but the mesh is 505
# The SEAMAP report for 2016 does not mention anything about shallow bongos or tows with 505 mesh.
# Link to SEAMAP report: gsmfc.org/publications/GSMFC%20Number%20261.pdf

# generate a list of the unique stations:
stations<- unique(seamap_bongos$P_STA_NO) # there are 119, as stated in the report.

# initialize a dataframe for the bluefin stations:
seamap_bluefin<- data.frame()

# run a loop to get all the positive and zero stations:
for (i in stations){
  # pull out the rows for this station
  subdata<- seamap_bongos[seamap_bongos$P_STA_NO==i,]
  # check if there are Thunnus thynnus
  j<- which(subdata$TAXON=="Thunnus thynnus")
  if (length(j)>0){
    newrow<- subdata[j,]
  } else if (length(j)==0){
    newrow<- subdata[1,]
  }
  seamap_bluefin<- rbind(seamap_bluefin, newrow)
}
# Calculate abundance as N per 10 sq. meters. DEPTH_EMAX is the maximum pressure sensor reading during the tow
seamap_bluefin$Abundance<- 10*seamap_bluefin$TOT_LARVAE/seamap_bluefin$VOL_FILT*seamap_bluefin$DEPTH_EMAX

