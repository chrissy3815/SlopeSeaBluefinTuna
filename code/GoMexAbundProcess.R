## SEAMAP data processing:

# read in the data that Glenn Zapfe sent me for 2016:
seamap_bongos<- read.csv(here('data','SPRING2016_TUNA_LARVAE_bongos_ENVIRO.csv'))
# seamap_shallow_bongos<- read.csv(here('data','SPRING2016_TUNA_LARVAE_shallow_bongos_ENVIRO.csv'))
# I'm not going to use the "Shallow Bongos" dataset because they seem to go to similar depths but the mesh is 505
# The SEAMAP report for 2016 does not mention anything about shallow bongos or tows with 505 mesh.
# Link to SEAMAP report: gsmfc.org/publications/GSMFC%20Number%20261.pdf

# generate a list of the unique stations:
stations<- unique(seamap_bongos$P_STA_NO)


seamap_zeros<- seamap_bongos[seamap_bongos$TAXON=="NO TUNA LARVAE CAUGH",]
seamap_bluefin<- seamap_bongos[seamap_bongos$TAXON=="Thunnus thynnus",]
seamap_bongos<- rbind(seamap_zeros, seamap_bluefin)
seamap_bongos$Abundance<- 10*seamap_bongos$TOT_LARVAE/seamap_bongos$VOL_FILT*seamap_bongos$DEPTH_EMAX
# remove the repeated rows:
I<- which(seamap_bongos$P_STA_NO==63)
seamap_bongos<- seamap_bongos[-I[1],]
I<- which(seamap_bongos$P_STA_NO==66)
seamap_bongos<- seamap_bongos[-I[1],]