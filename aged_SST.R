## A bit of code for calculating the average collection temperature of aged larvae

# run the Slope Sea processing script:
source('SlopeSeaOtoProcess.R')
# run the Gulf of Mexico processing script:
source('GoMexOtoProcess.R')

mean(SS_oto_data$SST)

# add temperatures to the GOM data:
seamap_stdbongos<- read.csv('data/SPRING2016_TUNA_LARVAE_bongos_ENVIRO.csv')
seamap_shlwbongos<- read.csv('data/SPRING2016_TUNA_LARVAE_shallow_bongos_ENVIRO.csv')
seamap_bongoSST<- unique(rbind(seamap_stdbongos[,c("SAMPLE_NO", "TEMPSURF")], 
                        seamap_shlwbongos[,c("SAMPLE_NO", "TEMPSURF")]))

GOM_oto_data$SST<- NA
for (i in 1:length(GOM_oto_data$Fish)){
  sampleid<- GOM_oto_data$Sample_No[i]
  J<- which(seamap_bongoSST$SAMPLE_NO==sampleid)
  if (length(J)>0){
    GOM_oto_data$SST<- seamap_bongoSST$TEMPSURF[J]
  }
}

mean(GOM_oto_data$SST)
