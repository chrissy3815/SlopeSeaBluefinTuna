## This is for making tables for my chapter/manuscript on bluefin tuna:
library(readxl)
library(here)

source(here("code",'SlopeSeaAbundProcess.R'))

# table to write to a csv:
table_SS<- tunacounts_SS[,c("Date", "Cruise", "Station", "Operation", 
                            "Latitude", "Longitude", "BottomDepth", "SST", 
                            "Nbluefin", "Abundance")]
write.csv(table_SS, file='results/SlopeSeaSampling_Table.csv')

## Make a quick plot of abundance vs. SST:
setEPS()
postscript('results/SS2016_SST_vs_Abundance.eps', height=5, width=6)
plot(table_SS$SST, table_SS$Abundance, pch=19,
     xlab="Sea Surface Temperature", ylab="Abundance")
dev.off()
