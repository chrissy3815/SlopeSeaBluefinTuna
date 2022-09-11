############################################################################
##
## Conceptual map figure for proposals and talks
##
############################################################################

# required packages:
library(oce)
library(ocedata)
library(ncdf4)
library(here)
library(marmap)
library(sp)

# Load in the coastlines
data("coastlineWorldFine", package="ocedata")
# load bathymetry from marmap
b <- as.topo(getNOAA.bathy(-40, -120, 15, 55, keep=TRUE))
# colormap for water:
blue.col <- colorRampPalette(c("steelblue", "lightseagreen"))

setEPS()
postscript(here("doc","ConceptualBaseMap.eps"), height=5.5, width=6)
# png(here("doc","ConceptualBaseMap.png"), height=5.5, width=5.5, units="in", res=300)
mapPlot(coastlineWorldFine, projection="+proj=aea +lat_1=30 +lat_2=50 +lon_0=-70",
        longitudelim = c(-62, -99),
        latitudelim = c(25, 50), col='grey', 
        grid=FALSE, axes = FALSE)
mapImage(b, col=oceColorsGebco, breaks=seq(-1600, 0, 200))
mapPolygon(coastlineWorldFine, col='light grey', lwd=0.5)
mapGrid(dlongitude = 10, dlatitude = 10, lwd=0.25)
mapAxis(cex.axis = 0.5)
dev.off()
