### Chunks of rejected code:


png(filename='results/SS_GOM_2016_agelength.png', 
    height=6.5, width=7.5, units= 'in', res=300, pointsize=12)
plot(GOM_oto_data$Increments, GOM_oto_data$SL_mm_EtOH, pch=19, col='grey', 
     xlab='Daily Increments',  ylab='Length (mm)', main='Slope Sea', 
     xlim=c(0,14), ylim=c(1, 9), cex=1.25, cex.lab=1.5, cex.axis=1.5)
exes<- 0:15
whys<- exes*summary(GOM_agelength)$coefficients[2,1]+summary(GOM_agelength)$coefficients[1,1]
lines(exes, whys, col='grey25')
text(1, 7.5, "SL=2.85+0.37*DI", pos=4, col='grey25')
# line for the sub8inc data:
exes<- 0:8
whys<- summary(GOM_agelength_sub8inc)$coefficients[1,1]+summary(GOM_agelength_sub8inc)$coefficients[2,1]*exes
lines(exes, whys, col='grey25')
text(1, 6.75, "SL=2.47+0.46*DI", pos=4, cex=1.2, col='grey25')
# add SS data:
points(SS_oto_data$Increments, SS_oto_data$Length, pch=1)
whys<- exes*summary(SS_agelength)$coefficients[2,1]+summary(SS_agelength)$coefficients[1,1]
# add the line for the linear fit
lines(exes, whys)
# add text for the best-fit line for Slope Sea:
text(8, 3, 'SL = 3.07 + 0.37*DI')
dev.off()

png(filename='results/SS_GOM_2016_RadiusAtAge.png', 
    height=6.5, width=7.5, units= 'in', res=300, pointsize=12)
plot(GOM_oto_data$Increments, GOM_oto_data$Radius, pch=19, col='grey', 
     xlab='Daily Increments',  ylab='Radius (um)', 
     xlim=c(0,14), ylim=c(10,55), cex.lab=1.5, cex.axis=1.5)
exes<- 0:15
whys<- exes*summary(GOM_radatage)$coefficients[2,1]+summary(GOM_radatage)$coefficients[1,1]
lines(exes, whys, col='grey25')
text(1, 55, "R=9.44+3.16*DI", pos=4, col='grey25')
# line for the sub8inc data:
exes<- 0:8
whys<- summary(GOM_radatage_sub8inc)$coefficients[1,1]+summary(GOM_radatage_sub8inc)$coefficients[2,1]*exes
lines(exes, whys, col='grey25')
text(1, 50, "Radius=9.81+3.05*DI", pos=4, cex=1.2, col='grey25')
# add SS data:
points(SS_oto_data$Increments, SS_oto_data$Radius, pch=1)
whys<- exes*summary(SS_radatage)$coefficients[2,1]+summary(SS_radatage)$coefficients[1,1]
# add the line for the linear fit
lines(exes, whys)
# add text for the best-fit line for Slope Sea:
text(1, 35, "Radius=10.88+2.64*DI", pos=4)
dev.off()