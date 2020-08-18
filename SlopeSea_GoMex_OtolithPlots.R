#### Run otolith processing steps and make figures

# run the Slope Sea processing script:
source('SlopeSeaOtoProcess.R')
# run the Gulf of Mexico processing script:
source('GoMexOtoProcess.R')

## Construct age-length relationships for the Slope Sea reads:
# make a linear model of age and length
SS_agelength<- lm(Length~Increments, data=SS_oto_data)
summary(SS_agelength)
# make a second model for only individuals with 0-4 increments
SS_agelength_sub4inc<- lm(Length~Increments, data=SS_oto_data[SS_oto_data$Increments<=4,])
summary(SS_agelength_sub4inc)

# build linear models for the GoMex:
GOM_agelength<- lm(SL_mm_EtOH~Increments, data=GOM_oto_data)
summary(GOM_agelength)
GOM_agelength_sub8inc<- lm(SL_mm_EtOH~Increments, data=GOM_oto_data[GOM_oto_data$Increments<9,])
summary(GOM_agelength_sub8inc)
GOM_agelength_sub4inc<- lm(SL_mm_EtOH~Increments, data=GOM_oto_data[GOM_oto_data$Increments<=4,])
summary(GOM_agelength_sub4inc)

## Slope Sea plots:
png(filename='results/SlopeSea2016_SizeAtAge2.png', height=6.5, width=7.5, 
    units= 'in', res=300)
plot(SS_oto_data$Increments, SS_oto_data$Length, xlab='Daily Increments', 
     ylab='Standard Length (mm)', pch=1, cex=1.25, cex.lab=1.5, cex.axis=1.5,
     xlim=c(0,13), ylim=c(2, 8))
exes<- 0:10
whys<- exes*summary(SS_agelength)$coefficients[2,1]+summary(SS_agelength)$coefficients[1,1]
# add the line for the overall SS linear fit
lines(exes, whys, lwd=1.5)
# add the line for the 0-8 inc GoMex fit:
whys<- exes*summary(GOM_agelength_sub8inc)$coefficients[2,1]+summary(GOM_agelength_sub8inc)$coefficients[1,1]
lines(exes, whys, col='grey', lwd=1.5)
# add the line for the 0-4 inc SS fit:
whys<- exes*summary(SS_agelength_sub4inc)$coefficients[2,1]+summary(SS_agelength_sub4inc)$coefficients[1,1]
lines(exes, whys, lty=2, lwd=1.5)
# add the line for the 0-4 inc GoMex fit:
whys<- exes*summary(GOM_agelength_sub4inc)$coefficients[2,1]+summary(GOM_agelength_sub4inc)$coefficients[1,1]
lines(exes, whys, lty=2, lwd=1.5, col='grey')
legend("topleft", legend=c("SS 0-8 inc", "SS 0-4 inc", "GOM 0-8 inc", "GOM 0-4 inc"), 
       lty=c(1,2,1,2), lwd=1.5, col=c("black", "black", "grey", "grey"))
dev.off()

# radius at age:
png(filename='results/SlopeSea2016_RadiusAtAge.png', height=6.5, width=7.5, 
    units= 'in', res=300)
plot(SS_oto_data$Increments, SS_oto_data$Radius, xlab='Daily Increments', 
     ylab='Otolith Radius (um)', pch=19)
SS_radatage<- lm(Radius~Increments, data=SS_oto_data)
summary(SS_radatage)
exes<- 0:14
whys<- summary(SS_radatage)$coefficients[1,1]+summary(SS_radatage)$coefficients[2,1]*exes
lines(exes, whys)
text(1, 35, "Radius=10.88+2.64*DI", pos=4, cex=1.5)
dev.off()


## Gulf of Mexico otolith figures:
# plot of size at age (to file)
png(filename='results/GoMex2016_SizeAtAge_2.png', height=6.5, width=7.5, 
    units= 'in', res=300)
plot(GOM_oto_data$Increments, GOM_oto_data$SL_mm_EtOH, xlab='Daily Increments', 
     ylab='Standard Length (mm)', pch=1, cex=1.25, cex.lab=1.5, cex.axis=1.5,
     xlim=c(0,13), ylim=c(2, 8))
exes<- 0:14
whys<- summary(GOM_agelength)$coefficients[1,1]+summary(GOM_agelength)$coefficients[2,1]*exes
lines(exes, whys, lwd=1.5, lty=3)
# text(1, 6.75, "SL=2.85+0.37*DI", pos=4)
# line for fish up to 8 increments, to match with the SS data:
whys2<- summary(GOM_agelength_sub8inc)$coefficients[1,1]+summary(GOM_agelength_sub8inc)$coefficients[2,1]*exes
lines(exes, whys2, lwd=1.5)
# text(1, 7.5, "SL=2.47+0.46*DI", pos=4, cex=1.5)
# add the line for the 0-4 inc GoMex fit:
whys<- exes*summary(GOM_agelength_sub4inc)$coefficients[2,1]+summary(GOM_agelength_sub4inc)$coefficients[1,1]
lines(exes, whys, lty=2, lwd=1.5)
# add the line for the overall SS linear fit
whys<- exes*summary(SS_agelength)$coefficients[2,1]+summary(SS_agelength)$coefficients[1,1]
lines(exes, whys, lwd=1.5, col='grey')
# add the line for the 0-4 inc SS fit:
whys<- exes*summary(SS_agelength_sub4inc)$coefficients[2,1]+summary(SS_agelength_sub4inc)$coefficients[1,1]
lines(exes, whys, lty=2, lwd=1.5, col='grey')
legend("topleft", legend=c("GOM 0-13 inc", "GOM 0-8 inc", "GOM 0-4 inc", "SS 0-8 inc", "SS 0-4 inc"), 
       lty=c(3,1,2,1,2), lwd=1.5, col=c("black", "black", "black", "grey", "grey"))
dev.off()

# radius at age (to get rid of shrinkage differences)
png(filename='results/GoM2016_RadiusAtAge.png', height=6.5, width=7.5, 
    units= 'in', res=300)
plot(GOM_oto_data$Increments, GOM_oto_data$Radius, xlab='Daily Increments', 
     ylab='Otolith Radius (um)', pch=19)
GOM_radatage<- lm(Radius~Increments, data=GOM_oto_data)
summary(GOM_radatage)
exes<- 0:14
whys<- summary(GOM_radatage)$coefficients[1,1]+summary(GOM_radatage)$coefficients[2,1]*exes
lines(exes, whys)
text(1, 55, "Radius=9.44+3.16*DI", pos=4, cex=1.5)
# line for fish up to 8 increments, to match with the SS data:
GOM_radatage_sub8inc<- lm(Radius~Increments, data=GOM_oto_data[GOM_oto_data$Increments<9,])
summary(GOM_radatage_sub8inc)
exes2<- 1:8
whys2<- summary(GOM_radatage_sub8inc)$coefficients[1,1]+summary(GOM_radatage_sub8inc)$coefficients[2,1]*exes2
lines(exes2, whys2, lty=2)
text(1, 50, "Radius=9.81+3.05*DI", pos=4)
dev.off()


## Increment width figure:
# Daily increment widths: only for increments with n=3 or more
counts<- apply(incwidthSS, 2, function(x){sum(!is.na(x))})
I<- which(counts>=3)
subincSS<- incwidthSS[,I]
means<- apply(subincSS[,-1], 2, FUN=mean, na.rm=TRUE)
SE<- function(x){sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x)))}
serror<- apply(subincSS[,-1], 2, FUN=SE)
maxes<- means+serror

# do the same for the GoM, but only fish with 8 increments or less
counts<- apply(incwidthGOM_sub8inc, 2, function(x){sum(!is.na(x))})
I<- which(counts>=3)
subincGOM<- incwidthGOM_sub8inc[,I]
meansGOM<- apply(subincGOM[,-1], 2, FUN=mean, na.rm=TRUE)
serrorGOM<- apply(subincGOM[,-1], 2, FUN=SE)
maxesGOM<- meansGOM+serrorGOM
xlimz<- c(0.5,8.5)
ylimz<- c(0, max(c(maxes, maxesGOM), na.rm=TRUE))

png(filename='results/SS_GOM_2016_IncWidth.png', height=6.5, width=7.5, units= 'in', res=300)
exes<- seq(from=0.9, by=1, length.out = length(means))
plot(exes, means, pch=19, xlab='Increment', ylab='Mean Increment Width (um)',
     xlim=xlimz, ylim=ylimz, cex=1.25, cex.lab=1.5, cex.axis=1.5, 
     xaxp=c(1,8,7))
arrows(exes, means-serror, exes, means+serror, length=0.05, angle=90, code=3)
exes<- seq(from=1.1, by=1, length.out = length(meansGOM))
points(exes, meansGOM, pch=2, cex=1.25)
arrows(exes, meansGOM-serrorGOM, exes, meansGOM+serrorGOM, length=0.05, angle=90, code=3)
legend('topleft', legend=c('Slope Sea', 'Gulf of Mexico'), pch=c(19, 2), cex=1.5)
dev.off()


## ANCOVA to test for significantly different slopes in the growth lines: 
# fish with fewer than 8 increments from GoMex and SS
# get the columns we want for SS:
subSS<- SS_oto_data[SS_oto_data$Increments<9,c("Increments", "Length")]
subSS$Region<- "SS"
# get the columns we want for GoMex:
subGOM<- GOM_oto_data[GOM_oto_data$Increments<9, c("Increments", "SL_mm_EtOH")]
names(subGOM)<- c("Increments", "Length")
subGOM$Region<- "GoMex"
foraov<- rbind(subSS, subGOM)
# run the ancova:
aovmodel<- aov(Length~Increments*Region, data=foraov)
summary(aovmodel)
aovmodel2<- aov(Length~Increments+Region, data=foraov)
summary(aovmodel2)