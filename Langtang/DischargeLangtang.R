################################################################################
# DischargeLangtang - prepare all discharge time series from raw data and rating curves
# 
# DischargeLangtang.R
#
# ReadMe: 
# 
# (1) Kyanjing
# (2) Langshisha
# (3) Lirung

# Output:

#
# Created:          2019/08/15
# Latest Revision:  2019/08/22
#
# Jakob F Steiner| PhD candidate | Faculty of Geosciences | Universiteit Utrecht | Princetonlaan 8a, 3584 CB Utrecht 
# Vening Meinesz building, room 4.30 | P.O. Box 80.115, 3508 TC Utrecht | j.f.steiner@uu.nl | www.uu.nl/staff/jfsteiner | www.mountainhydrology.org 
################################################################################
# clear entire workspace (excl. packages)
rm(list = ls())
gc()

# define &-sign for pasting string-elements
'&' <- function(...) UseMethod('&')
'&.default' <- .Primitive('&')
'&.character' <- function(...) paste(...,sep='')

# install necessary packages if not available yet via install.packages()
library(pacman)
library(RColorBrewer)
require('stats')
require('ggplot2')
require('lubridate')

#source("F:\\PHD\\BasicCode\\timeseriesAnalysis\\TSaggregate.R")

path_RawData_Kyanjin <- 'F:\\PhD\\Research\\CatchmentHydrology\\Discharge\\ProcessingRawData\\RAWData\\Kyangjin'
path_RawData_Langshisha <- 'F:\\PhD\\Research\\CatchmentHydrology\\Discharge\\ProcessingRawData\\RAWData\\Langshisha'
path_RawData_Lirung <- 'F:\\PhD\\Research\\CatchmentHydrology\\Discharge\\ProcessingRawData\\RAWData\\Lirung'


##########################
# Get Discharge data Kyanjing
##########################

# Convert Pressure to dh
DataStation <- read.csv(path_RawData_Kyanjin&'\\Kyangjin_HS_ICIMOD.csv')
DateStation <- as.POSIXlt(paste(DataStation$DATE,DataStation$TIME,sep=' '),format='%Y-%m-%d %H:%M:%S')

# Convert data where only presure is available
noNA <- which(!is.na(DataStation$P_HOBO))
DataStation$ATMPRESS[62318:62388] <- 64.2
pCorr <- DataStation$P_HOBO - (DataStation$ATMPRESS)
pCorr[pCorr <= 6] <- NA

wHCorr <- pCorr / 1000 * 1000 / 9.81

# Correct time series 
DataStation$TS1 <- DataStation$SENSORDEPTH
DataStation$TS2 <- DataStation$STGHT

shift <- DataStation$TS1 - DataStation$TS2
DataStation$TS1[shift >=0.3] <- NA # The pressure tranducer shifted


DataStation$TS_final <- DataStation$TS2
DataStation$TS_final[which(!is.na(DataStation$TS1)&is.na(DataStation$TS2))] <- DataStation$TS1[which(!is.na(DataStation$TS1)&is.na(DataStation$TS2))] 
DataStation$TS_final[!is.na(wHCorr)] <- wHCorr[!is.na(wHCorr)]

# Get rating curve
DataAll <- read.csv(path_RawData_Kyanjin&'\\Kyangjin_RatingCurve.csv')
# Remove data with quality flag (9999)
faultySH <- which(DataAll$QUAL == 9999)
DataAll <- DataAll[-faultySH,]

date <- as.POSIXlt(DataAll$DATE,format='%m/%d/%Y')
datetime <- as.POSIXlt(paste(DataAll$DATE,DataAll$TIME,sep=' '),format='%m/%d/%Y %H:%M:%S')
hRating <- hour(datetime) / 24
mRating <- month(datetime)
mRating[is.na(mRating)] <- 1

# Transitions from high to low flow rating curve
cutoff_low_14 <- 1.05 # FINAL VALUE 1.05
cutoff_high_14 <- 1.2 # FINAL VALUE 1.2
cutoff_mean_14 <- 1.15 # FINAL VALUE 1.15

cutoff_low_19 <- 1.2  # FINAL VALUE 1.2
cutoff_high_19 <- 1.4 # FINAL VALUE 1.4
cutoff_mean_19 <- 1.4 # FINAL VALUE 1.4

#ID_man_low <- which(!is.na(DataAll$MANSTAGHT) & is.na(DataAll$PRESSSTGHT) & DataAll$MANSTAGHT<=cutoff)
ID_salt <- which(DataAll$TECH==1)

ID_pres <- which(!is.na(DataAll$PRESSSTGHT))
ID_man <- which(!is.na(DataAll$MANSTAGHT) & is.na(DataAll$PRESSSTGHT))
DataAll$finalStageHeight <- DataAll$PRESSSTGHT
DataAll$finalStageHeight[ID_man] <- DataAll$MANSTAGHT[ID_man]
DataAll$finalStageHeight <- DataAll$finalStageHeight
equ <- function(stage, a,C, n)   numericDeriv(quote(a+C*(stage)^n),c("a","C", "n"), parent.frame())


###############
# Rating Curve 2013 to 2019
###############
tim_2 <- which(year(datetime)<=2019 & year(datetime)>=2013)
datafile_2 <- DataAll[tim_2,]

cutoff_low <- cutoff_low_19
cutoff_high <- cutoff_high_19
cutoff_mean <- cutoff_mean_19
######## 
# Extract only low flows
########

lowID <- which(datafile_2$finalStageHeight<cutoff_high)

# Writing data into an R dataframe:
gaugings_low <- data.frame(stage=c(datafile_2$finalStageHeight[lowID])) # Original Data in [m]
discharge_low <- c(datafile_2$Q[lowID] / 1000) # Original Data in [m3/s]

# Fitting the power law
# Note that the start argument requires a list with initial (rough!) estimates of the coefficients to be estimated
power.nls<-nls(discharge_low~equ(stage,a, C, n),  data = gaugings_low,
               start=list(a=20,C=2.73,  n=1.5),
               trace = T,
               lower=list(a = 0,C=0.1,  n = 0),
               upper=list(a = 50,C=20,  n=10),
               algorithm="port",
               control=nls.control(maxiter = 100))
# Viewing the model summary and accessing estimated constants
summary(power.nls)
C_low<-coef(power.nls)["C"]
a_low<-coef(power.nls)["a"]
n_low<-coef(power.nls)["n"]
cor(discharge_low,predict(power.nls))

# Compute standard errors for predicted values
x<-seq(0, cutoff_mean, 0.003) # range for visualization
se.fit <- sqrt(apply(attr(predict(power.nls,list(stage = x)),"gradient"),1, function(x) sum(vcov(power.nls)*outer(x,x))))

conf_low <- cbind(x[x<cutoff_high],se.fit[x<cutoff_high])

######## 
# Extract only high flows
######## 

highID <- which(datafile_2$finalStageHeight>cutoff_low)


# Writing data into an R dataframe:
gaugings_high <- data.frame(stage=c(datafile_2$finalStageHeight[highID])) # Original Data in [m]
discharge_high <- c(datafile_2$Q[highID] / 1000) # Original Data in [m3/s]

# Fitting the power law
# Note that the start argument requires a list with initial (rough!) estimates of the coefficients to be estimated
power.nls2<-nls(discharge_high~equ(stage,a, C, n),  data = gaugings_high,
                start=list(a=20,C=2.73,  n=1.5),
                trace = T,
                lower=list(a = 0,C=0.1,  n = 0),
                upper=list(a = 50,C=20,  n=10),
                algorithm="port", control=nls.control(maxiter = 100))
# Viewing the model summary and accessing estimated constants
summary(power.nls)
C_high<-coef(power.nls2)["C"]
a_high<-coef(power.nls2)["a"]
n_high<-coef(power.nls2)["n"]
cor(discharge_high,predict(power.nls2))
confint(power.nls2,  level = 0.95)

# Compute standard errors for predicted values
x<-seq(0, 2.5, 0.003) # range for visualization
se.fit2 <- sqrt(apply(attr(predict(power.nls2,list(stage = x)),"gradient"),1, function(x) sum(vcov(power.nls2)*outer(x,x))))

conf_high <- cbind(x[x>=cutoff_low],se.fit2[x>=cutoff_low])


# Visualize Rating curve 
#png(file=paste0(path_code,'/',FinalFigureName,'_Raw.png'), res = 160,width=1200,height=1200)
# Plotting data and fitted curve
oldpar<-par(mar=c(3.5,3.75,1.5,0.5), mgp=c(2.35,1,0))
x <- x<-seq(0, cutoff_high, 0.003)
x1 <- x
matplot(x, predict(power.nls, list(stage=x))+outer(se.fit[1:length(x)], qnorm(c(0.5, 0.025, 0.975))), type="l", lty=c(1,2,2), col=c("black", "red", "red"), xlab="stage (m)", ylab=expression(paste("discharge (", m^3,s^-1, ")")), ylim=c(0, 30), xlim=c(0, 2.5), main="Rating curve 2013 - 2019")
with(gaugings_low, points(stage, discharge_low))
x <- x<-seq(cutoff_high, 2.5, 0.003)
matplot(x, predict(power.nls2, list(stage=x))+outer(se.fit2[length(x1)+1:length(x)], qnorm(c(0.5, 0.025, 0.975))), type="l", lty=c(1,2,2), col=c("black", "blue", "blue"), add=T,xlim=c(1.2,3))
with(gaugings_high, points(stage, discharge_high),xlim=c(1.2,3))
##points(gaugings$stage[1:4], discharge[1:4], col="red", pch=20, cex=1.5)
##legend("topleft", legend="Highflow discharges where 1m is added to the stage", pch=20, col="red", cex=)
grid()
par(oldpar)
grid()

par_1619 <- rbind(cbind(a_low,C_low,n_low),cbind(a_high,C_high,n_high))
colnames(par_1619) <- c('a','C','n')
rownames(par_1619) <- c('low','high')
conf_range_1619 <- rbind(conf_low,conf_high)
colnames(conf_range_1619) <- cbind('stage height','covariance')

write.csv(par_1619, file = path_RawData_Kyanjin&'\\Kyangjin_Params_1319.csv')
write.csv(conf_range_1619, file = path_RawData_Kyanjin&'\\Kyangjin_confidence_1319.csv')

# Produce final Discharge time series

Qfinal <- DataStation$TS_final * 0 
Qfinal[which(DataStation$TS_final<=cutoff_high_19)] <- par_1619[1,1] + par_1619[1,2]*DataStation$TS_final[which(DataStation$TS_final<=cutoff_high_19)]^par_1619[1,3]
Qfinal[which(DataStation$TS_final>cutoff_high_19)] <- par_1619[2,1] + par_1619[2,2]*DataStation$TS_final[which(DataStation$TS_final>cutoff_high_19)]^par_1619[2,3]

# Water Temperature merged
DataStation$TEMP[!is.na(DataStation$TEMP2)] <- DataStation$TEMP2[!is.na(DataStation$TEMP2)]
finDF <- data.frame(cbind(DataStation$DATE,DataStation$TIME,Qfinal,DataStation$TEMP))
colnames(finDF) <- c('DATE','TIME','Q [m3/s]','water temperature [C]')

write.csv(finDF, file = path_RawData_Kyanjin&'\\Kyanjing_Q_20132019.csv',row.names=F)


##########################
# Get Discharge data Langshisha old station (2013 to 2016)
##########################

# For the Langshisha Station we have relatively few data points and hence the rating curve 
# has a larger uncetainty. I suggest to play around with the range (right now it is set to
# a very tight 2.5%/97.5% percentile)
# Convert Pressure to dh
DataStation <- read.csv(path_RawData_Langshisha&'\\Langshisha_HS_ICIMOD.csv')
DateStation <- as.POSIXlt(paste(DataStation$DATE,DataStation$TIME,sep=' '),format='%Y-%m-%d %H:%M:%S')

# Get rating curve
DataAll <- read.csv(path_RawData_Langshisha&'\\LangshishaRatingCurve_old.csv')

date <- as.POSIXlt(DataAll$DATE,format='%m/%d/%Y')
datetime <- as.POSIXlt(paste(DataAll$DATE,DataAll$TIME,sep=' '),format='%m/%d/%Y %H:%M:%S')
hRating <- hour(datetime) / 24
mRating <- month(datetime)
mRating[is.na(mRating)] <- 1


equ <- function(stage, a,C, n)   numericDeriv(quote(a+C*(stage)^n),c("a","C", "n"), parent.frame())

# Writing data into an R dataframe:
gaugings_high <- data.frame(stage=c(DataAll$PRESSSTGHT)) # Original Data in [m]
discharge_high <- c(DataAll$Q / 1000) # Original Data in [m3/s]

# Fitting the power law
# Note that the start argument requires a list with initial (rough!) estimates of the coefficients to be estimated
power.nls2<-nls(discharge_high~equ(stage,a, C, n),  data = gaugings_high,
                start=list(a=20,C=2.73,  n=1.5),
                trace = T,
                lower=list(a = 0,C=0.1,  n = 0),
                upper=list(a = 50,C=20,  n=10),
                algorithm="port", control=nls.control(maxiter = 100))
# Viewing the model summary and accessing estimated constants
C_high<-coef(power.nls2)["C"]
a_high<-coef(power.nls2)["a"]
n_high<-coef(power.nls2)["n"]
cor(discharge_high,predict(power.nls2))
confint(power.nls2,  level = 0.95)

# Compute standard errors for predicted values
x<-seq(0, 0.5, 0.003) # range for visualization
se.fit2 <- sqrt(apply(attr(predict(power.nls2,list(stage = x)),"gradient"),1, function(x) sum(vcov(power.nls2)*outer(x,x))))

# Visualize the final rating curve
#png(file=paste0(path_code,'/',FinalFigureName,'_Raw.png'), res = 160,width=1200,height=1200)
# Plotting data and fitted curve
oldpar<-par(mar=c(3.5,3.75,1.5,0.5), mgp=c(2.35,1,0))
x<-seq(0, 0.5, 0.003)
x1 <- x
matplot(x, predict(power.nls2, list(stage=x))+outer(se.fit[1:length(x)], qnorm(c(0.5, 0.025, 0.975))), type="l", lty=c(1,2,2), col=c("black", "red", "red"), xlab="stage (m)", ylab=expression(paste("discharge (", m^3,s^-1, ")")), ylim=c(0, 2), xlim=c(0, 0.5))
with(gaugings_high, points(stage, discharge_high))
##points(gaugings$stage[1:4], discharge[1:4], col="red", pch=20, cex=1.5)
##legend("topleft", legend="Highflow discharges where 1m is added to the stage", pch=20, col="red", cex=)
grid()
par(oldpar)
grid()

conf_range <- cbind(x,se.fit2[1:length(x)])

par_langshold <- rbind(cbind(a_high,C_high,n_high))
colnames(par_langshold) <- c('a','C','n')

write.csv(par_langshold, file = path_RawData_Langshisha&'\\Langshisha_old_Params.csv')
write.csv(conf_range, file = path_RawData_Langshisha&'\\Langshisha_old_confidence.csv')


Qfinal_Langshisha <- DataStation$SHGT * 0 
Qfinal_Langshisha <- par_langshold[1,1] + par_langshold[1,2]*(DataStation$SHGT/100/9.81)^par_langshold[1,3]

finDF <- data.frame(cbind(DataStation$DATE,DataStation$TIME,Qfinal_Langshisha))
colnames(finDF) <- c('DATE','TIME','Q [m3/s]')

write.csv(finDF, file = path_RawData_Langshisha&'\\Langshisha_Q_20132016.csv',row.names=F,overwrite=T)

##########################
# Get Discharge data Langshisha new station (2017)
##########################

# Convert Pressure to dh
DataStation <- read.csv(path_RawData_Langshisha&'\\201704_LangshishaHS.csv')
DateStation <- as.POSIXlt(paste(DataStation$DATE,DataStation$TIME,sep=' '),format='%Y-%m-%d %H:%M:%S')

DataStation$STGHT <- (DataStation$WPRESS - DataStation$APRESS) / 9.81
DataStation$STGHT[DataStation$STGHT<0] <- NA

# Get rating curve
DataAll <- read.csv(path_RawData_Langshisha&'\\LangshishaRatingCurve_new.csv')

date <- as.POSIXlt(DataAll$DATE,format='%m/%d/%Y')
datetime <- as.POSIXlt(paste(DataAll$DATE,DataAll$TIME,sep=' '),format='%m/%d/%Y %H:%M:%S')
hRating <- hour(datetime) / 24
mRating <- month(datetime)
mRating[is.na(mRating)] <- 1

equ <- function(stage, a,C, n)   numericDeriv(quote(a+C*(stage)^n),c("a","C", "n"), parent.frame())

# Writing data into an R dataframe:
gaugings_high <- data.frame(stage=c(DataAll$PRESSSTGHT)) # Original Data in [m]
discharge_high <- c(DataAll$Q / 1000) # Original Data in [m3/s]

# Fitting the power law
# Note that the start argument requires a list with initial (rough!) estimates of the coefficients to be estimated
power.nls2<-nls(discharge_high~equ(stage,a, C, n),  data = gaugings_high,
                start=list(a=20,C=2.73,  n=1.5),
                trace = T,
                lower=list(a = 0,C=0.1,  n = 0),
                upper=list(a = 50,C=20,  n=10),
                algorithm="port", control=nls.control(maxiter = 100))
# Viewing the model summary and accessing estimated constants
C_high<-coef(power.nls2)["C"]
a_high<-coef(power.nls2)["a"]
n_high<-coef(power.nls2)["n"]
cor(discharge_high,predict(power.nls2))
confint(power.nls2,  level = 0.95)

# Compute standard errors for predicted values
x<-seq(0, 0.7, 0.003) # range for visualization
se.fit2 <- sqrt(apply(attr(predict(power.nls2,list(stage = x)),"gradient"),1, function(x) sum(vcov(power.nls2)*outer(x,x))))

#png(file=paste0(path_code,'/',FinalFigureName,'_Raw.png'), res = 160,width=1200,height=1200)
# Plotting data and fitted curve
oldpar<-par(mar=c(3.5,3.75,1.5,0.5), mgp=c(2.35,1,0))
x <- x<-seq(0, 0.7, 0.003)
x1 <- x
matplot(x, predict(power.nls2, list(stage=x))+outer(se.fit[1:length(x)], qnorm(c(0.5, 0.025, 0.975))), type="l", lty=c(1,2,2), col=c("black", "red", "red"), xlab="stage (m)", ylab=expression(paste("discharge (", m^3,s^-1, ")")), ylim=c(0, 3), xlim=c(0, 0.75))
with(gaugings_high, points(stage, discharge_high))
##points(gaugings$stage[1:4], discharge[1:4], col="red", pch=20, cex=1.5)
##legend("topleft", legend="Highflow discharges where 1m is added to the stage", pch=20, col="red", cex=)
grid()
par(oldpar)
grid()

conf_range <- cbind(x,se.fit2[1:length(x)])

par_langshnew <- rbind(cbind(a_high,C_high,n_high))
colnames(par_langshnew) <- c('a','C','n')

write.csv(par_langshnew, file = path_RawData_Langshisha&'\\Langshisha_new_Params.csv')
write.csv(conf_range, file = path_RawData_Langshisha&'\\Langshisha_new_confidence.csv')

Qfinal_Langshisha <- DataStation$STGHT * 0 
Qfinal_Langshisha <- par_langshnew[1,1] + par_langshnew[1,2]*(DataStation$STGHT)^par_langshnew[1,3]
Qfinal_Langshisha[Qfinal_Langshisha>25] <- NA

finDF <- data.frame(cbind(DataStation$DATE,DataStation$TIME,Qfinal_Langshisha,DataStation$WT,DataStation$AT))
colnames(finDF) <- c('DATE','TIME','Q [m3/s]','water temperature [C]','air temperature [C]')

write.csv(finDF, file = path_RawData_Langshisha&'\\Langshisha_Q_20172018.csv',row.names=F)

##########################
# Get Discharge data Lirung Station
##########################

# Convert Pressure to dh
DataStation <- read.csv(path_RawData_Lirung&'\\Lirung_HS_ICIMOD.csv')
DateStation <- as.POSIXlt(paste(DataStation$DATE,DataStation$TIME,sep=' '),format='%Y-%m-%d %H:%M:%S')


# Get rating curve
DataAll <- read.csv(path_RawData_Lirung&'\\Lirung_HS_ICIMOD_RatingCurve.csv')

date <- as.POSIXlt(DataAll$DATE,format='%Y-%m-%d')
datetime <- as.POSIXlt(paste(DataAll$DATE,DataAll$HOUR,sep=' '),format='%Y-%m-%d %H:%M:%S')
hRating <- hour(datetime) / 24
mRating <- month(datetime)
mRating[is.na(mRating)] <- 1

equ <- function(stage, a,C, n)   numericDeriv(quote(a+C*(stage)^n),c("a","C", "n"), parent.frame())

# Writing data into an R dataframe:

Rlevel <- DataAll$RLEV
Rlevel[-5]
Qmeas <- DataAll$DISC
Qmeas[-5]

gaugings_high <- data.frame(stage=c(Rlevel)) # Original Data in [m]
discharge_high <- c(Qmeas) # Original Data in [m3/s]


# Fitting the power law
# Note that the start argument requires a list with initial (rough!) estimates of the coefficients to be estimated
power.nls2<-nls(discharge_high~equ(stage,a, C, n),  data = gaugings_high,
                start=list(a=20,C=2.73,  n=1.5),
                trace = T,
                lower=list(a = 0,C=0.1,  n = 0),
                upper=list(a = 50,C=20,  n=10),
                algorithm="port", control=nls.control(maxiter = 100))
# Viewing the model summary and accessing estimated constants
C_high<-coef(power.nls2)["C"]
a_high<-coef(power.nls2)["a"]
n_high<-coef(power.nls2)["n"]
cor(discharge_high,predict(power.nls2))
confint(power.nls2,  level = 0.95)

# Compute standard errors for predicted values
x<-seq(0, 1, 0.003) # range for visualization
se.fit2 <- sqrt(apply(attr(predict(power.nls2,list(stage = x)),"gradient"),1, function(x) sum(vcov(power.nls2)*outer(x,x))))

#png(file=paste0(path_code,'/',FinalFigureName,'_Raw.png'), res = 160,width=1200,height=1200)
# Plotting data and fitted curve
oldpar<-par(mar=c(3.5,3.75,1.5,0.5), mgp=c(2.35,1,0))
x1 <- x
matplot(x, predict(power.nls2, list(stage=x))+outer(se.fit[1:length(x)], qnorm(c(0.5, 0.025, 0.975))), type="l", lty=c(1,2,2), col=c("black", "red", "red"), xlab="stage (m)", ylab=expression(paste("discharge (", m^3,s^-1, ")")), ylim=c(0, 3.3), xlim=c(0, 1))
with(gaugings_high, points(stage, discharge_high))
##points(gaugings$stage[1:4], discharge[1:4], col="red", pch=20, cex=1.5)
##legend("topleft", legend="Highflow discharges where 1m is added to the stage", pch=20, col="red", cex=)
grid()
par(oldpar)
grid()

par_lirung <- rbind(cbind(a_high,C_high,n_high))
colnames(par_lirung) <- c('a','C','n')

conf_range <- cbind(x,se.fit2[1:length(x)])

write.csv(par_lirung, file = path_RawData_Lirung&'\\Lirung_Params.csv')
write.csv(conf_range, file = path_RawData_Lirung&'\\Lirung_confidence.csv')

Qfinal_Lirung <- DataStation$SGHT * 0 
Qfinal_Lirung <- par_lirung[1,1] + par_lirung[1,2]*(DataStation$SGHT)^par_lirung[1,3]

finDF <- data.frame(cbind(DataStation$DATE,DataStation$TIME,Qfinal_Lirung,DataStation$WTEMP))
colnames(finDF) <- c('DATE','TIME','Q [m3/s]','water temperature [C]')

write.csv(finDF, file = path_RawData_Lirung&'\\Lirung_Q_all.csv',row.names=F)