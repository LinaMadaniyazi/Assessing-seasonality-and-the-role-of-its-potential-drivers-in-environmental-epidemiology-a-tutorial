# SUPPLEMENTARY MATERIAL:
# "QUANTIFYING SEASONALITY CHANGES IN ENVIRONMENTAL EPIDEMIOLOGY STUDIES "
# Lina Madaniyazi, Aurelio Tobias, Yoonhee Kim, Yeonseung Chung, Ben Armstrong, Masahiro Hashizume
# Journal 2021
#
# Jan 2022
##############################################################################

library(mgcv);library(tsModel);library(dlnm);library(splines);library(mondate);library(ggplot2); library(devtools)

# LOAD THE DATA
df<-read.csv(url('https://raw.githubusercontent.com/gasparrini/2014_gasparrini_BMCmrm_Rcodedata/master/london.csv'))
##############################################################################

# YEAR
df$year <- as.factor(as.character(df$year))
years <- unique(df$year)

#DAY OF WEEK
df$dow <- as.factor(as.character(df$dow))

#STRATUM WITH YEAR AND DOW TO CONTROL FOR LONG TERM TREND AND EFFECT OF DAY OF WEEK
df$stratum <- as.factor(df$year:df$dow)

# CREATE SEASONALITY INDICATOR (DAY OF YEAR)
for(t in years){
  tempind <- df$year == t
  temp <- as.numeric( strftime( df$date[tempind] , format = "%j") )
  if( length(temp) == 365 ){ temp[60:365] <- 61:366 }
  else if( length(temp) == 364 ){ temp[60:364] <- 61:365 }

  df$seasonal[tempind] <- temp
}

# DERIVE THE CROSS-BASIS FOR TEMPERATURE
cb.temp <- crossbasis(df$tmean, lag=21,
                      argvar=list(fun="ns",knots=quantile(df$tmean,c(.25,.50,.75),na.rm=T)),
                      arglag=list(fun="ns",knots=logknots(21,3)))

# DERIVE THE CYCLIC SPLINE FOR SEASONALITY
source("cyclic.R")
spline.season <- cyclic( df$seasonal, df = 4  )
##############################################################################
# RUN THE MODEL AND OBTAIN COEFFICIENTS

# RUN THE MODEL
# TEMPERATURE UNADJUSTED
fit1<- glm(death ~ spline.season + factor(stratum), data=df, family=quasipoisson(),na.action="na.exclude")
# TEMPERATURE ADJUSTED
fit2<- glm(death ~ spline.season + factor(stratum) +cb.temp,  data=df, family=quasipoisson(),na.action="na.exclude")

# COEFFICIENTS FROM CYCLIC SPLINE
coef1<- coef(fit1)[ 2:5 ]
vcov1<- vcov(fit1)[ 2:5, 2:5 ]
coef2<- coef(fit2)[ 2:5 ]
vcov2<- vcov(fit2)[ 2:5, 2:5 ]

# COEFFICIENTS FOR EACH DAY OF YEAR (LOGRR, SE)
logrr1 <- matrix(NA, 1, 366 )
logrr_se1<- matrix(NA, 1, 366 )
logrr2 <- matrix(NA, 1, 366 )
logrr_se2 <- matrix(NA, 1, 366 )

Trange <- 1:366

bvar<- cyclic( df$seasonal, df = 4  )

# LOGRR
l <- length(Trange)
for(j in 1:l){ logrr1[1,j] <-  t(bvar[j,]) %*% as.numeric( coef1)  }
for(j in 1:l){ logrr2[1,j] <-  t(bvar[j,]) %*% as.numeric( coef2)  }

# SE OF LOGRR
for(j in 1:l){ logrr_se1[1,j] <- sqrt( as.numeric( t(bvar[j,]) %*% vcov1 %*% (bvar[j,]) ) ) }
for(j in 1:l){ logrr_se2[1,j] <- sqrt( as.numeric( t(bvar[j,]) %*% vcov2%*% (bvar[j,]) ) ) }

# CENTERING LOGRR AT TROUGH
maxcen<- apply(logrr1, 1, which.min )
maxcen2<- apply(logrr2, 1, which.min )

bvar_cen<- bvar[maxcen,]
bvar_cen2<- bvar[maxcen2,]

# CENTERED LOGRR
l <- length(Trange)
for(j in 1:l){ logrr1[1,j] <- as.numeric( t(bvar[j,]-bvar_cen) %*% as.numeric( coef1) ) }
for(j in 1:l){ logrr2[1,j] <- as.numeric( t(bvar[j,]-bvar_cen2) %*% as.numeric( coef2) )}

# SD FOR CENTERED LOGRR
for(j in 1:l){ logrr_se1[1,j] <- sqrt( as.numeric( t(bvar[j,]-bvar_cen) %*% vcov1%*% (bvar[j,]-bvar_cen) ) ) }
for(j in 1:l){ logrr_se2[1,j] <- sqrt( as.numeric( t(bvar[j,]-bvar_cen2) %*% vcov2%*% (bvar[j,]-bvar_cen2) ) ) }

##############################################################################
# SUMMARY AND COMPARISON OF KEY FEATURES

# TROUGH (95%CI)
source("findmin.R")

trough1_est <- findmin(spline.season,fit1)
trough1<-findmin(spline.season,fit1,sim=T)

trough2_est <- findmin(spline.season,fit2)
trough2<- findmin(spline.season,fit2,sim=T)

trough<-rbind(trough1,trough2)
rownames(trough)<-c("Unadjusted","Adjusted")

# PEAK (95%CI)
source("findmax.R")

peak1_est <- findmax(spline.season,fit1)
peak1<-findmax(spline.season,fit1,sim=T)

peak2_est <- findmax(spline.season,fit2)
peak2<-findmax(spline.season,fit2,sim=T)

peak<-rbind(peak1,peak2)
rownames(peak)<-c("Unadjusted","Adjusted")

# TIMINGS
timings<-cbind(peak,trough)
print(timings)

# PEAK-TO-TROUGH (95%CI)
max1<-apply(logrr1, 1, which.max )
logrr1_max<-apply(logrr1, 1, max,na.rm=TRUE)
se1<- logrr_se1[cbind(seq_along(max1), max1)]
ptr1 <-cbind(peak=exp(logrr1_max), peak.low=exp(logrr1_max-1.96*se1),peak.high=exp( logrr1_max+1.96*se1))

max2<-apply(logrr2, 1, which.max )
logrr2_max <-apply(logrr2, 1, max,na.rm=TRUE)
se2 <- logrr_se2 [cbind(seq_along(max2), max2)]
ptr2 <-cbind(peak=exp(logrr2_max), peak.low=exp( logrr2_max-1.96*se2),peak.high=exp( logrr2_max+1.96*se2))

ptr<-rbind(Unadjusted=ptr1,Adjusted=ptr2)

# DIFFERENCE BETWEEN PTR BEFORE AND AFTER ADJUSTMENT
# ABSOLUTE DIFFERENCE
delta<-logrr2_max-logrr1_max
delta_se<-sqrt(se1^2+se2^2)
ad_ptr<-c(change=delta,change.low=delta-1.96*delta_se,change.high=delta+1.96*delta_se)
#RELATIVE DIFFERENCE
rd_ptr<-exp(ad_ptr)
print(rd_ptr)

# ATTRIBUTABLE FRACTION
source("attrs.R")
af1_est <- attrs(df$seasonal,spline.season,df,fit1,type="af",tot=T)
simaf1 <- attrs (df$seasonal,spline.season,df,fit1,type="af",tot=T,sim = T)
af1.low <- quantile(simaf1,c(2.5)/100)
af1.up <- quantile(simaf1,c(97.5)/100)
af1<-c(af=af1_est,af1.low,af1.up)

af2_est <- attrs (df$seasonal,spline.season,df,fit2,type="af",tot=T) #un-adjusted
simaf2 <- attrs (df$seasonal,spline.season,df,fit2,type="af",tot=T,sim = T)
af2.low <- quantile(simaf2,c(2.5)/100)
af2.up <- quantile(simaf2,c(97.5)/100)
af2<-c(af=af2_est,af2.low,af2.up)

af<-rbind(Unadjusted=af1,Adjusted=af2)
colnames(af)[2:3]<-c("af.low","af.up")

# DIFFERENCE BETWEEN AF BEFORE AND AFTER ADJUSTMENT
delta<-af2_est-af1_est
delta_se<-sqrt(((af1.up-af1.low)/1.96)^2+((af2.up-af2.low)/1.96)^2)
delta_af<-c(change=delta,change.low=delta-1.96*delta_se,change.high=delta+1.96*delta_se)

#SUMMARY AND COMPARISON OF KEY FEATURES (TABLE 2)
table2<-rbind(cbind(timings,ptr,af),Changes=c(peak2_est-peak1_est,NA,NA,trough2_est-trough1_est,NA,NA,rd_ptr,delta_af))
table2<-round(table2,digits = 3)
print(table2)


# FIGURE 2 (WITHOUT POINTS AND CI FOR PEAK AND TROUGH)#
# CREATE DAY OF YEAR FOR X AXIS
df$data<-as.Date(df$date,"%Y-/%m-/%d")
DY <-as.Date( df$date[df$year == 2000])
dateY <- format(DY, format="%b %d")
M <- mondate("1-1-2000")
FDM <- as.Date( rev(M - 1:12) )
FDM <- c( FDM, as.Date( DY[length(DY)] ) )
dateM <- format(FDM, format="%b %d")
monIND <- which( as.character(dateY) %in% as.character(dateM) )

plot(Trange, exp(logrr1), type='l', xlab="Day of year", ylab="RR (95% CI)",  ylim=c(0.9,1.4), lwd=2, xaxt='n',col="black" )
axis(1, at=Trange[monIND], labels=FALSE)
text(x=Trange[monIND], y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), labels=dateM, srt=45, adj=1, xpd=TRUE, cex=0.8)
polygon(c(Trange, rev(Trange)), c(exp(logrr1-1.96*logrr_se1), rev(exp(logrr1+1.96*logrr_se1))),
        col = adjustcolor("black", alpha.f = 0.10) ,border = NA)
# 95%CI FOR RR
lines(Trange, exp(logrr1-1.96*logrr_se1),lty=2,col="black" )
lines(Trange, exp(logrr1+1.96*logrr_se1),lty=2,col="black" )

lines(Trange, exp(logrr2),lwd=2,col="red")
polygon(c(Trange, rev(Trange)), c(exp(logrr2-1.96*logrr_se2), rev(exp(logrr2+1.96*logrr_se2))),
        col = adjustcolor("red", alpha.f = 0.10) ,border = NA)
# 95%CI FOR RR
lines(Trange, exp(logrr2-1.96*logrr_se2),lty=2,col="red")
lines(Trange, exp(logrr2+1.96*logrr_se2),lty=2,col="red")

abline(h=1, col="black")

##############################################################################
# SENSITIVITY ANALYSIS (COSINOR FUNCTION) #

# Run the model (Cosinor)
fit1_f<- glm(death ~  sin(2*pi*seasonal/366) + cos(2*pi*seasonal/366) + factor(dow)+ factor(stratum),data=df, family=quasipoisson(),na.action="na.exclude")                   # temperature unadjusted
fit2_f<- glm(death ~ sin(2*pi*seasonal/366) + cos(2*pi*seasonal/366) + factor(dow)+ factor(stratum) +cb.temp,  data=df,family=quasipoisson(),na.action="na.exclude")          # temperature adjusted

# DERIVE RR FOR EACH DAY-OF-YEAR; PTR; AF
# RR
b01 <- coef(fit1_f)[1]
alpha1 <- coef(fit1_f)[2]
beta1 <- coef(fit1_f)[3]

b02 <- coef(fit2_f)[1]
alpha2 <- coef(fit2_f)[2]
beta2 <- coef(fit2_f)[3]

doy<-c(1:366)
f0_1<-beta1*cos(2*pi*(doy/366))+ alpha1*sin(2*pi*(doy/366))
f1<-exp(f0_1)
f0_2<-beta2*cos(2*pi*(doy/366))+alpha2*sin(2*pi*(doy/366))
f2<-exp(f0_2)

# PTR
ptr1_f<-max(f1)/min(f1)
ptr2_f<-max(f2)/min(f2)

# AF
cases<-matrix(by(df$death,df$seasonal,FUN = mean))
af1_f <- 1-exp(-log(f1/min(f1)))
an1_f <- af1_f*cases[,1]
isna <- is.na(an1_f)
an1_f <- sum(an1_f,na.rm=T)
an1_f <- an1_f/sum(cases[!isna])

af2_f <- 1-exp(-log(f2/min(f2)))
an2_f <- af2_f*cases[,1]
isna <- is.na(an2_f)
an2_f <- sum(an2_f,na.rm=T)
an2_f <- an2_f/sum(cases[!isna])

# Function for computing the Q-AIC
QAIC <- function(model) {
  phi <- summary(model)$dispersion
  loglik <- sum(dpois(model$y, model$fitted.values, log=TRUE))
  return(-2*loglik + 2*summary(model)$df[3]*phi)
}

# Table S1 (Reporting QAIC)
# Cyclic Spline
QAIC(fit1);QAIC(fit2)
# Cosinor
QAIC(fit1_f); QAIC(fit2_f)

# Figure S1
plot(Trange, exp(logrr1), type='l', xlab="Day of year", ylab="RR (95% CI)",  ylim=c(0.98,1.35), lwd=4, xaxt='n',col="orangered",cex.lab=1.2 ) # Fitter curve from cyclic spline before temperature adjustment
axis(1, at=Trange[monIND], labels=FALSE)
text(x=Trange[monIND], y=par()$usr[3]-0.02*(par()$usr[4]-par()$usr[3]), labels=dateM, srt=45, adj=1, xpd=TRUE, cex=1)
lines(Trange,exp(logrr2), lwd=4, xaxt='n',col="orangered",lty=2) # Fitter curve from cyclic spline after temperature adjustment

lines(Trange,f1/min(f1),col="blue", lwd=4) # Fitter curve from cosinor before temperature adjustment
lines(Trange,f2/min(f2),lwd=4,col="blue",lty=2) # Fitter curve from cosinor after temperature adjustment
abline(h=1, col="black")

legend(x = "topright",          # Position
       inset = c(-0.02, 0),
       bty = "n", # Removes the legend box
       legend = c("Cyclic spline,temperature unadjustment", "Cyclic spline, temperature adjustment","Cosinor, temperature unadjustment","Cosnior, temperature adjustment"),  # Legend texts
       lty = c(1, 2,1,2),           # Line types
       col = c("orangered", "orangered","blue","blue"),           # Line colors
       cex = 1, # Change legend size
       lwd = 2,                 # Line width
       xpd = TRUE,y.intersp=2)











