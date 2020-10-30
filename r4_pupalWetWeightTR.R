
library(plyr)        
library(reshape)
library(ggplot2)
library(nlme)
library(lmeSplines)
library(AICcmodavg)
library(MuMIn)
library(bbmle)

source("r1_queries.R")
#*******************************************Data*******************************************************
pupae$date_of_emergence_min <- as.Date(pupae$date_of_emergence_min,format="%d/%m/%Y")
pupae$larviposition_date <- as.Date(pupae$larviposition_date,format="%d/%m/%Y")
pupae$date_of_death <- as.Date(pupae$date_of_death,format="%d/%m/%Y")
pupae$abortion <- as.integer(pupae$abortion)
pupae$name <- as.factor(pupae$name)
pupae$mAgeDays <- pupae$larviposition_date - pupae$date_of_emergence_min
pupae$mAgeDays <- as.numeric(pupae$mAgeDays)
levels(pupae$name) <- c("Control","Mating delay","Nutritional stress")

pupWeight <- pupae[pupae$abortion %in% 0,] # take out abortions
pupWeight <- pupWeight[!pupWeight$wet_weight%in%NA,]     # remove those - 31 of 1239 for which no wet weight data - either lost or crushed

nuts <- pupWeight[pupWeight$name %in% "Nutritional stress",]
ctrl <- pupWeight[pupWeight$name %in% "Control",]
mate <- pupWeight[pupWeight$name %in% "Mating delay",]

#*******************assuming 1 break point**********************
nllFunc<-function(logp,dat=ctrl){
  p <- exp(logp)
dat$time1 <- p   # create dummy variable for the times before the breakpoint
dat$time1[dat$mAgeDays<p] <- dat$mAgeDays[dat$mAgeDays <p] # here assume breakpoint is 60 days
dat$time2 <- p  # create dummy variable for the times after the breakpoint
dat$time2[dat$mAgeDays>=p] <- dat$mAgeDays[dat$mAgeDays>=p]

lmeCtrl <- lmeControl(opt='optim')
 mod <- lme(wet_weight ~ time1 + time2
             ,random= ~1+time1 + time2| adults_id
             , data=dat,method="ML",control=lmeCtrl)

nll <- as.numeric(-logLik(mod))
return(nll)
}

fitCtrl <- mle2(nllFunc,start=list(logp=log(60)),method="Brent"
              ,data=ctrl
              ,parameters=list(logp)
              ,lower=log(20),upper=log(100))

exp(coef(fitCtrl))

fitMate <- mle2(nllFunc,start=list(logp=log(60)),method="Brent"
                ,data=mate
                ,parameters=list(logp)
                ,lower=log(20),upper=log(100))

exp(coef(fitMate))

fitNuts <- mle2(nllFunc,start=list(logp=log(60)),method="Brent"
                ,data=nuts
                ,parameters=list(logp)
                ,lower=log(20),upper=log(100))

exp(coef(fitNuts))


p1 <- 50
p2 <- 75


#*******************assuming 2 break points**********************
nllFunc<-function(logp1,logp2,dat=ctrl){
  p1 <- exp(logp1)
  p2 <- exp(logp2)
  if( (p1 < 45)  | (p1 >= 55) | (p2 < 55)  | (p2 >= 85)){ nll <- 1000000  }else{
  dat$time1 <- p1   # create dummy variable for the times before the breakpoint
  dat$time1[dat$mAgeDays<p1] <- dat$mAgeDays[dat$mAgeDays <p1] # here assume breakpoint is 60 days
  dat$time2 <- dat$mAgeDays
  dat$time2 <- p1  # create dummy variable for the times after the breakpoint
  dat$time2[(dat$mAgeDays>=p1)&(dat$mAgeDays<p2)] <- dat$mAgeDays[(dat$mAgeDays>=p1)&(dat$mAgeDays<p2)]
  dat$time2[(dat$mAgeDays>p2)] <- p2
  dat$time3 <- p2  # create dummy variable for the times after the breakpoint
  dat$time3[dat$mAgeDays>=p2] <- dat$mAgeDays[dat$mAgeDays>=p2]
  
  lmeCtrl <- lmeControl(opt='optim')
  mod <- lme(wet_weight ~ time1 + time2 + time3
            ,random= ~1+time1 + time2 + time3 | adults_id
             , data=dat,method="ML",control=lmeCtrl)
  
  nll <- as.numeric(-logLik(mod))
 }
  return(nll)
}

fitCtrl <- mle2(nllFunc,start=list(logp1=log(50),logp2=log(75)),method="Nelder-Mead"
                ,data=ctrl
                ,parameters=list(logp1,logp2)
                ,skip.hessian=T)

exp(coef(fitCtrl))


fitMate <- mle2(nllFunc,start=list(logp1=log(50),logp2=log(60)),method="Nelder-Mead"
                ,data=mate
                ,parameters=list(logp1,logp2))

exp(coef(fitMate))


fitNuts <- mle2(nllFunc,start=list(logp=log(60)),method="Brent"
                ,data=nuts
                ,parameters=list(logp)
                ,lower=log(20),upper=log(100))

exp(coef(fitNuts))

