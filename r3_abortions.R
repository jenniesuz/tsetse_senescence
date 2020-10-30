# This script contains the analysis for the abortion data 
source("r1_queries.R")

library(ggplot2)
library(binom)
library(geepack)
library(lme4)
library(nlme)
library(AICcmodavg)
library(MuMIn)
library(reshape)

#*****************************Data***************************
mother.larvipositions$date_of_emergence_min <- as.Date(mother.larvipositions$date_of_emergence_min,format="%d/%m/%Y")
mother.larvipositions$larviposition_date <- as.Date(mother.larvipositions$larviposition_date,format="%d/%m/%Y")
mother.larvipositions$date_of_death <- as.Date(mother.larvipositions$date_of_death,format="%d/%m/%Y")
mother.larvipositions$abortion <- as.integer(mother.larvipositions$abortion)
mother.larvipositions$name <- as.factor(mother.larvipositions$name)
mother.larvipositions$mAge <- mother.larvipositions$larviposition_date - mother.larvipositions$date_of_emergence_min
mother.larvipositions$mAge <- as.numeric(mother.larvipositions$mAge)
levels(mother.larvipositions$name) <- c("Control","Mating delay","Nutritional stress")

abortion.data <- mother.larvipositions[mother.larvipositions$abortion %in% 1,]
abortTreat <- ddply(abortion.data,.(name,adults_id),summarise,l=length(abortion))
#******************************************************************

ctrl <- mother.larvipositions[mother.larvipositions$name %in% "Control",]
mate <- mother.larvipositions[mother.larvipositions$name %in% "Mating delay",]
nuts <- mother.larvipositions[mother.larvipositions$name %in% "Nutritional stress",]
#************************************************************************************


#****************************Function to fit models*******************************
fitModFunc <- function(form,modType,dat,whichAIC){
  if(modType=="me"){
   mod <- glmer(form
        ,family=binomial
        ,data=dat)
  }else{
   mod <- glm(form
        ,family=binomial
        ,data=dat)
  }
  
  out <- summary(mod)
  coefs <- coefficients(out)
  coefsci <- confint(mod)
  
  ll <- logLik(mod)
  parms <- attributes(ll)$df 
  
  if (whichAIC =="aic"){
    aic <- AIC(mod)
  }else{aic <- AICc(mod)}
  
  return(list(coefs,coefsci,parms,ll,aic))
}

#*************************************Model formula*******************************
m1 <- abortion ~ mAge + (1+mAge| adults_id)
m2 <- abortion ~ mAge + (1| adults_id)
m3 <- abortion ~ mAge 
m4 <- abortion ~ 1
#********************************************************************************

#******************************Function to summarise outputs*********************
modelSummaryFunc <- function(modelFitOutputs,form){
  
  model <-form

  k <- modelFitOutputs[[3]]
  ll <- as.numeric(modelFitOutputs[[4]])
  aic <- modelFitOutputs[[5]]
  
  modelSummary <- c(model,k,ll,aic)
  return(modelSummary)
  
}
#********************************************************************************

#***************************Function to return coefficients**********************
coefSummaryFunc <- function(modelFit,model=1){
  
  coefs <- modelFit[[1]]
  est <- coefs[,1]
  ci <- modelFit[[2]]
  parms<-row.names(coefs)
  coefsS <- cbind.data.frame(modelNumber=rep(model,length(names(est))),parameter=parms,lower=ci[names(est),1],est,upper=ci[names(est),2])
  row.names(coefsS) <- c()
  return(coefsS)
}
#*****************************************************************************

#************Control treatment*************

#fitm1C <- fitModFunc(form=m1,dat=ctrl,whichAIC="aic",modType="me") # singular fit

fitm2C <- fitModFunc(form=m2,dat=ctrl,whichAIC="aic",modType="me")
fitm3C <- fitModFunc(form=m3,dat=ctrl,whichAIC="aic",modType="fe")
fitm4C <- fitModFunc(form=m4,dat=ctrl,whichAIC="aic",modType="fe")


summarym2C <- modelSummaryFunc(fitm2C,form=m2)
summarym3C <- modelSummaryFunc(fitm3C,form=m3)
summarym4C <- modelSummaryFunc(fitm4C,form=m4)

list <- list(summarym2C,summarym3C,summarym4C)
model <- as.character(sapply(list,'[[',1))
k <- sapply(list,'[[',2)
ll <- sapply(list,'[[',3)
aic <- sapply(list,'[[',4)

modelSummary <- cbind.data.frame(model,k,ll,aic)
modelSummary <- modelSummary[order(modelSummary$aic),]

modelSummary$deltaAIC <- round(modelSummary$aic - min(modelSummary$aic,na.rm=T),3)
modelSummary$weights <- round(Weights(modelSummary$aic),3)
modelSummary$modelNumber <- row.names(modelSummary)
modelSummary <- modelSummary[,c(7,1,2,3,4,5,6)]

row.names(modelSummary) <- c()

saveRDS(modelSummary, file = "modelSummaryAbortionControl.rds")

coefm2C <- coefSummaryFunc(fitm2C,model=2)
coefm3C <- coefSummaryFunc(fitm3C,model=3)

coefsTable <- rbind.data.frame(coefm3C,coefm2C)
coefsTable <- coefsTable[order(coefsTable$parameter),]
row.names(coefsTable) <- c()
saveRDS(coefsTable,file="coefsAbortionControl.rds")

#***************Mating delay treatment**********************
#fitm1M <- fitModFunc(form=m1,dat=mate,whichAIC="aic",modType="me") # singular fit

fitm2M <- fitModFunc(form=m2,dat=mate,whichAIC="aic",modType="me")
fitm3M <- fitModFunc(form=m3,dat=mate,whichAIC="aic",modType="fe")
fitm4M <- fitModFunc(form=m4,dat=mate,whichAIC="aic",modType="fe")



summarym2M <- modelSummaryFunc(fitm2M,form=m2)
summarym3M <- modelSummaryFunc(fitm3M,form=m3)
summarym4M <- modelSummaryFunc(fitm4M,form=m4)

list <- list(summarym2M,summarym3M,summarym4M)
model <- as.character(sapply(list,'[[',1))
k <- sapply(list,'[[',2)
ll <- sapply(list,'[[',3)
aic <- sapply(list,'[[',4)

modelSummary <- cbind.data.frame(model,k,ll,aic)
modelSummary <- modelSummary[order(modelSummary$aic),]

modelSummary$deltaAIC <- round(modelSummary$aic - min(modelSummary$aic,na.rm=T),3)
modelSummary$weights <- round(Weights(modelSummary$aic),3)
modelSummary$modelNumber <- row.names(modelSummary)
modelSummary <- modelSummary[,c(7,1,2,3,4,5,6)]
row.names(modelSummary) <- c()
saveRDS(modelSummary, file = "modelSummaryAbortionMate.rds")


coefm2M <- coefSummaryFunc(fitm2M,model=2)
coefm3M <- coefSummaryFunc(fitm3M,model=3)
coefsTable <- rbind.data.frame(coefm3M,coefm2M)
coefsTable <- coefsTable[order(coefsTable$parameter),]
row.names(coefsTable) <- c()
saveRDS(coefsTable,file="coefsAbortionMate.rds")

#**********Nutritional stress treatment******
#fitm1N <- fitModFunc(form=m1,dat=nuts,whichAIC="aic",modType="me") # singular fit

fitm2N <- fitModFunc(form=m2,dat=nuts,whichAIC="aic",modType="me")
fitm3N <- fitModFunc(form=m3,dat=nuts,whichAIC="aic",modType="fe")
fitm4N <- fitModFunc(form=m4,dat=nuts,whichAIC="aic",modType="fe")

summarym2N <- modelSummaryFunc(fitm2N,form=m2)
summarym3N <- modelSummaryFunc(fitm3N,form=m3)
summarym4N <- modelSummaryFunc(fitm4N,form=m4)

list <- list(summarym2N,summarym3N,summarym4N)
model <- as.character(sapply(list,'[[',1))
k <- sapply(list,'[[',2)
ll <- sapply(list,'[[',3)
aic <- sapply(list,'[[',4)

modelSummary <- cbind.data.frame(model,k,ll,aic)
modelSummary <- modelSummary[order(modelSummary$aic),]

modelSummary$deltaAIC <- round(modelSummary$aic - min(modelSummary$aic,na.rm=T),3)
modelSummary$weights <- round(Weights(modelSummary$aic),3)
modelSummary$modelNumber <- row.names(modelSummary)
modelSummary <- modelSummary[,c(7,1,2,3,4,5,6)]
row.names(modelSummary) <- c()
saveRDS(modelSummary, file = "modelSummaryAbortionNuts.rds")


coefm2N <- coefSummaryFunc(fitm2N,model=2)
coefm3N <- coefSummaryFunc(fitm3N,model=3)
coefsTable <- rbind.data.frame(coefm2N,coefm3N)
coefsTable <- coefsTable[order(coefsTable$parameter),]
row.names(coefsTable) <- c()
saveRDS(coefsTable,file="coefsAbortionNuts.rds")


#*****************************Predict*************************
cMod <- glm(m3
            ,family=binomial
            ,data=ctrl)
mMod <-  glm(m3
             ,family=binomial
             ,data=mate)
nMod <- glmer(m2
              ,family=binomial
              ,data=nuts)
  
#****************Control*******************
ilink <- family(cMod)$linkinv
ctrlP <- with(ctrl,
                data.frame(mAge = seq(min(mAge), max(mAge),
                                            length = 100)))
ctrlP <- cbind(ctrlP, predict(cMod, ctrlP, type = "link", se.fit = TRUE)[1:2])
ctrlP <- transform(ctrlP, Fitted = ilink(fit), Upper = ilink(fit + (1.96 * se.fit)),
                     Lower = ilink(fit - (1.96 * se.fit)))
ctrlP$name <- "Control"
#************Mating delay******************
ilink <- family(mMod)$linkinv
mateP <- with(mate,
              data.frame(mAge = seq(min(mAge), max(mAge),
                                    length = 100)))
mateP <- cbind(mateP, predict(mMod, mateP, type = "link", se.fit = TRUE)[1:2])
mateP <- transform(mateP, Fitted = ilink(fit), Upper = ilink(fit + (1.96 * se.fit)),
                   Lower = ilink(fit - (1.96 * se.fit)))
mateP$name <- "Mating delay"
#***********Nutritional stress************
ilink <- family(nMod)$linkinv
nutsP <- with(nuts,
              data.frame(mAge = seq(min(mAge), max(mAge),
                                    length = 100)))
nutsP$Fitted <- predict(nMod,re.form=NA,newdat=nutsP,type="response")
nutsP$fit <- as.numeric(nutsP$Fitted)
merBootN <- bootMer(nMod, function(x) predict(x, newdata = nutsP,re.form=NA,type="response"), nsim = 100, re.form = NA)

nutsP$Lower <- apply(merBootN$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE)))
nutsP$Upper <- apply(merBootN$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))


nutsP$name <- "Nutritional stress"


newdat <- rbind.data.frame(ctrlP[,c("mAge","Fitted","Upper","Lower","name")],mateP[,c("mAge","Fitted","Upper","Lower","name")],nutsP[,c("mAge","Fitted","Upper","Lower","name")])

newdat$name <- as.factor(newdat$name)
levels(newdat$name) <- c("Control","Mating delay","Nutritional stress")

#**********************Plot*******************************************

tiff("Fig2_abortions.tiff", height = 5, width = 5, units = 'in', compression="lzw", res=400)
ggplot(mother.larvipositions,aes(as.numeric(mAge),abortion)) +
  scale_y_continuous( breaks=c(0,0.25,0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1)) +
  #geom_jitter(width = 0.1, height = 0.1,col="darkgrey",size=0.5) +
  geom_ribbon(data = newdat, aes(ymin = Lower, ymax = Upper, x = mAge,fill=name)
              , alpha = 0.5, inherit.aes = FALSE) +
  geom_line(data=newdat,aes(as.numeric(mAge),Fitted,col=name)) +
  labs(  x="Maternal age (days)"
         ,y="Probability of abortion"
         #,title="b)"
  ) +
  scale_color_manual(values = c("#875777","#eab051","#c0b9ac")) + 
  scale_fill_manual(values = c("#875777","#eab051","#c0b9ac")) +
  theme_set(theme_bw()) +
  theme(axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=8)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=9)
        ,legend.position =c(0.2,0.8)
        ,legend.title = element_blank()
        ,strip.background = element_rect(colour="white", fill="white")
        ,panel.border = element_blank()
  ) 
dev.off()




tiff("FigRawAbortions.tiff", height = 3, width = 6, units = 'in', compression="lzw", res=400)
ggplot(mother.larvipositions) +
  geom_jitter(aes(x=mAge,y=abortion)
              ,position = position_jitter(height = .02)
              ,alpha=0.4
              ,size=0.5) +
  ylab("Abortion") +
  xlab("Maternal age (days)") +
  scale_y_continuous(breaks=c(0,1),labels=c("0","1")) +
  theme_set(theme_bw()) +
  theme(axis.line = element_line(color = 'black')
        ,text=element_text(size=9)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=7)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,legend.position =c(0.2,0.9)
        ,legend.title = element_blank()
        ,strip.background = element_rect(colour="white", fill="white")
        ,panel.border = element_blank()
  ) +
  facet_wrap(~name)
dev.off()



#**********Nutritional stress treatment with mothers still alive only******
deadNuts
nuts <- nuts[!nuts$adults_id %in% deadNuts,]

fitm2N <- fitModFunc(form=m2,dat=nuts,whichAIC="aic",modType="me")
fitm3N <- fitModFunc(form=m3,dat=nuts,whichAIC="aic",modType="fe")
fitm4N <- fitModFunc(form=m4,dat=nuts,whichAIC="aic",modType="fe")

summarym2N <- modelSummaryFunc(fitm2N,form=m2)
summarym3N <- modelSummaryFunc(fitm3N,form=m3)
summarym4N <- modelSummaryFunc(fitm4N,form=m4)

list <- list(summarym2N,summarym3N,summarym4N)
model <- as.character(sapply(list,'[[',1))
k <- sapply(list,'[[',2)
ll <- sapply(list,'[[',3)
aic <- sapply(list,'[[',4)

modelSummary <- cbind.data.frame(model,k,ll,aic)
modelSummary <- modelSummary[order(modelSummary$aic),]

modelSummary$deltaAIC <- round(modelSummary$aic - min(modelSummary$aic,na.rm=T),3)
modelSummary$weights <- round(Weights(modelSummary$aic),3)
modelSummary$modelNumber <- row.names(modelSummary)
modelSummary <- modelSummary[,c(7,1,2,3,4,5,6)]
row.names(modelSummary) <- c()
saveRDS(modelSummary, file = "modelSummaryAbortionNutsAlive.rds")


coefm2N <- coefSummaryFunc(fitm2N,model=2)
coefm3N <- coefSummaryFunc(fitm3N,model=3)
coefsTable <- rbind.data.frame(coefm2N,coefm3N)
coefsTable <- coefsTable[order(coefsTable$parameter),]
row.names(coefsTable) <- c()
saveRDS(coefsTable,file="coefsAbortionNutsAlive.rds")
