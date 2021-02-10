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



mother.larvipositions$date_of_emergence_min <- as.Date(mother.larvipositions$date_of_emergence_min
                                              , format="%d/%m/%Y") # convert date columns to date format
mother.larvipositions$date_of_death <- as.Date(mother.larvipositions$date_of_death,format="%d/%m/%Y")

mother.larvipositions$name <- as.factor(mother.larvipositions$name)  # change treatment to a factor
# days to death
mother.larvipositions$motherAgeatDeath <- mother.larvipositions$date_of_death - mother.larvipositions$date_of_emergence_min
mother.larvipositions$motherAgeatDeath <- as.numeric(mother.larvipositions$motherAgeatDeath)
mother.larvipositions$motherAgeatDeath[mother.larvipositions$motherAgeatDeath %in% NA] <- 100

#****************************Function to fit models*******************************
fitModFunc <- function(form,modType,dat,whichAIC){
  if(modType=="me"){
   mod <- glmer(form
        ,family=binomial
        ,data=dat
       # , glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5))
        )
  }else{
   mod <- glm(form
        ,family=binomial
        ,data=dat)
  }
  
  out <- summary(mod)
 # coefs <- coefficients(out)
 # coefsci <- confint(mod,devtol = 1e-04)
  
  ll <- logLik(mod)
  parms <- attributes(ll)$df 
  
  if (whichAIC =="aic"){
    aic <- AIC(mod)
  }else{aic <- AICc(mod)}
  
  #return(list(coefs,coefsci,parms,ll,aic))
  list(parms,ll,aic)
}

#*************************************Model formula*******************************
m1 <- abortion ~ mAge*name + motherAgeatDeath*name + (1| adults_id) #r
m1a <- abortion ~ mAge*name + motherAgeatDeath + (1| adults_id) #r
m1b <- abortion ~ mAge + motherAgeatDeath*name + (1| adults_id) #r
m1c <- abortion ~ mAge + motherAgeatDeath + name + (1| adults_id) #r
m2 <- abortion ~ mAge + (1| adults_id )   #r
m2a <- abortion ~ mAge*name + (1| adults_id)   #r
m2b <- abortion ~ mAge + name + (1| adults_id)   #r
m3 <- abortion ~ mAge*name + motherAgeatDeath*name 
m3a <- abortion ~ mAge*name + motherAgeatDeath
m3b <- abortion ~ mAge + motherAgeatDeath*name
m3c <- abortion ~ mAge + motherAgeatDeath + name
m3d <- abortion ~ mAge + motherAgeatDeath
m3e <- abortion ~ mAge*name
m3f <- abortion ~ mAge + name
m4 <- abortion ~ mAge
m5 <- abortion ~ 1
#********************************************************************************

#******************************Function to summarise outputs*********************
modelSummaryFunc <- function(modelFitOutputs,form){
  
  model <-form

  k <- modelFitOutputs[[1]]
  ll <- as.numeric(modelFitOutputs[[2]])
  aic <- modelFitOutputs[[3]]
  
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
fitm1 <- fitModFunc(form=m1,dat=mother.larvipositions,whichAIC="aic",modType="me") #  failed to converge
fitm1a <- fitModFunc(form=m1a,dat=mother.larvipositions,whichAIC="aic",modType="me") # failed to converge
fitm1b <- fitModFunc(form=m1b,dat=mother.larvipositions,whichAIC="aic",modType="me") # failed to converge
fitm1c <- fitModFunc(form=m1c,dat=mother.larvipositions,whichAIC="aic",modType="me") # rescale variables?
fitm2 <- fitModFunc(form=m2,dat=mother.larvipositions,whichAIC="aic",modType="me") 
fitm2a <- fitModFunc(form=m2a,dat=mother.larvipositions,whichAIC="aic",modType="me") # rescale variables, very large eigenvalue ratio?
fitm2b <- fitModFunc(form=m2b,dat=mother.larvipositions,whichAIC="aic",modType="me") 
fitm3 <- fitModFunc(form=m3,dat=mother.larvipositions,whichAIC="aic",modType="fe") 
fitm3a <- fitModFunc(form=m3a,dat=mother.larvipositions,whichAIC="aic",modType="fe") 
fitm3b <- fitModFunc(form=m3b,dat=mother.larvipositions,whichAIC="aic",modType="fe") 
fitm3c <- fitModFunc(form=m3c,dat=mother.larvipositions,whichAIC="aic",modType="fe") 
fitm3d <- fitModFunc(form=m3d,dat=mother.larvipositions,whichAIC="aic",modType="fe") 
fitm3e <- fitModFunc(form=m3e,dat=mother.larvipositions,whichAIC="aic",modType="fe") 
fitm3f <- fitModFunc(form=m4,dat=mother.larvipositions,whichAIC="aic",modType="fe") 
fitm4 <- fitModFunc(form=m4,dat=mother.larvipositions,whichAIC="aic",modType="fe") 
fitm5 <-fitModFunc(form=m5,dat=mother.larvipositions,whichAIC="aic",modType="fe") 

summarym2 <- modelSummaryFunc(fitm2,form=m2)
summarym2b <- modelSummaryFunc(fitm2b,form=m2b)
summarym3 <- modelSummaryFunc(fitm3,form=m3) 
summarym3a <- modelSummaryFunc(fitm3a,form=m3a)
summarym3b <- modelSummaryFunc(fitm3b,form=m3b)
summarym3c <- modelSummaryFunc(fitm3c,form=m3c)
summarym3d <- modelSummaryFunc(fitm3d,form=m3d)
summarym3e <- modelSummaryFunc(fitm3e,form=m3e)
summarym3f <- modelSummaryFunc(fitm3f,form=m3f)
summarym4 <- modelSummaryFunc(fitm4,form=m4)
summarym5 <- modelSummaryFunc(fitm5,form=m5)

list <- list(summarym2,summarym2b,summarym3,summarym3a,summarym3b,summarym3c,summarym3d,summarym3e,summarym3f,summarym4,summarym5)
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

saveRDS(modelSummary, file = "modelSummaryAbortionInclTreatment.rds")

#*****************************Predict*************************

nMod <- glmer(m2b
              ,family=binomial
              ,data=mother.larvipositions)

ilink <- family(nMod)$linkinv

mother.larvipositions$Fitted <- predict(nMod
                                        ,re.form=NA
                                        ,newdat=mother.larvipositions
                                        ,type="response")

mother.larivpositions$fit <- as.numeric(mother.larvipositions$Fitted)
merBootN <- bootMer(nMod, function(x) predict(x, newdata = mother.larvipositions,re.form=NA,type="response"), nsim = 100, re.form = NA)

mother.larvipositions$Lower <- apply(merBootN$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE)))
mother.larvipositions$Upper <- apply(merBootN$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))

mother.larvipositions$name <- as.factor(mother.larvipositions$name)
levels(mother.larvipositions$name) <- c("Control","Mating delay","Nutritional stress")

#**********************Plot*******************************************


#***********calculate mean weight and confidence intervals for females per 10 days***********

mother.larvipositions$mAgeBins <- cut(as.numeric(mother.larvipositions$mAge)
                                      ,c(17,27,37,47,57,67,77,87,97)
                                      ,labels=c("22","32","42","52","62","72","82","92"))
ctrl <- mother.larvipositions[mother.larvipositions$name %in% "Control",]
mate <- mother.larvipositions[mother.larvipositions$name %in% "Mating delay",]
nuts <- mother.larvipositions[mother.larvipositions$name %in% "Nutritional stress",]

#************************************************************************************
ctrl.means <- lapply(unique(ctrl$mAgeBins),function(x){
  temp <- ctrl[ctrl$mAgeBins %in% x,]
  m <- binom.confint(sum(temp$abortion),length(temp$abortion),method="exact")
  m$mAgeBin <- x
  return(m)
})
ctrl.means <- do.call(rbind.data.frame,ctrl.means)

mate.means <- lapply(unique(mate$mAgeBins),function(x){
  temp <- mate[mate$mAgeBins %in% x,]
  m <- binom.confint(sum(temp$abortion),length(temp$abortion),method="exact")
  m$mAgeBin <- x
  return(m)
})
mate.means <- do.call(rbind.data.frame,mate.means)

nuts.means <- lapply(unique(nuts$mAgeBins),function(x){
  temp <- nuts[nuts$mAgeBins %in% x,]
  m <- binom.confint(sum(temp$abortion),length(temp$abortion),method="exact")
  m$mAgeBin <- x
  return(m)
})
nuts.means <- do.call(rbind.data.frame,nuts.means)

nuts.means$name <- "Nutritonal stress"
ctrl.means$name <- "Control"
mate.means$name <- "Mating delay"
abort.means <- rbind.data.frame(nuts.means,ctrl.means,mate.means)

abort.means$mAgeBin <- as.numeric(as.character(abort.means$mAgeBin))
abort.means$name <- as.factor(abort.means$name)
levels(abort.means$name) <- c("Control","Mating delay","Nutritional stress")

names(abort.means)[4] <- "abortion"
names(abort.means)[7] <- "mAge"


tiff("Fig2_abortionsInclTreatment.tiff", height = 5, width = 5, units = 'in', compression="lzw", res=400)
ggplot(mother.larvipositions,aes(as.numeric(mAge),abortion)) +
  scale_y_continuous( breaks=c(0,0.25,0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1)) +
  #geom_jitter(width = 0.1, height = 0.1,col="darkgrey",size=0.5) +
  geom_point(data=abort.means,aes(mAge,abortion,col=name)) +
  geom_errorbar(data=abort.means,aes(x=mAge,ymin=lower,ymax=upper,col=name),width=.1) +
  
  geom_ribbon(data = mother.larvipositions, aes(ymin = Lower, ymax = Upper, x = mAge,fill=name)
              , alpha = 0.5, inherit.aes = FALSE) +
  geom_line(data=mother.larvipositions,aes(as.numeric(mAge),Fitted,col=name)) +
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

