# This script analyses pupal wet weight and produces Figure 4 of the manuscript
source("r1_queries.R")

library(plyr)        
library(reshape)
library(ggplot2)
library(nlme)
library(AICcmodavg)
library(MuMIn)

#*******************************************Data*******************************************************
pupae$date_of_emergence_min <- as.Date(pupae$date_of_emergence_min,format="%d/%m/%Y")
pupae$larviposition_date <- as.Date(pupae$larviposition_date,format="%d/%m/%Y")
pupae$date_of_death <- as.Date(pupae$date_of_death,format="%d/%m/%Y")
pupae$abortion <- as.integer(pupae$abortion)
pupae$name <- as.factor(pupae$name)
pupae$mAgeDays <- pupae$larviposition_date - pupae$date_of_emergence_min
pupae$mAgeDays <- as.numeric(pupae$mAgeDays)
levels(pupae$name) <- c("Control","Mating delay","Nutritional stress")


adult.deaths$date_of_emergence_min <- as.Date(adult.deaths$date_of_emergence_min
                                              , format="%d/%m/%Y") # convert date columns to date format
adult.deaths$date_of_death <- as.Date(adult.deaths$date_of_death,format="%d/%m/%Y")

adult.deaths$name <- as.factor(adult.deaths$name)  # change treatment to a factor
# days to death
adult.deaths$mAgeatLastObs <- adult.deaths$date_of_death - adult.deaths$date_of_emergence_min
adult.deaths$mAgeatLastObs <- as.numeric(adult.deaths$mAgeatLastObs)

pupae$mAgeatLastObs <- NA

mad <- lapply(unique(pupae$adults_id),function(x){
  temp <- pupae[pupae$adults_id %in% x,] 
  temp$mAgeatLastObs <- adult.deaths$mAgeatLastObs[adult.deaths$adults_id %in% x]
  return(temp)
})


pupae <- do.call(rbind,mad)
pupae$mAgeatLastObs[pupae$mAgeatLastObs %in% NA] <- 100

names(pupae)
names(pupae)[24] <- "mAgeatLastObs"

pupWeight <- pupae[pupae$abortion %in% 0,] # take out abortions
ddply(pupWeight,.(name),summarise,l=length(wet_weight))

numWeighed <- ddply(pupWeight,.(adults_id),summarise,num=length(wet_weight))
mean(numWeighed$num) # 4.8
range(numWeighed$num) # 1 - 8

numWeighed <- ddply(pupWeight,.(name,adults_id),summarise,num=length(wet_weight))
ddply(numWeighed,.(name),summarise,mean(num))
ddply(numWeighed,.(name),summarise,range(num))

pupWeight <- pupWeight[!pupWeight$wet_weight%in%NA,]     # remove those - 31 of 1239 for which no wet weight data - either lost or crushed

#**********************wet weight vs fat************

pupWeight$fat <- pupWeight$dry_weight - pupWeight$residual_dry_weight
fat <- pupWeight[pupWeight$fat > 1,]
fat <- fat[!fat$fat %in% NA,]

corr <- lm(fat ~ wet_weight, data=fat)
summary(corr)


#tiff("FigWeightFat.tiff", height = 5, width = 5, units = 'in', compression="lzw", res=400)
ggplot(fat) +
  geom_point(aes(x=wet_weight,y=fat,col=name),size=0.8) +
  labs(x="Wet weight (mg)"
       ,y="Fat (mg)"
       ,col="Treatment") +
  scale_color_manual(values=c("#875777","#eab051","#c0b9ac")) +
  theme_set(theme_bw()) +
  theme(axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=8)
        ,legend.key.size = unit(0.5,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=7)
        ,legend.position =c(0.9,0.2)
        ,legend.title = element_text(size=7)
        ,strip.background = element_rect(colour="white", fill="white")
        ,panel.border = element_blank()
        #,strip.text.x = element_blank()
  )
#dev.off()



#***********calculate mean weight and confidence intervals for females per 10 days***********
pupWeight$mAgeBins <- cut(as.numeric(pupWeight$mAgeDays)
                          ,c(17,27,37,47,57,67,77,87,97)
                          ,labels=c("22","32","42","52","62","72","82","92"))


weight.means <- ddply(pupWeight,.(name,mAgeBins),summarise,mean.weight = mean(wet_weight,na.rm=T))
weight.sd <- ddply(pupWeight,.(name,mAgeBins),summarise,sd.weight =sd(wet_weight,na.rm=T))
weight.size <- ddply(pupWeight,.(name,mAgeBins),summarise,size=length(wet_weight))
weight.ci <- 1.96*(weight.sd$sd.weight/sqrt(weight.size$size))
weight.means$upper <- weight.means$mean.weight+weight.ci
weight.means$lower <- weight.means$mean.weight-weight.ci
weight.means

weight.means$mAgeBins <- as.numeric(as.character(weight.means$mAgeBins))
weight.means$name <- as.factor(weight.means$name)
levels(weight.means$name) <- c("Control","Mating delay","Nutritional stress")

weight.means <- weight.means[-10,]

#******************************************************************************************
nuts <- pupWeight[pupWeight$name %in% "Nutritional stress",]
ctrl <- pupWeight[pupWeight$name %in% "Control",]
mate <- pupWeight[pupWeight$name %in% "Mating delay",]

# 
mean(mate$wet_weight,na.rm=T) # 32.57
mean(ctrl$wet_weight,na.rm=T) # 33.54
mean(nuts$wet_weight,na.rm=T) # 28.02


#*********************Fit model func*****************
fitModRand <- function(fixedEffects,randomEffects,dat,whichAIC){
  lmeCtrl <- lmeControl(opt='optim')
 if(is.na(randomEffects)){
   mod <- lm(fixedEffects,data=dat)
   coefsci <- confint(mod)
   coefs <- coefficients(mod)
   coefs <- cbind.data.frame(lower=coefsci[,1],est.=coefs,upper=coefsci[,2])
   rand <- c(NA,NA)
 }else{
  mod <- lme(fixedEffects
            ,random = randomEffects
            ,data=dat,method="ML",control=lmeCtrl)
  coefs <- intervals(mod,which="fixed")
  coefs <- data.frame(coefs$fixed)
  names(coefs) <- c("lower","est.","upper")
  rand <- VarCorr(mod)
  rand <- data.frame(rand[1:length(rand[1,]),1:2])
  }
  out <- summary(mod)

  ll <- logLik(mod)
  parms <- attributes(ll)$df - 1

  if (whichAIC =="aic"){
    aic <- AIC(mod)
  }else{aic <- AICc(mod)}#,nobs=length(dat$name))}

  return(list(coefs,parms,ll,aic,rand))
}

#************************************************************************
bartlett.test(wet_weight ~ name,data=pupWeight)
#*********************
lmeCtrl <- lmeControl(opt='optim')

mod <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2)+mAgeatLastObs+name
           ,random = ~1+mAgeDays+I(mAgeDays^2)| adults_id
           ,data=pupWeight,method="ML",control=lmeCtrl)
plot(mod)
pupWeight$resids <- residuals(mod)
plot(pupWeight$resids,pupWeight$mAge)
boxplot(pupWeight$resids~pupWeight$name)


bartlett.test(resids~ name,data=pupWeight)


#********************************************************************

mod <- gls(wet_weight ~ mAgeDays+I(mAgeDays^2)+mAgeatLastObs+name
           ,data=pupWeight)
plot(mod)
pupWeight$resids <- residuals(mod)
plot(pupWeight$resids,pupWeight$mAge)
boxplot(pupWeight$resids~pupWeight$name)
AIC(mod)

vn <- varIdent(form= ~ 1 | name)

mod <- gls(wet_weight ~ mAgeDays+I(mAgeDays^2)+mAgeatLastObs+name
           ,data=pupWeight,weights=vn)

summary(mod)
pupWeight$resids <- residuals(mod)
boxplot(pupWeight$resids~pupWeight$name)
AIC(mod)









#********************************************************************


# formulas incl. random effects
fe1 <- wet_weight ~ mAgeDays+I(mAgeDays^2)+I(mAgeDays^3)+mAgeatLastObs
re1 <- ~1+mAgeDays+I(mAgeDays^2)+I(mAgeDays^3)| adults_id
re1a <- ~1| adults_id

fe2 <-  wet_weight ~ mAgeDays+I(mAgeDays^2)+mAgeatLastObs
re2 <-  ~1+mAgeDays+I(mAgeDays^2)| adults_id
re2a <- ~1| adults_id

fe3 <- wet_weight ~ log(mAgeDays)+mAgeatLastObs
re3 <- ~1+log(mAgeDays)| adults_id
re3a <- ~1| adults_id

fe4 <- wet_weight ~ mAgeDays+mAgeatLastObs
re4 <- ~1+mAgeDays| adults_id
re4a <- ~1| adults_id

fe5 <- wet_weight ~ mAgeatLastObs
re5 <- ~1| adults_id



combinations <- list(c(fe1,re1)
                     ,c(fe1,re1a)
                     ,c(fe1,NA)
                     ,c(fe2,re2)
                     ,c(fe2,re2a)
                     ,c(fe2,NA)
                     ,c(fe3,re3)
                     ,c(fe3,re3a)
                     ,c(fe3,NA)
                     ,c(fe4,re4)
                     ,c(fe4,re4a)
                     ,c(fe4,NA)
                     ,c(fe5,re5)
                     ,c(fe5,NA))


modelSummaryFunc <- function(modelFitOutputs,name= "modelSummaryWetWeightControl.rds"){
  
  fixedEffects <- as.character(sapply(modelFitOutputs,'[[',1))
  randomEffects <- as.character(sapply(modelFitOutputs,'[[',2))
  
  outPuts <- lapply(modelFitOutputs,'[[',3)
  k <- sapply(outPuts,'[[',2)
  ll <- sapply(outPuts,'[[',3)
  aic <- sapply(outPuts,'[[',4)
  
  modelSummary <- cbind.data.frame(fixedEffects,randomEffects,k,ll,aic)
  modelSummary$modelNumber <- row.names(modelSummary)
  modelSummary <- modelSummary[order(modelSummary$aic),]
  modelSummary <- modelSummary[,c(6,1,2,3,4,5)]
  
  modelSummary$deltaAIC <- round(modelSummary$aic - min(modelSummary$aic),3)
  modelSummary$weights <- round(Weights(modelSummary$aic),3)
  row.names(modelSummary) <- c()
  saveRDS(modelSummary, file = name)
  return(modelSummary)
  
}

coefSummaryFunc<- function(modelFitOutputs,name= "coefsWetWeightControl.rds",mods=modelsC){

  outPuts <- lapply(modelFitOutputs,'[[',3)
  aics <- as.numeric(lapply(outPuts,'[[',4))
  cfs <- lapply(outPuts,'[[',1)
  
  cfs2 <- lapply(1:length(cfs),function(x){
    temp <- cfs[[x]]
    temp$modelNumber <- x
    temp$parameter <- row.names(temp)
    row.names(temp) <- c()
    temp$aic <- aics[[x]]
    return(temp[,c(4,5,1,2,3,6)])
  })
  
  cfs3 <- do.call(rbind,cfs2)
  cfs3 <- cfs3[order(cfs3$aic),]
  cfs3 <- cfs3[,-6]
  cfs3 <- cfs3[order(cfs3$parameter),]
  modelsWeightZero <- mods$modelNumber[mods$weights == 0]
  cfs3 <- cfs3[!cfs3$modelNumber %in% as.numeric(modelsWeightZero),]
  row.names(cfs3) <- c()
  saveRDS(cfs3, file = name)
  
  op <- lapply(modelFitOutputs,'[[',3)
  randomEffects <- lapply(op,'[[',5)
  re2 <- lapply(1:length(randomEffects),function(x){
    temp <- randomEffects[[x]]
    if(!is.na(temp)){
    temp$modelNumber <- x
    temp$parameter <- row.names(temp)
    row.names(temp) <- c()
    temp$aic <- aics[[x]]
    return(temp[,c(3,4,1,2,5)])
    }
  })
  
  re3 <- do.call(rbind,re2)
  re3 <- re3[order(re3$aic),]
  re3 <- re3[,-5]
  re3[,3] <- round(as.numeric(as.character(re3[,3])),3)
  re3[,4] <- round(as.numeric(as.character(re3[,4])),3)
  re3 <- re3[order(re3$parameter),]
  re3 <- re3[!re3$modelNumber %in% as.numeric(modelsWeightZero),]
  row.names(re3) <- c()
  saveRDS(re3,file=paste("RE",name))
  
  return(cfs3)
}
#***********************************Control treatment***************************************


fitModsC <- lapply(4:length(combinations),function(x){
  temp <- combinations[[x]]
  fe <- as.formula(as.character(temp[1]))
  re <-  NA
  if(!is.na(temp[2])){re <- as.formula(as.character(temp[2]))}
  output <- fitModRand(fixedEffects=fe,randomEffects=re,dat=ctrl,whichAIC="aicc")
  return(list(temp[1],temp[2],output))
})

modelsC <- modelSummaryFunc(fitModsC,name="modelSummaryWetWeightControl.rds")

coefSummaryFunc(modelFitOutputs=fitModsC,name="coefsWetWeightControl.rds",mods=modelsC)



#**************Predict**************************
lmeCtrl <- lmeControl(opt='optim')

modCtrl <- lme(fe2
           ,random = re2
           ,data=ctrl,method="ML",control=lmeCtrl)

modMate <- lme(fe2
               ,random = re2
               ,data=mate,method="ML",control=lmeCtrl)

modNuts <- lme(fe2
               ,random = re2
               ,data=nuts,method="ML",control=lmeCtrl)

summary(modCtrl)
summary(modMate)
summary(modNuts)

ctrl$pred1 <- predict(modCtrl,newdata=ctrl,level=0)
ctrl$pred2 <- predict(modCtrl,newdata=ctrl,level=1)
mate$pred1 <- predict(modMate,newdata=mate,level=0)
mate$pred2 <- predict(modMate,newdata=mate,level=1)
nuts$pred1 <- predict(modNuts,newdata=nuts,level=0)
nuts$pred2 <- predict(modNuts,newdata=nuts,level=1)

Designmat <- model.matrix(eval(eval(modCtrl$call$fixed)[-2]), ctrl[-ncol(ctrl)])
predvar <- diag(Designmat %*% modCtrl$varFix %*% t(Designmat))
ctrl$SE <- sqrt(predvar) 
ctrl$CI <- 1.96*sqrt(predvar+modCtrl$sigma^2)

Designmat <- model.matrix(eval(eval(modMate$call$fixed)[-2]), mate[-ncol(mate)])
predvar <- diag(Designmat %*% modMate$varFix %*% t(Designmat))
mate$SE <- sqrt(predvar) 
mate$CI <- 1.96*sqrt(predvar+modMate$sigma^2)#

Designmat <- model.matrix(eval(eval(modNuts$call$fixed)[-2]), nuts[-ncol(nuts)])
predvar <- diag(Designmat %*% modNuts$varFix %*% t(Designmat))
nuts$SE <- sqrt(predvar) 
nuts$CI <- 1.96*sqrt(predvar+modNuts$sigma^2)



all <- rbind.data.frame(ctrl,mate,nuts)
all$name <- as.factor(all$name)
levels(all$name) <- c("Control","Mating delay","Nutritional stress")


#*******************************Plot************************************************
tiff("Fig3_wet_weight.tiff", height = 3, width = 6, units = 'in', compression="lzw", res=400)
ggplot(all,aes()) +
  scale_color_manual(values=c("#875777","#eab051","#c0b9ac")) +
  geom_line(data=all,aes(x=mAgeDays,y=pred2,group=adults_id,col=name),alpha=0.5,lwd=.2) + 
  geom_line(data=all,aes(x=mAgeDays,y=pred1,col=name)) + 
  xlim(15,100) +
  geom_point(data=weight.means,aes(x=mAgeBins,y=mean.weight,col=name),size=0.8) +
  geom_errorbar(data=weight.means,aes(x=mAgeBins,ymin=lower, ymax=upper,col=name), width=.1) +
  labs(  x="Maternal age (days)"
         ,y="Offspring wet weight (mg)"
  )+ 
  theme_set(theme_bw()) +
  theme(
    axis.line = element_line(color = 'black')
    ,text=element_text(size=10)
    ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
    ,axis.text=element_text(size=8)
    ,legend.key.size = unit(0.8,"line")
    ,legend.background = element_blank()
    ,legend.text=element_text(size=9)
    ,legend.position ="none"
    ,legend.title = element_blank()
    ,strip.background = element_rect(colour="white", fill="white")
    ,panel.border = element_blank()
  ) +
  facet_wrap(~name,ncol=3,scales="free_x")
dev.off()

#********plot population prediction and prediction intervals**************
#*
#*
ggplot(all,aes()) +
  scale_color_manual(values=c("#875777","#eab051","#c0b9ac")) +
 # geom_line(data=all,aes(x=mAgeDays,y=pred2,group=adults_id,col=name),alpha=0.5,lwd=.2) + 
  geom_line(data=all,aes(x=mAgeDays,y=pred1,col=name)) + 
  xlim(15,100) +
  geom_point(data=weight.means,aes(x=mAgeBins,y=mean.weight,col=name),size=0.8) +
  geom_errorbar(data=weight.means,aes(x=mAgeBins,ymin=lower, ymax=upper,col=name), width=.1) +
  labs(  x="Maternal age (days)"
         ,y="Offspring wet weight (mg)"
  )+ 
  geom_ribbon(data = all, aes(ymin = pred1-CI, ymax = pred1+CI, x =mAgeDays,fill=name)
              , alpha = 0.5, inherit.aes = FALSE) +
  
  theme_set(theme_bw()) +
  theme(
    axis.line = element_line(color = 'black')
    ,text=element_text(size=10)
    ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
    ,axis.text=element_text(size=8)
    ,legend.key.size = unit(0.8,"line")
    ,legend.background = element_blank()
    ,legend.text=element_text(size=9)
    ,legend.position ="none"
    ,legend.title = element_blank()
    ,strip.background = element_rect(colour="white", fill="white")
    ,panel.border = element_blank()
  ) 



#*************************raw data plot**********************

tiff("FigRawWetweight.tiff", height = 3, width = 6, units = 'in', compression="lzw", res=400)
ggplot(pupWeight) +
  geom_point(aes(x=mAgeDays,y=wet_weight)
             ,size=0.5) +
  ylab("Offspring wet weight (mg)") +
  xlab("Maternal age (days)") +
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


#************residuals against explanatory variable***********
ggplot(all,aes()) +
  geom_point(data=all,aes(x=mAgeDays,y=wet_weight-pred2,group=adults_id),alpha=0.5,lwd=.2) + 
  xlim(15,100) +
  labs(  x="Mother age (days)"
         ,y="Offspring wet weight (mg)"
  )+ 
  theme_set(theme_bw()) +
  theme(
    axis.line = element_line(color = 'black')
    ,text=element_text(size=10)
    ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
    ,axis.text=element_text(size=8)
    ,legend.key.size = unit(0.8,"line")
    ,legend.background = element_blank()
    ,legend.text=element_text(size=9)
    ,legend.position ="none"
    ,legend.title = element_blank()
    ,strip.background = element_rect(colour="white", fill="white")
    ,panel.border = element_blank()
  ) +
  facet_wrap(~name,ncol=3,scales="free_x")


pupae[!pupae$wet_weight %in% NA,]
lmeCtrl <- lmeControl(opt='optim')

mod <- lme(wet_weight ~ mAgeDays*name+I(mAgeDays^2)*name+mAgeatLastObs*name
           ,random = ~1+mAgeDays+I(mAgeDays^2)| adults_id
           ,data=pupWeight,method="ML",control=lmeCtrl)
