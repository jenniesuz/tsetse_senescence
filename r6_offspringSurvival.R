# This script analyses offspring survival and produces figure 5 of the manuscript

library(plyr)          # packages required
library(reshape)
library(ggplot2)
library(gridExtra)
library(AICcmodavg)
library(MuMIn)
library(lme4)
library(nlme)
source("r1_queries.R")    # query database 

#******************************************Data*********************************************
pupae$date_emerged <- as.Date(pupae$date_emerged,format=c("%d/%m/%Y"))
pupae$date_of_death <- as.Date(pupae$date_of_death,format=c("%d/%m/%Y"))
pupae$date_of_emergence_min <- as.Date(pupae$date_of_emergence_min, format="%d/%m/%Y")  
pupae$larviposition_date <- as.Date(pupae$larviposition_date, format="%d/%m/%Y")        
pupae$mAgeDays <- pupae$larviposition_date - pupae$date_of_emergence_min
pupae$mAgeDays <- as.numeric(pupae$mAgeDays)

pupEmerged <- pupae[pupae$emerged %in% "1.0",]            
pupEmerged <- pupEmerged[!pupEmerged$wet_weight %in% NA,]
pupEmerged$daysSurv <- pupEmerged$date_of_death - pupEmerged$date_emerged 
pupEmerged$dead <- 1 

pupEmerged$name <- as.factor(pupEmerged$name)

levels(pupEmerged$name) <- c("Control","Mating delay","Nutritional \n stress")


tiff("fig_offspringSurvHist.tiff", height = 3, width = 6, units = 'in', compression="lzw", res=400)
ggplot(pupEmerged) +
  geom_histogram(aes(daysSurv)) +
  labs(  x="Days surviving"
         ,y="Number of flies"
  ) +
  theme_set(theme_bw()) +
  theme(axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=8)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=9)
        ,legend.position =c(0.2,0.9)
        ,legend.title = element_blank()
        ,strip.background = element_rect(colour="white", fill="white")
        ,panel.border = element_blank()
  ) +
  facet_wrap(~name) 
dev.off()

weights <- quantile(pupEmerged$wet_weight)
pupEmerged$sex <- as.factor(pupEmerged$sex)
levels(pupEmerged$sex) <- c("F","M")

#tiff("fig_offspringWeightSex.tiff", height = 5, width = 4, units = 'in', compression="lzw", res=400)
weightSexPlot <- ggplot(pupEmerged) +
  geom_boxplot(aes(y=wet_weight,x=sex,fill=name)) +
  scale_fill_manual(values=c("#875777","#eab051","#c0b9ac")) +
  ylab("Wet weight") +
  xlab("Sex") +
  labs(title="a)") +
  theme_set(theme_bw()) +
  theme(axis.line = element_line(color = 'black')
        ,text=element_text(size=7)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=7)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=7)
        ,legend.position ="none"
        ,legend.title = element_blank()
        ,strip.background = element_rect(colour="white", fill="white")
        ,panel.border = element_blank()
  ) +
  facet_wrap(~name) 
#dev.off()

pupEmerged$daysSurv <- as.numeric(pupEmerged$daysSurv)

ctrl <- pupEmerged[pupEmerged$name %in% "Control",]
mate <- pupEmerged[pupEmerged$name %in% "Mating delay",]
nuts <- pupEmerged[pupEmerged$name %in% "Nutritional \n stress",]

#*****************************************control*********************************************
fitModRand <- function(fixedEffects,randomEffects,dat,whichAIC){
  lmeCtrl <- lmeControl(opt='optim',msMaxIter=100)
  
  
  
  if(is.na(randomEffects)){
    mod <- tryCatch(
      
      lm(fixedEffects,data=dat)
      ,error = function(cond) {
        message("Here's the original error message:")
        message(cond)
        return(NA)
      },
      
      finally = {
        message(paste("Processed"))
      }
      
    )
    
    if(!is.na(mod)){
    coefsci <- confint(mod)
    coefs <- coefficients(mod)
    coefs <- cbind.data.frame(lower=coefsci[,1],est.=coefs,upper=coefsci[,2])
    rand <- c(NA,NA) 
    }else{
      coefs <- c(NA,NA,NA)
      rand <- c(NA,NA)
    }

  }else{
    mod <- tryCatch( lme(fixedEffects
               ,random = randomEffects
               ,data=dat,method="ML",control=lmeCtrl)
               
               ,error = function(cond) {
                 message("Here's the original error message:")
                 message(cond)
                 return(NA)
               },
               
               finally = {
                 message(paste("Processed"))
               }
    )
    
    if(!is.na(mod)){
    coefs <- intervals(mod,which="fixed")
    coefs <- data.frame(coefs$fixed)
    names(coefs) <- c("lower","est.","upper")
    rand <- VarCorr(mod)
    rand <- data.frame(rand[1:length(rand[1,]),1:2])
    }else{
      coefs <- c(NA,NA,NA)
      rand <- c(NA,NA)
    }
  }
  
  if(!is.na(mod)){out <- summary(mod)
  
  ll <- logLik(mod)
  parms <- attributes(ll)$df - 1

  if (whichAIC =="aic"){
    aic <- AIC(mod)
  }else{aic <- AICc(mod)}
  
  }else{
    ll <- NA
    parms <- NA
    aic <- NA
  }
  
  
  return(list(coefs,parms,ll,aic,rand))
}



# formulas incl. random effects
fe1 <- daysSurv ~ wet_weight + sex + mAgeDays+I(mAgeDays^2)+I(mAgeDays^3)
fe1a <- daysSurv ~ wet_weight + mAgeDays+I(mAgeDays^2)+I(mAgeDays^3)
re1 <- ~1+mAgeDays+I(mAgeDays^2)+I(mAgeDays^3)| adults_id
re1a <- ~1| adults_id

fe2 <-  daysSurv ~ wet_weight + sex + mAgeDays+I(mAgeDays^2)
fe2a <- daysSurv ~ wet_weight + mAgeDays+I(mAgeDays^2)
re2 <-  ~1+mAgeDays+I(mAgeDays^2)| adults_id
re2a <- ~1| adults_id

fe3 <- daysSurv ~ wet_weight + sex + log(mAgeDays)
fe3a <-  daysSurv ~ wet_weight + log(mAgeDays)
re3 <- ~1+log(mAgeDays)| adults_id
re3a <- ~1| adults_id

fe4 <- daysSurv ~ wet_weight + sex + mAgeDays
fe4a <- daysSurv ~ wet_weight + mAgeDays
re4 <- ~1+mAgeDays| adults_id
re4a <- ~1| adults_id

fe5 <- daysSurv ~ wet_weight + sex
re5 <- ~1| adults_id

fe6 <- daysSurv ~ wet_weight
re6 <- ~1| adults_id

fe7 <- daysSurv ~ 1
re7 <- ~1| adults_id

combinations <- list(c(fe1,re1)
                     ,c(fe1a,re1)
                     ,c(fe1,re1a)
                     ,c(fe1a,re1a)
                     ,c(fe1,NA)
                     ,c(fe1a,NA)
                     ,c(fe2,re2)
                     ,c(fe2a,re2)
                     ,c(fe2,re2a)
                     ,c(fe2a,re2a)
                     ,c(fe2,NA)
                     ,c(fe2a,NA)
                     ,c(fe3,re3)
                     ,c(fe3a,re3)
                     ,c(fe3,re3a)
                     ,c(fe3a,re3a)
                     ,c(fe3,NA)
                     ,c(fe3a,NA)
                     ,c(fe4,re4)
                     ,c(fe4a,re4)
                     ,c(fe4,re4a)
                     ,c(fe4a,re4a)
                     ,c(fe4,NA)
                     ,c(fe4a,NA)
                     ,c(fe5,re5)
                     ,c(fe5,NA)
                     ,c(fe6,re6)
                     ,c(fe6,NA)
                     ,c(fe7,re7)
                     ,c(fe7,NA))


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
  
  
  modelSummary$deltaAIC <- round(modelSummary$aic - min(modelSummary$aic,na.rm=T),3)
  modelSummary$weights <- NA
  modelSummary$weights[!modelSummary$aic %in% NA] <- round(Weights(modelSummary$aic[!modelSummary$aic %in% NA]),3)
  row.names(modelSummary) <- c()
  saveRDS(modelSummary, file = name)
  return(modelSummary)
  
}

coefSummaryFunc<- function(modelFitOutputs,name= "coefsWetWeightControl.rds",mods=modelsC){
  
  outPuts <- lapply(modelFitOutputs,'[[',3)
  aics <- as.numeric(lapply(outPuts,'[[',4))
  cfs <- lapply(outPuts,'[[',1)
  cfs2 <- lapply(1:length(cfs),function(x){
    temp <- data.frame(cfs[[x]])
    if(!is.na(temp)){
    temp$modelNumber <- x
    temp$parameter <- row.names(temp)
    row.names(temp) <- c()
    temp$aic <- aics[[x]]
    return(temp[,c(4,5,1,2,3,6)])
    }
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

#*********************************CONTROL*************************


fitModsC <- lapply(1:length(combinations),function(x){
  temp <- combinations[[x]]
  fe <- as.formula(as.character(temp[1]))
  re <- NA
  if(!is.na(temp[2])){re <- as.formula(as.character(temp[2]))}
  output <- fitModRand(fixedEffects=fe,randomEffects=re,dat=ctrl,whichAIC="aicc")
  return(list(temp[1],temp[2],output))
})


modelsC <- modelSummaryFunc(fitModsC,name="modelSummaryStarvControl.rds")

coefSummaryFunc(modelFitOutputs=fitModsC,name="coefsStarvControl.rds",mods=modelsC)

#*********************************MATE DELAY*************************


fitModsM <- lapply(1:length(combinations),function(x){
  temp <- combinations[[x]]
  fe <- as.formula(as.character(temp[1]))
  re <- NA
  if(!is.na(temp[2])){re <- as.formula(as.character(temp[2]))}
  output <- fitModRand(fixedEffects=fe,randomEffects=re,dat=mate,whichAIC="aicc")
  return(list(temp[1],temp[2],output))
})


modelsM <- modelSummaryFunc(fitModsM,name="modelSummaryStarvMate.rds")

coefSummaryFunc(modelFitOutputs=fitModsM,name="coefsStarvMate.rds",mods=modelsM)

#*********************************NUTRITIONAL STRESS*************************

fitModsN <- lapply(1:length(combinations),function(x){
  temp <- combinations[[x]]
  fe <- as.formula(as.character(temp[1]))
  re <- NA
  if(!is.na(temp[2])){re <- as.formula(as.character(temp[2]))}
  output <- fitModRand(fixedEffects=fe,randomEffects=re,dat=nuts,whichAIC="aicc")
  return(list(temp[1],temp[2],output))
})


modelsN <- modelSummaryFunc(fitModsN,name="modelSummaryStarvNuts.rds")

coefSummaryFunc(modelFitOutputs=fitModsN,name="coefsStarvNuts.rds",mods=modelsN)



#*****************Predict effect of age*******************
 lmeCtrl <- lmeControl(opt='optim',msMaxIter=100)

ctrlMod <- lm(fe2 ,data=ctrl)
mateMod <- lm(fe2, data=mate)
nutsMod <- lmer(daysSurv ~ wet_weight + sex + mAgeDays + I(mAgeDays^2) + I(mAgeDays^3) + (1  | adults_id)
               ,data=nuts,REML=F)

#*************Control**********
predFuncC <- function(weightDat=as.numeric(weights[1])) {
  ctrlAgedat <- cbind.data.frame(wet_weight=rep(rep(weightDat,100),2),mAgeDays=rep(c(1:100),2),sex=c(rep("M",100),rep("F",100)))
  ctrlPred <- data.frame(predict(ctrlMod
                       ,newdata=ctrlAgedat
                       ,interval="confidence"))
  ctrlAgedat$pred <- ctrlPred$fit
  ctrlAgedat$lwr <- ctrlPred$lwr
  ctrlAgedat$upr <- ctrlPred$upr
  return(ctrlAgedat)
}

ctrlAge1 <- predFuncC(as.numeric(weights[1]))
ctrlAge1$q <- "Quartile 1"
ctrlAge2 <- predFuncC(as.numeric(weights[2]))
ctrlAge2$q <- "Quartile 2"
ctrlAge3 <- predFuncC(as.numeric(weights[3]))
ctrlAge3$q <- "Quartile 3"
ctrlAge4 <- predFuncC(as.numeric(weights[4]))
ctrlAge4$q <- "Quartile 4"

ctrlAge <- rbind.data.frame(ctrlAge1,ctrlAge2,ctrlAge3,ctrlAge4)

#*******************mating delay*******************
predFuncM <- function(weightDat=as.numeric(weights[1])) {
  mateAgedat <- cbind.data.frame(wet_weight=rep(rep(weightDat,100),2),mAgeDays=rep(c(1:100),2),sex=c(rep("M",100),rep("F",100)))
  matePred <- data.frame(predict(mateMod
                                 ,newdata=mateAgedat
                                 ,interval="confidence"))
  mateAgedat$pred <- matePred$fit
  mateAgedat$lwr <- matePred$lwr
  mateAgedat$upr <- matePred$upr
  return(mateAgedat)
}

mateAge1 <- predFuncM(as.numeric(weights[1]))
mateAge1$q <- "Quartile 1"
mateAge2 <- predFuncM(as.numeric(weights[2]))
mateAge2$q <- "Quartile 2"
mateAge3 <- predFuncM(as.numeric(weights[3]))
mateAge3$q <- "Quartile 3"
mateAge4 <- predFuncM(as.numeric(weights[4]))
mateAge4$q <- "Quartile 4"

mateAge <- rbind.data.frame(mateAge1,mateAge2,mateAge3,mateAge4)

#****************Nutritional stress************************
#*
predFuncN <- function(weightDat=as.numeric(weights[1])) {
  nutsAgedat <- cbind.data.frame(wet_weight=rep(rep(weightDat,100),2),mAgeDays=rep(c(1:100),2),sex=c(rep("M",100),rep("F",100)),adults_id=rep(rep(113,100),2))
  nutsPred <- data.frame(predict(nutsMod
                                 ,newdata=nutsAgedat))
  nutsAgedat$pred <- nutsPred[,1]
  merBootN <- bootMer(nutsMod, function(x) predict(x, newdata = nutsAgedat), nsim = 100)
  
  nutsAgedat$lwr <- apply(merBootN$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE)))
  nutsAgedat$upr <- apply(merBootN$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))
  
  return(nutsAgedat)
}

nutsAge1 <- predFuncN(as.numeric(weights[1]))
nutsAge1$q <- "Quartile 1"
nutsAge2 <- predFuncN(as.numeric(weights[2]))
nutsAge2$q <- "Quartile 2"
nutsAge3 <- predFuncN(as.numeric(weights[3]))
nutsAge3$q <- "Quartile 3"
nutsAge4 <- predFuncN(as.numeric(weights[4]))
nutsAge4$q <- "Quartile 4"

nutsAge <- rbind.data.frame(nutsAge1,nutsAge2,nutsAge3,nutsAge4)

ctrlAge$name <- "Control"
mateAge$name <- "Mating delay"
nutsAge$name <- "Nutrtional stress"

nutsAge <- nutsAge[,-4]

ageEffect <- rbind.data.frame(ctrlAge,mateAge,nutsAge)
ageEffect$name <- as.factor(ageEffect$name)
levels(ageEffect$name) <- c("Control","Mating delay","Nutritional stress")
levels(ageEffect$sex) <- c("F","M")

cols <- c("#cccccc"
          ,"#969696"
          ,"#636363"
          ,"#252525")

survAge.plot <- ggplot(ageEffect) +
  geom_line(aes(x=mAgeDays,y=pred,col=name,linetype=q)) +
  ylim(2.5,12) +
  scale_linetype_manual(values=c(1,2,3,4)) +
  scale_color_manual(values=c("#875777","#eab051","#c0b9ac"),guide=F) +
  ylab("Predicted days surviving") +
  xlab("Maternal age (days)") +
  labs(title="b)"
    #   ,color="Treatment"
       ,linetype="Wet weight") +
  theme_set(theme_bw()) +
  theme(axis.line = element_line(color = 'black')
        ,text=element_text(size=7)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=7)
        ,legend.key.size = unit(0.5,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=5)
        ,legend.margin=margin(0,0,0,0)
        ,legend.box.margin=margin(-10,-2,-10,-10)
        # ,legend.position =c(0.9,0.05)
          ,legend.title = element_text(size=5)
        ,strip.background = element_rect(colour="white", fill="white")
        ,panel.border = element_blank()
        #,strip.text.x = element_blank()
  ) +
  facet_wrap(~name+sex,ncol=2,labeller = label_wrap_gen(multi_line=FALSE))


pupEmerged$sex <- as.factor(pupEmerged$sex)
levels(pupEmerged$sex) <- c("F","M")

tiff("FigRawSexWeight.tiff", height = 3.5, width = 6, units = 'in', compression="lzw", res=400)
ggplot(pupEmerged) +
  geom_point(aes(x=mAgeDays,y=daysSurv,col=wet_weight)
             ,size=0.5) +
  ylab("Offspring survival (days)") +
  xlab("Maternal age (days)") +
  scale_color_continuous(low="#D3D3D3",high="black") +
  labs(colour="Wet weight") +
  theme_set(theme_bw()) +
  theme(axis.line = element_line(color = 'black')
        ,text=element_text(size=7)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=7)
      #  ,legend.position =c(0.95,0.5)
       # ,legend.title = element_text(size=5)
        ,strip.background = element_rect(colour="white", fill="white")
        ,panel.border = element_blank()
  ) +
  facet_wrap(~name+sex,labeller = label_wrap_gen(multi_line=FALSE))
dev.off()


tiff("Fig4_offspringSurv.tiff", height = 4, width =5.5, units = 'in', compression="lzw", res=400)
grid.arrange(weightSexPlot,survAge.plot,ncol=2,widths = 1:2)
dev.off()


#***************nutritional stress with alive females only***********
deadNuts
nuts <- nuts[!nuts$adults_id %in% deadNuts,]


fitModsN <- lapply(1:length(combinations),function(x){
  temp <- combinations[[x]]
  fe <- as.formula(as.character(temp[1]))
  re <- NA
  if(!is.na(temp[2])){re <- as.formula(as.character(temp[2]))}
  output <- fitModRand(fixedEffects=fe,randomEffects=re,dat=nuts,whichAIC="aicc")
  return(list(temp[1],temp[2],output))
})


modelsN <- modelSummaryFunc(fitModsN,name="modelSummaryStarvNutsAlive.rds")

coefSummaryFunc(modelFitOutputs=fitModsN,name="coefsStarvNutsAlive.rds",mods=modelsN)

