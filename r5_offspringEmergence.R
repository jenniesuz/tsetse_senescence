
# This script analyses the probabiltiy of offspring emergence and produces figure 5a of the manuscript
library(plyr)         
library(reshape)
library(ggplot2)
library(binom)
library(gridExtra)
library(lme4)
library(AICcmodavg)
library(MuMIn)



source("r1_queries.R")    # query database 


#*********************************************Data********************************************
pupae$date_emerged <- as.Date(pupae$date_emerged,format=c("%d/%m/%Y"))
pupae$date_of_death <- as.Date(pupae$date_of_death,format=c("%d/%m/%Y"))
pupae$date_of_emergence_min <- as.Date(pupae$date_of_emergence_min, format="%d/%m/%Y")  
pupae$larviposition_date <- as.Date(pupae$larviposition_date, format="%d/%m/%Y")        
pupae$mAgeDays <- pupae$larviposition_date - pupae$date_of_emergence_min
pupae$mAgeDays <- as.numeric(pupae$mAgeDays)

pupae <- pupae[pupae$abortion %in% "0.0",]
pupae <- pupae[!pupae$treatment %in% "fat analysis",]

pupae$emerged <- as.integer(pupae$emerged)

pupae <- pupae[!pupae$wet_weight %in% NA,]

emergedByTreat <- ddply(pupae,.(name),summarise,em=sum(emerged),totem=length(emerged))
#***********************************************************************************************

ctrl <- pupae[pupae$name %in% "control",]
mate <- pupae[pupae$name %in% "mate_delay",]
nuts <- pupae[pupae$name %in% "nutrition",]

sum(1-ctrl$emerged) #15
sum(1-mate$emerged) #14
sum(1-nuts$emerged) #20


#********************************************Analysis***************************************

fitModFunc <- function(form,modType,dat,whichAIC){
  #glmeCtrl <- glmerControl(optimizer='optim')
  
  mod <- tryCatch( 
  
  if(modType=="me"){
    glmer(form
                 ,family=binomial
                 ,data=dat)
  }else{
    glm(form
               ,family=binomial
               ,data=dat)
  }
  ,
  warning = function(w) {
    message("Here's the original error message:")
    message(w)
    return(NA)
  },
  
  error = function(cond) {
    message("Here's the original error message:")
    message(cond)
    return(NA)
  },
  
  finally = {
    message(paste("Processed"))
  }
  
  )
  
  sing<-F
  
  if(is.na(mod)==F){
    
    if(modType=="me"){
       if(isSingular(mod)==T){sing<-T}
    }
      
      if(sing == F){
        out <- summary(mod)

        coefs <- coefficients(out)
        
     #   coefsci <-   confint(mod)
       
        ll <- logLik(mod)
         parms <- attributes(ll)$df 
  
       if (whichAIC =="aic"){
         aic <- AIC(mod)
       }else{aic <- AICc(mod)}
  
       # return(list(coefs,coefsci,parms,ll,aic))
         return(list(coefs,parms,ll,aic))
         
    }else{
        return(list(-999,-999,-999,-999,-999))
    }
  }else{return(list(NA,NA,NA,NA,NA))}
  
  }
  
modelSummaryFunc <- function(modelFitOutputs,form){
  
  model <-form
  
  k <- modelFitOutputs[[3]]
  ll <- as.numeric(modelFitOutputs[[4]])
  aic <- modelFitOutputs[[5]]
  
  modelSummary <- c(model,k,ll,aic)
  return(modelSummary)
  
}




m1 <- emerged ~ wet_weight + mAgeDays + I(mAgeDays^2) + I(mAgeDays^3) + (1 + wet_weight  |adults_id ) # incl. additional random effects likely unidentifiable
m2 <- emerged ~ wet_weight + mAgeDays + I(mAgeDays^2) + I(mAgeDays^3) + (1 |adults_id )
m3 <- emerged ~ wet_weight + mAgeDays + I(mAgeDays^2) + I(mAgeDays^3) 

m4 <- emerged ~ wet_weight + mAgeDays + I(mAgeDays^2) + (1 + wet_weight |adults_id )
m5 <- emerged ~ wet_weight + mAgeDays + I(mAgeDays^2) + (1 |adults_id )
m6 <- emerged ~ wet_weight + mAgeDays + I(mAgeDays^2)

m7 <- emerged ~ wet_weight + log(mAgeDays) + (1 + wet_weight |adults_id )
m8 <- emerged ~ wet_weight + log(mAgeDays) + (1 |adults_id )
m9 <- emerged ~ wet_weight + log(mAgeDays) 

m10 <- emerged ~ wet_weight + mAgeDays + (1 + wet_weight |adults_id )
m11 <- emerged ~ wet_weight + mAgeDays + (1 |adults_id )
m12 <- emerged ~ wet_weight + mAgeDays 

m13 <- emerged ~ wet_weight + (1 + wet_weight |adults_id )
m14 <- emerged ~ wet_weight + (1 |adults_id )
m15 <- emerged ~ wet_weight 
m16 <- emerged ~ 1

mT <- c("me","me","fe","me","me","fe","me","me","fe","me","me","fe","me","me","fe","fe")
  

combinations <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16)
#************************Control treatment*************************
fitModsC <- lapply(1:length(combinations),function(x){
  form <- formula(as.character(combinations[x]))
  modelFit <- fitModFunc(form=form,modType=mT[x],dat=ctrl,whichAIC="aicc")
  return(modelFit)
})


modSum <- lapply(1:length(combinations),function(x){
  f <- combinations[x]
  modO <- fitModsC[[x]]
  modSum <- modelSummaryFunc(modO,f)
  return(modSum)
})

model <- as.character(sapply(modSum,'[[',1))
k <- sapply(modSum,'[[',2)
ll <- sapply(modSum,'[[',3)
aic <- sapply(modSum,'[[',4)

modelSummary <- cbind.data.frame(model,k,ll,aic)
modelSummary$k[modelSummary$k %in% -999] <- NA
modelSummary$ll[modelSummary$ll %in% -999] <- NA
modelSummary$comment <- NA
modelSummary$comment[modelSummary$aic %in% -999] <- "Singular fit"
modelSummary$comment[modelSummary$aic %in% NA] <- "Convergence failure"
modelSummary$aic[modelSummary$aic %in% -999] <- NA 

modelSummary <- modelSummary[order(modelSummary$aic),]

modelSummary$deltaAIC <- round(modelSummary$aic - min(modelSummary$aic,na.rm=T),3)
modelSummary$weights[is.na(modelSummary$aic)==F] <- round(Weights(modelSummary$aic[is.na(modelSummary$aic)==F]),3)
saveRDS(modelSummary, file = "modelSummaryEmergenceControl.rds")

#***********************Mating delay***********************
fitModsM <- lapply(1:length(combinations),function(x){
  form <- formula(as.character(combinations[x]))
  modelFit <- fitModFunc(form=form,modType=mT[x],dat=mate,whichAIC="aicc")
  return(modelFit)
})


modSum <- lapply(1:length(combinations),function(x){
  f <- combinations[x]
  modO <- fitModsM[[x]]
  modSum <- modelSummaryFunc(modO,f)
  return(modSum)
})

model <- as.character(sapply(modSum,'[[',1))
k <- sapply(modSum,'[[',2)
ll <- sapply(modSum,'[[',3)
aic <- sapply(modSum,'[[',4)

modelSummary <- cbind.data.frame(model,k,ll,aic)
modelSummary$k[modelSummary$k %in% -999] <- NA
modelSummary$ll[modelSummary$ll %in% -999] <- NA
modelSummary$comment <- NA
modelSummary$comment[modelSummary$aic %in% -999] <- "Singular fit"
modelSummary$comment[modelSummary$aic %in% NA] <- "Convergence failure"
modelSummary$aic[modelSummary$aic %in% -999] <- NA 

modelSummary <- modelSummary[order(modelSummary$aic),]

modelSummary$deltaAIC <- round(modelSummary$aic - min(modelSummary$aic,na.rm=T),3)
modelSummary$weights[is.na(modelSummary$aic)==F] <- round(Weights(modelSummary$aic[is.na(modelSummary$aic)==F]),3)
saveRDS(modelSummary, file = "modelSummaryEmergenceMate.rds")



#**********************Nutritional stress***********************
nuts <- pupae[pupae$name %in% "nutrition",]
nuts$notEmerged <- 1-nuts$emerged

# incl. additional random effects likely unidentifiable


nuts$mAgeDays <- scale(nuts$mAgeDays)
nuts$wet_weight <- scale(nuts$wet_weight)

fitModsN <- lapply(1:length(combinations),function(x){
  form <- formula(as.character(combinations[x]))
  modelFit <- fitModFunc(form=form,modType=mT[x],dat=nuts,whichAIC="aicc")
  return(modelFit)
})


library("brglm2")
library("blme")

test <- glmer(m1,family="binomial",data=nuts) # problem is complete separation
# many of the individual adults had all their offspring emerge

summary(test2 <- glm(emerged ~ wet_weight + mAgeDays + I(mAgeDays^2) + I(mAgeDays^3),data=nuts,family="binomial"))


test <- glmer(notEmerged ~ wet_weight + mAgeDays + I(mAgeDays^2) + I(mAgeDays^3) + (1 + wet_weight  |adults_id ) 
    , data =nuts, family = binomial, method="detect_separation")



cmod_blme_L2 <- bglmer(m1,family="binomial",data=nuts,
                       fixef.prior = normal(cov = diag(9,4)))





modSum <- lapply(1:length(combinations),function(x){
  f <- combinations[x]
  modO <- fitModsN[[x]]
  modSum <- modelSummaryFunc(modO,f)
  return(modSum)
})

model <- as.character(sapply(modSum,'[[',1))
k <- sapply(modSum,'[[',2)
ll <- sapply(modSum,'[[',3)
aic <- sapply(modSum,'[[',4)

modelSummary <- cbind.data.frame(model,k,ll,aic)
modelSummary$k[modelSummary$k %in% -999] <- NA
modelSummary$ll[modelSummary$ll %in% -999] <- NA
modelSummary$comment <- NA
modelSummary$comment[modelSummary$aic %in% -999] <- "Singular fit"
modelSummary$comment[modelSummary$aic %in% NA] <- "Convergence failure"
modelSummary$aic[modelSummary$aic %in% -999] <- NA 

modelSummary <- modelSummary[order(modelSummary$aic),]

modelSummary$deltaAIC <- round(modelSummary$aic - min(modelSummary$aic,na.rm=T),3)
modelSummary$weights[is.na(modelSummary$aic)==F] <- round(Weights(modelSummary$aic[is.na(modelSummary$aic)==F]),3)
saveRDS(modelSummary, file = "modelSummaryEmergenceNuts.rds")



mod <- tryCatch( 
  
  glmer(emerged ~ wet_weight + mAgeDays + I(mAgeDays^2) + (1 | adults_id)
        ,data=nuts
        ,family="binomial")
  
  ,
  warning = function(w) {
    message("Here's the original error message:")
    message(w)
    return(NA)
  },
  
  error = function(cond) {
    message("Here's the original error message:")
    message(cond)
    return(NA)
  },
  
  finally = {
    message(paste("Processed"))
  }
  
)



car::vif(nuts2)






#************predict***********
weights <- quantile(nuts$wet_weight)

ilink <- family(nuts2)$linkinv


predictFunc <- function(weightDat=as.numeric(weights[1])) {
  ageP <- with(nuts,
               data.frame(mAgeDays = seq(min(mAgeDays), max(mAgeDays),
                                         length = 100)
                          ,wet_weight=rep(weightDat,100)))
  ageP <- cbind(ageP, predict(nuts2, ageP, type = "link", se.fit = TRUE)[1:2])
  ageP <- transform(ageP, Fitted = ilink(fit), Upper = ilink(fit + (1.96 * se.fit)),
                    Lower = ilink(fit - (1.96 * se.fit)))
  
  return(ageP)
}

ageP1 <- predictFunc(as.numeric(weights[1]))
ageP1$q <- "Quartile 1"
ageP2 <- predictFunc(as.numeric(weights[2]))
ageP2$q <- "Quartile 2"
ageP3 <- predictFunc(as.numeric(weights[3]))
ageP3$q <- "Quartile 3"
ageP4 <- predictFunc(as.numeric(weights[4]))
ageP4$q <- "Quartile 4"

ageP <- rbind.data.frame(ageP1,ageP2,ageP3,ageP4)


cols <- c("#cccccc"
          ,"#969696"
          ,"#636363"
          ,"#252525")

agePlot <- ggplot(ageP) +
  geom_line(aes(x=mAgeDays,y=Fitted,col=q,linetype=q)) +
  ylim(0.2,1) +
  ylab("Predicted probability of emergence") +
  xlab("Mother age (days)") +
  #labs(col="Wet weight quartile") +
  scale_linetype_manual("", values=c(1,2,3,4)) +
  scale_color_manual("",values=cols) +
  theme_set(theme_bw()) +
  theme(axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=8)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=9)
        ,legend.position =c(0.1,0.2)
        ,strip.background = element_rect(colour="white", fill="white")
        ,panel.border = element_blank()
  )


tiff("fig_nutEmergence.tiff", height = 3, width = 6, units = 'in', compression="lzw", res=400)
agePlot
dev.off()
