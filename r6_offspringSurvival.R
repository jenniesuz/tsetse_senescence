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

levels(pupEmerged$name) <- c("Control","Mating delay","Nutritional stress")


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


tiff("fig_offspringWeightSex.tiff", height = 5, width = 4, units = 'in', compression="lzw", res=400)
ggplot(pupEmerged) +
  geom_boxplot(aes(y=wet_weight,x=sex)) +
  ylab("Wet weight") +
  xlab("Sex") +
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


#*****************************************control*********************************************
ctrl <- pupEmerged[pupEmerged$name %in% "Control",]
#*******************random effects***********
y1 <- lme(as.numeric(daysSurv) ~ wet_weight + sex + mAgeDays + I(mAgeDays^2), random = ~ 1 + mAgeDays + I(mAgeDays^2)| adults_id,data=ctrl,method="ML")
y2 <- lme(as.numeric(daysSurv) ~ wet_weight + sex + mAgeDays + I(mAgeDays^2), random = ~ 1 | adults_id,data=ctrl,method="ML")
y3 <- lm(as.numeric(daysSurv) ~ wet_weight + sex + mAgeDays + I(mAgeDays^2), data=ctrl)
AICc(y1) 
AICc(y2)
AICc(y3)
aictab(list(y2))
aictab(list(y3))
Weights(c(AICc(y2),AICc(y3)))
#****************fixed effects***************
y4 <- lm(as.numeric(daysSurv) ~ wet_weight + mAgeDays + I(mAgeDays^2), data=ctrl)
y5 <- lm(as.numeric(daysSurv) ~ wet_weight + sex + mAgeDays, data=ctrl)
y6 <- lm(as.numeric(daysSurv) ~ wet_weight + mAgeDays, data=ctrl)
y7 <- lm(as.numeric(daysSurv) ~ wet_weight + sex, data=ctrl)
y8 <- lm(as.numeric(daysSurv) ~ wet_weight, data=ctrl)
y9 <- lm(as.numeric(daysSurv) ~ 1,data=ctrl)
AICc(y3) 
AICc(y4)
AICc(y5) 
AICc(y6)
AICc(y7)
AICc(y8)
AICc(y9)
aictab(list(y3,y4,y5,y6,y7,y8,y9))
Weights(c(AICc(y3),AICc(y4),AICc(y5),AICc(y6),AICc(y7),AICc(y8),AICc(y9)))
#***********************************************
summary(y3)
#car::vif(y4)

#*****************Predict effect of age*******************
predFuncC <- function(weightDat=as.numeric(weights[1])) {
  ctrlAgedat <- cbind.data.frame(wet_weight=rep(rep(weightDat,100),2),mAgeDays=rep(c(1:100),2),sex=c(rep("m",100),rep("f",100)))
  ctrlPred <- data.frame(predict(y3
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



#****************************************mating delay*************************************
mate <- pupEmerged[pupEmerged$name %in% "Mating delay",]
#*******************random effects**************
y1 <- lme(as.numeric(daysSurv) ~ wet_weight + sex + mAgeDays + I(mAgeDays^2), random = ~ 1 + mAgeDays + I(mAgeDays^2)| adults_id,data=mate,method="ML")
y2 <- lme(as.numeric(daysSurv) ~ wet_weight + sex + mAgeDays + I(mAgeDays^2), random = ~ 1 | adults_id,data=mate,method="ML")
y3 <- lm(as.numeric(daysSurv) ~ wet_weight + sex + mAgeDays + I(mAgeDays^2), data=mate)
AICc(y2) 
AICc(y3)
aictab(list(y2))
aictab(list(y3))
Weights(c(AICc(y2),AICc(y3)))
#*******************fixed effects******************
y4 <- lm(as.numeric(daysSurv) ~ wet_weight + mAgeDays + I(mAgeDays^2), data=mate)
y5 <- lm(as.numeric(daysSurv) ~ wet_weight + sex + mAgeDays, data=mate)
y6 <- lm(as.numeric(daysSurv) ~ wet_weight + sex ,data=mate)
y7 <- lm(as.numeric(daysSurv) ~ wet_weight + mAgeDays,data=mate)
y8 <- lm(as.numeric(daysSurv) ~ wet_weight,data=mate)
y9 <- lm(as.numeric(daysSurv) ~ 1,data=mate)
AICc(y3) 
AICc(y4) 
AICc(y5) 
AICc(y6) 
AICc(y7)
AICc(y8)
AICc(y9)

aictab(list(y3,y4,y5,y6,y7,y8,y9))
Weights(c(AICc(y3),AICc(y4),AICc(y5),AICc(y6),AICc(y7),AICc(y8),AICc(y9)))
summary(y3)

#************************predict************************
#***************age**************

predFuncM <- function(weightDat=as.numeric(weights[1])) {
  mateAgedat <- cbind.data.frame(wet_weight=rep(rep(weightDat,100),2),mAgeDays=rep(c(1:100),2),sex=c(rep("m",100),rep("f",100)))
  matePred <- data.frame(predict(y3
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

#***********************************nutritional stress***************************
nuts <- pupEmerged[pupEmerged$name %in% "Nutritional stress",]
#*********************random effects***********************
y1 <- lmer(as.numeric(daysSurv) ~ wet_weight + sex + mAgeDays + I(mAgeDays^2) + (1 + mAgeDays + I(mAgeDays^2)| adults_id),data=nuts,REML=F)
y2 <- lmer(as.numeric(daysSurv) ~ wet_weight + sex + mAgeDays + I(mAgeDays^2) + (1 | adults_id),data=nuts,REML=F)
y3 <- lm(as.numeric(daysSurv) ~ wet_weight + sex + mAgeDays + I(mAgeDays^2) , data=nuts)
AICc(y2) 
AICc(y3) 
aictab(list(y2))
aictab(list(y3))

summary(y2)
Weights(c(AICc(y2),AICc(y3)))
#********************fixed effects****************
y4 <- lmer(as.numeric(daysSurv) ~ wet_weight  + mAgeDays + I(mAgeDays^2) + ( 1 | adults_id),data=nuts,REML=F)
y5 <- lmer(as.numeric(daysSurv) ~ wet_weight + sex + mAgeDays + ( 1 | adults_id),data=nuts,REML=F)
y6 <- lmer(as.numeric(daysSurv) ~ wet_weight + sex + ( 1 | adults_id),data=nuts,REML=F)
y7 <- lmer(as.numeric(daysSurv) ~ wet_weight + mAgeDays + ( 1 | adults_id),data=nuts,REML=F)
y8 <- lmer(as.numeric(daysSurv) ~ wet_weight + ( 1 | adults_id),data=nuts,REML=F)
y9 <- lmer(as.numeric(daysSurv) ~ 1 + ( 1 | adults_id),data=nuts,REML=F)
aictab(list(y2,y4,y5,y6,y7,y8,y9))

Weights(c(AICc(y2),AICc(y4),AICc(y5),AICc(y6),AICc(y7),AICc(y8),AICc(y9)))

summary(y2)


#***********************predict**************
#***age****

predFuncN <- function(weightDat=as.numeric(weights[1])) {
  nutsAgedat <- cbind.data.frame(wet_weight=rep(rep(weightDat,100),2),mAgeDays=rep(c(1:100),2),sex=c(rep("m",100),rep("f",100)),adults_id=rep(rep(113,100),2))
  nutsPred <- data.frame(predict(y2
                                 ,newdata=nutsAgedat))
  nutsAgedat$pred <- nutsPred[,1]
  merBootN <- bootMer(y2, function(x) predict(x, newdata = nutsAgedat), nsim = 100)
  
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
levels(ageEffect$sex) <- c("Female","Male")

cols <- c("#cccccc"
          ,"#969696"
          ,"#636363"
          ,"#252525")

survAge.plot <- ggplot(ageEffect) +
  geom_line(aes(x=mAgeDays,y=pred,col=q,linetype=q)) +
  ylim(2.5,12) +
  scale_linetype_manual("", values=c(1,2,3,4)) +
  scale_color_manual("",values=cols) +
  ylab("Predicted days surviving") +
  xlab("Mother age (days)") +
  theme_set(theme_bw()) +
  theme(axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=8)
        ,legend.key.size = unit(0.5,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=6)
        ,legend.position =c(0.9,0.05)
        ,legend.title = element_blank()
        ,strip.background = element_rect(colour="white", fill="white")
        ,panel.border = element_blank()
        #,strip.text.x = element_blank()
  ) +
  facet_wrap(~name+sex,ncol=2,labeller = label_wrap_gen(multi_line=FALSE))


survAge.plot

tiff("Fig5_offspringSurv.tiff", height = 6, width = 4.5, units = 'in', compression="lzw", res=400)
survAge.plot
dev.off()

pupEmerged$sex <- as.factor(pupEmerged$sex)
levels(pupEmerged$sex) <- c("Female","Male")

tiff("FigRawSurv.tiff", height = 5, width = 6, units = 'in', compression="lzw", res=400)
ggplot(pupEmerged) +
  geom_point(aes(x=mAgeDays,y=daysSurv,col=wet_weight)
             ,size=0.5) +
  ylab("Offspring survival (days)") +
  xlab("Mother age (days)") +
  scale_color_continuous(low="#D3D3D3",high="black") +
  labs(colour="Wet weight") +
  theme_set(theme_bw()) +
  theme(axis.line = element_line(color = 'black')
        ,text=element_text(size=8)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=7)
      #  ,legend.position =c(0.95,0.5)
       # ,legend.title = element_blank()
        ,strip.background = element_rect(colour="white", fill="white")
        ,panel.border = element_blank()
  ) +
  facet_wrap(~name+sex,labeller = label_wrap_gen(multi_line=FALSE))
dev.off()

