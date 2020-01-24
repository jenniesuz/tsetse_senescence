# This script analyses the abortion data and produces Figure 3 of the manuscript
source("r1_queries.R")

library(ggplot2)
library(binom)
library(nlme)
library(geepack)
library(lme4)
library(AICcmodavg)
library(MuMIn)

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


#****************************Models****************************
#************control group*************
ctrl <- mother.larvipositions[mother.larvipositions$name %in% "Control",]
length(ctrl$name) # 629

cMod1 <- glmer(abortion ~ mAge  + (1+mAge| adults_id)
               ,family=binomial
               ,data=ctrl)  # singular fit

cMod2 <- glmer(abortion ~ mAge + (1| adults_id)
             ,family=binomial
             ,data=ctrl) 

cMod3 <- glm(abortion ~ mAge 
             ,family=binomial
             ,data=ctrl)

cMod4 <- glm(abortion ~ 1
             ,family=binomial
             ,data=ctrl)

AIC(cMod2) 
AIC(cMod3) 
AIC(cMod4) 

aictab(list(cMod2), c("cMod2"))
aictab(list(cMod3,cMod4), c("cMod3","cMod4"))

Weights(c(AIC(cMod2),AIC(cMod3),AIC(cMod4)))

summary(cMod3)
#*********************************************


#**********nutritional stress group******
nuts <- mother.larvipositions[mother.larvipositions$name %in% "Nutritional stress",]

nMod1 <- glmer(abortion ~ mAge  + (1+mAge| adults_id)
               ,family=binomial
               ,data=nuts)  # singular fit

nMod2 <- glmer(abortion ~ mAge  + (1| adults_id)
               ,family=binomial
               ,data=nuts) 

nMod3 <- glm(abortion ~ mAge 
             ,family=binomial
             ,data=nuts)

nMod4 <- glm(abortion ~ 1
             ,family=binomial
             ,data=nuts)

AIC(nMod1) 
AIC(nMod2) 
AIC(nMod3) 
AIC(nMod4) 
aictab(list(nMod2))
aictab(list(nMod3,nMod4), c("nMod3","nMod4"))

Weights(c(AIC(nMod2),AIC(nMod3),AIC(nMod4)))

summary(nMod2)
#***************************************************************

#*********************mating delay group**********************
mate <- mother.larvipositions[mother.larvipositions$name %in% "Mating delay",]

mMod1 <- glmer(abortion ~ mAge + (1+mAge| adults_id)
               ,family=binomial
               ,data=mate)  # singular fit

mMod2 <- glmer(abortion ~ mAge + (1| adults_id)
               ,family=binomial
               ,data=mate) 

mMod3 <- glm(abortion ~ mAge 
             ,family=binomial
             ,data=mate)


mMod4 <- glm(abortion ~ 1
             ,family=binomial
             ,data=mate)

AIC(mMod2) 
AIC(mMod3) 
AIC(mMod4) 
aictab(list(mMod2))
aictab(list(mMod3,mMod4), c("mMod3","mMod4"))

Weights(c(AIC(mMod2),AIC(mMod3),AIC(mMod4)))

summary(mMod3)
#****************************************************************

#**************************Predict****************

#****************Control*******************
ilink <- family(cMod3)$linkinv
ctrlP <- with(ctrl,
                data.frame(mAge = seq(min(mAge), max(mAge),
                                            length = 100)))
ctrlP <- cbind(ctrlP, predict(cMod3, ctrlP, type = "link", se.fit = TRUE)[1:2])
ctrlP <- transform(ctrlP, Fitted = ilink(fit), Upper = ilink(fit + (1.96 * se.fit)),
                     Lower = ilink(fit - (1.96 * se.fit)))
ctrlP$name <- "Control"
#************Mating delay******************
ilink <- family(mMod3)$linkinv
mateP <- with(mate,
              data.frame(mAge = seq(min(mAge), max(mAge),
                                    length = 100)))
mateP <- cbind(mateP, predict(mMod3, mateP, type = "link", se.fit = TRUE)[1:2])
mateP <- transform(mateP, Fitted = ilink(fit), Upper = ilink(fit + (1.96 * se.fit)),
                   Lower = ilink(fit - (1.96 * se.fit)))
mateP$name <- "Mating delay"
#***********Nutritional stress************
ilink <- family(nMod2)$linkinv
nutsP <- with(nuts,
              data.frame(mAge = seq(min(mAge), max(mAge),
                                    length = 100)))
nutsP$Fitted <- predict(nMod2,re.form=NA,newdat=nutsP,type="response")
nutsP$fit <- as.numeric(nutsP$Fitted)
merBootN <- bootMer(nMod2, function(x) predict(x, newdata = nutsP,re.form=NA,type="response"), nsim = 100, re.form = NA)

nutsP$Lower <- apply(merBootN$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE)))
nutsP$Upper <- apply(merBootN$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))


nutsP$name <- "Nutritional stress"


newdat <- rbind.data.frame(ctrlP[,c("mAge","Fitted","Upper","Lower","name")],mateP[,c("mAge","Fitted","Upper","Lower","name")],nutsP[,c("mAge","Fitted","Upper","Lower","name")])

newdat$name <- as.factor(newdat$name)
levels(newdat$name) <- c("Control","Mating delay","Nutritional stress")

#**********************Plot*******************************************

tiff("Fig3_abortions.tiff", height = 3, width = 6, units = 'in', compression="lzw", res=400)
ggplot(mother.larvipositions,aes(as.numeric(mAge),abortion)) +
  scale_y_continuous( breaks=c(0,0.25,0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1)) +
  geom_jitter(width = 0.1, height = 0.1,col="darkgrey",size=0.5) +
  geom_ribbon(data = newdat, aes(ymin = Lower, ymax = Upper, x = mAge),
              fill = "grey", alpha = 0.5, inherit.aes = FALSE) +
 geom_line(data=newdat,aes(as.numeric(mAge),Fitted)) +
   labs(  x="Mother age (days)"
       ,y="Probability of abortion"
       #,title="b)"
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
  facet_wrap(~name,ncol=3)
dev.off()


