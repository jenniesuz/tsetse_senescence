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



pupWeight <- pupae[pupae$abortion %in% 0,] # take out abortions
ddply(pupWeight,.(name),summarise,l=length(wet_weight))

numWeighed <- ddply(pupWeight,.(adults_id),summarise,num=length(wet_weight))
mean(numWeighed$num) # 4.8
range(numWeighed$num) # 1 - 8

numWeighed <- ddply(pupWeight,.(name,adults_id),summarise,num=length(wet_weight))
ddply(numWeighed,.(name),summarise,mean(num))
ddply(numWeighed,.(name),summarise,range(num))

pupWeight <- pupWeight[!pupWeight$wet_weight%in%NA,]     # remove those - 31 of 1239 for which no wet weight data - either lost or crushed



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


mean(mate$wet_weight,na.rm=T) # 32.57
mean(ctrl$wet_weight,na.rm=T) # 33.54
mean(nuts$wet_weight,na.rm=T) # 28.02




#*******************Control treatment****************
#************************RANDOM EFFECTS***************************
rand3 <- lm(wet_weight ~ mAgeDays+I(mAgeDays^2)+I(mAgeDays^3), data=ctrl)
rand2 <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2)+I(mAgeDays^3), random = ~1 | adults_id, data=ctrl,method="ML")
rand1 <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2)+I(mAgeDays^3), random = ~1+mAgeDays| adults_id, data=ctrl,method="ML")

aictab(list(rand1,rand2))
aictab(list(rand3))

AIC(rand1)
AIC(rand2)
AIC(rand3)

Weights(c(AIC(rand3),AIC(rand2),AIC(rand1)))

#************autocorrelation**************
a1 <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2)+I(mAgeDays^3) , random = ~1+mAgeDays| adults_id, data=ctrl,method="ML")
a2 <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2)+I(mAgeDays^3) , random = ~1+mAgeDays| adults_id, corr=corCompSymm(form =~mAgeDays),data=ctrl,method="ML")
a3 <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2)+I(mAgeDays^3) , random = ~1+mAgeDays| adults_id, corr=corAR1(form=~mAgeDays),data=ctrl,method="ML")
# #*******************************************
aictab(list(a1,a2,a3))
AIC(a1)
AIC(a2)
AIC(a3)
Weights(c(AIC(a1),AIC(a2),AIC(a3)))

#*************************Fixed effects****************
c1 <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2)+I(mAgeDays^3) , random = ~1+mAgeDays| adults_id,data=ctrl,method="ML")
c2 <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2) , random = ~1+mAgeDays| adults_id,data=ctrl,method="ML")
c3 <- lme(wet_weight ~ mAgeDays , random = ~1+mAgeDays| adults_id,data=ctrl,method="ML")
c4 <- lme(wet_weight ~ 1, random = ~1| adults_id,data=ctrl,method="ML")
aictab(list(c1,c2,c3,c4))
AIC(c1)
AIC(c2)
AIC(c3)
AIC(c4)
Weights(c(AIC(c1),AIC(c2),AIC(c3),AIC(c4)))

summary(c2)
modCtrl <- c2

plot(ctrl$wet_weight,residuals(c2))
plot(ctrl$mAgeDays
     ,residuals(c2)
     ,ylim=c(-20,20)
)

#*******************Mating delay treatment****************
#************************RANDOM EFFECTS***************************
rand3 <- lm(wet_weight ~ mAgeDays+I(mAgeDays^2)+I(mAgeDays^3), data=mate)
rand2 <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2)+I(mAgeDays^3), random = ~1 | adults_id, data=mate,method="ML")
rand1 <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2)+I(mAgeDays^3), random = ~1+mAgeDays| adults_id, data=mate,method="ML")

aictab(list(rand1,rand2))
aictab(list(rand3))
Weights(c(AIC(rand3),AIC(rand2),AIC(rand1)))


a1 <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2)+I(mAgeDays^3) , random = ~1+mAgeDays| adults_id, data=mate,method="ML")
a2 <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2)+I(mAgeDays^3) , random = ~1+mAgeDays| adults_id, corr=corCompSymm(form =~mAgeDays),data=mate,method="ML")
a3 <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2)+I(mAgeDays^3) , random = ~1+mAgeDays| adults_id, corr=corAR1(form=~mAgeDays),data=mate,method="ML")
# #*******************************************
aictab(list(a1,a2,a3))
AICc(a1)
AICc(a2)
AICc(a3)

Weights(c(AICc(a1),AICc(a2),AICc(a3)))


c1 <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2)+I(mAgeDays^3) , random = ~1+mAgeDays| adults_id,data=mate,method="ML")
c2 <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2) , random = ~1+mAgeDays| adults_id, data=mate,method="ML")
c3 <- lme(wet_weight ~ mAgeDays , random = ~1+mAgeDays| adults_id,data=mate,method="ML")
c4 <- lme(wet_weight ~ 1, random = ~1| adults_id, data=mate,method="ML")
aictab(list(c1,c2,c3,c4))
length(mate[,1])/9
AIC(c1)
AIC(c2)
AIC(c3)
AIC(c4)
Weights(c(AICc(c1),AICc(c2),AICc(c3),AICc(c4)))

summary(c2)
modMate <- c2

plot(mate$mAgeDays
     ,residuals(c2)
     ,ylim=c(-20,20)
     )


#*******************Nutritional stress treatment****************
#************************RANDOM EFFECTS***************************
rand3 <- lm(wet_weight ~ mAgeDays+I(mAgeDays^2)+I(mAgeDays^3), data=nuts)
rand2 <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2)+I(mAgeDays^3), random = ~1 | adults_id, data=nuts,method="ML")
rand1 <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2)+I(mAgeDays^3), random = ~1+mAgeDays| adults_id, data=nuts,method="ML")

aictab(list(rand1,rand2))
aictab(list(rand3))
Weights(c(AIC(rand3),AIC(rand2),AIC(rand1)))



a1 <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2)+I(mAgeDays^3) , random = ~1+mAgeDays| adults_id, data=nuts,method="ML")
a2 <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2)+I(mAgeDays^3) , random = ~1+mAgeDays| adults_id, corr=corCompSymm(form =~mAgeDays),data=nuts,method="ML")
a3 <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2)+I(mAgeDays^3) , random = ~1+mAgeDays| adults_id, corr=corAR1(form=~mAgeDays),data=nuts,method="ML")
# #*******************************************
aictab(list(a1,a2,a3))
AICc(a1)
AICc(a2)
AICc(a3)

Weights(c(AICc(a1),AICc(a2),AICc(a3)))


c1 <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2)+I(mAgeDays^3) , random = ~1+mAgeDays| adults_id,data=nuts,method="ML")
c2 <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2) , random = ~1+mAgeDays| adults_id, data=nuts,method="ML")
c3 <- lme(wet_weight ~ mAgeDays , random = ~1+mAgeDays| adults_id, data=nuts,method="ML")
c4 <- lme(wet_weight ~ 1, random = ~1| adults_id, data=nuts,method="ML")
aictab(list(c1,c2,c4))

AIC(c1)
AIC(c2)
AIC(c3)
AIC(c4)
Weights(c(AICc(c1),AICc(c2),AICc(c4)))

summary(c2)
modNuts <- c2

plot(nuts$mAgeDays,residuals(c2))
plot(fitted(c2),residuals(c2))
#***********************************************


#**************Predict**************************
modTreat1 <- lme(wet_weight ~ mAgeDays*name+I(mAgeDays^2)*name , random = ~1+mAgeDays| adults_id,data=pupWeight,method="ML")
AIC(modTreat1)
modTreat2 <- lme(wet_weight ~ mAgeDays+I(mAgeDays^2)+name , random = ~1+mAgeDays| adults_id,data=pupWeight,method="ML")
AIC(modTreat2)


summary(modTreat2)

summary(modCtrl)
summary(modMate)
summary(modNuts)

ctrl$pred1 <- predict(modCtrl,newdata=ctrl,level=0)
ctrl$pred2 <- predict(modCtrl,newdata=ctrl,level=1)
mate$pred1 <- predict(modMate,newdata=mate,level=0)
mate$pred2 <- predict(modMate,newdata=mate,level=1)
nuts$pred1 <- predict(modNuts,newdata=nuts,level=0)
nuts$pred2 <- predict(modNuts,newdata=nuts,level=1)


all <- rbind.data.frame(ctrl,mate,nuts)
all$name <- as.factor(all$name)
levels(all$name) <- c("Control","Mating delay","Nutritional stress")

#*******************************Plot************************************************
tiff("Fig4_wet_weight.tiff", height = 3, width = 6, units = 'in', compression="lzw", res=400)
ggplot(all,aes()) +
  scale_color_manual(values=c("#191919","#5E5E5E","#808080")) +
  geom_line(data=all,aes(x=mAgeDays,y=pred2,group=adults_id),col="darkgrey",lwd=.2) + 
  geom_line(data=all,aes(x=mAgeDays,y=pred1),col="black") + 
  xlim(15,100) +
  geom_point(data=weight.means,aes(x=mAgeBins,y=mean.weight),size=0.8) +
  geom_errorbar(data=weight.means,aes(x=mAgeBins,ymin=lower, ymax=upper), width=.1,alpha=0.5) +
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
    ,legend.position =c(0.8,0.2)
    ,legend.title = element_blank()
    ,strip.background = element_rect(colour="white", fill="white")
    ,panel.border = element_blank()
  ) +
  facet_wrap(~name,ncol=3,scales="free_x")
dev.off()
