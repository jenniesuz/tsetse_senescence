# This script analyses offspring survival and produces figure 5 of the manuscript

library(plyr)          # packages required
library(reshape)
library(ggplot2)
library(gridExtra)
library(AICcmodavg)
library(MuMIn)
library(lme4)
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


tiff("fig_offspringSurvHist.tiff", height = 3, width = 6, units = 'in', compression="lzw", res=400)
ggplot(pupEmerged) +
  geom_histogram(aes(daysSurv)) +
  labs(  x="Days surviving"
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

#*****************************************control*********************************************
ctrl <- pupEmerged[pupEmerged$name %in% "control",]
#*******************random effects***********
y1 <- lme(as.numeric(daysSurv) ~ wet_weight + mAgeDays + I(mAgeDays^2), random = ~ 1 + mAgeDays + I(mAgeDays^2)| adults_id,data=ctrl,method="ML")
y2 <- lme(as.numeric(daysSurv) ~ wet_weight + mAgeDays + I(mAgeDays^2), random = ~ 1 | adults_id,data=ctrl,method="ML")
y3 <- lm(as.numeric(daysSurv) ~ wet_weight + mAgeDays + I(mAgeDays^2), data=ctrl)
AICc(y1) 
AICc(y2)
AICc(y3)
aictab(list(y1,y2))
aictab(list(y3))
Weights(c(AICc(y1),c(AICc(y2),AICc(y3))))
#****************fixed effects***************
y4 <- lm(as.numeric(daysSurv) ~ wet_weight + mAgeDays, data=ctrl)
y5 <- lm(as.numeric(daysSurv) ~ wet_weight, data=ctrl)
y6 <- lm(as.numeric(daysSurv) ~ 1,data=ctrl)
AICc(y3) 
AICc(y4)
AICc(y5) 
AICc(y6) 
aictab(list(y3,y4,y5,y6))
Weights(c(AICc(y3),AICc(y4),AICc(y5),AICc(y6)))
#***********************************************
summary(y3)
#car::vif(y4)

#*****************Predict effect of age*******************
predFuncC <- function(weightDat=as.numeric(weights[1])) {
  ctrlAgedat <- cbind.data.frame(wet_weight=rep(weightDat,100),mAgeDays=c(1:100))
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
mate <- pupEmerged[pupEmerged$name %in% "mate_delay",]
#*******************random effects**************
y1 <- lme(as.numeric(daysSurv) ~ wet_weight + mAgeDays + I(mAgeDays^2), random = ~ 1 + mAgeDays + I(mAgeDays^2)| adults_id,data=mate,method="ML")
y2 <- lme(as.numeric(daysSurv) ~ wet_weight + mAgeDays + I(mAgeDays^2), random = ~ 1 | adults_id,data=mate,method="ML")
y3 <- lm(as.numeric(daysSurv) ~ wet_weight + mAgeDays + I(mAgeDays^2), data=mate)
AICc(y2) 
AICc(y3)
aictab(list(y2))
aictab(list(y3))
Weights(c(AICc(y2),AICc(y3)))
#*******************fixed effects******************
y4 <- lm(as.numeric(daysSurv) ~ wet_weight + mAgeDays , data=mate)
y5 <- lm(as.numeric(daysSurv) ~ wet_weight, data=mate)
y6 <- lm(as.numeric(daysSurv) ~ 1,data=mate)
AICc(y3) 
AICc(y4) 
AICc(y5) 
AICc(y6) 
aictab(list(y3,y4,y5,y6))
Weights(c(AICc(y3),AICc(y4),AICc(y5),AICc(y6)))
summary(y3)

#************************predict************************
#***************age**************

predFuncM <- function(weightDat=as.numeric(weights[1])) {
  mateAgedat <- cbind.data.frame(wet_weight=rep(weightDat,100),mAgeDays=c(1:100))
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
nuts <- pupEmerged[pupEmerged$name %in% "nutrition",]
#*********************random effects***********************
y1 <- lmer(as.numeric(daysSurv) ~ wet_weight + mAgeDays + I(mAgeDays^2) + (1 + mAgeDays + I(mAgeDays^2)| adults_id),data=nuts,REML=F)
y2 <- lmer(as.numeric(daysSurv) ~ wet_weight + mAgeDays + I(mAgeDays^2) + (1 | adults_id),data=nuts,REML=F)
y3 <- lm(as.numeric(daysSurv) ~ wet_weight + mAgeDays + I(mAgeDays^2) , data=nuts)
AICc(y2) 
AICc(y3) 
aictab(list(y2))
aictab(list(y3))

summary(y2)
Weights(c(AICc(y2),AICc(y3)))
#********************fixed effects****************
y4 <- lmer(as.numeric(daysSurv) ~ wet_weight+mAgeDays + ( 1 | adults_id),data=nuts,REML=F)
y5 <- lmer(as.numeric(daysSurv) ~ wet_weight + ( 1 | adults_id),data=nuts,REML=F)
y6 <- lmer(as.numeric(daysSurv) ~ 1 + ( 1 | adults_id),data=nuts,REML=F)

aictab(list(y2,y4,y5,y6))

Weights(c(AICc(y2),AICc(y4),AICc(y5),AICc(y6)))

y4a <- lme(as.numeric(daysSurv) ~ wet_weight+mAgeDays, random = ~ 1 | adults_id,data=nuts,method="ML")

summary(y2)


#***********************predict**************
#***age****

predFuncN <- function(weightDat=as.numeric(weights[1])) {
  nutsAgedat <- cbind.data.frame(wet_weight=rep(weightDat,100),mAgeDays=c(1:100),adults_id=rep(113,100))
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

nutsAge <- nutsAge[,-3]

ageEffect <- rbind.data.frame(ctrlAge,mateAge,nutsAge)
ageEffect$name <- as.factor(ageEffect$name)
levels(ageEffect$name) <- c("Control","Mating delay","Nutritional stress")

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
  xlab("Mother age") +
  labs(title="b)") +
  theme_set(theme_bw()) +
  theme(axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=8)
        ,legend.key.size = unit(0.5,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,legend.position ="none"
        ,legend.title = element_blank()
        ,strip.background = element_rect(colour="white", fill="white")
        ,panel.border = element_blank()
        #,strip.text.x = element_blank()
  ) +
  facet_wrap(~name)


tiff("Fig5_offspringEmSurv.tiff", height = 6, width = 6, units = 'in', compression="lzw", res=400)
grid.arrange(agePlot,survAge.plot ,nrow=2,ncol=2,widths=c(2,1)
             ,layout_matrix = rbind(c(1, NA),
                                   c(2,2)))
dev.off()

