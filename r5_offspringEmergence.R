
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



#********************************************Analysis***************************************

# ****************************Control group*****************************
ctrl <- pupae[pupae$name %in% "control",]
#********random effects********
rand1 <- glm(emerged ~ 1,data=ctrl,family=binomial)
rand2 <- glm(emerged ~ wet_weight + mAgeDays + I(mAgeDays^2) , data=ctrl,family=binomial)
rand3 <- glmer(emerged ~ wet_weight + mAgeDays + I(mAgeDays^2) + (1|adults_id), data=ctrl,family=binomial)
rand4 <- glmer(emerged ~ wet_weight + mAgeDays + I(mAgeDays^2) + (1 + mAgeDays + I(mAgeDays^2)|adults_id ))

summary(rand2)
###********fixed effects**********
ctrl1 <- glm(emerged ~ wet_weight + mAgeDays + I(mAgeDays^2),data=ctrl,family=binomial)
ctrl2 <- glm(emerged ~ wet_weight + mAgeDays,data=ctrl,family=binomial)
ctrl3  <- glm(emerged ~ wet_weight, data=ctrl,family=binomial)
ctrl4 <- glm(emerged ~ 1, data=ctrl,family=binomial)
AICc(ctrl1) 
AICc(ctrl2) 
AICc(ctrl3)
AICc(ctrl4) 
aictab(list(ctrl1,ctrl2,ctrl3,ctrl4))
Weights(c(AICc(ctrl1),AICc(ctrl2),AICc(ctrl3),AICc(ctrl4)))
#**************************************************

#**********mating delay group**********
mate <- pupae[pupae$name %in% "mate_delay",]
rand1 <- glm(emerged ~ 1,data=mate,family=binomial)
rand2 <- glm(emerged ~ wet_weight + mAgeDays + I(mAgeDays^2), data=mate,family=binomial)
rand3 <- glmer(emerged ~ wet_weight + mAgeDays + I(mAgeDays^2) + (1|adults_id), data=mate,family=binomial)
rand4 <- glmer(emerged ~ wet_weight + mAgeDays + I(mAgeDays^2) + (1 + mAgeDays + I(mAgeDays^2)|adults_id ),data=mate,family=binomial)

###fixed effects
mate1 <- glm(emerged ~ wet_weight + mAgeDays + I(mAgeDays^2) ,data=mate,family=binomial)
mate2 <- glm(emerged ~ wet_weight + mAgeDays ,data=mate,family=binomial)
mate3  <- glm(emerged ~ wet_weight , data=mate,family=binomial)
mate4 <- glm(emerged ~ 1, data=mate,family=binomial)
AICc(mate1)
AICc(mate2)
AICc(mate3) 
AICc(mate4) 
aictab(list(mate1,mate2,mate3,mate4))
Weights(c(AICc(mate1),AICc(mate2),AICc(mate3),AICc(mate4)))

#************************************************


# ********************************Nutritional stress group*********************** 
#******random effects*******
nuts <- pupae[pupae$name %in% "nutrition",]
rand1 <- glm(emerged ~ 1,data=nuts,family=binomial)
rand2 <- glm(emerged ~ wet_weight+mAgeDays + I(mAgeDays^2) , data=nuts,family=binomial)
rand3 <- glmer(emerged ~ wet_weight+mAgeDays + I(mAgeDays^2)  + (1|adults_id), data=nuts,family=binomial)
rand4 <- glmer(emerged ~ wet_weight+mAgeDays + I(mAgeDays^2)  + (1+mAgeDays + I(mAgeDays^2)|adults_id), data=nuts,family=binomial)
# convergence failure/ singular fit

#*******fixed effects*******
nuts1  <- glm(emerged ~ wet_weight + mAgeDays + I(mAgeDays^2) , data=nuts,family=binomial)
nuts2 <- glm(emerged ~ wet_weight + mAgeDays,data=nuts,family=binomial)
nuts3 <- glm(emerged ~ wet_weight ,data=nuts,family=binomial)
nuts4 <- glm(emerged ~ 1 , data=nuts,family=binomial)
AICc(nuts1) 
AICc(nuts2) 
AICc(nuts3) 
AICc(nuts4) 
aictab(list(nuts1,nuts2,nuts3,nuts4))
Weights(c(AICc(nuts1),AICc(nuts2),AICc(nuts3),AICc(nuts4)))

summary(nuts2)
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
