

library("plyr")          # packages required
library("reshape")
library("ggplot2")
library("survival")
library("survminer")
library("binom")
library("nlme")
library("gridExtra")
library("coxme")
library("AICcmodavg")
require(lme4)
require("MuMIn")



source("r3_queries.R")    # query database 
pupae$date_emerged <- as.Date(pupae$date_emerged,format=c("%d/%m/%Y"))
pupae$date_of_death <- as.Date(pupae$date_of_death,format=c("%d/%m/%Y"))
pupae$date_of_emergence_min <- as.Date(pupae$date_of_emergence_min, format="%d/%m/%Y")  # mother emergence date
pupae$larviposition_date <- as.Date(pupae$larviposition_date, format="%d/%m/%Y")        # date of larviposition

# add column for mother age in days
pupae$mAgeDays <- pupae$larviposition_date - pupae$date_of_emergence_min
pupae$mAgeDays <- as.numeric(pupae$mAgeDays)


pupae <- pupae[pupae$abortion %in% "0.0",]
ddply(pupae,.(name),summarise,length(wet_weight))

#1451 pupae
length(pupae[pupae$abortion %in% "0.0",1]) #1240



ddply(pupae,.(name),summarise,l=length(wet_weight))
# 555, 297, 357

length(pupae[pupae$treatment %in% "fat analysis",1]) # 437
length(pupae[!pupae$treatment %in% "fat analysis",1]) # 803


pupae <- pupae[!pupae$treatment %in% "fat analysis",]

length(pupae$emerged[pupae$emerged %in% "0.0"]) # 201
length(pupae$emerged[pupae$emerged %in% "1.0"]) # 736
boxplot(pupae$wet_weight ~ pupae$emerged)
boxplot(pupae$mAgeDays ~ pupae$emerged)




class(pupae$emerged)
pupae$emerged <- as.integer(pupae$emerged)

pupae <- pupae[!pupae$wet_weight %in% NA,]


emergedByTreat <- ddply(pupae,.(name),summarise,em=sum(emerged),totem=length(emerged))
#****************Analysis*****************

# ****************************control group*****************************
ctrl <- pupae[pupae$name %in% "control",]
#********random effects********
rand1 <- glm(emerged ~ 1,data=ctrl,family=binomial)
rand2 <- glm(emerged ~ wet_weight + mAgeDays + I(mAgeDays^2) , data=ctrl,family=binomial)
rand3 <- glmer(emerged ~ wet_weight + mAgeDays + I(mAgeDays^2) + (1|adults_id), data=ctrl,family=binomial)
rand4 <- glmer(emerged ~ wet_weight + mAgeDays + I(mAgeDays^2) + (1 + mAgeDays + I(mAgeDays^2)|adults_id ))
# can't fit including repeated measures
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
# singular fit. convergence failure

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


# ********************************nutrition group*********************** 
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

#************predict***********
weights <- quantile(nuts$wet_weight)

ilink <- family(nuts2)$linkinv

# weightP <- with(nuts,
#               data.frame(wet_weight = seq(min(wet_weight), max(wet_weight),
#                                           length = 100)
#               ,mAgeDays=rep(45,100)))
# weightP <- cbind(weightP, predict(nuts2, weightP, type = "link", se.fit = TRUE)[1:2])
# weightP <- transform(weightP, Fitted = ilink(fit), Upper = ilink(fit + (1.96 * se.fit)),
#                    Lower = ilink(fit - (1.96 * se.fit)))
# 


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
# weightPlot <- ggplot(weightP) +
#   geom_line(aes(x=wet_weight,y=Fitted)) +
#   geom_ribbon(data = weightP, aes(ymin = Lower, ymax = Upper, x = wet_weight),
#               fill = "grey", alpha = 0.5, inherit.aes = FALSE) +
#   xlab("Wet weight (mg)") + 
#   ylab("Predicted probability of emergence") +
#   theme_set(theme_bw()) +
#   theme(axis.line = element_line(color = 'black')
#         ,text=element_text(size=10)
#         ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
#         ,axis.text=element_text(size=8)
#         ,legend.key.size = unit(0.8,"line")
#         ,legend.background = element_blank()
#         ,legend.text=element_text(size=9)
#         ,legend.position =c(0.2,0.9)
#         ,legend.title = element_blank()
#         ,strip.background = element_rect(colour="white", fill="white")
#         ,panel.border = element_blank()
#   )

cols <- c("#cccccc"
  ,"#969696"
  ,"#636363"
  ,"#252525")

agePlot <- ggplot(ageP) +
  geom_line(aes(x=mAgeDays,y=Fitted,col=q,linetype=q)) +
  ylim(0.2,1) +
  ylab("Predicted probability of emergence") +
  xlab("Mother age (days)") +
  labs (title="a)") +
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
       # ,legend.position =c(0.1,0.2)
        ,strip.background = element_rect(colour="white", fill="white")
        ,panel.border = element_blank()
  )

  #geom_ribbon(data = ageP, aes(ymin = Lower, ymax = Upper, x = mAgeDays,group=q)
  #            ,fill = "grey", alpha = 0.5, inherit.aes = FALSE) +
  
tiff("fig_nutEmergence.tiff", height = 3, width = 6, units = 'in', compression="lzw", res=400)
grid.arrange(weightPlot,agePlot,ncol=2)
dev.off()


ggplot(mate, aes(x = wet_weight, colour = emreged)) +
  geom_line(stat = "density") 
