
# This script produces Figure 2 of the manuscript: Kaplan-Meier & Cox proportional hazards model for female survival
# also parameteric survival analysis

source("r1_queries.R")   # query database

library(plyr)         
library(reshape)
library(ggplot2)
library(grid)
library(gridExtra)
library(survival)
library(survminer)
library(nlme)
library(muhaz)
library(AICcmodavg)
library(MuMIn)

#*********************Total deaths by treatment group**********************
sapply(unique(adult.deaths$name),function(x){               
  temp <- adult.deaths[adult.deaths$name %in% x,]
  return(length(temp$date_of_death[!temp$date_of_death %in% NA]))
})
#**************************************************************************

#*********************Format data*****************************************
adult.deaths$date_of_emergence_min <- as.Date(adult.deaths$date_of_emergence_min
                                              , format="%d/%m/%Y") # convert date columns to date format

adult.deaths$date_of_death <- as.Date(adult.deaths$date_of_death,format="%d/%m/%Y")

adult.deaths$name <- as.factor(adult.deaths$name)  # change treatment to a factor

# calculate days to death
adult.deaths$dayOfDeath <- adult.deaths$date_of_death - adult.deaths$date_of_emergence_min
adult.deaths$dayOfDeath <- as.numeric(adult.deaths$dayOfDeath)

# create a column for event 'dead' 1 or 'alive' 0 (at last obeservation)
adult.deaths$dead <- NA
adult.deaths$dead[!(adult.deaths$date_of_death %in% NA)] <- 1
adult.deaths$dead[adult.deaths$date_of_death %in% NA] <- 0
adult.deaths$dayOfDeath[adult.deaths$date_of_death %in% NA] <- 100
#**************************************************************************

#************************KM survival***************************************
survObj <- Surv(time=adult.deaths$dayOfDeath,event=adult.deaths$dead)

survFit <- survfit(survObj ~ name, data=adult.deaths,conf.type="log-log")

summary(survFit)   
#************************************************************************

#**********************Cumulative count of deaths***********************
dates <- seq.Date(adult.deaths$date_of_emergence_min[1]
                  ,adult.deaths$date_of_emergence_min[1]+max(survFit$time)
                  ,by=1) # create string of days since experiment start

# want to generate a cumulative count over time by treatment
deathsByDateTreat <- lapply(unique(adult.deaths$name),function(y){
  temp <- adult.deaths[adult.deaths$name %in% y,]
  deathsByDate <- sapply(1:length(dates),function(x){
    cum.dates <- dates[1:x]
    return(length(temp$date_of_death[temp$date_of_death %in% cum.dates]))
  })
  return(deathsByDate)
})

deathsByDateTreat <- do.call(cbind.data.frame,deathsByDateTreat)
names(deathsByDateTreat) <- unique(adult.deaths$name)
deathsByDateTreat$date <- dates

deathsByDateTreat <- melt(deathsByDateTreat,id.vars=c("date"))

deathsByDateTreat$starting.size <- numeric(length(deathsByDateTreat$date))

deathsByDateTreat$starting.size[deathsByDateTreat$variable %in% "control"] <- treatment$starting_size[treatment$name %in% "control"]
deathsByDateTreat$starting.size[deathsByDateTreat$variable %in% "nutrition"] <- treatment$starting_size[treatment$name %in% "nutrition"]
deathsByDateTreat$starting.size[deathsByDateTreat$variable %in% "mate_delay"] <- treatment$starting_size[treatment$name %in% "mate_delay"]

deathsByDateTreat$day <- rep(seq(1,length(deathsByDateTreat[deathsByDateTreat$variable %in% "control",1]),1),3)
names(deathsByDateTreat)[3] <- "dead"
deathsByDateTreat$surv <- (deathsByDateTreat$starting.size - deathsByDateTreat$dead)/deathsByDateTreat$starting.size
names(deathsByDateTreat)[2] <- "name"
#*************************************************************************

#**************************hazard function*********************************
hazCtrl <- muhaz(adult.deaths$dayOfDeath[adult.deaths$name %in% "control"]
                 ,adult.deaths$dead[adult.deaths$name %in% "control"]
                 ,max.time=100
                 ,bw.grid=30
                 , bw.method="global",b.cor="none")


hazMate <- muhaz(adult.deaths$dayOfDeath[adult.deaths$name %in% "mate_delay"]
                 ,adult.deaths$dead[adult.deaths$name %in% "mate_delay"]
                 ,max.time=100
                 ,bw.grid=30
                 , bw.method="global",b.cor="none")

hazNut <- muhaz(adult.deaths$dayOfDeath[adult.deaths$name %in% "nutrition"]
             ,adult.deaths$dead[adult.deaths$name %in% "nutrition"]
             ,max.time=100
             ,bw.grid=30
            , bw.method="global",b.cor="none")

haz <- c(hazCtrl$haz.est,hazMate$haz.est,hazNut$haz.est)
haz <- cbind.data.frame(haz=haz,name=c(rep("Control",101)
                                       ,rep("Mating delay",101)
                                       ,rep("Nutritional stress",101)
                                       )
                        )
#*****************************************************************************

#****************************hazard plot***************************************
hazardPlot <- ggplot(haz) +
  scale_color_manual(values=c("#191919","#5E5E5E","#808080"))  +
  geom_line(aes(x=rep(0:100,3),y=haz,col=name,linetype=name)) +
  theme_set(theme_bw()) +
  ggtitle("Hazard") +
  xlab("Time since emergence (days)") +
  ylab("Smoothed hazard") +
  theme(axis.line = element_line(color = 'black')
        ,text=element_text(size=7)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=7)
        ,plot.title=element_text(size=7)
        ,legend.key.size = unit(0.6,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=6)
        ,legend.position =c(0.3,0.9)
        ,legend.title = element_blank()
        ,strip.background = element_rect(colour="white", fill="white")
        ,panel.border = element_blank()
  )

#***********************Fig 2****************************************
tiff("Fig2_mother_survival.tiff", height = 3, width = 5, units = 'in', compression="lzw", res=400)
ggsv <- ggsurvplot(survFit, data=adult.deaths
           ,linetype="strata"
           ,palette=c("#191919","#5E5E5E","#808080")
           ,ylim=c(0.5,1)
           ,xlab="Time since emergence (days)"
           ,title="Survival"
           ,font.x=c(7)
           ,font.y=c(7)
           ,font.title=c(7)
           ,font.tickslab=c(6)
           ,font.legend=c(5)
           ,font.caption=c(5) 
           ,conf.int = TRUE
           ,legend =c(0.3,0.3)
           ,risk.table=F
           ,risk.table.y.text.col = T # colour risk table text annotations.
           ,risk.table.y.text = FALSE
           ,risk.table.font=c(3)
           ,risk.table.title="Number alive"
           ,legend.title = "Treatment:"
           ,legend.labs = c("Control", "Mating delay","Nutritional stress")) #+
grid.arrange(ggsv$plot,hazardPlot,nrow=1)
dev.off()

#**************************************Cox PHM***********************************

fit.coxph <- coxph(survObj ~ 1,data=adult.deaths)
fit.coxph1 <- coxph(survObj ~ name , 
                   data = adult.deaths)


aictab(list(fit.coxph,fit.coxph1))

AIC(fit.coxph)
AIC(fit.coxph1)

#****************************************************************************


#*********************Parametric analysis*************************************
nuts <- adult.deaths[adult.deaths$name %in% "nutrition",]
ctrl <- adult.deaths[adult.deaths$name %in% "control",]
mate <- adult.deaths[adult.deaths$name %in% "mate_delay",]

survFitexpNut <- survreg(Surv(time=nuts$dayOfDeath,event=nuts$dead) ~ 1, dist="exponential" )
survFitweibNut <- survreg(Surv(time=nuts$dayOfDeath,even=nuts$dead) ~ 1,dist="weibull")
AIC(survFitexpNut) 
AIC(survFitweibNut) 
Weights(c(AIC(survFitexpNut),AIC(survFitweibNut)))

survFitexpCtrl <- survreg(Surv(time=ctrl$dayOfDeath,event=ctrl$dead) ~ 1, dist="exponential" )
survFitweibCtrl <- survreg(Surv(time=ctrl$dayOfDeath,even=ctrl$dead) ~ 1,dist="weibull")
AIC(survFitexpCtrl) 
AIC(survFitweibCtrl) 
Weights(c(AIC(survFitexpCtrl),AIC(survFitweibCtrl)))


survFitexpMate <- survreg(Surv(time=mate$dayOfDeath,event=mate$dead) ~ 1, dist="exponential" )
survFitweibMate <- survreg(Surv(time=mate$dayOfDeath,even=mate$dead) ~ 1,dist="weibull")
AIC(survFitexpMate) 
AIC(survFitweibMate) 
Weights(c(AIC(survFitexpMate),AIC(survFitweibMate)))

#*******************************************************************************











