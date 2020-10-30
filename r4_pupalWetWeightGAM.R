# This script analyses pupal wet weight and produces Figure 4 of the manuscript
source("r4_pupalWetWeight.R")

library(plyr)        
library(reshape)
library(ggplot2)
library(AICcmodavg)
library(MuMIn)
library(mgcv)

#***** GAMM mgcv - cross validation can determine optimal amount of smoothing
#*******************Control treatment****************
ctrl$adults_idF <- as.factor(ctrl$adults_id)

#******************choose knots**********************
ctrlGAMfits <- lapply(3:10,function(x){
  g1Gam <- gam(wet_weight ~ s(mAgeDays,bs="cr",k=x)+s(adults_idF,bs="re")
               ,family=gaussian,data=ctrl,se=T,method="ML")
  return(c(x,AICc(g1Gam,nobs=length(ctrl$name))))
  
})
ctrlGAMfits <- do.call(rbind.data.frame,ctrlGAMfits)
names(ctrlGAMfits) <- c("Number of knots","AICc")
saveRDS(ctrlGAMfits,file="GAMfitsctrl.rds")
#*****************************************************

g1Gam <- gam(wet_weight ~ s(mAgeDays,bs="cr",k=3)+s(adults_idF,bs="re")
             ,family=gaussian,data=ctrl,se=T,method="ML")

g1M2 <- predict(g1Gam
                ,level=1
                ,exclude = 's(adults_idF)'
                ,type="response",se=T,unconditional=T)

ctrl$pred <- g1M2$fit
ctrl$se <- g1M2$se.fit

#***************************Mating delay*******************************

mate$adults_idF <- as.factor(mate$adults_id)

#*****************choose knots****************
mateGAMfits <- lapply(3:10,function(x){
  g1Gam <- gam(wet_weight ~ s(mAgeDays,bs="cr",k=x)+s(adults_idF,bs="re")
               ,family=gaussian,data=mate,se=T,method="ML")
  return(c(x,AICc(g1Gam,nobs=length(mate$name))))
  
})
mateGAMfits <- do.call(rbind.data.frame,mateGAMfits)
names(mateGAMfits) <- c("Number of knots","AICc")
saveRDS(mateGAMfits,file="GAMfitsmate.rds")
#********************************************

g1Gam <- gam(wet_weight ~ s(mAgeDays,k=4,bs="cr")+s(adults_idF,bs="re")
             ,family=gaussian,data=mate,se=T)
g1M2 <- predict(g1Gam
                ,level=1
                ,exclude = 's(adults_idF)'
                ,type="response",se=T,unconditional=T)

mate$pred <- g1M2$fit
mate$se <- g1M2$se.fit


#*******************Nutritional stress treatment****************
nuts$adults_idF <- as.factor(nuts$adults_id)

#**********************choose knots***************
nutsGAMfits <- lapply(3:10,function(x){
  g1Gam <- gam(wet_weight ~ s(mAgeDays,bs="cr",k=x)+s(adults_idF,bs="re")
               ,family=gaussian,data=nuts,se=T,method="ML")
  return(c(x,AICc(g1Gam,nobs=length(nuts$name))))
  
})
nutsGAMfits <- do.call(rbind.data.frame,nutsGAMfits)
names(nutsGAMfits) <- c("Number of knots","AICc")
saveRDS(nutsGAMfits,file="GAMfitsnuts.rds")
#************************************************

g1Gam <- gam(wet_weight ~ s(mAgeDays,k=4,bs="cr")+s(adults_idF,bs="re")
             ,family=gaussian,data=nuts,se=T,method="ML")
g1M2 <- predict(g1Gam
                ,level=1
                ,exclude = 's(adults_idF)'
                ,type="response",se=T,unconditional=T)

nuts$pred <- g1M2$fit
nuts$se <- g1M2$se.fit


#*******************************Plot************************************************
all <- rbind.data.frame(ctrl,mate,nuts)
all$name <- as.factor(all$name)
levels(all$name) <- c("Control","Mating delay","Nutritional stress")


#tiff("FigSX_wet_weightGam.tiff", height = 3, width = 6, units = 'in', compression="lzw", res=400)
plotWWgam <- ggplot(all,aes()) +
  scale_color_manual(values=c("#875777","#eab051","#c0b9ac")) +
  geom_line(data=all,aes(x=mAgeDays,y=pred1),linetype=2) + 
  geom_line(data=all,aes(x=mAgeDays,y=pred,col=name),size=0.8) +
  geom_line(aes(x=mAgeDays,y=pred-2*se,col=name)) +
  geom_line(aes(x=mAgeDays,y=pred+2*se,col=name)) +
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
#dev.off()

saveRDS(plotWWgam,file="FigSuppGAMWeight.rds")

