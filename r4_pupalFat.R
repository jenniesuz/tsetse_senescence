# This script plots pupal fat against pupal wet weight for the supplementary
# material
source("r1_queries.R")

library(plyr)        
library(reshape)
library(ggplot2)


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
pupWeight <- pupWeight[!pupWeight$wet_weight%in%NA,]     # remove those - 31 of 1239 for which no wet weight data - either lost or crushed
pupWeight <- pupWeight[!pupWeight$residual_dry_weight%in%NA,]
pupWeight$fat <- pupWeight$dry_weight - pupWeight$residual_dry_weight

cor(pupWeight$wet_weight,pupWeight$fat,method="pearson")


#*******************************Plot************************************************
tiff("fig_fatWetWeight.tiff", height = 4, width = 4, units = 'in', compression="lzw", res=400)
ggplot(pupWeight,aes()) +
  scale_color_manual(values=c("#191919","#5E5E5E","#808080")) +
  geom_point(aes(x=wet_weight,y=fat,col=name)) + 
  labs(  x="Offspring wet weight (mg)"
         ,y="Offspring fat (mg)"
  ) + 
  theme_set(theme_bw()) +
  theme(
    axis.line = element_line(color = 'black')
    ,text=element_text(size=8)
    ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
    ,axis.text=element_text(size=7)
    ,legend.key.size = unit(0.8,"line")
    ,legend.background = element_blank()
    ,legend.text=element_text(size=7)
    ,legend.position =c(0.2,0.9)
    ,legend.title = element_blank()
    ,strip.background = element_rect(colour="white", fill="white")
    ,panel.border = element_blank()
  ) 

dev.off()
