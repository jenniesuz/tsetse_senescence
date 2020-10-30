### PRELIMINARY ANALYSIS OF COLONY DATA ###
library(ggplot2)
library(plyr)
library(lme4)
library(survival)


## first analyses...

dat0<-read.csv("dat_colonyData2019.csv")
dat0$rowID <- row.names(dat0)

dat<-dat0
table(duplicated(dat))
dat[duplicated(dat),]
dat<-dat[!duplicated(dat),] # remove fully duplicated rows

# convert tray (from 12 to 1, youngest to oldest) to 'maternal age' (from 1 to 12 weeks of age)
tray.ages <- 12:1
mum.ages <- 1:12

dat$MatAge <- 0
for(i in 1:12)
{
	dat$MatAge[dat$Tray==tray.ages[i]]<-mum.ages[i]
	
}

# quickly remove obvious errors
dat$Mg_avg[dat$Mg_avg==0]<-NA; 

dat.nna <- na.omit(dat)
plot(dat.nna$rowID,dat.nna$Mg_avg)

# remove outliers: 4sd from mean? 

points(dat.nna$rowID[abs(dat.nna$Mg_avg-mean(dat.nna$Mg_avg))>3*sd(dat.nna$Mg_avg)],dat.nna$Mg_avg[abs(dat.nna$Mg_avg-mean(dat.nna$Mg_avg))>3*sd(dat.nna$Mg_avg)],pch=16,col="red")

# remove outliers more than 3sd from mean: 
dat$Mg_avg[abs(dat$Mg_avg-mean(na.omit(dat$Mg_avg)))>4*sd(na.omit(dat$Mg_avg))]=NA

plot(dat$Mg_avg)


#***********************Plot for supp mat******************
#*********************************************************

# get average and sd weight for females from different trays
means <- tapply(dat.nna$Mg_avg,dat.nna$MatAge,mean)
sd <- tapply(dat.nna$Mg_avg,dat.nna$MatAge,sd)
n <- tapply(dat.nna$Mg_avg,dat.nna$MatAge,length)
se <- sd/sqrt(n)


tiff("FigColonyData.tiff", height = 3, width = 6, units = 'in', compression="lzw", res=400)

par(mfrow=c(1,2))
plot(2:12,means[-1],pch=19,ylim=c(25,30),cex=0.5,cex.axis=0.5,cex.lab=0.8
     ,xlim=c(1,12),type="p"
     ,bty="n",xlab="Mother age (weeks)",ylab="Average pupa weight (mg)")
arrows(2:12,means[-1]-1.96*se[-1],2:12,means[-1]+1.96*se[-1],angle=90,length=0)

# add predicted quadratic effect:
# m1 <- lm(Mg_avg~MatAge+I(MatAge^2),data=dat)
# newdf<-c()
# newdf$MatAge <- seq(1,12,1)
# newdf$pred <- predict(m1,newdf)
# lines(newdf$MatAge,newdf$pred,col="grey")

mtext("(a)",side=3,adj=0,line=1)


## also plot the average number of abortions (per total pupa)

means <- lapply(unique(dat.nna$MatAge),function(x){
  temp <- dat.nna[dat.nna$MatAge %in% x,]
  temp$AbortDenom <- temp$P_count + temp$U_count + temp$A_count
  temp$AbortNum <- temp$U_count + temp$A_count
  meanci <- binom::binom.confint(x=temp$AbortNum
                                 ,n=temp$AbortDenom
                                 ,method="exact")
  mean <- mean(meanci$mean)
  lower <- mean(meanci$lower)
  upper <- mean(meanci$upper)
  
  return(c("MatAge"=x,"mean"=mean,"lower"=lower,"upper"=upper))
})

means <- as.data.frame(do.call(rbind,means))
means <- means[order(means$MatAge),]

# means <- tapply(dat.nna$Abort,dat.nna$MatAge,mean)
# sd <- tapply(dat.nna$Abort,dat.nna$MatAge,sd)
# n <- tapply(dat.nna$Abort,dat.nna$MatAge,length)
# se <- sd/sqrt(n)

plot(2:12,means$mean[-1],xlim=c(1,12),cex=0.5,cex.axis=0.5,cex.lab=0.8,ylim=c(0,0.5),pch=19,type="p",xlab="Mother age (weeks)",ylab="Proportion aborted young")
arrows(2:12,means$lower[-1],ylim=c(0,0.5),2:12,means$upper[-1],angle=90,length=0)






# add predicted quadratic effect:
# m2 <- lm(Abort~MatAge+I(MatAge^2)+I(MatAge^3),data=dat.nna)
# newdf<-c()
# newdf$MatAge <- seq(1,12,1)
# newdf$pred <- predict(m2,newdf)
# lines(newdf$MatAge,newdf$pred,col="grey")

mtext("(b)",side=3,adj=0,line=1)

dev.off()

## add cohort column
dat$Date <- as.Date(dat$Date,format="%d/%m/%y")

dat <- dat[order(dat$Date,dat$Tray),]
dat$Weekday <- weekdays(dat$Date) 

dat$Week <- strftime(as.POSIXlt(dat$Date),format="%Y/%W")
dat$WeekTray <- paste(dat$Week,dat$Tray,sep="/")

# give all the Tray 12s a unique code

dat1 <- unique(dat[,c("Tray","Week","WeekTray")])
dat1 <- dat1[order(dat1$Week,-rank(dat1$Tray)),]

dat1$CohortID <- NA
dat1$CohortID[1:12] <- 12:1

length(dat1[dat1$Tray==12,1]) # 165

dat1$CohortID[dat1$Tray==12] <- 12:176

for(i in 1:length(dat1[,1]))
{
	if(is.na(dat1$CohortID[i])) dat1$CohortID[i] = dat1$CohortID[i-13]
}


table(table(dat$Week))
table(table(dat1$CohortID))

dat <- merge(dat,dat1[,c("WeekTray","CohortID")],by="WeekTray",all.x=TRUE)

length(unique(dat$CohortID))

# 176 unique Cohort IDs

dat2 <- dat[order(dat$CohortID,dat$Date),]

# save as data file
# write.csv(dat2,"ColonyData_Cohort.csv",quote=F,row.names=F)

# start with Cohort 12 for this graph
dat2 <- dat2[dat2$CohortID>11,]

all_cohorts <- unique(dat2$CohortID)


# add cohort age, assuming 1 day old when brought into colony
# note this is not the same as maternal age (tray) because accounts for Mon, Wed, Fri collections
# may be simpler to analyze on weekly scale but did using separate days-in-week for a finer time resolution... 

dat2$Age = NA

for(i in 1:length(all_cohorts))
{
	dat2$Age[dat2$CohortID==all_cohorts[i]] <- dat2$Date[dat2$CohortID==all_cohorts[i]] - min(dat2$Date[dat2$CohortID==all_cohorts[i]])
	
}

par(mfrow=c(1,1))
cols <- gray.colors(176)
plot(1:85,rep(0,85),ylim=c(17,37),xlab="Age",ylab="Pupal weight")

for(i in 1:length(all_cohorts))
{
	plotdat <- dat2[dat2$CohortID==all_cohorts[i],]	
	points(plotdat$Age,plotdat$Mg_avg,pch=16,type="o",cex=0.5,col=rgb(0.5,0.5,0.5,0.1))
		
}

# run analysis to look at effect of age and cohort on pupal mass

dat2$zAge <- mean

mod1 <- lmer(Mg_avg ~ Age + I(Age^2) + (1|CohortID), data=dat2)
mod2 <- lmer(Mg_avg ~ Age + I(Age^2) + (Age|CohortID), data=dat2)

plot(dat$Mg_avg)


# suggests re-scale because different scales of predictors. Mean center and standardize age? 

dat2$zAge <- (dat2$Age - mean(dat2$Age))/sd(dat2$Age)

mod3 <- lmer(Mg_avg ~ zAge + I(zAge^2) + (1|CohortID), data=dat2)
mod4 <- lmer(Mg_avg ~ zAge + I(zAge^2) + (zAge|CohortID), data=dat2)

anova(mod3,mod4) # random slope model significantly better? 

# use Andrew's predict.model function to get predictions from lmer
 predict.model = function(model, data, use.random=TRUE) {

    MM = model.matrix(terms(model),data)
    fixed.effects = fixef(model)

    #make initial (fixed effects) predictions
    predictions = MM%*%fixed.effects

    random.effects = ranef(model)

    if(use.random) {
        #GO THROUGH EACH RANDOM EFFECT
        for(random.effect in names(random.effects)) {

            effect.coef = random.effects[[random.effect]]    #the coefficients for this random effect
            effect.MM = as.data.frame(MM[,colnames(MM)%in%colnames(effect.coef)])         #the data relevant to this random effect
            colnames(effect.MM) = colnames(MM)[colnames(MM)%in%colnames(effect.coef)]

            #...AND EACH LEVEL OF THAT RANDOM EFFECT
            for(effect.level in rownames(effect.coef)) {

                level.rows = data[,random.effect]==effect.level

                level.MM = as.matrix(effect.MM[level.rows,])
                #colnames(level.MM) = colnames(effect.MM)
                level.coef = as.matrix(t(effect.coef[rownames(effect.coef)==effect.level,]))
                #colnames(level.coef) = colnames(effect.coef)

                predictions[level.rows] = predictions[level.rows] + level.MM%*%level.coef
            }
        }
    }

    return(predictions)
    }


dat2$pred <- predict.model(mod3, dat2)
dat2$predRS <- predict.model(mod4, dat2)

plotdat <- na.omit(dat2[,c("zAge","Mg_avg","pred","predRS","CohortID")])

mod3 <- lmer(Mg_avg ~ zAge + I(zAge^2) + (1|CohortID), data=plotdat)
mod4 <- lmer(Mg_avg ~ zAge + I(zAge^2) + (zAge|CohortID), data=plotdat)


g <- ggplot(plotdat, aes(x=zAge, y=Mg_avg, colour=CohortID)) + geom_line()#color="gray60")
g #Plot the observed points
g + geom_line(aes(y=pred), color="tomato") #Plot the modeled points


chkdat <- fortify(mod4)

p1.0 <- ggplot(chkdat, aes(x=zAge, y=Mg_avg)) +
  geom_point() +
  geom_line(aes(x=zAge, y=.fitted), size=.75) +
  facet_wrap(~CohortID, nrow=5)
p1.0


## TO DO: May be easier to run analysis on weekly averages (lose some power, but also reduces noise?) especially as have mortality data on this sampling schedule too.

dat2$G_total[dat2$G_total==0]=NA
dat2$P_count[dat2$P_count==0]=NA

dat2$Aborted <- dat2$U_count + dat2$A_count
dat2.nna <- na.omit(dat2)

# weekly averages for each cohort:
wkdat <- ddply(dat2.nna,.(Week,CohortID),summarize,
	trayNo = min(Tray),
	matAge = min(MatAge),
	dateStart = min(Date),
	totalPupaNo = sum(P_count),
	totalPupaMass = sum(G_total),
	meanPupaMass = mean(Mg_avg),
	totalAborted = sum(Aborted),
	noDays = length(Mg_avg))

# limit graphs to when have 3 days per cohort for that week
wkdat1 <- wkdat[wkdat$noDays==3,]
wkdat1$totalAborted[wkdat1$totalAborted>100]=NA

par(mfrow=c(1,1))
plot(wkdat1$dateStart,rep(0,length(wkdat1$dateStart)),ylim=c(20,440),xlab="Date",ylab="Total offspring number",type="n")
for(i in sample(all_cohorts,60))
{
	plotdat <- wkdat1[wkdat1$CohortID==all_cohorts[i],]
	points(plotdat$dateStart,plotdat$totalPupaNo,pch=16,type="o",cex=0.5,col=cols[i])
		
}

summary(wkdat1)

n = 165
cols = rep(rgb(0.5,0.5,0.5,0.2),n) #rainbow(n) #rep(rgb(0.5,0.5,0.5,0.2),n)
plot_cohorts = sample(all_cohorts, n)

par(mfrow=c(2,2))

plot(wkdat1$matAge,rep(0,length(wkdat1$matAge)),ylim=c(0,13),xlab="Maternal age",ylab="Total offspring mass",type="n")
for(i in 1:length(plot_cohorts))
{
	plotdat <- wkdat1[wkdat1$CohortID==plot_cohorts[i],]
	points(plotdat$matAge,plotdat$totalPupaMass,pch=16,type="o",cex=0.5,col=cols[i])		
}

plot(wkdat1$matAge,rep(0,length(wkdat1$matAge)),ylim=c(20,440),xlab="Maternal age",ylab="Total offspring number",type="n")
for(i in 1:length(plot_cohorts))
{
	plotdat <- wkdat1[wkdat1$CohortID==plot_cohorts[i],]
	points(plotdat$matAge,plotdat$totalPupaNo,pch=16,type="o",cex=0.5,col=cols[i])		
}

plot(wkdat1$matAge,rep(0,length(wkdat1$matAge)),ylim=c(24,33),xlab="Maternal age",ylab="Mean offspring mass",type="n")
for(i in 1:length(plot_cohorts))
{
	plotdat <- wkdat1[wkdat1$CohortID==plot_cohorts[i],]
	points(plotdat$matAge,plotdat$meanPupaMass,pch=16,type="o",cex=0.5,col=cols[i])		
}


plot(wkdat1$matAge,rep(0,length(wkdat1$matAge)),ylim=c(0,60),xlab="Maternal age",ylab="Total aborted",type="n")
for(i in 1:length(plot_cohorts))
{
	plotdat <- wkdat1[wkdat1$CohortID==plot_cohorts[i],]
	points(plotdat$matAge,plotdat$totalAborted,pch=16,type="o",cex=0.5,col=cols[i])		
}

write.csv(wkdat1,"ColonyDat_Weekly_July2016.csv",quote=F,row.names=F)


##### try to get age-specific fecundity and survival (*check other sources => put several papers in dropbox)

# age-specific fecundity

wkdat<-read.csv("ColonyDat_Weekly_July2016.csv")

survdat <- read.csv("AgeSpecificMortality_LRH.csv")

# add survival object

# survdat$SurvObj <- with(survdat, Surv(FlyAge, status==2)) => doesn't work, not sure how to deal with this properly

# use scatter smooth to predict the proportion of females alive at each day

m1 <- loess(PropFemAlive ~ FlyAge, data=survdat)

survdat$predProp <- predict(m1)
roundUp <- function(x){ifelse(round(abs(x-trunc(x)),1) == 0.5, trunc(x+0.5),round(x))}

survdat$noFem = survdat$predProp*12*48
survdat$predFem<- roundUp(survdat$noFem)

survdat$trayNo <- survdat$Tray

wkdat <- merge(wkdat, survdat[,c("trayNo","predFem")],all.x=TRUE, by = "trayNo")

wkdat$Fecundity <- (wkdat$totalPupaNo / wkdat$predFem)*0.5

# make data-frame of age-specific survival and fecundity (no. females produced = fecundity / 2)

x <- data.frame(with(wkdat, tapply(Fecundity, trayNo, mean)))
names(x) <- "meanFecundity"
x$trayNo <- rownames(x)

x <- x[,c(2,1)]

head(survdat) 

lifetable <- merge(x, survdat[,c("trayNo","predProp")], by = "trayNo", all.x=T)

all.ages <- ((12:1)-1)*7
all.ages <- all.ages[1:10]
lifetable$Age <- 0
wkdat$Age <- 0

for(i in 1:10) lifetable$Age[lifetable$trayNo==i] = all.ages[i]
for(i in 1:10) wkdat$Age[wkdat$trayNo==i] = all.ages[i]

lifetable <- lifetable[order(lifetable$Age),]

lifetable <- lifetable[,c(1,4,2,3)]
colnames(lifetable) <- c("trayNo", "Age", "meanFecundityFem", "propSurviving")

#write.csv(lifetable,"lifeTable_byWeek.csv",quote=F,row.names=F)

head(wkdat)

mortadd <- survdat[,c("trayNo","predProp")]
colnames(mortadd) <- c("trayNo", "Mortality")

wkdat <- merge(wkdat, mortadd, by="trayNo", all.x=T)

wkdat <- wkdat[order(wkdat$CohortID,wkdat$matAge),]
#write.csv(wkdat, "lifeTable_byCohort.csv",quote=F,row.names=F)

### check -- why no tray no = 10? 

# first look at total pupa mass

mm1 <- lmer(totalPupaMass ~ matAge + I(matAge^2) + (1|CohortID), data=wkdat)
mm1.1 <- lmer(totalPupaMass ~ matAge + I(matAge^2) + (matAge|CohortID), data=wkdat)
mm1.2 <- lmer(totalPupaMass ~ matAge + I(matAge^2) + I(matAge^3) + (matAge|CohortID), data=wkdat)
m1 <- lm(totalPupaMass ~ matAge + I(matAge^2), data=wkdat)

anova(mm1, mm2)
AIC(mm1.1);AIC(mm1.2)

# compare AIC - better to have random effect (& esp if slope)
# random effect accounts for 0.60/1.95 = 30% variance

wkdat$predTotalMass = predict.model(mm1.2, wkdat, use.random=TRUE)
wkdat$PRTotalMass = wkdat$predTotalMass - wkdat$totalPupaMass


preddat <- data.frame(matAge = seq(3,12,0.1))
preddat$Fecundity = seq(0.04, 0.4, length.out=91)
preddat$meanPupaMass = seq(23,32.5,length.out=91)
preddat$propAborted = seq(0,0.76,length.out=91)
preddat$predTotalMass <- 0
preddat$totalPupaMass <- 0

preddat$predTotalMass = predict.model(mm1.2,preddat,use.random=FALSE)


means <- tapply(wkdat$totalPupaMass, wkdat$matAge,mean)
sd <- tapply(wkdat$totalPupaMass, wkdat$matAge,sd)
se <- tapply(wkdat$totalPupaMass, wkdat$matAge,sd)/sqrt(tapply(wkdat$totalPupaMass, wkdat$matAge,length))

plot(preddat$matAge, preddat$predTotalMass, type="l", ylim=c(2,10), ylab="Total offspring mass (Mg)", xlab= "Maternal age (weeks)")
points(3:12,means,pch=16)
arrows(3:12,means-se,3:12,means+se,angle=90,length=0)



# now look at fecundity (rate/proportion) => here use the number of surviving females in measure

mm2 <- lmer(Fecundity ~ matAge + I(matAge^2) + (1|CohortID), data=wkdat)
mm2.1 <- lmer(Fecundity ~ matAge + I(matAge^2) + (matAge|CohortID), data=wkdat)
mm2.2 <- lmer(Fecundity ~ matAge + I(matAge^2) + I(matAge^3) + (matAge|CohortID), data=wkdat)

m2 <- lm(Fecundity ~ matAge + I(matAge^2), data=wkdat)
 
AIC(m2);AIC(mm2);AIC(mm2.1)
summary(mm2) # variance 0.0008797/0.0018832 = 47%

wkdat$predFecundity = predict.model(mm2.1, wkdat, use.random=TRUE)
preddat$predFecundity = predict.model(mm2.1, preddat, use.random=FALSE)

means <- tapply(wkdat$Fecundity, wkdat$matAge,mean)
sd <- tapply(wkdat$Fecundity, wkdat$matAge,sd)
se <- tapply(wkdat$Fecundity, wkdat$matAge,sd)/sqrt(tapply(wkdat$Fecundity, wkdat$matAge,length))

plot(preddat$matAge, preddat$predFecundity, type="l", ylim=c(0.13, 0.3), ylab="Fecundity (offspring per female)", xlab= "Maternal age (weeks)")
points(3:12,means,pch=16)
arrows(3:12,means-se,3:12,means+se,angle=90,length=0)



# now look at average offspring mass

mm3 <- lmer(meanPupaMass ~ matAge + I(matAge^2) + (1|CohortID), data=wkdat)
mm3.1 <- lmer(meanPupaMass ~ matAge + I(matAge^2) + (matAge|CohortID), data=wkdat)
mm3.2 <- lmer(meanPupaMass ~ matAge + I(matAge^2) + I(matAge^3) + (1|CohortID), data=wkdat)
m3 <- lm(meanPupaMass ~ matAge + I(matAge^2), data=wkdat)
 
AIC(m3);AIC(mm3);AIC(mm3.1);AIC(mm3.2)
summary(mm3) # variance = 0.71/1.06 = 71%

wkdat$predMeanMass = predict.model(mm3.2, wkdat, use.random=TRUE)
preddat$predMeanMass = predict.model(mm3.2, preddat, use.random=FALSE)

means <- tapply(wkdat$meanPupaMass, wkdat$matAge,mean)
sd <- tapply(wkdat$meanPupaMass, wkdat$matAge,sd)
se <- tapply(wkdat$meanPupaMass, wkdat$matAge,sd)/sqrt(tapply(wkdat$meanPupaMass, wkdat$matAge,length))

plot(preddat$matAge, preddat$predMeanMass, type="l", ylim=c(26.5,29.2), ylab=expression("Mean offspring mass ("~mu~"g)"), xlab= "Maternal age (weeks)")
points(3:12,means,pch=16)
arrows(3:12,means-se,3:12,means+se,angle=90,length=0)


# now abortion rate - but how to plot? 
# almost looks gaussian = range from 0 through to 60

mm4 <- lmer(totalAborted ~ matAge + I(matAge^2) + (1|CohortID), data=wkdat)
m4 <- lm(totalAborted ~ matAge + I(matAge^2), data=wkdat)
mm4.1 <- lmer(totalAborted ~ matAge + (1|CohortID), data=wkdat)

AIC(m4);AIC(mm4)
summary(mm4) # variance = 11/51 = 21%
anova(mm4,mm5) # linear fit is better


wkdat$propAborted = wkdat$totalAborted/wkdat$totalPupaNo
wkdat$logPropAborted = log(wkdat$propAborted + 0.0001)
preddat$logPropAborted = log(preddat$propAborted + 0.0001)

mm5 <- lmer(log(propAborted+0.0001) ~ matAge + I(matAge^2) + (1|CohortID), data=wkdat)
mm5.1 <- lmer(log(propAborted+0.0001) ~ matAge + I(matAge^2) + (matAge|CohortID), data=wkdat)
mm5.2 <- lmer(log(propAborted+0.0001) ~ matAge + I(matAge^2) + I(matAge^3) + (matAge|CohortID), data=wkdat)

AIC(mm5);AIC(mm5.1);AIC(mm5.2)

wkdat$predAborted = predict.model(mm5.1, wkdat, use.random=TRUE)
preddat$predAborted = predict.model(mm5.1, preddat, use.random=FALSE)

wkdat1 <- na.omit(wkdat)
means <- tapply(wkdat1$propAborted, wkdat1$matAge,mean)
sd <- tapply(wkdat1$propAborted, wkdat1$matAge,sd)
se <- tapply(wkdat1$propAborted, wkdat1$matAge,sd)/sqrt(tapply(wkdat1$propAborted, wkdat1$matAge,length))

plot(preddat$matAge, exp(preddat$predAborted), type="l", ylim=c(0.02,0.3), ylab="Proportion aborted", xlab= "Maternal age (weeks)")
points(3:12,means,pch=16)
arrows(3:12,means-se,3:12,means+se,angle=90,length=0)






wkdat$predAborted = predict.model(mm5, wkdat, use.random=TRUE)


length(unique(wkdat$CohortID))

par(mfrow=c(2,2))

# total mass
plot(wkdat$matAge, wkdat$predTotalMass, type="n", xlab="Maternal age (weeks)", ylab="Total offspring mass (Mg)")
for(i in 1:length(unique(wkdat$CohortID)))
{
	plotdat <- wkdat[wkdat$CohortID == unique(wkdat$CohortID)[i],]
	lines(plotdat$matAge, plotdat$predTotalMass, col=rgb(0.5,0.5,0.5,0.2),lwd=1)
}

# fecundity
plot(wkdat$matAge, wkdat$predFecundity, type="n", xlab="Maternal age (weeks)", ylab="Fecundity (offspring per female)")
for(i in 1:length(unique(wkdat$CohortID)))
{
	plotdat <- wkdat[wkdat$CohortID == unique(wkdat$CohortID)[i],]
	lines(plotdat$matAge, plotdat$predFecundity, col=rgb(0.5,0.5,0.5,0.2),lwd=1)
}

# mean mass
plot(wkdat$matAge, wkdat$predMeanMass, type="n", xlab="Maternal age (weeks)", ylab="Mean offspring mass (Mg)")
for(i in 1:length(unique(wkdat$CohortID)))
{
	plotdat <- wkdat[wkdat$CohortID == unique(wkdat$CohortID)[i],]
	lines(plotdat$matAge, plotdat$predMeanMass, col=rgb(0.5,0.5,0.5,0.2),lwd=1)
}

# total aborted
plot(wkdat$matAge, wkdat$predAborted, type="n", xlab="Maternal age (weeks)", ylab="Number offspring aborted")
for(i in 1:length(unique(wkdat$CohortID)))
{
	plotdat <- wkdat[wkdat$CohortID == unique(wkdat$CohortID)[i],]
	lines(plotdat$matAge, plotdat$predAborted, col=rgb(0.5,0.5,0.5,0.2),lwd=1)
}


# check resid 
par(mfrow=c(2,2)) # don't look too bad so think model fits are OK
boxplot(wkdat$predTotalMass - wkdat$totalPupaMass ~ wkdat$matAge)
boxplot(wkdat$predFecundity - wkdat$Fecundity ~ wkdat$matAge)
boxplot(wkdat$predMeanMass - wkdat$meanPupaMass ~ wkdat$matAge)
boxplot(wkdat$predAborted - wkdat$totalAborted ~ wkdat$matAge)


## calculate r for each week using uniroot

eulerlotka <- function(r,x,L,m) sum(L*m*exp(-r*x))-1

# first test with leslie ranson data:
x=seq(8,72,length=9)
m=c( 0.6504, 2.3939, 2.9727, 2.4662, 1.7043, 1.0815, 0.6683, 0.4286, 0.3000 )
L=c( 0.83349, 0.73132, 0.58809, 0.43343, 0.29277, 0.18126, 0.10285, 0.05348, 0.02549 )
r.range<- c(0, 10) 
eulerlotka <- function(r) sum(L * m * exp(-r * x)) - 1 
res <- uniroot(f = eulerlotka, interval = r.range, tol = 1e-8) 
res$root 

EL_dat <- c()

for(i in 1:length(unique(wkdat$CohortID)))
{
	solvdat <- wkdat[wkdat$CohortID == unique(wkdat$CohortID)[i],]
	
	if(length(solvdat[,1])<5) next
	else 
	{
		x = solvdat$Age
		m = solvdat$Fecundity
		L = solvdat$Mortality
		r = NA
		eulerlotka <- function(r) sum(L * m * exp(-r * x)) - 1 
		try(res <- uniroot(f = eulerlotka, interval = r.range, tol = 1e-8)) 
		r = res$root
			
		newline = c(i,r)
		
		EL_dat <- rbind(EL_dat, newline)
	}
		
}

EL_dat = data.frame(EL_dat)
names(EL_dat) = c("CohortID","r")

write.csv(EL_dat, "r_byCohort.csv", quote=F, row.names=F)

hist(EL_dat$r, breaks=50)

# now plot cohort data in ascending order of r... 

plot_Cohorts = EL_dat$CohortID[order(EL_dat$r)]

chkdat <- merge(wkdat,EL_dat,by="CohortID")
chkdat <- chkdat[order(chkdat$r),]

mod4 <- lmer(totalPupaMass ~ matAge + I(matAge^2) + (1|CohortID), data=chkdat)

chkdat <- fortify(mod4)

p1.0 <- ggplot(chkdat, aes(x= matAge, y= totalPupaMass)) +
  geom_point() +
  geom_line(aes(x= matAge, y=.fitted), size=.75) +
  facet_wrap(~CohortID, nrow=5)
p1.0

slopedf <- data.frame(as.matrix(ranef(mod4)))
names(slopedf) <- "InterceptVal"
slopedf$CohortID <- row.names(slopedf)



# Lyon talk 
wkdat$matAge2 = round(wkdat$Age/7)
with(wkdat[wkdat$matAge2>1,],boxplot(meanPupaMass ~ matAge2, col="darkred", outline=F,xlab="Maternal age (weeks)", ylab="Offspring weight (mg)") )
table(wkdat$matAge)


