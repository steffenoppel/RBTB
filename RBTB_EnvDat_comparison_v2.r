############################################################################################################
#######  TROPICBIRD FORAGING ENVIRONMENT COMPARISON ST HELENA - MADELAINE  ###################################################
############################################################################################################

## collaboration with Jacob Gonzalez Solis, Ngone Diop, Laura Zango and Anneale Beard
## analysis written by steffen.oppel@rspb.org.uk on 9 Jan 2017

#### REVISION TO INCLUDE MORE ENVIRONMENTAL VARIABLES STARTED ON 9 JULY 2018



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD PACKAGES AND PREVIOUSLY SAVED DATA FROM script "RBTB_EnvData_data_prep_v3.R"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(spdep)
library(lme4)
require(sp)
require(rgdal)
require(rgeos)
library(randomForest)
#library(randomForestSRC)     ## tried once on 24 July and immediately crashed R, so gave up
#library(ggRandomForests)
library(verification)
library(tidyr)
library(tidyverse)
library(dplyr)
library(EMbC)
library(data.table)
library(lubridate)
library(scales)
source("C:\\STEFFEN\\RSPB\\Statistics\\RF_partial_plot.r")
setwd("C:\\STEFFEN\\MANUSCRIPTS\\submitted\\RBTB_StHelena")

load("C:\\STEFFEN\\MANUSCRIPTS\\submitted\\RBTB_StHelena\\REV1_env_data_background.RData")
#load("A:\\MANUSCRIPTS\\submitted\\RBTB_StHelena\\RBTB_environmental_analysis_input.RData")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SUMMARISE AND COMPARE BACKGROUND VS TRACK LOCATIONS 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(RBTB_ENV)
head(BGRD_ENV)

dim(RBTB_ENV %>% filter(loc=="StHelena") %>% filter(state_EMBC=="ext_search"))


### set up loop over month, loc and variable
summary_t_test<-data.frame()

for (l in unique(BGRD_ENV$loc)){
  
  locTRACK<- RBTB_ENV %>% filter(loc==l)
  locBGR<- BGRD_ENV %>% filter(loc==l)
  months<- unique(locTRACK$month)
  
  for (m in months){
    
    ## subset the environmental variables to the month in question
    monthBG<- locBGR %>% filter(month==m)
    
    ## subset the bird tracking data to the month in question
    monthTRACK<- locTRACK %>% filter(month==m)
    
    for(v in unique(monthBG$Var)){
      
      ## subset the environmental variables to the variable in question
      BGx<- monthBG %>% filter(Var==v) %>% filter(!is.na(tmp.vec.long))
      
      ## subset the bird tracking data to the variable in question
      TRACKx<- monthTRACK %>% dplyr::select(v) %>% filter(!is.na(v))
      names(TRACKx)<-c("tmp.vec.long")
      
      ## combine data and perform simple t-test
      testdat<- t.test(x=TRACKx$tmp.vec.long, y=BGx$tmp.vec.long, alternative="two.sided")
      out_t_test<-data.frame(loc=l, month=m, Var=v,mean_BGR=testdat$estimate[2],mean_TRACK=testdat$estimate[1], p_val=testdat$p.value,diffLCL=testdat$conf.int[1],diffUCL=testdat$conf.int[2])
      summary_t_test<-rbind(summary_t_test, out_t_test)     
      
    }
    
  }
  
}



#fwrite(summary_t_test,"Env_data_comparison.csv")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Question 1:  Do environmental conditions differ between foraging and travelling locations 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

head(RBTB_ENV)

## SELECT ONLY FORAGING AND COMMUTING LOCATIONS ##
q1df<-RBTB_ENV %>% filter(state_EMBC %in% c("relocation","ext_search")) %>%   ### changed on 24 July 2018 to adjust for new data provided by Laura
  mutate(speed=sqrt(Wind_U^2+Wind_V^2)) %>%
  gather(variable, value, -ID_GPS_deployment, -loc, -DateTime, -Longitude, -Latitude, -Sexe, -Breeding_statut, -state_EMBC, -month) %>%
  mutate(forage=ifelse(state_EMBC=="ext_search",1,0))  %>%
  mutate(behav=ifelse(state_EMBC=="relocation","travelling","foraging")) %>%
  filter(!is.na(value))     ## remove all NA values

names(q1df)[c(1,2,6,7)]<-c("trip_id","Colony","Sex","breeding_status")

head(q1df)
dim(q1df)


## CREATE OUTPUT SUMMARY
Table3<-q1df %>% group_by(variable,Colony,behav) %>%
  summarise(n=length(value),mean=mean(value), sd=sd(value)) %>%
  mutate(p_var=0,p_int=0,beta=0,se=0,R2=0)






##### RUN ANALYSIS OVER ALL ENVIRONMENTAL VARIABLES ###########

variables<-unique(q1df$variable)

pdf("RBTB_env_histograms_forage_travel.pdf", width=6, height=8)
for(v in variables) {

m0<-glmer(forage~Colony+breeding_status+(1|trip_id), data=q1df[q1df$variable==v,], family=binomial)			## no effect of SST
m1<-glmer(forage~Colony+breeding_status+value+(1|trip_id), data=q1df[q1df$variable==v,], family=binomial)		## SST differs between forage and travel
m2<-glmer(forage~Colony*value+breeding_status+(1|trip_id), data=q1df[q1df$variable==v,], family=binomial)		## difference in SST between forage and travel varies by colony
pvals<-anova(m2,m1,m0)
Table3$p_var[Table3$variable==v]<-as.data.frame(pvals)[2,8]
Table3$p_int[Table3$variable==v]<-as.data.frame(pvals)[3,8]

modout<-summary(m1)
Table3$beta[Table3$variable==v]<-modout$coefficients[4,1]
Table3$se[Table3$variable==v]<-modout$coefficients[4,2]


# PLOT HISTOGRAM

histplotV<-q1df %>% filter(variable==v) %>%
  ggplot() + geom_histogram(aes(x=value, fill=behav), alpha=0.4, position="identity") +
  facet_wrap("Colony", ncol=1, scales = "fixed")+
  xlab(v) +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=16, color="black"),
        axis.title=element_text(size=18),
        strip.background=element_rect(fill="white", colour="black"),
        strip.text.x=element_text(size=18, color="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())
print(histplotV)

# Calculation of the variance in fitted values
VarF <- var(as.vector(fixef(m2) %*% t(m2@pp$X)))

# R2GLMM - conditional R2GLMM for full model
Table3$R2[Table3$variable==v]<-(VarF + VarCorr(m2)$trip_id[1])/(VarF + VarCorr(m2)$trip_id[1] + pi^2/3)


}
dev.off()


### FORMAT TABLE 3 FOR MANUSCRIPT ###
head(Table3)

Table3<-Table3 %>% arrange(variable, Colony, behav)

#fwrite(Table3,"Table3.csv")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Question 2:  Do environmental conditions differ between foraging locations and general background?
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

head(RBTB_ENV)
head(BGRD_ENV)
dim(BGRD_ENV)
dim(BGRD_ENV[is.na(BGRD_ENV$tmp.vec.long),])

BGRD_ENV$state_EMBC<-"background"

### FORMAT THE DATA FOR MODELLING INPUT ###

bgrd<- BGRD_ENV %>% dplyr::select(state_EMBC,lat, long, loc, month, Var, tmp.vec.long, date) %>%
  spread(key=Var, value=tmp.vec.long) %>%
  mutate(Colony=as.factor(loc)) %>%
  mutate(month=as.factor(month)) %>%
  mutate(speed=sqrt(Wind_U^2+Wind_V^2))
bgrd<- bgrd[complete.cases(bgrd),]
dim(bgrd)
head(bgrd)
plot(lat~long, bgrd)

### SAMPLE RANDOM BACKGROUND LOCATIONS TO SELECT ONLY 3-times AS MANY AS FORAGING LOCATIONS
## no longer necessary as complete cases are only 2870...
#q2df<-rbind(q2df, BGRD_ENV[sample(1:dim(BGRD_ENV)[1],size=dim(q2df)[1]*3),])



### NEED TO COMBINE THE TWO DATA FRAMES ###

q2df<-RBTB_ENV %>% filter(state_EMBC=="ext_search") %>%   ### changed on 24 July 2018 to adjust for new data provided by Laura
  gather(key="Var", value="tmp.vec.long", -ID_GPS_deployment, -loc, -DateTime, -Longitude, -Latitude, -Sexe, -Breeding_statut, -state_EMBC, -month) %>%
  dplyr::select(state_EMBC,Latitude, Longitude, loc, month, Var, tmp.vec.long, DateTime) %>%
  spread(key=Var, value=tmp.vec.long) %>%
  mutate(Colony=as.factor(loc)) %>%
  mutate(month=as.factor(month)) %>%
  mutate(speed=sqrt(Wind_U^2+Wind_V^2)) 
names(q2df)<-names(bgrd)
dim(q2df)
#q2df<- q2df[complete.cases(q2df),]         ### this eliminates ~500 locations from St Helena, so we use rfImpute instead
dim(q2df)

head(q2df)
q2df<-rbind(q2df,bgrd)
head(q2df)
dim(q2df)

plot(lat~long, q2df, col=as.numeric(state_EMBC))





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FITTING RANDOM FOREST MODELS FOR EACH STUDY AREA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##### SUBSET DATA FOR ST HELENA ##############

StHel<-q2df %>% filter(loc=="StHelena")
names(StHel)
StHel$state_EMBC<-droplevels(StHel$state_EMBC)
table(StHel$state_EMBC) ## quick assessment how balanced the data are

## IMPUTE MISSING DATA FOR ST HELENA ###
dim(StHel[complete.cases(StHel),])  ## we lose a lot of data, so we impute values rather than discard them
StHel<-rfImpute(state_EMBC~lat+long+BATHY+botTemp+CHLA+COLONY+EKE+MixLayDep+NPP+Salin+SEAMOUNT+slope+SSH+SST+SSTnight+TUNA+speed+month, data=StHel, iter=10,ntree=1500, mtry=3)
StHel$FORAGE<-ifelse(StHel$state_EMBC=="ext_search",1,0)  ## for numeric response/regressionn setting
table(StHel$FORAGE)


## CORRELATION BETWEEN PREDICTOR VARIABLES TO EXCLUDE HIGHLY CORRELATED VARIABLES
head(StHel)
CORMAT<-cor(StHel[,c(4:18)])
as.data.frame(CORMAT) %>% mutate(Var1=row.names(CORMAT)) %>%
  gather(key="variable", value="correlation",-Var1) %>%
  filter(abs(correlation)>0.7) %>%
  filter(correlation<1)



## FIT AND EVALUATE RANDOM FOREST MODEL IN CLASSIFICATION MODE
RF_STH<-randomForest(state_EMBC~BATHY+botTemp+COLONY+EKE+MixLayDep+CHLA+Salin+SEAMOUNT+slope+SSH+SST+TUNA+speed+month, data=StHel, ntree=1500, mtry=6, importance=T)
RF_STH
varImpPlot(RF_STH)
# partialPlot(RF_STH,StHel,"TUNA")
# partialPlot(RF_STH,StHel,"speed")
# partialPlot(RF_STH,StHel,"SST")


##   ASSESS RESIDUAL AUTOCORRELATION
### this is not possible with a classification response, hence we fit the same model with a regression response
REG_STH<-randomForest(FORAGE~BATHY+botTemp+COLONY+EKE+MixLayDep+CHLA+Salin+SEAMOUNT+slope+SSH+SST+TUNA+speed+month, data=StHel, ntree=1500, mtry=6, importance=T)
StHel$PRED<-predict(REG_STH, OOB=T)

# calculate residual Moran's I for each month (all months combined is nonsense, as same locations have similar values in all months)
MORANIres<-data.frame()
for (m in unique(StHel$month)){
  StHelm<-StHel[StHel$month==m,]
  #REG_STH<-randomForest(FORAGE~BATHY+botTemp+COLONY+EKE+MixLayDep+CHLA+Salin+SEAMOUNT+slope+SSH+SST+TUNA+speed+month, data=StHelm, ntree=1500, mtry=6, importance=T)
  #StHelm$PRED<-predict(REG_STH, OOB=T)
  SPECIES.coords <- as.matrix(cbind(StHelm$long, StHelm$lat))
  SPECIES.dists<-knn2nb(knearneigh(SPECIES.coords, k= 50, longlat=TRUE))		### adjust k= to the number of nearest neighbours
  coord.list <-make.sym.nb(SPECIES.dists) 
  coord.list <- nb2listw(coord.list,glist=NULL,style="W",zero.policy=FALSE) 
  mi<-moran.test(StHelm$FORAGE-StHelm$PRED, coord.list, randomisation=TRUE, zero.policy=TRUE, na.action=na.omit, spChk=NULL, adjust.n=TRUE)
  out<-data.frame(month=m, MoranI=mi$estimate[1],Moran_p=mi$p.value, rsq=mean(REG_STH$rsq))
  MORANIres<-rbind(MORANIres,out)
}

MORANIres




#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
##### SUBSET DATA FOR SENEGAL ##############

Sene<-q2df %>% filter(loc=="Senegal")
Sene$state_EMBC<-droplevels(Sene$state_EMBC)
table(Sene$state_EMBC)

## IMPUTE MISSING DATA ###
dim(Sene[complete.cases(Sene),])
Sene<-rfImpute(state_EMBC~lat+long+BATHY+botTemp+CHLA+COLONY+EKE+MixLayDep+NPP+Salin+SEAMOUNT+slope+SSH+SST+SSTnight+TUNA+speed+month, data=Sene, iter=10,ntree=1500, mtry=3)

Sene$FORAGE<-ifelse(Sene$state_EMBC=="ext_search",1,0)
table(Sene$FORAGE)


## CORRELATION BETWEEN PREDICTOR VARIABLES TO EXCLUDE 

CORMAT<-cor(Sene[,c(4:18)])
as.data.frame(CORMAT) %>% mutate(Var1=row.names(CORMAT)) %>%
  gather(key="variable", value="correlation",-Var1) %>%
  filter(abs(correlation)>0.7) %>%
  filter(correlation<1)


## FIT AND EVALUATE RANDOM FOREST MODEL IN CLASSIFICATION MODE
RF_SEN<-randomForest(state_EMBC~BATHY+CHLA+COLONY+EKE+MixLayDep+NPP+Salin+slope+SSH+SST+TUNA+speed+month, data=Sene, ntree=1500, mtry=8, importance=T,na.action=na.omit)
RF_SEN
varImpPlot(RF_SEN)
# partialPlot(RF_SEN,Sene,"TUNA")
# partialPlot(RF_SEN,Sene,"CHLA")
# partialPlot(RF_SEN,Sene,"COLONY")


##   ASSESS RESIDUAL AUTOCORRELATION
### this is not possible with a classification response, hence we fit the same model with a regression response
REG_SEN<-randomForest(FORAGE~BATHY+CHLA+COLONY+EKE+MixLayDep+NPP+Salin+slope+SSH+SST+TUNA+speed+month, data=Sene, ntree=1500, mtry=8, importance=T,na.action=na.omit)
Sene$PRED<-predict(REG_SEN, OOB=T)

# calculate residual Moran's I for each month (all months combined is nonsense, as same locations have similar values in all months)

MORANIres2<-data.frame()
for (m in unique(Sene$month)){
  StHelm<-Sene[Sene$month==m,]
  SPECIES.coords <- as.matrix(cbind(StHelm$long, StHelm$lat))
  SPECIES.dists<-knn2nb(knearneigh(SPECIES.coords, k= 50, longlat=TRUE))		### adjust k= to the number of nearest neighbours
  coord.list <-make.sym.nb(SPECIES.dists) 
  coord.list <- nb2listw(coord.list,glist=NULL,style="W",zero.policy=FALSE) 
  mi<-moran.test(StHelm$FORAGE-StHelm$PRED, coord.list, randomisation=TRUE, zero.policy=TRUE, na.action=na.omit, spChk=NULL, adjust.n=TRUE)
  out<-data.frame(month=m, MoranI=mi$estimate[1],Moran_p=mi$p.value, rsq=mean(REG_SEN$rsq))
  MORANIres2<-rbind(MORANIres2,out)
}

### CONCLUDE THAT SPATIAL AUTOCORRELATION IS NO PROBLEM BECAUSE IT HAS BEEN RESOLVED BY THE MODEL (=no spatial autocorrelation in residuals)
MORANIres2
MORANIres



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### TESTING WHETHER PREDICTIONS WORK ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


VALDAT<- Sene %>% mutate(SHpred=predict(RF_STH,newdata=Sene)) %>%
  dplyr::select(state_EMBC,FORAGE,SHpred) %>%
  filter(!is.na(SHpred))
1-(table(VALDAT$state_EMBC, VALDAT$SHpred)[1,1]/sum(table(VALDAT$state_EMBC, VALDAT$SHpred)[1,]))
# roc.area(VALDAT$FORAGE, VALDAT$SHpred)


VALDAT<- StHel %>% mutate(SENpred=predict(RF_SEN,newdata=StHel)) %>%
  dplyr::select(state_EMBC,FORAGE,SENpred) %>%
  filter(!is.na(SENpred))

table(VALDAT$state_EMBC, VALDAT$SENpred)
1-(table(VALDAT$state_EMBC, VALDAT$SENpred)[1,1]/sum(table(VALDAT$state_EMBC, VALDAT$SENpred)[1,]))

# roc.area(VALDAT$FORAGE, VALDAT$SENpred)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE AND PLOT VARIABLE IMPORTANCE  ##############
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##### St Helena ##########

VAR<-importance(RF_STH, type=1)			
IMP<-data.frame(variable=row.names(VAR), IncMSE=VAR[,1])
IMP<-IMP[order(IMP$IncMSE, decreasing=T),]  ## SORTED BY node homogeneity
IMP$rel_imp<-round((IMP$IncMSE/IMP$IncMSE[1])*100,2)

#write.table(IMP,"clipboard", row.names=F, sep="\t")


##### Senegal ##########

VARsen<-importance(RF_SEN, type=1)			
IMPsen<-data.frame(variable=row.names(VARsen), IncMSE=VARsen[,1])
IMPsen<-IMPsen[order(IMPsen$IncMSE, decreasing=T),]  ## SORTED BY node homogeneity
IMPsen$rel_imp<-round((IMPsen$IncMSE/IMPsen$IncMSE[1])*100,2)




pdf("RBTB_Variable_Importance_RF.pdf", width=9, height=13)
par(mfrow=c(2,1),mar=c(5,6,2,1))
barplot(IMP$rel_imp[10:1], horiz=T, names.arg=row.names(IMP)[10:1], xlim=c(0,100), las=1,xlab="Relative importance (%)", col='lightgray',main="St Helena")
barplot(IMPsen$rel_imp[10:1], horiz=T, names.arg=row.names(IMPsen)[10:1], xlim=c(0,100), las=1,xlab="Relative importance (%)", col='lightgray',main="Senegal")
dev.off()



TABLE4<- merge(IMP,IMPsen,by="variable", all=T)
#write.table(TABLE4,"clipboard", row.names=F, sep="\t")




# PLOT HISTOGRAM FOR ONE VARIABLE
hist(StHel$speed[StHel$forage==1],100, col="red", main="", xlab="wind speed (m/s)", ylim=c(0,200))
hist(StHel$speed[StHel$forage==0],100, main="", axes=F, add=T, ylim=c(0,200))








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##### PLOT VARIABLE RELATIONSHIPS  ##############
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
source("C:\\STEFFEN\\RSPB\\Statistics\\RF_partial_plot.r")

## for 2 variables per colony ##

SENOUT<-RF_partialPlot(data=Sene, RF=RF_SEN, variables=c("CHLA","COLONY"))
STHOUT<-RF_partialPlot(data=StHel, RF=RF_STH, variables=c("TUNA","Salin"))


#pdf("C:\\STEFFEN\\MANUSCRIPTS\\submitted\\RBTB_StHelena\\Fig6.pdf", width=13, height=13)



### COMBINE OUTPUT OF THE TWO IN ONE PLOT ###

plotdat<-(SENOUT[[1]] %>% mutate(Col="Madelaine"))
plotdat<-rbind(plotdat,(STHOUT[[1]] %>% mutate(Col="StHelena")))
head(plotdat)

## relabel the variables
plotdat$VARIABLE<-ifelse(plotdat$VARIABLE=="COLONY","distance to colony",plotdat$VARIABLE)
plotdat$VARIABLE<-ifelse(plotdat$VARIABLE=="TUNA","n Scombridae species",plotdat$VARIABLE)
plotdat$VARIABLE<-ifelse(plotdat$VARIABLE=="CHLA","chlorophyll a conc.",plotdat$VARIABLE)
plotdat$VARIABLE<-ifelse(plotdat$VARIABLE=="Salin","salinity",plotdat$VARIABLE)


out<- plotdat %>% mutate(pred.num=ifelse(pred=="ext_search",1,0)) %>%
  group_by(VARIABLE, VALUE,Col) %>%
  summarise(mean=mean(pred.num, na.rm=T), sd=sd(pred.num, na.rm=T)) %>%
  mutate(ucl=mean+0.5*sd, lcl=mean-0.5*sd)


## CREATE PLOT FOR ALL VARIABLES

ggplot(out) + geom_line(aes(x=VALUE, y=mean), linetype=1)+
  geom_line(aes(x=VALUE, y=lcl), linetype=2)+
  geom_line(aes(x=VALUE, y=ucl), linetype=2)+
  facet_wrap(~Col+VARIABLE, ncol=2, scales="free_x") +
  
  ## format axis ticks
  xlab("Range of variable value")+
  ylab("predicted proportion of foraging locations") +
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=18, color="black"),
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"),
        strip.text.y=element_text(size=18, color="black"),
        axis.title.y=element_text(margin=margin(0,20,0,0)), 
        strip.background=element_rect(fill="white", colour="black"))







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT PREDICTIONS FOR ST HELENA AND FOR MADELAINE BASED ON MODEL FROM OTHER COLONY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### loading maps in ggplot2 fails when tidyverse is loaded due to conflict with purrr::map


library(maps)
library(mapdata)
WorldData<- world_map <- maps::map("world", ".", exact = FALSE, plot = FALSE, fill = TRUE) %>% fortify()

RF_STH<-randomForest(FORAGE~BATHY+botTemp+COLONY+EKE+MixLayDep+CHLA+Salin+SEAMOUNT+slope+SSH+SST+TUNA+speed, data=StHel, ntree=1500, mtry=6, importance=T)
Sene$pred<-predict(RF_STH,newdata=Sene)

RF_MAD<-randomForest(FORAGE~BATHY+CHLA+COLONY+EKE+MixLayDep+NPP+Salin+slope+SSH+SST+TUNA+speed, data=Sene, ntree=1500, mtry=8, importance=T,na.action=na.omit)
StHel$pred<-predict(RF_MAD,newdata=StHel)
#SHplot<-STATS_INPUT[STATS_INPUT$Colony=="StHelena",]
#MADplot<-STATS_INPUT[STATS_INPUT$Colony=="Madelaine",]



ggplot()+geom_rect(data=Sene[Sene$state_EMBC =="background",], aes(xmin=long-0.4,ymin=lat-0.4,xmax=long+0.4,ymax=lat+0.4, fill = pred)) +
	scale_fill_gradient(name = 'Predicted prob. \n of foraging', low="white", high="red", na.value = "white", guide = "colourbar", limits=c(0, 1))+
	scale_x_continuous(limits=range(Sene$long))+
	scale_y_continuous(limits=range(Sene$lat))+
	theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())+
	theme(strip.text.y = element_text(size = 10, colour = "black"), strip.text.x = element_text(size = 10, colour = "black"))+
	geom_point(data=Sene[Sene$state_EMBC =="ext_search",], aes(x=long, y=lat))+
	geom_map(data=WorldData, map=WorldData, aes(map_id=region, x=long, y=lat), fill="darkgreen", colour="black", size=0.25)+
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
	  legend.text=element_text(size=16),
	  legend.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank()) 



ggplot()+geom_rect(data=StHel[StHel$state_EMBC =="background",], aes(xmin=long-0.4,ymin=lat-0.4,xmax=long+0.4,ymax=lat+0.4, fill = pred)) +
  scale_fill_gradient(name = 'Predicted prob. \n of foraging', low="white", high="red", na.value = "white", guide = "colourbar", limits=c(0, 1))+
  scale_x_continuous(limits=range(StHel$long))+
  scale_y_continuous(limits=range(StHel$lat))+
  theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())+
  theme(strip.text.y = element_text(size = 10, colour = "black"), strip.text.x = element_text(size = 10, colour = "black"))+
  geom_point(data=StHel[StHel$state_EMBC =="ext_search",], aes(x=long, y=lat))+
  geom_map(data=WorldData, map=WorldData, aes(map_id=region, x=long, y=lat), fill="darkgreen", colour="black", size=0.25)+
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank()) 



############ ASSESS RESIDUAL SPATIAL AUTOCORRELATION ########



#########################  TEST FOR AUTOCORRELATION USING MORAN's I AND GEARY's C##########################
################ taken from http://www.ats.ucla.edu/stat/r/faq/morans_i.htm   ################

MoranI<-data.frame(season=c(1:13))
MoranI$statistic<-0
MoranI$p<-0
for (i in 1:13){								### loop to calculate Morans I for each season separately
  b<-subset(TRAIN, Period==i)
  b<-subset(b, COUNT>0)
  SPECIES.coords <- as.matrix(cbind(b$LONG, b$LAT))
  SPECIES.dists<-knn2nb(knearneigh(SPECIES.coords, k= 50, longlat=TRUE))		### adjust k= to the number of nearest neighbours
  coord.list <-make.sym.nb(SPECIES.dists) 
  coord.list <- nb2listw(coord.list,glist=NULL,style="W",zero.policy=FALSE) 
  out<-moran.test(b$COUNT, coord.list, randomisation=TRUE, zero.policy=TRUE, na.action=na.omit, spChk=NULL, adjust.n=TRUE)
  out2<-geary.test(b$BASH, coord.list, randomisation=TRUE, zero.policy=TRUE, spChk=NULL, adjust.n=TRUE)
  MoranI[i,2]<-as.numeric(out$estimate[1])
  MoranI[i,3]<-out$p
  MoranI[i,4]<-as.numeric(out2$estimate[1])
  MoranI[i,5]<-out2$p
}
MoranI
write.table(MoranI, "Spatial_autocorrelation.csv", sep=",", row.names=F)







######################################## NOT CHANGED DURING REVISION #############################################






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT ENVIRONMENTAL VARIABLES FOR DIFFERENT BEHAVIOURAL STATES WITH RESPECT TO BATHYMETRY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

head(STATS_INPUT)

q3df<-STATS_INPUT[STATS_INPUT$Colony=="Madelaine" & STATS_INPUT$variable=='BATHY',]
q3df$Behaviour<-as.factor(c("rest","forage","travel","search")[q3df$STATE])

pdf("RBTB_Senegal_bathymetry_histogram.pdf", width=8, height=13)
ggplot(q3df, aes(x=value,y=(..count..)/sum(..count..))) +   facet_wrap(~breeding_status, ncol=1)+

    geom_histogram(aes(fill = Behaviour), data = q3df[q3df$STATE==2,], alpha = 0.4,binwidth = 100) + 
    geom_histogram(aes(fill = Behaviour), data = q3df[q3df$STATE==3,], alpha = 0.4,binwidth = 100) +
    geom_histogram(aes(fill = Behaviour), data = q3df[q3df$STATE==4,], alpha = 0.4,binwidth = 100) +
    xlab("Depth of ocean (m)") +
    ylab("Proportion of locations in behavioural state") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=14, color="black"),
        legend.text=element_text(size=14, color="black"),
        legend.title=element_text(size=18, color="black"),  
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())
dev.off()


ddply(q3df,c("breeding_status","Behaviour"), summarise,mean=mean(value, na.rm=T), sd=sd(value, na.rm=T))














#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT SST FOR ALL AFRICAN TRACKS BY MONTH
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

head(STATS_INPUT)

plotdf<-STATS_INPUT %>%
	filter(Colony=="Madelaine") %>%
	filter(variable=="SST") %>%
	mutate(Month=format(DateTime, format="%b"))

range(plotdf$value)

africa <- map_data("world")

ggplot(plotdf, aes(x=Longitude,y=Latitude, colour=value)) +   facet_wrap(~Month, ncol=3)+
	geom_point(pch=16, size=2)+

    	xlab("Longitude") +
    	ylab("Latitude") +

	scale_colour_gradient(name = 'SST (C)', low="white", high="red", guide = "colourbar", limits=c(19.75, 28))+

  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=14, color="black"),
        legend.text=element_text(size=14, color="black"),
        legend.title=element_text(size=18, color="black"),  
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())









#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT BACKGROUND 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

head(q2df)

plotdf<-q2df %>%
	filter(behav=="background") %>%
	filter(Colony=="Madelaine") %>%
	select(Longitude,Latitude,SST,Month) %>%							##night,NPP,SST,CHLA,BATHY,Month,speed
	gather(key=variable, value=value,-Latitude, -Longitude,-Month)

plotdf$value<-as.numeric(plotdf$value)
#head(plotdf)
#dim(plotdf)
#str(plotdf)
#unique(plotdf$Month)

ggplot()+geom_rect(data=plotdf, aes(xmin=Longitude-0.2,ymin=Latitude-0.2,xmax=Longitude+0.2,ymax=Latitude+0.2, fill = value)) +
	facet_wrap(~Month, ncol=3)+
	scale_fill_gradient(name = 'SST (C)', low="white", high="red", na.value = "white", guide = "colourbar", limits=c(20, 30))+
	theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())+
	theme(strip.text.y = element_text(size = 10, colour = "black"), strip.text.x = element_text(size = 10, colour = "black"))+

  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
	  legend.text=element_text(size=16),
	  legend.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank()) #+





plotdf<-q2df %>%
	filter(behav=="background") %>%
	filter(Colony=="StHelena") %>%
	select(Longitude,Latitude,SST,Month) %>%							##night,NPP,SST,CHLA,BATHY,Month,speed
	gather(key=variable, value=value,-Latitude, -Longitude,-Month)

plotdf$value<-as.numeric(plotdf$value)
range(plotdf$value, na.rm=T)
#head(plotdf)
#dim(plotdf)
#str(plotdf)
#unique(plotdf$Month)

ggplot()+geom_rect(data=plotdf, aes(xmin=Longitude-0.2,ymin=Latitude-0.2,xmax=Longitude+0.2,ymax=Latitude+0.2, fill = value)) +
	facet_wrap(~Month, ncol=2)+
	scale_fill_gradient(name = 'SST (C)', low="white", high="red", na.value = "white", guide = "colourbar", limits=c(16, 24))+
	theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())+
	theme(strip.text.y = element_text(size = 10, colour = "black"), strip.text.x = element_text(size = 10, colour = "black"))+

  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
	  legend.text=element_text(size=16),
	  legend.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank()) #+






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Question 4:  DOES FORAGING AND RESTING PROPORTION DIFFER BETWEEN DAY AND NIGHT?
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### MERGE THE DATA THAT WE NEED ###

names(tracks2)
q4df1<- tracks2 %>%
      select(c(1,3,4,6,8,9,20:22)) %>%
      mutate(Colony="StHelena") %>%
      mutate(ID=as.character(ID)) %>%
      mutate(ID=as.character(trip_id))

names(tracksSEN)
q4df2<- tracksSEN %>%
  select(c(1,3,4,7,6,2,19:21)) %>%
  mutate(Colony="Madelaine") %>%
  mutate(ID=as.character(ID)) %>%
  mutate(ID=as.character(trip_id))
  
names(q4df2)<-names(q4df1)

q4df<-bind_rows(q4df2,q4df1)
head(q4df)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE sunrise/sunset FOR EACH LOCATION AND DIVIDE INTO DAY AND NIGHT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sunr <- function(lon, lat, ts, dir) {
    sunriset(as.matrix(data.frame(lon, lat)), ts, POSIXct.out=TRUE, direction=dir)$time
  }

  q4df<-q4df %>% mutate(sunrise = sunr(Longitude, Latitude, DateTime, dir='sunrise')) %>%
  mutate(sunset = sunr(Longitude, Latitude, DateTime, dir='sunset')) %>%
  mutate(night=ifelse(DateTime>sunrise & DateTime<sunset,0,1)) %>%
  #mutate(moon=moon.illumination(DateTime)) %>%
  mutate(hour=hour(q4df$DateTime))
  
  q4df$forage<-ifelse(q4df$STATE==3,0,1)
  q4df$Behaviour<-ifelse(q4df$STATE==3,"travelling",ifelse(q4df$STATE==1,"resting","foraging"))

head(q4df)



########### PLOTTING HISTOGRAM FOR BEHAVIOUR #######

q4sum<-q4df %>%
  mutate(count=1) %>%
  group_by(Colony,night,Behaviour) %>%
  summarise(locs=sum(count))



ggplot(q4sum) + geom_bar(aes(y = locs, x = night, fill = Behaviour), position = "fill", stat="identity")+
  facet_wrap("Colony", ncol=1, scales = "fixed")+
  scale_y_continuous(labels = percent_format()) +
  scale_x_continuous(name="",breaks=c(0,1),labels = c('day','night')) +
  ylab("Proportion of locations") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=14, color="black", angle=45, vjust=0.5), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())



########### PLOTTING HISTOGRAM FOR BEHAVIOUR SPLIT BY BREEDING STATUS #######

q4sum<-q4df %>%
  mutate(count=1) %>%
  group_by(Colony,night,Behaviour,breeding_status) %>%
  summarise(locs=sum(count))



ggplot(q4sum) + geom_bar(aes(y = locs, x = night, fill = Behaviour), position = "fill", stat="identity")+
  facet_wrap(~Colony+breeding_status, ncol=2, scales = "fixed")+
  scale_y_continuous(labels = percent_format()) +
  scale_x_continuous(name="",breaks=c(0,1),labels = c('day','night')) +
  ylab("Proportion of locations") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=14, color="black", angle=45, vjust=0.5), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())
