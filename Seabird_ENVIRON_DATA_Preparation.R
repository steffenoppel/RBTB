############################################################################################################
#######  SEABIRD ENVIRONMENTAL DATA PREPARATION FOR TRACKING DATA  #########################################
############################################################################################################
#### based on ST HELENA TROPICBIRD MARINE DISTRIBUTION ANALYSIS
#### forked from RBTB_EnvData_summary2016.r
#### REVISED IN JULY 2018 TO INCLUDE ADDITIONAL ENVIRONMENTAL VARIABLES
#### changed approach to first summarise background environment
#### added bird locations to extract environmental data manually, as not all variables available on Movebank

#### downloaded data from 
## http://marine.copernicus.eu/services-portfolio/access-to-products/?option=com_csw&view=details&product_id=GLOBAL_REANALYSIS_PHY_001_025
## higher resolution than the data from Movebank

### UPDATE 16 JULY 2018: loaded all environmental rasters and added tuna species richness
### FINALISED 23 JULY 2018: included newly classified tracking data from Laura
### completely overhauled spatial arrangement as many grid cells ended up with NA for some variables



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD PACKAGES AND CUSTOM SCRIPTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(maps)
require(mapdata)
require(adehabitatHR)
require(adehabitatLT)
require(gpclib)
require(foreign)
require(maptools)
require(geosphere)
require(sp)
require(rgdal)
require(rgeos)
library(raster)
library(ncdf4)
library(data.table)
library(tidyverse)
library(lubridate)
library(rworldmap)
library(rworldxtra)  ## for high resolution world map

### SPECIFY DATES TO DOWNLOAD TIME SERIES FROM COPERNICUS SERVER
## 2015-08-16 13:00:00,2015-09-16 01:00:00,2015-10-16 13:00:00,2015-11-16 00:00:00,2015-12-16 12:00:00,2016-01-16 12:00:00,2016-02-15 12:00:00,2016-03-16 12:00:00

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~ 			LOAD ENVIRONMENTAL VARIABLES FROM COPERNICUS				   ~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

### added during revision on 9 July 2018
### script from Ascension adapted to St Helena and Senegal


setwd("C:\\STEFFEN\\MANUSCRIPTS\\submitted\\RBTB_StHelena")

ncdfs<-list.files(pattern=".nc")


envdat2<-data.frame()

for (f in c(2,3)){
  
  # open a NetCDF file
  ncin <- nc_open(ncdfs[f])
  #print(ncin)
  
  #lon <- ncvar_get(ncin, "lon")-360			### longitudes are in degrees EAST, so we subtract 360 to get negative values
  #lat <- ncvar_get(ncin, "lat", verbose = F)
  
  lon <- ncvar_get(ncin, "longitude")			
  lat <- ncvar_get(ncin, "latitude", verbose = F)
  
  t <- ncvar_get(ncin, "time")
  tunits <- ncatt_get(ncin, "time", "units")
  start<-as.POSIXct("1950-01-01 00:00:00")		### convert times to dates - NEEDS MANUAL STAT DATE DEFINITION
  dates<-start+(t*3600)
  nt <- dim(t)
  
  ## CREATE A FLAT FILE WITH MONTH AS COLUMN
  lonlatDF <- expand.grid(lon, lat, dates)
  names(lonlatDF)<-c('long','lat','date')
  ## rasterFromXYZ(lonlatDF)
  
  #### FOR EACH OF THE VARIABLES EXTRACT DATA ####
  ### different files and variable names for 2015 and 2016
  CMEMS_DATA<-data.frame()
    varnames<-c("zos","mlotst","so","bottomT","uo","vo")
    colnames<-c("SSH","MixLayDep","Salin","botTemp","CurVel_U","CurVel_V")
  
  # if (f %in% c(2,4)){
  # varnames<-c("zos","mlotst","so","bottomT","uo","vo")
  # colnames<-c("SSH","MixLayDep","Salin","botTemp","CurVel_U","CurVel_V")
  # }else{
  # varnames<-c("v","u","bottomT","salinity","mlp","ssh")
  # colnames<-c("CurVel_V","CurVel_U","botTemp","Salin","MixLayDep","SSH")
  # }
  
  for (v in 1:length(varnames)){
    dname<-varnames[v]
    #dname<-"sokaraml"
    
    tmp.array <- ncvar_get(ncin, dname)
    dlname <- ncatt_get(ncin, dname, "long_name")
    dunits <- ncatt_get(ncin, dname, "units")
    fillvalue <- ncatt_get(ncin, dname, "_FillValue")		### extract the FillValue (=NA)
    tmp.array[tmp.array == fillvalue$value] <- NA			### replaces random fill value with NA
    tmp.vec.long <- as.vector(tmp.array)
    if(length(dim(tmp.array))>3){tmp.vec.long <- as.vector(tmp.array[,,1,])}
    tmp.df02 <- data.frame(cbind(lonlatDF, tmp.vec.long))
    head(tmp.df02)
    
    # SUMMARISE ENVIRONMENTAL VARIABLES FOR EACH MONTH
    # monthly data - not averaged
    out <- tmp.df02 ##%>% group_by(long,lat) %>%
    #  summarise(mean=mean(tmp.vec.long, na.rm=T))
    out$Var<-colnames[v]
    out$loc<-ifelse(f<3,"Senegal","StHelena")
    CMEMS_DATA<-rbind(CMEMS_DATA,as.data.frame(out))
    
    
  }	# close loop over each variable
  
  head(CMEMS_DATA)
  dim(CMEMS_DATA)
  
  nc_close(ncin)
  
  envdat2<-rbind(envdat2,CMEMS_DATA)
  
} # close loop over each ncdf file




### Calculate Eddy Kinetic Energy

envdat2 <- envdat2 %>% spread(key=Var, value=tmp.vec.long) %>%
  mutate(EKE=(CurVel_U^2 + CurVel_V^2)/2) %>%
  gather(key="Var", value="tmp.vec.long", -date,-lat,-long,-loc)

head(envdat2)
str(envdat2)



### bbox for each study area ###

envdat2 %>% group_by(loc) %>% summarise(latl=min(lat), latu=max(lat), longl=min(long), longu=max(long))


### MAKE TIME PERIODS MATCH FOR DIFFERENT DATA SOURCES

envdat2<-envdat2 %>%
  mutate(date=as.Date(date)) %>%
  mutate(date=if_else(day(date)>15,date-1,date)) %>%
  mutate(month=month.abb[month(date)]) %>%
  arrange(date)

unique(envdat2$date)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD ENVIRONMENTAL VARIABLES FROM MOVEBANK
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

envdatBGR<-read.table("MOVEBANK_ENV_DATA.csv", header=T, sep=",")
names(envdatBGR)[12:20]<-c("SSTv","SSTnight","BATHY","NPP","SSTnightV","SST","Wind_U","CHLA","Wind_V")
head(envdatBGR)


### FORMAT MOVEBANK DATA TO SAME SHAPE AS OTHER ENV DATA

envdatBGR<-envdatBGR %>%
    filter(tag.local.identifier=="99999") %>%
    dplyr::select(location.long, location.lat,timestamp,SSTnight,NPP,SST,Wind_U,CHLA,Wind_V, BATHY) %>%
    gather(key="Var", value="tmp.vec.long", -timestamp,-location.lat,-location.long) %>%
    mutate(timestamp=ymd_hms(timestamp)) %>%
    mutate(loc=ifelse(location.lat>0,"Senegal","StHelena")) %>%
  mutate(date=as.Date(timestamp)) %>%
  mutate(date=if_else(day(date)>15,date-1,date)) %>%
  mutate(month=month.abb[month(date)]) %>%
  dplyr::select(location.long, location.lat,date,loc, Var,tmp.vec.long,month) %>%
  arrange(date)

unique(envdatBGR$date)
names(envdatBGR)<-names(envdat2)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MERGE ALL ENV DATA INTO A SINGLE DATA FRAME 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### COMBINE DATA ###
envdatBGR<-rbind(envdatBGR,envdat2)
head(envdatBGR)
dim(envdatBGR)

bgrd<-envdatBGR %>% dplyr::select(lat, long, loc, month, Var, tmp.vec.long, date) %>%
  spread(key=Var, value=tmp.vec.long)
head(bgrd)

unique(envdatBGR$long)
unique(envdatBGR$lat)

plot(lat~long, BGRD_ENV, pch=18, cex=0.1)			### looks like a fairly regular grid!


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE BATHYMETRIC GRADIENT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### CREATE A RASTER OBJECT FOR ST HELENA ###
## create a static raster for bathymetry
bathgrid<- envdatBGR %>% filter(Var=="BATHY") %>% 
  filter(loc=="StHelena") %>% dplyr::select(long,lat,tmp.vec.long) %>%
  group_by(long,lat) %>%
  summarise(BATHY=min(tmp.vec.long, na.rm=T))
head(bathgrid)

#bathras<-rasterFromXYZ(bathgrid, crs=CRS("+proj=longlat +ellps=WGS84"), digits=5)    ## only works for regular sized cells
bathsp <- SpatialPoints(bathgrid[,1:2], proj4string=CRS("+proj=longlat + datum=wgs84"))		### these are the raw location data
bathpix <- SpatialPixelsDataFrame(bathsp, data=bathgrid[,3], tolerance=0.01)
#delmepol <- as(delmepx,"SpatialPolygonsDataFrame")
SHbathras<- raster(bathpix)

### calculate slope
slopegridSH<-terrain(SHbathras, opt='slope', unit='degrees', neighbors=4)
slopegridSH




### CREATE A RASTER OBJECT FOR SENEGAL ###
## create a static raster for bathymetry
bathgrid<- envdatBGR %>% filter(Var=="BATHY") %>% 
  filter(loc=="Senegal") %>% dplyr::select(long,lat,tmp.vec.long) %>%
  group_by(long,lat) %>%
  summarise(BATHY=min(tmp.vec.long, na.rm=T))
head(bathgrid)

#bathras<-rasterFromXYZ(bathgrid, crs=CRS("+proj=longlat +ellps=WGS84"), digits=5)    ## only works for regular sized cells
bathsp <- SpatialPoints(bathgrid[,1:2], proj4string=CRS("+proj=longlat + datum=wgs84"))		### these are the raw location data
bathpix <- SpatialPixelsDataFrame(bathsp, data=bathgrid[,3], tolerance=0.01)
#delmepol <- as(delmepx,"SpatialPolygonsDataFrame")
SENbathras<- raster(bathpix)

### calculate slope
slopegridSEN<-terrain(SENbathras, opt='slope', unit='degrees', neighbors=4)
slopegridSEN





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE DISTANCE TO SEAMOUNTS AND COLONY AND CREATE RASTER OBJECTS 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seamounts<-fread("StHelena_seamounts.csv")
colony<-data.frame(loc=unique(envdatBGR$loc), long=c(-17.47163,-5.725067), lat=c(14.65414,-15.928676))


############ CALCULATE DISTANCE TO NEAREST SEAMOUNT #############

gridlocs<- envdatBGR %>% filter(Var=="BATHY") %>% 
  group_by(long,lat, loc) %>%
  summarise(SEAMOUNT=min(tmp.vec.long, na.rm=T))
seamountDIST<-spDists(as.matrix(gridlocs[,1:2]),as.matrix(seamounts[,4:3]), longlat=T)
gridlocs$SEAMOUNT<- apply(seamountDIST,1,min)

## CREATE A RASTER OBJECT FOR ST HELENA SEAMOUNT DISTANCE ##
tempsp <- SpatialPoints(gridlocs[gridlocs$loc=="StHelena",1:2], proj4string=CRS("+proj=longlat + datum=wgs84"))		### these are the raw location data
temppix <- SpatialPixelsDataFrame(tempsp, data=gridlocs[gridlocs$loc=="StHelena",4], tolerance=0.01)
SHseamras<- raster(temppix)

## CREATE A RASTER OBJECT FOR SENEGAL SEAMOUNT DISTANCE ##
tempsp <- SpatialPoints(gridlocs[gridlocs$loc=="Senegal",1:2], proj4string=CRS("+proj=longlat + datum=wgs84"))		### these are the raw location data
temppix <- SpatialPixelsDataFrame(tempsp, data=gridlocs[gridlocs$loc=="Senegal",4], tolerance=0.01)
SENseamras<- raster(temppix)


############ CALCULATE DISTANCE TO COLONY #############

gridlocs<- envdatBGR %>% filter(Var=="BATHY") %>% 
  group_by(long,lat, loc) %>%
  summarise(COLONY=min(tmp.vec.long, na.rm=T))
colDIST<-spDists(as.matrix(gridlocs[,1:2]),as.matrix(colony[,2:3]), longlat=T)
gridlocs$COLONY<- apply(colDIST,1,min)

## CREATE A RASTER OBJECT FOR ST HELENA COLONY DISTANCE ##
tempsp <- SpatialPoints(gridlocs[gridlocs$loc=="StHelena",1:2], proj4string=CRS("+proj=longlat + datum=wgs84"))		### these are the raw location data
temppix <- SpatialPixelsDataFrame(tempsp, data=gridlocs[gridlocs$loc=="StHelena",4], tolerance=0.01)
SHcolras<- raster(temppix)

## CREATE A RASTER OBJECT FOR SENEGAL COLONY DISTANCE ##
tempsp <- SpatialPoints(gridlocs[gridlocs$loc=="Senegal",1:2], proj4string=CRS("+proj=longlat + datum=wgs84"))		### these are the raw location data
temppix <- SpatialPixelsDataFrame(tempsp, data=gridlocs[gridlocs$loc=="Senegal",4], tolerance=0.01)
SENcolras<- raster(temppix)







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE RASTER OF TUNA SPECIES RICHNESS 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## downloaded on 13 JULY 2018 from
## https://www.aquamaps.org/CreateCSV.php?download_option=rich&taxon=Scombridae&user_session=11
# Species Richness Mapping Parameters for Scombridae (with a probability of occurrence > 0.5)
# Citation: Species Richness Map for Scombridae. 
# www.aquamaps.org, version of Aug. 2016. Web. Accessed 13 Jul. 2018.
# Cite AquaMaps itself as: Kaschner, K., K. Kesner-Reyes, C. Garilao, J. Rius-Barile, T. Rees, and R. Froese. 2016. AquaMaps: Predicted range maps for aquatic species. World wide web electronic publication, www.aquamaps.org, Version 08/2016. 

tuna<-fread("Tuna_richness.csv")
head(tuna)
names(tuna)[3]<-"TUNAsp"

## CREATE A RASTER OBJECT FOR TUNA SPECIES RICHNESS ##
tempsp <- SpatialPoints(tuna[,2:1], proj4string=CRS("+proj=longlat + datum=wgs84"))		### these are the raw location data
temppix <- SpatialPixelsDataFrame(tempsp, data=tuna[,3], tolerance=0.01)
TUNAras<- raster(temppix)

## RESCALE TUNA RASTER FOR ST HELENA ##
tempsp <- SpatialPoints(gridlocs[gridlocs$loc=="StHelena",1:2], proj4string=CRS("+proj=longlat + datum=wgs84"))		### these are the raw location data
temppix <- SpatialPixelsDataFrame(tempsp, data=data.frame(TUNA=raster::extract(TUNAras,tempsp, method='bilinear')), tolerance=0.01)
SHtunaras<- raster(temppix)

##  RESCALE TUNA RASTER FOR SENEGAL ##
tempsp <- SpatialPoints(gridlocs[gridlocs$loc=="Senegal",1:2], proj4string=CRS("+proj=longlat + datum=wgs84"))		### these are the raw location data
temppix <- SpatialPixelsDataFrame(tempsp, data=data.frame(TUNA=raster::extract(TUNAras,tempsp, method='bilinear')), tolerance=0.01)
SENtunaras<- raster(temppix)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD BIRD TRACKING DATA (with environmental data from Movebank)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### revised 23 July after Laura sent reclassified tracking data

setwd("C:\\STEFFEN\\MANUSCRIPTS\\submitted\\RBTB_StHelena")
RAW<-read.table("RBTB_tracks_classified_EmBC.csv", header=T, sep=",")

RBTB<- RAW %>% 
  mutate(DateTime=ymd_hms(Date_Time)) %>%
  mutate(year=as.numeric(year(DateTime))) %>%
  mutate(loc=ifelse(Colony=="Madeleine","Senegal","StHelena")) %>%
  dplyr::select(ID_GPS_deployment,loc, DateTime,Longitude, Latitude, Sexe, Breeding_statut, state_EMBC)


RBTB$month<-month.abb[month(RBTB$DateTime)]
RBTB$month<-ifelse(RBTB$month=="Aug","Sep",RBTB$month)  ## omit the few days at the end of august on St Helena
RBTB$month[RBTB$month=="Nov"]<-ifelse(day(RBTB$DateTime[RBTB$month=="Nov"])<15,"Oct","Dec")  ## omit the few days in November on St Helena



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# OVERLAY BIRD TRACKING DATA WITH BACKGROUND 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### THE TWO DATA FRAMES THAT WE WANT TO OVERLAY
head(RBTB)
head(envdatBGR)


### THE STATIONARY RASTERS THAT WE HAVE ARE
# bathymetry created above: SENbathras; SHbathras
# bathy slope created above: slopegridSEN; slopegridSH
# dist seamount: SENseamras, SHseamras
# dist colony: SENcolras, SHcolras
# tuna richness: TUNAras





###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
### CREATE RASTER STACK FROM FLAT DATA FRAME FOR EACH TIME INTERVAL #######
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########

### check that unique time periods exist in both data frames
month.abb[month(unique(envdatBGR$date))]
unique(RBTB$month)

RBTB_ENV<-data.frame()
BGRD_ENV<-data.frame()

for (l in unique(envdatBGR$loc)){
  
  locTRACK<- RBTB %>% filter(loc==l)
  locBGR<- envdatBGR %>% filter(loc==l)
  months<- unique(locTRACK$month)
  
  for (m in months){
  
  ## subset the environmental variables to the month in question
  monthBG<- locBGR %>% filter(month==m) %>% filter(Var!="BATHY") %>% select(long,lat,Var,tmp.vec.long) %>%
    spread(key=Var, value=tmp.vec.long)
  #head(monthBG) 
  
  ## create a raster brick for dynamic variables
  #monthras<-rasterFromXYZ(monthBG, crs=CRS("+proj=longlat +ellps=WGS84"), digits=5)
  tempsp <- SpatialPoints(monthBG[,1:2], proj4string=CRS("+proj=longlat + datum=wgs84"))		### these are the raw location data
  temppix <- SpatialPixelsDataFrame(tempsp, data=monthBG[,c(3,7)], tolerance=0.01)       ### need a second column, doesn't matter which one
  monthras<- raster(temppix)
  for (col in 4:15){
    temppix <- SpatialPixelsDataFrame(tempsp, data=monthBG[,c(col,7)], tolerance=0.01)
    monthras<- addLayer(monthras,raster(temppix))
  }
  
  ## subset the bird tracking data to the month in question
  monthTRACK<- locTRACK %>% filter(month==m)
  monthSPDF<-SpatialPointsDataFrame(coords=monthTRACK[,4:5], proj4string=CRS("+proj=longlat + datum=wgs84"), data = monthTRACK)
  #dim(monthSPDF)
  
  ## OVERLAY BIRD DATA WITH BACKGROUND GRID
  EXT1<-raster::extract(monthras,monthSPDF, method='bilinear')
  
  ## FOR STATIC RASTERS USE THE LOCATION-SPECIFIC RASTER
  if(l=="StHelena"){
    EXT2<-raster::extract(SHbathras,monthSPDF, method='bilinear')
    EXT3<-raster::extract(slopegridSH,monthSPDF, method='bilinear')
    EXT4<-raster::extract(SHcolras,monthSPDF, method='bilinear')
    EXT5<-raster::extract(SHseamras,monthSPDF, method='bilinear')
    monthras<- addLayer(monthras,SHbathras,slopegridSH,SHcolras,SHseamras,SHtunaras)
  }else{
    EXT2<-raster::extract(SENbathras,monthSPDF, method='bilinear')
    EXT3<-raster::extract(slopegridSEN,monthSPDF, method='bilinear')
    EXT4<-raster::extract(SENcolras,monthSPDF, method='bilinear')
    EXT5<-raster::extract(SENseamras,monthSPDF, method='bilinear')
    monthras<- addLayer(monthras,SENbathras,slopegridSEN,SENcolras,SENseamras,SENtunaras)
  }
  
  ## OVERLAY BIRD DATA WITH TUNA GRID
  EXT6<-raster::extract(TUNAras,monthSPDF, method='bilinear')
  
  ## WRITE OUT BIRD ENVIRONMENT DATA 
  RBTB_ENV<-rbind(RBTB_ENV,cbind(monthTRACK,EXT1,EXT2, EXT3, EXT4,EXT5, EXT6))


  ## WRITE OUT THE ENVIRONMENTAL DATA
  out<-as.data.frame(rasterToPoints(monthras))
  names(out)[1:2]<-c("long","lat")
  out<- out %>% gather(key="Var", value="tmp.vec.long", -lat,-long) %>%
   mutate(loc=l, month=m)
  out$date=(locBGR %>% filter(month==m) %>% select(date))[1,1]
  head(out)
  BGRD_ENV<-rbind(BGRD_ENV,out[,c(5,6,1,2,7,3,4)])
  
  }
  
}
names(RBTB_ENV)[23:27]<-c("BATHY","slope","COLONY","SEAMOUNT","TUNA")
head(RBTB_ENV)
dim(RBTB_ENV)








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REMOVE ON-LAND LOCATIONS AND SUMMARISE ENVIRONMENTAL VARIABLES FOR EACH STUDY AREA AND TIME PERIOD
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### USE semi_join FILTER FUNCTION


setwd("C:\\STEFFEN\\MANUSCRIPTS\\submitted\\RBTB_StHelena")
#### OVERLAY WITH LANDMASK BECAUSE OTHERWISE NONSENSE OCCURS 
dim(BGRD_ENV)
head(BGRD_ENV)
unique(BGRD_ENV$Var)
plot(lat~long, BGRD_ENV)


## get a world map
worldmap <- getMap(resolution = "high")
dim(worldmap)

## overlay and remove points that fall on land
uniquepoints<-BGRD_ENV %>% group_by(lat, long, loc) %>%
	summarise(bullshit=max(tmp.vec.long)) %>%
	select(lat, long, loc)


envSPDF<-SpatialPointsDataFrame(coords=uniquepoints[,2:1], proj4string=CRS("+proj=longlat +datum=WGS84"), data = uniquepoints)
proj4string(worldmap)<-proj4string(envSPDF)
envred<-envSPDF %over% worldmap
head(envred)
str(envred)
dim(BGRD_ENV)
INCLUDECOORDS<-as.data.frame(coordinates(envSPDF)[is.na(envred$ScaleRank),])
BGRD_ENV<-semi_join(BGRD_ENV,INCLUDECOORDS, by=c('lat','long'))
dim(BGRD_ENV)
plot(lat~long, BGRD_ENV)




######### PRODUCE SUMMARY OF BACKGROUND DATA AT SEA ##############

ENVmeans.sem <- BGRD_ENV %>% #filter(is.na(envred$ScaleRank)) %>%
  group_by(loc,month,Var) %>%
  summarise(mean=mean(tmp.vec.long, na.rm=T), sd=sd(tmp.vec.long, na.rm=T)) %>%
  mutate(lower=mean-(0.5*sd), upper=mean+(0.5*sd))

ENVmeans.sem 
unique(ENVmeans.sem$Var)
fwrite(ENVmeans.sem,"RBTB_environment_background.csv")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SUMMARISE ENVIRONMENTAL VARIABLES FOR BIRD TRACKING DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(RBTB_ENV)
#names(RBTB_ENV)[c(15,20,23:26)]<-c('CHLA1','NPP1','SST1','SSTnight1','Wind_U1','Wind_V1')

RBTB_ENV$SST


ENVmeans.track <- RBTB_ENV %>% dplyr::select(c(2,8:27)) %>%
  filter(!is.na(state_EMBC)) %>%
  gather(key="Var", value="tmp.vec.long",-loc,-month,-state_EMBC) %>%
  group_by(loc,month,Var,state_EMBC) %>%
  summarise(mean=mean(tmp.vec.long, na.rm=T), sd=sd(tmp.vec.long, na.rm=T)) %>%
  mutate(lower=mean-(0.5*sd), upper=mean+(0.5*sd))
ENVmeans.track 
unique(ENVmeans.track$Var)
fwrite(ENVmeans.track,"RBTB_environment_track_locs.csv")



save.image("REV1_env_data_background.RData")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TROUBLESHOOT WHY THERE ARE SO MANY NA IN BGRD_ENV DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(lat~long, BGRD_ENV)













