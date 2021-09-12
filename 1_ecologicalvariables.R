##This R script shows how we calculated the climatic and risk variables. Note that we are sharing this script for transparency purposes only--we have not supplied the Environmental_variables file because this file contains GPS points for each field site. In order to protect the privacy of participants, who for the most part come from small communities with few members, we have chosen not to publish these GPS points.


###install.packages("dismo")
###install.packages("raster")
##install.packages("dplyr")
##install.packages("usdm")
##install.packages("geosphere")
##install.packages("ncdf4")
##install.packages("tdyr")
library(dismo)
library(dplyr)
library(usdm)
library(geosphere)
library(ncdf4)
library(raster)
library(tidyverse)

########################
###CLIMATIC VARIABLES###
########################

##Import datasets with lat and long for each field site
Eco<-read.csv("Environmental_variables.csv",header=T)
colnames(Eco)[1]  <- "society"

###Import University of East Anglia climate data, and make each a dataset, from http://data.ceda.ac.uk/badc/cru/data/cru_ts/cru_ts_4.03/data
nc.pre <- nc_open("cru_ts4.03.1901.2018.pre.dat.nc")
nc.pre ##precipitation is in mm/month
pre <- brick("cru_ts4.03.1901.2018.pre.dat.nc", varname="pre")

nc.tmx<-nc_open("cru_ts4.03.1901.2018.tmx.dat.nc")
nc.tmx ##maximum monthly temp in degrees C
tmx<-brick("cru_ts4.03.1901.2018.tmx.dat.nc", varname="tmx")

nc.tmn<-nc_open("cru_ts4.03.1901.2018.tmn.dat.nc")
nc.tmn ##minimum monthly temp in degrees C
tmn<-brick("cru_ts4.03.1901.2018.tmn.dat.nc", varname="tmn")

##number of weather stations contributing to pre and temp data
stn.pre <- brick("cru_ts4.03.1901.2018.pre.dat.nc", varname="stn")
stn.tm<-brick("cru_ts4.03.1901.2018.tmx.dat.nc", varname="stn")

##Select the relevant GPS points, putting lon before lat
samples<-Eco%>%select(society,society_lon,society_lat)
samples<-samples%>%column_to_rownames(var="society")
par(mar=c(0,0,0,0))
plot(pre$X1901.01.16)
points(samples, pch=21)

##prepare renaming of year and month columns
years <- 1901:2018
month <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

##Using the extract function in the raster package, extract relevant data, rename variables, and make long
pre.sites <- data.frame(raster::extract(pre, samples, ncol=2))
names(pre.sites) <- paste(rep(years, each=12), rep(month, times=118), sep="_")
pre.sites$society <- Eco$society
pre_data<-gather(pre.sites,key,value,-society)
pre_data<-pre_data%>%rename(yearMonth=key,pre=value)
summary(pre_data)

tmx.sites <- data.frame(raster::extract(tmx, samples, ncol=2))
names(tmx.sites) <- paste(rep(years, each=12), rep(month, times=118), sep="_")
tmx.sites$society <- Eco$society
tmx_data<-gather(tmx.sites,key,value,-society)
tmx_data<-tmx_data%>%rename(yearMonth=key,tmx=value)
summary(tmx_data)

tmn.sites <- data.frame(raster::extract(tmn, samples, ncol=2))
names(tmn.sites) <- paste(rep(years, each=12), rep(month, times=118), sep="_")
tmn.sites$society <- Eco$society
tmn_data<-gather(tmn.sites,key,value,-society)
tmn_data<-tmn_data%>%rename(yearMonth=key,tmn=value)
summary(tmn_data)

stn.pre.sites <- data.frame(raster::extract(stn.pre, samples, ncol=2))
names(stn.pre.sites) <- paste(rep(years, each=12), rep(month, times=118), sep="_")
stn.pre.sites$society <- Eco$society
stn.pre_data<-gather(stn.pre.sites,key,value,-society)
stn.pre_data<-stn.pre_data%>%rename(yearMonth=key,stn.pre=value)
summary(stn.pre_data)

stn.tm.sites <- data.frame(raster::extract(stn.tm, samples, ncol=2))
names(stn.tm.sites) <- paste(rep(years, each=12), rep(month, times=118), sep="_")
stn.tm.sites$society <- Eco$society
stn.tm_data<-gather(stn.tm.sites,key,value,-society)
stn.tm_data<-stn.tm_data%>%rename(yearMonth=key,stn.tm=value)
summary(stn.tm_data)

##Merge the datasets together, and name them CRU_data
data<-list(stn.tm_data,stn.pre_data,tmn_data,tmx_data,pre_data)

CRU_data<- Reduce(
  function(x, y, ...) merge(x, y, by=c("society","yearMonth"),all=T,...), 
  data
)
head(CRU_data)
summary(CRU_data)

##make month and year separate columns
CRU_data<-CRU_data%>%separate(yearMonth,c("year","month"))

##add other variables
CRU_data<-merge(CRU_data,Eco,by="society")

##Split so that each society is its own dataframe within a list
Split<-split(CRU_data,CRU_data$society)

##Subset to include only years which preceded the end of data collection by 30 years inclusively
Split<-lapply(Split,subset,year<=endYear&year>=endYear-29)

##Calculate the number of weather stations contributing to each type of data, and the maximum absolute value for each field site
stn<-function(x){
  x%>%
    summarise(mean.stn.pre=mean(stn.pre),sd.stn.pre=sd(stn.pre),mean.stn.tm=mean(stn.tm),sd.stn.tm=sd(stn.tm))
}

mean.stn<-lapply(Split,stn)
mean.stn<-do.call(rbind.data.frame,mean.stn)
mean.stn <- tibble::rownames_to_column(mean.stn, "society")
mean.stn

##calculate month-by-month 30-year means, i.e. average all januaries, all februaries, etc.
means<-function(x){
  x%>%
    group_by(month)%>%
    summarise(pre=mean(pre),tmx=mean(tmx),tmn=mean(tmn))
}

Split<-lapply(Split,means)

##make the months in order
ordered<-function(x){
  x$month<-factor(x$month,levels=month)
  x<-x[order(x$month),]
  
}

Split<-lapply(Split,ordered)

##Calculate means and SD for Mikea because CV>100%
mean(Split[["Mikea"]]$pre) ##43.03
sd(Split[["Mikea"]]$pre) ##55.08

##Run function to extract the relevant biovars, the specifics of which can be found https://pubs.usgs.gov/ds/691/ds691.pdf
bio_results<-lapply(Split,function(x){
  prc<-x$pre
  tmn<-x$tmn
  tmx<-x$tmx
  bio<-biovars(prc,tmn,tmx)
  return(bio)
})

bio_results<-do.call(rbind.data.frame,bio_results)


##Rename columns to understand variables easier, based on https://cran.r-project.org/web/packages/dismo/dismo.pdf
bio<-bio_results%>%rename(meanAnnualTemp=bio1,meanDiurnalRange=bio2,isothermality=bio3,tempSeasonality=bio4,maxTempWarmestMonth=bio5,minTempColdestMonth=bio6,tempAnnualRange=bio7,meanTempWettestQ=bio8,meanTempDriestQ=bio9,meanTempWarmestQ=bio10,meanTempColdestQ=bio11,totalAnnualPrec=bio12,precWettestMonth=bio13,precDriestMonth=bio14,precSeasonality=bio15,precWettestQ=bio16,precDriestQ=bio17,precWarmestQ=bio18,precColdestQ=bio19)
bio <- tibble::rownames_to_column(bio, "society")

##Merge datasets, and rename columns for clarity
data<-list(Eco,mean.stn,bio)
Eco_total<- Reduce(
  function(x, y, ...) merge(x, y, by="society",all=T,...), 
  data
)

##Select variables of interest to check VIF for multicolinearity
Eco_test1<-Eco_total%>%select(tempSeasonality,meanAnnualTemp,totalAnnualPrec,precSeasonality,NPP)

##vif very high for temp seasonality and mean annual temperature (>10)
cor(Eco_test1)
vif(Eco_test1)

##without temp seasonality
Eco_test2<-Eco_total%>%select(meanAnnualTemp,totalAnnualPrec,precSeasonality,NPP)

##vif <3
cor(Eco_test2)
vif(Eco_test2)

###########################
###FOR DANGEROUS MAMMALS###
###########################

##Pantheria database  citation: Kate E. Jones, Jon Bielby, Marcel Cardillo, Susanne A. Fritz, Justin O'Dell, C. David L. Orme, Kamran Safi, Wes Sechrest, Elizabeth H. Boakes, Chris Carbone, Christina Connolly, Michael J. Cutts, Janine K. Foster, Richard Grenyer, Michael Habib, Christopher A. Plaster, Samantha A. Price, Elizabeth A. Rigby, Janna Rist, Amber Teacher, Olaf R. P. Bininda-Emonds, John L. Gittleman, Georgina M. Mace, and Andy Purvis. 2009. PanTHERIA: a species-level database of life history, ecology, and geography of extant and recently extinct mammals. Ecology 90:2648. https://ecologicaldata.org/wiki/pantheria

##This database includes species specific info on: max & min lat and long, n/km2 density, and size
##download database from https://figshare.com/articles/dataset/Full_Archive/3531875

pantheria <- read.csv("PanTHERIA_1-0_WR93_Aug2008.csv")

##subset pantheria to relevant variables
summary(pantheria)
pantheria<-pantheria %>% select(MSW93_Order, MSW93_Family ,MSW93_Genus ,MSW93_Species , X5.1_AdultBodyMass_g, X21.1_PopulationDensity_n.km2, X26.2_GR_MaxLat_dd, X26.3_GR_MinLat_dd, X26.5_GR_MaxLong_dd, X26.6_GR_MinLong_dd)

##make a separate dataframe for each society
society<-Eco%>%select(society,society_lon,society_lat)
society<-split(society, society$society)

##make a function and merge each society to the pantheria database
merge_func <- function(x,y){merge(x, y)}
society_pantheria<-lapply(society, merge_func, pantheria)

##make a function and to retain only all animals ranges that overlap with the field sites
between<-function(x){
x$lon_match <-ifelse(x$society_lon<=x$X26.5_GR_MaxLong_dd & x$society_lon>= x$X26.6_GR_MinLong_dd,1,0)
x$lat_match<-ifelse(x$society_lat<=x$X26.2_GR_MaxLat_dd & x$society_lat>=x$X26.3_GR_MinLat_dd,1,0)
x$range<-ifelse(x$lon_match==1&x$lat_match==1,1,0)
x
}

society_range<-lapply(society_pantheria,between)

##now subset each dataframe according to range
society_range<-lapply(society_range,subset,range==1)

##now subset according to our species of interest
unique(pantheria$MSW93_Order)

##animals of interest are large carnivores, hippos, and elephants
society_range<-lapply(society_range,subset,MSW93_Order=="Proboscidea" |MSW93_Order=="Carnivora"|MSW93_Order=="Artiodactyla"& MSW93_Family=="Hippopotamidae")

##Select animals over 50000g (i.e. 50kg)
society_range<-lapply(society_range,subset,X5.1_AdultBodyMass_g>=50000)

##export animal lists by field sites for ethnographers to check over
for(i in names(society_range)){
  write.csv(society_range[[i]], paste0(i,"_dangerousAnimals.csv"))
}

##Import fieldworker corrected mammal list, select relevant pantheria variables
DMC<-read.csv("DangerousMammalsCorrected.csv",header=T)
PS<-pantheria%>%select(MSW93_Genus ,MSW93_Species ,X21.1_PopulationDensity_n.km2)
DMC_c<-merge(DMC,PS,by=c("MSW93_Genus","MSW93_Species"))


##Make table for supplement
DMC_c$Name<-paste(DMC_c$MSW93_Genus,DMC_c$MSW93_Species,sep=" ")
table<-data.frame(unclass(table(DMC_c$Name,DMC_c$society)))
table <- tibble::rownames_to_column(table, "Name")
density<-DMC_c%>%group_by(Name)%>%summarise(density=unique(X21.1_PopulationDensity_n.km2))
table<-merge(density,table,by="Name")
write.csv(table,"MammalTable.csv")

##sum the densities by society
densenkm2<-DMC_c%>%
  group_by(society)%>%
  summarise(sum_density=sum(X21.1_PopulationDensity_n.km2))

##merge with eco_total and select variables
Eco_total<-merge(Eco_total,densenkm2,by="society", all=TRUE)
Eco_total[is.na(Eco_total)] <- 0

#########################
###FOR VENOMOUS SNAKES###
#########################

##Medically relevant snake species per each country were extracted from WHO database on venomous snakes (https://apps.who.int/bloodproducts/snakeantivenoms/database/) and compared with integrated list available at https://github.com/joshlongbottom/snakebite/blob/master/snake_list.csv from his paper https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(18)31224-8/fulltext
##Per each species, it was visually checked if the area range overlapped with the location of the fieldsites.
##The resulting list of species per fieldsite was double checked by the ethnographers. 
##The total number of medically relevant venomous species of snakes was then calculated.

##In few cases, it was reported only the indication of presence of a whole genus(e.g. Micrurus spp), so that the total number of species belonging to the genus is not known. 
##This concerned especially coral snakes, for which species recognition is particularly complicated.
##The entry was considered as a single species, following the intuition that if multiple species of the genus were abundant, they would have been identified separatly.

##data pulled per country can be found here
##snake<-read_csv(url("https://raw.githubusercontent.com/joshlongbottom/snakebite/master/snake_list.csv"))

##import corrected snake file
snakes<-read.csv("snakes.csv",header=T)

##Make Table
snake_present<-subset(snakes,present==TRUE)
table<-data.frame(unclass(table(snake_present$species,snake_present$society)))
write.csv(table,"SnakeTable.csv")

##calculate snake sums
snakes$bin<-ifelse(snakes$present==TRUE,1,0)
snake_count<-snakes%>%
  group_by(society)%>%
  summarise(snake_count=sum(na.omit(bin)))

##merge with dataset
Eco_total<-merge(Eco_total,snake_count,by="society")

##Select relevant variables
Eco_total<-Eco_total%>%select(society,NPP,meanAnnualTemp,tempSeasonality,totalAnnualPrec,precSeasonality,dens_nkm2,snake_count)
write.csv(Eco_total,"Eco_total.csv")