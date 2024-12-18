#### PART 1 Presence+Absence point generation: Annually updated estimates of ebolavirus spillover potential accounting for changes to forests and human populations  ####


#' README:
#' The code below outlines the process of data preparation and cleaning to conduct
#' the analysis. It requires publicly available datasets that 
#' are too large to house on the github repository. 
#' Thus it will not be possible to run this script without first downloading necessary 
#' data and saving on the local drive.
#' We have made the final postprocessed dataset for analysis available 
#' on the github repository that is generated through this data 
#' preparation process. 


library(raster)
library(tidyverse)
library(sf)
library(sp)
library(rgdal)
library(tmap)
library(dismo)
library(exactextractr)
library(rgeos)
library(gbm)
library(seegSDM)

getwd()
setwd("C:/Users/pwv0/OneDrive - CDC/Filo Defor Data")

#### Prep presence points ####
evdpresence <- read.csv("EVD Spillover Linelist Coordinates/EVD_spill_csv.csv")
evdpresence <- evdpresence %>%
  mutate(x = CONFIRMED.EXACT.LONG) %>%
  mutate(y = CONFIRMED.EXACT.LAT) %>%
  mutate(long=x) %>%
  mutate(lat=y) %>%
  filter(Include != "No") %>%
  dplyr::select(ID,long,lat,x,y,Year,Strain)

evdpresencecoords <- evdpresence


#### Prep absence points ####
#' read in population raster that i will use for weighting the random absence points
pd <- raster("Covariate Data/Final Covariates/PopDensityLandScan.tif")
logpd <- log(pd+1)+1

# generate absence points
bg_logpd <- bgSample(logpd,
                  n = 10000,
                  prob = T, # turn off if you dont want to weight
                  replace = F,
                  spatial = T)

bg_logpd <- as.data.frame(bg_logpd)

# randomly assign a year 2001-2021 to each absence point
years<-seq.int(2001,2021,1)
Year<-NULL
for (i in 1:10000) {
  y<-sample(years,1)
  Year<-rbind(Year,y)
}

bg_logpdyr <-cbind(bg_logpd,Year)
bg_logpdyr <- bg_logpdyr %>% mutate(long=x,lat=y)
bg_logpdyr$ID <- 1:10000

bg_logpdyr <- bg_logpdyr %>%
  mutate(Strain = "Absence_logpdw") %>%
  select(ID,long,lat,x,y,Year,Strain)


# MERGE Pres and Abs
fulldf <- rbind(telfordprescoords,bg_logpdyr)
fulldf$NumID <- 1:nrow(fulldf)





#### Prep prediction grid ####
#' create grid for extraction using coordinates from 1x1km resolution raster
#' extract covariate values at each grid point for 2021 and 2022
grid <- raster("Covariate Data/BioClim/OnexOneraster.tif")
grid <- rasterToPoints(grid)
grid <- as.data.frame(grid)

# export and make buffers in QGIS without dissolve
write.csv(grid,"Covariate Data/Prediction Grid Extracted Vars/predictiongridcoords.csv")




