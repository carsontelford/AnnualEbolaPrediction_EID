#### PART 3: Prediction over grid of equatorial Africa ####


#' README:
#' The code below outlines the process of using the 100 fitted ensemble models
#' to make 100 separate predictions of the odds of spillover in each grid cell
#' across equatorial Africa, and then consolidating all of the ensemble predictions
#' to generate a single estimate. 
#' 
#' This process requires large datasets to be saved, including the fitted results from 
#' the 100 models within the ensemble and cannot be saved on the public repository due 
#' to file size restrictions. Therefore, the script below outlines the process for making the predictions, but
#' the data to do so must be saved locally. 
#' 
#' At each of the grid cells where we make predictions, we create donut buffers
#' and extract the values of each predictor using google earth engine. We then import
#' the extracted values to this program to generate raster grid of covariate values.
#' 
#' Scripts for covariate extraction in GEE are provided in the public repository.





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
library(gbm)
library(dismo)
library(maptools)
library(raster)
require(fields)
require(parallel)
require(snowfall)
library(seegSDM)


#### Prediction grid ####
#' read in data extracted from google earth engine

##### 2022 pred grid #####
##### join all covars, create FL loss vars, convert to rasters #####
fccovarsin <- read.csv("Covariate Data/Prediction Grid Extracted Vars/fccovars.csv")
othercovarsin <- read.csv("Covariate Data/Prediction Grid Extracted Vars/othercovars.csv")  
pccovars2022 <- read.csv("Covariate Data/Prediction Grid Extracted Vars/pccovars2022.csv")
pccovars2021 <- read.csv("Covariate Data/Prediction Grid Extracted Vars/pccovars2021.csv")
pccovars2020 <- read.csv("Covariate Data/Prediction Grid Extracted Vars/pccovars2020.csv")
fragcovar <- read.csv("Covariate Data/Prediction Grid Extracted Vars/fragcovars2022.csv")  
fragcovar.8 <- read.csv("Covariate Data/Prediction Grid Extracted Vars/fragcovars80pct2021.csv")

allcovars <- fccovarsin %>%
  left_join(fragcovar,by= c("X","x","y")) %>%
  left_join(fragcovar.8,by= c("X","x","y")) %>%
  left_join(othercovarsin, by= c("X","x","y")) %>%
  left_join(pccovars2022, by = c("X","x","y")) %>%
  left_join(pccovars2021, by = c("X","x","y")) %>%
  left_join(pccovars2020, by = c("X","x","y"))

allcovars <- allcovars %>%
  mutate(
    frag10km.7 = fragsy.7_10km*100,
    frag25km.7 = fragsy.7_25km*100,
    frag50km.7 = fragsy.7_50km*100,
    frag100km.7 = fragsy.7_100km*100,
    frag150km.7 = fragsy.7_150km*100,
    
    frag10km.8 = fragsy.8_10km*100,
    frag25km.8 = fragsy.8_25km*100,
    frag50km.8 = fragsy.8_50km*100,
    frag100km.8 = fragsy.8_100km*100,
    frag150km.8 = fragsy.8_150km*100,
    
    flsy10km = (  fc2021_10km-fc2022_10km)*10000,
    flsy25km = (  fc2021_25km-fc2022_25km)*10000,
    flsy50km = (  fc2021_50km-fc2022_50km)*10000,
    flsy100km = ( fc2021_100km-fc2022_100km)*10000,
    flsy150km = ( fc2021_150km-fc2022_150km)*10000,
    
    fl1yp10km = (  fc2020_10km-fc2021_10km)*10000,
    fl1yp25km = (  fc2020_25km-fc2021_25km)*10000,
    fl1yp50km = (  fc2020_50km-fc2021_50km)*10000,
    fl1yp100km = ( fc2020_100km-fc2021_100km)*10000,
    fl1yp150km = ( fc2020_150km-fc2021_150km)*10000,
    
    fl2yp10km = (  fc2019_10km-fc2020_10km)*10000,
    fl2yp25km = (  fc2019_25km-fc2020_25km)*10000,
    fl2yp50km = (  fc2019_50km-fc2020_50km)*10000,
    fl2yp100km = ( fc2019_100km-fc2020_100km)*10000,
    fl2yp150km = ( fc2019_150km-fc2020_150km)*10000,
    
    logpcsy10km = log( pc2022_10km+1),
    logpcsy25km = log( pc2022_25km+1),
    logpcsy50km = log( pc2022_50km+1),
    logpcsy100km = log(pc2022_100km+1),
    logpcsy150km = log(pc2022_150km+1),
    
    logpcfc10km = logpcsy10km*fc2022_10km,
    logpcfc25km = logpcsy25km*fc2022_25km,
    logpcfc50km = logpcsy50km*fc2022_50km,
    logpcfc100km = logpcsy100km*fc2022_100km,
    logpcfc150km = logpcsy150km*fc2022_150km,
    
    fcsy10km = fc2022_10km*100,
    fcsy25km = fc2022_25km*100,
    fcsy50km = fc2022_50km*100,
    fcsy100km = fc2022_100km*100,
    fcsy150km = fc2022_150km*100,
    
    pet_10km = pet_10km/1000,
    tempseason_10km = tempseason_10km/10) %>%
  dplyr::rename(pet = pet_10km,
                tempseason = tempseason_10km,
                precipseason = precipseason_10km,
                elevation = elevation_10km,
                ntlst = ntlst_10km) 


# create individual dfs
fcsy10km <- allcovars %>% dplyr::select(x,y,fcsy10km) 
fcsy25km <- allcovars %>% dplyr::select(x,y,fcsy25km)
fcsy50km <- allcovars %>% dplyr::select(x,y,fcsy50km)
fcsy100km <- allcovars %>% dplyr::select(x,y,fcsy100km)
fcsy150km <- allcovars %>% dplyr::select(x,y,fcsy150km)

flsy10km <- allcovars %>% dplyr::select(x,y,flsy10km) 
flsy25km <- allcovars %>% dplyr::select(x,y,flsy25km)
flsy50km <- allcovars %>% dplyr::select(x,y,flsy50km)
flsy100km <- allcovars %>% dplyr::select(x,y,flsy100km)
flsy150km <- allcovars %>% dplyr::select(x,y,flsy150km)

fl1yp10km <- allcovars %>% dplyr::select(x,y,fl1yp10km)
fl1yp25km <- allcovars %>% dplyr::select(x,y,fl1yp25km)
fl1yp50km <- allcovars %>% dplyr::select(x,y,fl1yp50km)
fl1yp100km <- allcovars %>% dplyr::select(x,y,fl1yp100km)
fl1yp150km <- allcovars %>% dplyr::select(x,y,fl1yp150km)

fl2yp10km <- allcovars %>% dplyr::select(x,y,fl2yp10km)
fl2yp25km <- allcovars %>% dplyr::select(x,y,fl2yp25km)
fl2yp50km <- allcovars %>% dplyr::select(x,y,fl2yp50km)
fl2yp100km <- allcovars %>% dplyr::select(x,y,fl2yp100km)
fl2yp150km <- allcovars %>% dplyr::select(x,y,fl2yp150km)

frag10km.7 <- allcovars %>% dplyr::select(x,y,frag10km.7)
frag25km.7 <- allcovars %>% dplyr::select(x,y,frag25km.7)
frag50km.7 <- allcovars %>% dplyr::select(x,y,frag50km.7)
frag100km.7 <- allcovars %>% dplyr::select(x,y,frag100km.7)
frag150km.7 <- allcovars %>% dplyr::select(x,y,frag150km.7)

frag10km.8 <- allcovars %>% dplyr::select(x,y,frag10km.8)
frag25km.8 <- allcovars %>% dplyr::select(x,y,frag25km.8)
frag50km.8 <- allcovars %>% dplyr::select(x,y,frag50km.8)
frag100km.8 <- allcovars %>% dplyr::select(x,y,frag100km.8)
frag150km.8 <- allcovars %>% dplyr::select(x,y,frag150km.8)

logpcsy10km <- allcovars %>% dplyr::select(x,y,logpcsy10km)
logpcsy25km <- allcovars %>% dplyr::select(x,y,logpcsy25km)
logpcsy50km <- allcovars %>% dplyr::select(x,y,logpcsy50km)
logpcsy100km <- allcovars %>% dplyr::select(x,y,logpcsy100km)
logpcsy150km <- allcovars %>% dplyr::select(x,y,logpcsy150km)

elevation <- allcovars %>% dplyr::select(x,y,elevation)
ntlst <- allcovars %>% dplyr::select(x,y,ntlst)
tempseason <- allcovars %>% dplyr::select(x,y,tempseason)
precipseason <- allcovars %>% dplyr::select(x,y,precipseason)
pet <- allcovars %>% dplyr::select(x,y,pet)
logpcfc10km <- allcovars %>% dplyr::select(x,y,logpcfc10km)

# convert to rasters
coordinates(fcsy10km) <- ~x+y
coordinates(fcsy25km) <- ~x+y
coordinates(fcsy50km) <- ~x+y
coordinates(fcsy100km) <- ~x+y
coordinates(fcsy150km) <- ~x+y

coordinates(flsy10km) <- ~x+y
coordinates(flsy25km) <- ~x+y
coordinates(flsy50km) <- ~x+y
coordinates(flsy100km) <- ~x+y
coordinates(flsy150km) <- ~x+y

coordinates(fl1yp10km) <- ~x+y
coordinates(fl1yp25km) <- ~x+y
coordinates(fl1yp50km) <- ~x+y
coordinates(fl1yp100km) <- ~x+y
coordinates(fl1yp150km) <- ~x+y

coordinates(fl2yp10km) <- ~x+y
coordinates(fl2yp25km) <- ~x+y
coordinates(fl2yp50km) <- ~x+y
coordinates(fl2yp100km) <- ~x+y
coordinates(fl2yp150km) <- ~x+y

coordinates(frag10km.7) <- ~x+y
coordinates(frag25km.7) <- ~x+y
coordinates(frag50km.7) <- ~x+y
coordinates(frag100km.7) <- ~x+y
coordinates(frag150km.7) <- ~x+y

coordinates(frag10km.8) <- ~x+y
coordinates(frag25km.8) <- ~x+y
coordinates(frag50km.8) <- ~x+y
coordinates(frag100km.8) <- ~x+y
coordinates(frag150km.8) <- ~x+y

coordinates(logpcsy10km) <- ~x+y
coordinates(logpcsy25km) <- ~x+y
coordinates(logpcsy50km) <- ~x+y
coordinates(logpcsy100km) <- ~x+y
coordinates(logpcsy150km) <- ~x+y

coordinates(elevation) <- ~x+y
coordinates(ntlst) <- ~x+y
coordinates(tempseason) <- ~x+y
coordinates(precipseason) <- ~x+y
coordinates(pet) <- ~x+y
coordinates(logpcfc10km) <- ~x+y

# gridded
gridded(fcsy10km) <- TRUE
gridded(fcsy25km) <- TRUE
gridded(fcsy50km) <- TRUE
gridded(fcsy100km) <- TRUE
gridded(fcsy150km) <- TRUE

gridded(flsy10km) <- TRUE
gridded(flsy25km) <- TRUE
gridded(flsy50km) <- TRUE
gridded(flsy100km) <- TRUE
gridded(flsy150km) <- TRUE

gridded(fl1yp10km) <- TRUE
gridded(fl1yp25km) <- TRUE
gridded(fl1yp50km) <- TRUE
gridded(fl1yp100km) <- TRUE
gridded(fl1yp150km) <- TRUE

gridded(fl2yp10km) <- TRUE
gridded(fl2yp25km) <- TRUE
gridded(fl2yp50km) <- TRUE
gridded(fl2yp100km) <- TRUE
gridded(fl2yp150km) <- TRUE

gridded(frag10km.7) <- TRUE
gridded(frag25km.7) <- TRUE
gridded(frag50km.7) <- TRUE
gridded(frag100km.7) <- TRUE
gridded(frag150km.7) <- TRUE

gridded(frag10km.8) <- TRUE
gridded(frag25km.8) <- TRUE
gridded(frag50km.8) <- TRUE
gridded(frag100km.8) <- TRUE
gridded(frag150km.8) <- TRUE

gridded(logpcsy10km) <- TRUE
gridded(logpcsy25km) <- TRUE
gridded(logpcsy50km) <- TRUE
gridded(logpcsy100km) <- TRUE
gridded(logpcsy150km) <- TRUE

gridded(elevation) <- TRUE
gridded(ntlst) <- TRUE
gridded(tempseason) <- TRUE
gridded(precipseason) <- TRUE
gridded(pet) <- TRUE
gridded(logpcfc10km) <- TRUE


# raster
fcsy10km <- raster(fcsy10km)
fcsy25km <- raster(fcsy25km)
fcsy50km <- raster(fcsy50km)
fcsy100km <- raster(fcsy100km)
fcsy150km <- raster(fcsy150km)
plot(fcsy100km)

flsy10km <- raster(flsy10km)
flsy25km <- raster(flsy25km)
flsy50km <- raster(flsy50km)
flsy100km <- raster(flsy100km)
flsy150km <- raster(flsy150km)
plot(flsy100km)

fl1yp10km <- raster(fl1yp10km)
fl1yp25km <- raster(fl1yp25km)
fl1yp50km <- raster(fl1yp50km)
fl1yp100km <- raster(fl1yp100km)
fl1yp150km <- raster(fl1yp150km)
plot(fl1yp50km)

fl2yp10km <- raster(fl2yp10km)
fl2yp25km <- raster(fl2yp25km)
fl2yp50km <- raster(fl2yp50km)
fl2yp100km <- raster(fl2yp100km)
fl2yp150km <- raster(fl2yp150km)
plot(fl2yp150km)

frag10km.7 <- raster(frag10km.7)
frag25km.7 <- raster(frag25km.7)
frag50km.7 <- raster(frag50km.7)
frag100km.7 <- raster(frag100km.7)
frag150km.7 <- raster(frag150km.7)
plot(frag150km.7$frag150km.7)

frag10km.8 <- raster(frag10km.8)
frag25km.8 <- raster(frag25km.8)
frag50km.8 <- raster(frag50km.8)
frag100km.8 <- raster(frag100km.8)
frag150km.8 <- raster(frag150km.8)
plot(frag150km.8$frag150km.8)

logpcsy10km <- raster(logpcsy10km)
logpcsy25km <- raster(logpcsy25km)
logpcsy50km <- raster(logpcsy50km)
logpcsy100km <- raster(logpcsy100km)
logpcsy150km <- raster(logpcsy150km)
plot(logpcsy150km)

elevation <- raster(elevation)
ntlst <- raster(ntlst)
tempseason <- raster(tempseason)
precipseason <- raster(precipseason)
pet <- raster(pet)
logpcfc10km <- raster(logpcfc10km)
plot(logpcfc10km)

spillover_grids2022 <- stack(fcsy10km,fcsy25km,fcsy50km,fcsy100km,fcsy150km,
                             flsy10km,flsy25km,flsy50km,flsy100km,flsy150km,
                             fl1yp10km,fl1yp25km,fl1yp50km,fl1yp100km,fl1yp150km,
                             fl2yp10km,fl2yp25km,fl2yp50km,fl2yp100km,fl2yp150km,
                             frag10km.7,frag25km.7,frag50km.7,frag100km.7,frag150km.7,
                             frag10km.8,frag25km.8,frag50km.8,frag100km.8,frag150km.8,
                             logpcsy10km,logpcsy25km,logpcsy50km,logpcsy100km,logpcsy150km,
                             elevation,ntlst,tempseason,precipseason,
                             pet,logpcfc10km)


#### 2021 pred grid ####
##### join all covars, create FL loss vars, convert to rasters #####
fccovarsin <- read.csv("Covariate Data/Prediction Grid Extracted Vars/fccovars_Dec15.csv")
othercovarsin <- read.csv("Covariate Data/Prediction Grid Extracted Vars/othercovars_Dec15.csv")  
pccovars2021 <- read.csv("Covariate Data/Prediction Grid Extracted Vars/pccovars2021_Apr3.csv")
pccovars2020 <- read.csv("Covariate Data/Prediction Grid Extracted Vars/pccovars2020_Apr3.csv")
fragcovar <- read.csv("Covariate Data/Prediction Grid Extracted Vars/fragcovars_Dec15.csv")


allcovars <- fccovarsin %>%
  left_join(fragcovar,by= c("X","x","y")) %>%
  left_join(othercovarsin, by= c("X","x","y")) %>%
  left_join(pccovars2021, by = c("X","x","y")) %>%
  left_join(pccovars2020, by = c("X","x","y")) 

allcovars <- allcovars %>%
  mutate(
    frag10km.7 = fragsy.7_10km*100,
    frag25km.7 = fragsy.7_25km*100,
    frag50km.7 = fragsy.7_50km*100,
    frag100km.7 = fragsy.7_100km*100,
    frag150km.7 = fragsy.7_150km*100,
    
    flsy10km = (  fc2020_10km-fc2021_10km)*10000,
    flsy25km = (  fc2020_25km-fc2021_25km)*10000,
    flsy50km = (  fc2020_50km-fc2021_50km)*10000,
    flsy100km = ( fc2020_100km-fc2021_100km)*10000,
    flsy150km = ( fc2020_150km-fc2021_150km)*10000,
    
    fl1yp10km = (  fc2019_10km-fc2020_10km)*10000,
    fl1yp25km = (  fc2019_25km-fc2020_25km)*10000,
    fl1yp50km = (  fc2019_50km-fc2020_50km)*10000,
    fl1yp100km = ( fc2019_100km-fc2020_100km)*10000,
    fl1yp150km = ( fc2019_150km-fc2020_150km)*10000,
    
    fl2yp10km = (  fc2018_10km-fc2019_10km)*10000,
    fl2yp25km = (  fc2018_25km-fc2019_25km)*10000,
    fl2yp50km = (  fc2018_50km-fc2019_50km)*10000,
    fl2yp100km = ( fc2018_100km-fc2019_100km)*10000,
    fl2yp150km = ( fc2018_150km-fc2019_150km)*10000,
    
    logpcsy10km = log( pc2021_10km+1),
    logpcsy25km = log( pc2021_25km+1),
    logpcsy50km = log( pc2021_50km+1),
    logpcsy100km = log(pc2021_100km+1),
    logpcsy150km = log(pc2021_150km+1),
    
    logpcfc10km = logpcsy10km*fc2021_10km,
    logpcfc25km = logpcsy25km*fc2021_25km,
    logpcfc50km = logpcsy50km*fc2021_50km,
    logpcfc100km = logpcsy100km*fc2021_100km,
    logpcfc150km = logpcsy150km*fc2021_150km,
    
    fcsy10km = fc2021_10km*100,
    fcsy25km = fc2021_25km*100,
    fcsy50km = fc2021_50km*100,
    fcsy100km = fc2021_100km*100,
    fcsy150km = fc2021_150km*100,
    
    pet_10km = pet_10km/1000,
    tempseason_10km = tempseason_10km/10) %>%
  rename(pet = pet_10km,
         tempseason = tempseason_10km,
         precipseason = precipseason_10km,
         elevation = elevation_10km,
         ntlst = ntlst_10km) 

# create individual dfs
fcsy10km <- allcovars %>% dplyr::select(x,y,fcsy10km) 
fcsy25km <- allcovars %>% dplyr::select(x,y,fcsy25km)
fcsy50km <- allcovars %>% dplyr::select(x,y,fcsy50km)
fcsy100km <- allcovars %>% dplyr::select(x,y,fcsy100km)
fcsy150km <- allcovars %>% dplyr::select(x,y,fcsy150km)

flsy10km <- allcovars %>% dplyr::select(x,y,flsy10km) 
flsy25km <- allcovars %>% dplyr::select(x,y,flsy25km)
flsy50km <- allcovars %>% dplyr::select(x,y,flsy50km)
flsy100km <- allcovars %>% dplyr::select(x,y,flsy100km)
flsy150km <- allcovars %>% dplyr::select(x,y,flsy150km)

fl1yp10km <- allcovars %>% dplyr::select(x,y,fl1yp10km)
fl1yp25km <- allcovars %>% dplyr::select(x,y,fl1yp25km)
fl1yp50km <- allcovars %>% dplyr::select(x,y,fl1yp50km)
fl1yp100km <- allcovars %>% dplyr::select(x,y,fl1yp100km)
fl1yp150km <- allcovars %>% dplyr::select(x,y,fl1yp150km)

fl2yp10km <- allcovars %>% dplyr::select(x,y,fl2yp10km)
fl2yp25km <- allcovars %>% dplyr::select(x,y,fl2yp25km)
fl2yp50km <- allcovars %>% dplyr::select(x,y,fl2yp50km)
fl2yp100km <- allcovars %>% dplyr::select(x,y,fl2yp100km)
fl2yp150km <- allcovars %>% dplyr::select(x,y,fl2yp150km)

frag10km.7 <- allcovars %>% dplyr::select(x,y,frag10km.7)
frag25km.7 <- allcovars %>% dplyr::select(x,y,frag25km.7)
frag50km.7 <- allcovars %>% dplyr::select(x,y,frag50km.7)
frag100km.7 <- allcovars %>% dplyr::select(x,y,frag100km.7)
frag150km.7 <- allcovars %>% dplyr::select(x,y,frag150km.7)

logpcsy10km <- allcovars %>% dplyr::select(x,y,logpcsy10km)
logpcsy25km <- allcovars %>% dplyr::select(x,y,logpcsy25km)
logpcsy50km <- allcovars %>% dplyr::select(x,y,logpcsy50km)
logpcsy100km <- allcovars %>% dplyr::select(x,y,logpcsy100km)
logpcsy150km <- allcovars %>% dplyr::select(x,y,logpcsy150km)

elevation <- allcovars %>% dplyr::select(x,y,elevation)
ntlst <- allcovars %>% dplyr::select(x,y,ntlst)
tempseason <- allcovars %>% dplyr::select(x,y,tempseason)
precipseason <- allcovars %>% dplyr::select(x,y,precipseason)
pet <- allcovars %>% dplyr::select(x,y,pet)
logpcfc10km <- allcovars %>% dplyr::select(x,y,logpcfc10km)

# convert to rasters
coordinates(fcsy10km) <- ~x+y
coordinates(fcsy25km) <- ~x+y
coordinates(fcsy50km) <- ~x+y
coordinates(fcsy100km) <- ~x+y
coordinates(fcsy150km) <- ~x+y

coordinates(flsy10km) <- ~x+y
coordinates(flsy25km) <- ~x+y
coordinates(flsy50km) <- ~x+y
coordinates(flsy100km) <- ~x+y
coordinates(flsy150km) <- ~x+y

coordinates(fl1yp10km) <- ~x+y
coordinates(fl1yp25km) <- ~x+y
coordinates(fl1yp50km) <- ~x+y
coordinates(fl1yp100km) <- ~x+y
coordinates(fl1yp150km) <- ~x+y

coordinates(fl2yp10km) <- ~x+y
coordinates(fl2yp25km) <- ~x+y
coordinates(fl2yp50km) <- ~x+y
coordinates(fl2yp100km) <- ~x+y
coordinates(fl2yp150km) <- ~x+y

coordinates(frag10km.7) <- ~x+y
coordinates(frag25km.7) <- ~x+y
coordinates(frag50km.7) <- ~x+y
coordinates(frag100km.7) <- ~x+y
coordinates(frag150km.7) <- ~x+y

coordinates(logpcsy10km) <- ~x+y
coordinates(logpcsy25km) <- ~x+y
coordinates(logpcsy50km) <- ~x+y
coordinates(logpcsy100km) <- ~x+y
coordinates(logpcsy150km) <- ~x+y

coordinates(elevation) <- ~x+y
coordinates(ntlst) <- ~x+y
coordinates(tempseason) <- ~x+y
coordinates(precipseason) <- ~x+y
coordinates(pet) <- ~x+y
coordinates(logpcfc10km) <- ~x+y


# gridded
gridded(fcsy10km) <- TRUE
gridded(fcsy25km) <- TRUE
gridded(fcsy50km) <- TRUE
gridded(fcsy100km) <- TRUE
gridded(fcsy150km) <- TRUE

gridded(flsy10km) <- TRUE
gridded(flsy25km) <- TRUE
gridded(flsy50km) <- TRUE
gridded(flsy100km) <- TRUE
gridded(flsy150km) <- TRUE

gridded(fl1yp10km) <- TRUE
gridded(fl1yp25km) <- TRUE
gridded(fl1yp50km) <- TRUE
gridded(fl1yp100km) <- TRUE
gridded(fl1yp150km) <- TRUE

gridded(fl2yp10km) <- TRUE
gridded(fl2yp25km) <- TRUE
gridded(fl2yp50km) <- TRUE
gridded(fl2yp100km) <- TRUE
gridded(fl2yp150km) <- TRUE

gridded(frag10km.7) <- TRUE
gridded(frag25km.7) <- TRUE
gridded(frag50km.7) <- TRUE
gridded(frag100km.7) <- TRUE
gridded(frag150km.7) <- TRUE

gridded(logpcsy10km) <- TRUE
gridded(logpcsy25km) <- TRUE
gridded(logpcsy50km) <- TRUE
gridded(logpcsy100km) <- TRUE
gridded(logpcsy150km) <- TRUE

gridded(elevation) <- TRUE
gridded(ntlst) <- TRUE
gridded(tempseason) <- TRUE
gridded(precipseason) <- TRUE
gridded(pet) <- TRUE
gridded(logpcfc10km) <- TRUE


# raster
fcsy10km <- raster(fcsy10km)
fcsy25km <- raster(fcsy25km)
fcsy50km <- raster(fcsy50km)
fcsy100km <- raster(fcsy100km)
fcsy150km <- raster(fcsy150km)
plot(fcsy100km)

flsy10km <- raster(flsy10km)
flsy25km <- raster(flsy25km)
flsy50km <- raster(flsy50km)
flsy100km <- raster(flsy100km)
flsy150km <- raster(flsy150km)
plot(flsy100km)

fl1yp10km <- raster(fl1yp10km)
fl1yp25km <- raster(fl1yp25km)
fl1yp50km <- raster(fl1yp50km)
fl1yp100km <- raster(fl1yp100km)
fl1yp150km <- raster(fl1yp150km)
plot(fl1yp50km)

fl2yp10km <- raster(fl2yp10km)
fl2yp25km <- raster(fl2yp25km)
fl2yp50km <- raster(fl2yp50km)
fl2yp100km <- raster(fl2yp100km)
fl2yp150km <- raster(fl2yp150km)
plot(fl2yp150km)

frag10km.7 <- raster(frag10km.7)
frag25km.7 <- raster(frag25km.7)
frag50km.7 <- raster(frag50km.7)
frag100km.7 <- raster(frag100km.7)
frag150km.7 <- raster(frag150km.7)
plot(frag150km.7$frag150km.7)

logpcsy10km <- raster(logpcsy10km)
logpcsy25km <- raster(logpcsy25km)
logpcsy50km <- raster(logpcsy50km)
logpcsy100km <- raster(logpcsy100km)
logpcsy150km <- raster(logpcsy150km)
plot(logpcsy150km)

elevation <- raster(elevation)
ntlst <- raster(ntlst)
tempseason <- raster(tempseason)
precipseason <- raster(precipseason)
pet <- raster(pet)
logpcfc10km <- raster(logpcfc10km)
plot(logpcfc10km)

spillover_grids2021 <- stack(fcsy10km,fcsy25km,fcsy50km,fcsy100km,fcsy150km,
                             flsy10km,flsy25km,flsy50km,flsy100km,flsy150km,
                             fl1yp10km,fl1yp25km,fl1yp50km,fl1yp100km,fl1yp150km,
                             fl2yp10km,fl2yp25km,fl2yp50km,fl2yp100km,fl2yp150km,
                             frag10km.7,frag25km.7,frag50km.7,frag100km.7,frag150km.7,
                             logpcsy10km,logpcsy25km,logpcsy50km,logpcsy100km,logpcsy150km,
                             elevation,ntlst,tempseason,precipseason,
                             pet,logpcfc10km)

#### PREDICT ####
# load All-species Full Model list
# this is 100 individually fitted models
load("Model_list/AllSpecies_FullModel.RData") 

##### 100mod preds #####
# separate each model
m1 <- model_list[[1]]
m2 <- model_list[[2]]
m3 <- model_list[[3]]
m4 <- model_list[[4]]
m5 <- model_list[[5]]
m6 <- model_list[[6]]
m7 <- model_list[[7]]
m8 <- model_list[[8]]
m9 <- model_list[[9]]
m10 <- model_list[[10]]
m11 <- model_list[[11]]
m12 <- model_list[[12]]
m13 <- model_list[[13]]
m14 <- model_list[[14]]
m15 <- model_list[[15]]
m16 <- model_list[[16]]
m17 <- model_list[[17]]
m18 <- model_list[[18]]
m19 <- model_list[[19]]
m20 <- model_list[[20]]
m21 <- model_list[[21]]
m22 <- model_list[[22]]
m23 <- model_list[[23]]
m24 <- model_list[[24]]
m25 <- model_list[[25]]
m26 <- model_list[[26]]
m27 <- model_list[[27]]
m28 <- model_list[[28]]
m29 <- model_list[[29]]
m30 <- model_list[[30]]
m31 <- model_list[[31]]
m32 <- model_list[[32]]
m33 <- model_list[[33]]
m34 <- model_list[[34]]
m35 <- model_list[[35]]
m36 <- model_list[[36]]
m37 <- model_list[[37]]
m38 <- model_list[[38]]
m39 <- model_list[[39]]
m40 <- model_list[[40]]
m41 <- model_list[[41]]
m42 <- model_list[[42]]
m43 <- model_list[[43]]
m44 <- model_list[[44]]
m45 <- model_list[[45]]
m46 <- model_list[[46]]
m47 <- model_list[[47]]
m48 <- model_list[[48]]
m49 <- model_list[[49]]
m50 <- model_list[[50]]
m51 <- model_list[[51]]
m52 <- model_list[[52]]
m53 <- model_list[[53]]
m54 <- model_list[[54]]
m55 <- model_list[[55]]
m56 <- model_list[[56]]
m57 <- model_list[[57]]
m58 <- model_list[[58]]
m59 <- model_list[[59]]
m60 <- model_list[[60]]
m61 <- model_list[[61]]
m62 <- model_list[[62]]
m63 <- model_list[[63]]
m64 <- model_list[[64]]
m65 <- model_list[[65]]
m66 <- model_list[[66]]
m67 <- model_list[[67]]
m68 <- model_list[[68]]
m69 <- model_list[[69]]
m70 <- model_list[[70]]
m71 <- model_list[[71]]
m72 <- model_list[[72]]
m73 <- model_list[[73]]
m74 <- model_list[[74]]
m75 <- model_list[[75]]
m76 <- model_list[[76]]
m77 <- model_list[[77]]
m78 <- model_list[[78]]
m79 <- model_list[[79]]
m80 <- model_list[[80]]
m81 <- model_list[[81]]
m82 <- model_list[[82]]
m83 <- model_list[[83]]
m84 <- model_list[[84]]
m85 <- model_list[[85]]
m86 <- model_list[[86]]
m87 <- model_list[[87]]
m88 <- model_list[[88]]
m89 <- model_list[[89]]
m90 <- model_list[[90]]
m91 <- model_list[[91]]
m92 <- model_list[[92]]
m93 <- model_list[[93]]
m94 <- model_list[[94]]
m95 <- model_list[[95]]
m96 <- model_list[[96]]
m97 <- model_list[[97]]
m98 <- model_list[[98]]
m99 <- model_list[[99]]
m100 <- model_list[[100]]

# predict with each of the 100 models
pred1 <- makePreds(m1,spillover_grids2022)
pred2 <- makePreds(m2,spillover_grids2022)
pred3 <- makePreds(m3,spillover_grids2022)
pred4 <- makePreds(m4,spillover_grids2022)
pred5 <- makePreds(m5,spillover_grids2022)
pred6 <- makePreds(m6,spillover_grids2022)
pred7 <- makePreds(m7,spillover_grids2022)
pred8 <- makePreds(m8,spillover_grids2022)
pred9 <- makePreds(m9,spillover_grids2022)
pred10 <- makePreds(m10,spillover_grids2022)
pred11 <- makePreds(m11,spillover_grids2022)
pred12 <- makePreds(m12,spillover_grids2022)
pred13 <- makePreds(m13,spillover_grids2022)
pred14 <- makePreds(m14,spillover_grids2022)
pred15 <- makePreds(m15,spillover_grids2022)
pred16 <- makePreds(m16,spillover_grids2022)
pred17 <- makePreds(m17,spillover_grids2022)
pred18 <- makePreds(m18,spillover_grids2022)
pred19 <- makePreds(m19,spillover_grids2022)
pred20 <- makePreds(m20,spillover_grids2022)
pred21 <- makePreds(m21,spillover_grids2022)
pred22 <- makePreds(m22,spillover_grids2022)
pred23 <- makePreds(m23,spillover_grids2022)
pred24 <- makePreds(m24,spillover_grids2022)
pred25 <- makePreds(m25,spillover_grids2022)
pred26 <- makePreds(m26,spillover_grids2022)
pred27 <- makePreds(m27,spillover_grids2022)
pred28 <- makePreds(m28,spillover_grids2022)
pred29 <- makePreds(m29,spillover_grids2022)
pred30 <- makePreds(m30,spillover_grids2022)
pred31 <- makePreds(m31,spillover_grids2022)
pred32 <- makePreds(m32,spillover_grids2022)
pred33 <- makePreds(m33,spillover_grids2022)
pred34 <- makePreds(m34,spillover_grids2022)
pred35 <- makePreds(m35,spillover_grids2022)
pred36 <- makePreds(m36,spillover_grids2022)
pred37 <- makePreds(m37,spillover_grids2022)
pred38 <- makePreds(m38,spillover_grids2022)
pred39 <- makePreds(m39,spillover_grids2022)
pred40 <- makePreds(m40,spillover_grids2022)
pred41 <- makePreds(m41,spillover_grids2022)
pred42 <- makePreds(m42,spillover_grids2022)
pred43 <- makePreds(m43,spillover_grids2022)
pred44 <- makePreds(m44,spillover_grids2022)
pred45 <- makePreds(m45,spillover_grids2022)
pred46 <- makePreds(m46,spillover_grids2022)
pred47 <- makePreds(m47,spillover_grids2022)
pred48 <- makePreds(m48,spillover_grids2022)
pred49 <- makePreds(m49,spillover_grids2022)
pred50 <- makePreds(m50,spillover_grids2022)
pred51 <- makePreds(m51,spillover_grids2022)
pred52 <- makePreds(m52,spillover_grids2022)
pred53 <- makePreds(m53,spillover_grids2022)
pred54 <- makePreds(m54,spillover_grids2022)
pred55 <- makePreds(m55,spillover_grids2022)
pred56 <- makePreds(m56,spillover_grids2022)
pred57 <- makePreds(m57,spillover_grids2022)
pred58 <- makePreds(m58,spillover_grids2022)
pred59 <- makePreds(m59,spillover_grids2022)
pred60 <- makePreds(m60,spillover_grids2022)
pred61 <- makePreds(m61,spillover_grids2022)
pred62 <- makePreds(m62,spillover_grids2022)
pred63 <- makePreds(m63,spillover_grids2022)
pred64 <- makePreds(m64,spillover_grids2022)
pred65 <- makePreds(m65,spillover_grids2022)
pred66 <- makePreds(m66,spillover_grids2022)
pred67 <- makePreds(m67,spillover_grids2022)
pred68 <- makePreds(m68,spillover_grids2022)
pred69 <- makePreds(m69,spillover_grids2022)
pred70 <- makePreds(m70,spillover_grids2022)
pred71 <- makePreds(m71,spillover_grids2022)
pred72 <- makePreds(m72,spillover_grids2022)
pred73 <- makePreds(m73,spillover_grids2022)
pred74 <- makePreds(m74,spillover_grids2022)
pred75 <- makePreds(m75,spillover_grids2022)
pred76 <- makePreds(m76,spillover_grids2022)
pred77 <- makePreds(m77,spillover_grids2022)
pred78 <- makePreds(m78,spillover_grids2022)
pred79 <- makePreds(m79,spillover_grids2022)
pred80 <- makePreds(m80,spillover_grids2022)
pred81 <- makePreds(m81,spillover_grids2022)
pred82 <- makePreds(m82,spillover_grids2022)
pred83 <- makePreds(m83,spillover_grids2022)
pred84 <- makePreds(m84,spillover_grids2022)
pred85 <- makePreds(m85,spillover_grids2022)
pred86 <- makePreds(m86,spillover_grids2022)
pred87 <- makePreds(m87,spillover_grids2022)
pred88 <- makePreds(m88,spillover_grids2022)
pred89 <- makePreds(m89,spillover_grids2022)
pred90 <- makePreds(m90,spillover_grids2022)
pred91 <- makePreds(m91,spillover_grids2022)
pred92 <- makePreds(m92,spillover_grids2022)
pred93 <- makePreds(m93,spillover_grids2022)
pred94 <- makePreds(m94,spillover_grids2022)
pred95 <- makePreds(m95,spillover_grids2022)
pred96 <- makePreds(m96,spillover_grids2022)
pred97 <- makePreds(m97,spillover_grids2022)
pred98 <- makePreds(m98,spillover_grids2022)
pred99 <- makePreds(m99,spillover_grids2022)
pred100 <- makePreds(m100,spillover_grids2022)
predstack100 <- stack(pred1,pred2,pred3,pred4,pred5,pred6,pred7,pred8,pred9,
                      pred10,pred11,pred12,pred13,pred14,pred15,pred16,pred17,pred18,pred19,
                      pred20,pred21,pred22,pred23,pred24,pred25,pred26,pred27,pred28,pred29,
                      pred30,pred31,pred32,pred33,pred34,pred35,pred36,pred37,pred38,pred39,
                      pred40,pred41,pred42,pred43,pred44,pred45,pred46,pred47,pred48,pred49,
                      pred50,pred51,pred52,pred53,pred54,pred55,pred56,pred57,pred58,pred59,
                      pred60,pred61,pred62,pred63,pred64,pred65,pred66,pred67,pred68,pred69,
                      pred70,pred71,pred72,pred73,pred74,pred75,pred76,pred77,pred78,pred79,
                      pred80,pred81,pred82,pred83,pred84,pred85,pred86,pred87,pred88,pred89,
                      pred90,pred91,pred92,pred93,pred94,pred95,pred96,pred97,pred98,pred99,
                      pred100)
plot(predstack100$layer.11)


##### Raw mean Odds Calc #####
#' calculate mean odds across the 100 models
rawoddsstack100 <- predstack100

predstack100 <- rawoddsstack100[[1]]*0
for (i in 1:100) {
  xyz <- rawoddsstack100[[i]]
  xyz2 <- as.data.frame(xyz)
  overallmean <- mean(xyz2[,1], na.rm=T)
  predrelodds1 <- xyz/overallmean
  predstack100 <- stack(predstack100,predrelodds1)
}
predstack100 <- dropLayer(predstack100,1)
reloddsmean <- mean(predstack100)
plot(reloddsmean)


