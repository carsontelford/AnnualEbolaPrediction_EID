
#### PART 2: Presence+Absence point generation: Annually updated estimates of ebolavirus spillover potential accounting for changes to forests and human populations  ####


#' README:
#' The code below is run with the post processed clean dataset of presence and 
#' absence points. At each location/year, covariate values have been extracted
#' using buffers of varying sizes. 


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


#### Descriptive ####
setwd("C:/Users/pwv0/OneDrive - CDC/Filo Defor Data")
tdf <- read.csv("R Code/EID files/tdf_clean.csv")

tdf2 <- tdf %>%
  mutate(flsy10km = fc1yp10km-fcsy10km,
         flsy25km = fc1yp25km-fcsy25km,
         flsy50km = fc1yp50km-fcsy50km,
         flsy100km = fc1yp100km-fcsy100km,
         flsy150km = fc1yp150km-fcsy150km,
         
         fl1yp10km = fc2yp10km-fc1yp10km,
         fl1yp25km = fc2yp25km-fc1yp25km,
         fl1yp50km = fc2yp50km-fc1yp50km,
         fl1yp100km = fc2yp100km-fc1yp100km,
         fl1yp150km = fc2yp150km-fc1yp150km,
         
         fl2yp10km =  fc3yp10km-fc2yp10km,
         fl2yp25km =  fc3yp25km-fc2yp25km,
         fl2yp50km =  fc3yp50km-fc2yp50km,
         fl2yp100km = fc3yp100km-fc2yp100km,
         fl2yp150km = fc3yp150km-fc2yp150km,
         
         logpcsy10km = log(pcsy10km+1),
         logpcsy25km = log(pcsy25km+1),
         logpcsy50km = log(pcsy50km+1),
         logpcsy100km = log(pcsy100km+1),
         logpcsy150km = log(pcsy150km+1),
         
         logpcfc10km = logpcsy10km*fcsy10km) %>%
  
  mutate(fcsy10km = fcsy10km*100,
         fcsy25km = fcsy25km*100,
         fcsy50km = fcsy50km*100,
         fcsy100km=fcsy100km*100,
         fcsy150km=fcsy150km*100,
         
         flsy10km = flsy10km*10000,
         flsy25km = flsy25km*10000,
         flsy50km = flsy50km*10000,
         flsy100km = flsy100km*10000,
         flsy150km = flsy150km*10000,
         fl1yp10km = fl1yp10km*10000,
         fl1yp25km = fl1yp25km*10000,
         fl1yp50km = fl1yp50km*10000,
         fl1yp100km = fl1yp100km*10000,
         fl1yp150km = fl1yp150km*10000,
         fl2yp10km = fl2yp10km*10000,
         fl2yp25km = fl2yp25km*10000,
         fl2yp50km = fl2yp50km*10000,
         fl2yp100km = fl2yp100km*10000,
         fl2yp150km = fl2yp150km*10000,
         
         frag10km.7 = frag.7allarea10km*100,
         frag25km.7 = frag.7allarea25km*100,
         frag50km.7 = frag.7allarea50km*100,
         frag100km.7 = frag.7allarea100km*100,
         frag150km.7 = frag.7allarea150km*100) %>%
  
  mutate(pet = pet/1000,
         tempseason = tempseason/10) %>%
  
  filter(Year>2000 & Year <2022) %>%
  dplyr::select(-c(fc1yp10km,fc1yp25km,fc1yp50km,fc1yp100km,fc1yp150km,
                   fc2yp10km,fc2yp25km,fc2yp50km,fc2yp100km,fc2yp150km,
                   fc3yp10km,fc3yp25km,fc3yp50km,fc3yp100km,fc3yp150km,
                   pcsy10km,pcsy25km,pcsy50km,pcsy100km,pcsy150km,
                   frag.7allarea10km,frag.7allarea25km,frag.7allarea50km,frag.7allarea100km,frag.7allarea150km))

table(tdf2$Strain)
tdf2$spillover <- ifelse(tdf2$Strain %in% c("Absence_logpdw"),0,1)



#### MODEL ####
tdf3 <- tdf2
brtdf <- tdf3 


#### CV AUC 3 FOLD ####
tdf4 <- brtdf %>%
  filter(Strain %in% c("Sudan","Bundibugyo", "Zaire","Absence_logpdw")) %>%
  arrange(spillover,Year)


# Define folds by Year
folds <- unique(tdf4$Year)  # Unique years as folds
n_folds <- length(folds)    # Number of folds

# Initialize a dataframe to store results
combinedfit <- data.frame()

# Leave-one-year-out cross-validation
for (fold in folds) {
  # Split data into training and testing
  train_data <- tdf4 %>% filter(Year != fold)
  test_data <- tdf4 %>% filter(Year == fold)
  
  presnum <- sum(train_data$spillover == 1)
  abspoints <- which(train_data$spillover == 0)
  prespoints <- which(train_data$spillover == 1)
  
  pab.reps <- 100
  dsl <- list()
  
  for (j in 1:pab.reps) {
    abs.samp <- sample(abspoints, presnum * 50, replace = TRUE)
    dse <- rbind(train_data[abs.samp, ], train_data[prespoints, ])
    dsl[[j]] <- dse
  }
  
  # Train ensemble models
  sfInit(cpus = 4, parallel = TRUE)
  sfLibrary(seegSDM)
  model_list <- sfLapply(
    dsl,
    gbm.step,
    gbm.x = c(7:42),
    gbm.y = 43,
    gbm.coords = 3:4,
    family = "bernoulli",
    tree.complexity = 3,
    learning.rate = 0.001,
    bag.fraction = 0.5
  )
  sfStop()
  
  # Predict on the test data
  preds1test <- data.frame(row = 1:nrow(test_data))
  for (j in 1:100) {
    preds <- predict.gbm(
      model_list[[j]],
      test_data,
      n.trees = model_list[[j]]$gbm.call$best.trees,
      type = "response"
    )
    preds <- as.data.frame(preds)
    preds$row <- 1:nrow(test_data)
    preds1test <- preds1test %>%
      left_join(preds, by = "row")
  }
  
  # Average predictions across models
  preds1test_2 <- preds1test %>% select(-row)
  fold_fit <- rowMeans(preds1test_2)
  
  # Combine fitted values with test data
  fold_fitted <- cbind(test_data, foldfit = fold_fit)
  combinedfit <- rbind(combinedfit, fold_fitted)
}

# Calculate AUC
library(pROC)
par(pty = "s")
roc(combinedfit$spillover,
    combinedfit$foldfit,
    plot = TRUE, legacy.axes = TRUE, print.auc = TRUE,
    main = "Leave-One-Year-Out Cross Validation")



#### AUCSensSpec  ####
# odds cutoff All species full
library(ROCR)
ROCR_pred_test <- prediction(combinedfit$foldfit,
                             combinedfit$spillover)
perf <- performance(ROCR_pred_test,'sens','spec')
perf
plot(perf)

cutoffs <- data.frame(cut=perf@alpha.values[[1]], spec=perf@x.values[[1]], 
                      sens=perf@y.values[[1]])
cutoffs$optim <- cutoffs$spec*cutoffs$sens
cutoffs <- cutoffs %>%
  arrange(desc(optim))
head(cutoffs)
# best cutoff is .008200968; 
# Sensitivity: .91
# Specificity: 0.7269



##### ROC plots #####
library(gridExtra)
library(grid)
# ROC Allspeices full model
mydatatable <- data.frame(
  Prediction = c("Spillover","No Spillover"),
  Spillover = c(19, 3),
  No_Spillover = c(2731, 7269)
)
colnames(mydatatable)[colnames(mydatatable) == "No_Spillover"] <- "No Spillover"
table_theme <- ttheme_default(
  core = list(fg_params=    list(cex = .4,lineheight=.1, hjust=.5,vjust=.5),
              padding = unit(c(2,2),'mm')),
  colhead = list(fg_params= list(cex = .35,lineheight=.1, hjust=.5,vjust=.5),
                 padding = unit(c(2,2),'mm'))
)
table_grob <- tableGrob(mydatatable,rows=NULL,theme=table_theme)


par(pty="s")
ASfull <- roc(combinedfit$spillover,
              combinedfit$foldfit,
              plot=T,print.auc=T)
theme_set(theme_bw())
p1 <- ggroc(ASfull, size = .9) +
  ggtitle("Multi-species Analysis (AUC = 0.88)") +
  theme(plot.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7)) +
  labs(x = "Specificity", y = "Sensitivity") +
  geom_hline(yintercept=.91, linetype="dashed",alpha=.6, 
             color = "red", size=.5) +
  geom_vline(xintercept=.72, linetype="dashed",alpha=.6, 
             color = "red", size=.5) +
  annotation_custom(
    grob = table_grob,
    xmin = 0, xmax = -.6,
    ymin = 0, ymax = .3
  )
p1



#### Final Model fit with all data ####
colnames(brtdf)
table(brtdf$Strain)

tdf4 <- brtdf
tdfexcl <- tdf4

presnum<-sum(tdfexcl$spillover==1)
abspoints <-which(tdfexcl$spillover==0)
prespoints <- which(tdfexcl$spillover==1)

pab.reps <- 100
# create data.frames
dsl <- {}
for(j in 1:pab.reps){
  abs.samp<-c(sample(abspoints,(presnum*50),replace=F))
  dse<-rbind(tdfexcl[abs.samp,],tdfexcl[prespoints,])
  # assign(paste("dse",j,sep=""),dse)
  dsl[[j]] <- dse
}

##### Run it #####

library(maptools)
require(fields)
require(parallel)

colnames(dsl[[1]])
sfInit(cpus = 4, parallel = TRUE)
sfLibrary(seegSDM)
model_list <- sfLapply(dsl,
                       runBRT,
                       gbm.x = c(7:42),
                       gbm.y = 43,
                       gbm.coords = 3:4,
                       family = "bernoulli",
                       tree.complexity = 3,
                       learning.rate = 0.0009,
                       n.folds=3,
                       bag.fraction = 0.5)
sfStop()
save(model_list,file="FittedModelEnsembles/AllspeciesFullModel_WeightedAbsence.RData")


#### Load Model Lists ####
load("C:/Users/pwv0/OneDrive - CDC/Filo Defor Data/Manuscripts/Model_list/AllSpecies_FullModelJuly24.RData")


#### Stats All-species stats ####
max_values <- sapply(model_list, function(model) median(model$model$fitted))
average_max_value <- mean(max_values)
average_max_value

pab.reps <- 100
list.stats<-matrix(nrow=pab.reps,ncol=8)
for(i in 1:pab.reps){
  list.stats[i,]<-as.numeric(as.character(model_list[[i]]$model$cv.statistics[c(1:6,9:10)]))
}
out.stats<-as.data.frame(list.stats)
names(out.stats)<-c("deviance.mean","deviance.se","correlation.mean","correlation.se",
                    "discrimination.mean","discrimination.se","cv.threshold","cv.threshold.se")
head(out.stats)

# calc overall auc 
mean_auc <- mean(out.stats$discrimination.mean)
meanauc_CI <- quantile(out.stats$discrimination.mean, probs = c(0.025,0.975))
seauc_CI <- quantile(out.stats$discrimination.se, probs = c(0.025,0.975))

aucall <- c(mean_auc,meanauc_CI)
aucall
aucseall <- c(mean_se,seauc_CI)
aucseall

##### Relaive influence plots ####
pab.reps<-100
nvars<-36
nres<-5
top.vars<-12

relinf<-matrix(NA,nrow=nvars,ncol=1+pab.reps)
ebola_relinf.ds<-as.data.frame(relinf)
names(dsl[[1]])
ebola_relinf.ds[,1]<-names(dsl[[1]][,c(7:42)])
names(ebola_relinf.ds)<-c("Var",1:pab.reps)
for(k in 1:pab.reps){
  tmp.sum<-summary(model_list[[k]]$model,plotit=F)
  for(i in 1:length(ebola_relinf.ds[,1])){
    ebola_relinf.ds[i,k+1]<-tmp.sum$rel.inf[tmp.sum$var==ebola_relinf.ds$Var[i]]
  }
}

nridf<-matrix(nrow=pab.reps,ncol=nvars)
for(i in 1:pab.reps){
  for(k in 1:nvars){
    nridf[i,k]<-ebola_relinf.ds[k,i+1]
  }
}

ebola_relinf.ds$varname<-ebola_relinf.ds$Var
ebola_relinf.ds$median<-NA
for(i in 1:length(ebola_relinf.ds[,1])){
  ebola_relinf.ds$median[i]<-median(as.numeric(ebola_relinf.ds[i,2:(pab.reps+1)]))
}
ebola_ord.ri<-ebola_relinf.ds[order(-ebola_relinf.ds$median),]
ebola_rim<-t(ebola_ord.ri[1:top.vars,2:101])

top12 <- arrange(ebola_relinf.ds, desc(median))
top12 <- top12[1:12,1]
top12
nameevdfull <- c("1. LogPC 100km","2. LogPC 10km","3. LogPC*FC 10km","4. Frag 150km",
                 "5. LogPC 150km","6. FL SY 10km","7. Precip. Ssn.","8. FC 10km",
                 "9. Temp. Ssn.","10. PET","11. FL 2YP 10km","12. FL SY 100km")

varname <- nameevdfull
type<-c("pd","pd","pd","pd",
        "pd","pd","pd","pd",
        "pd","pd","pd","pd")

ebola_relinf.ds$median<-NA
for(i in 1:length(ebola_relinf.ds[,1])){
  ebola_relinf.ds$median[i]<-median(as.numeric(ebola_relinf.ds[i,2:(pab.reps+1)]))
}
ebola_ord.ri<-ebola_relinf.ds[order(-ebola_relinf.ds$median),]
ebola_rim<-t(ebola_ord.ri[1:top.vars,2:101])

par(mar=c(10.8,6,3,3))
par(mfrow=c(1,1))
cv<-brewer.pal(4,'Dark2')
col.vec<- ifelse(type=='FL',cv[4],
                 ifelse(type=='pd',cv[1],
                        ifelse(type=='temp',cv[4],cv[4])))
boxplot(ebola_rim,names=varname,las=2,cex.axis=1.3,ylim=c(0,20),cex.lab=2,col=col.vec,ylab='',alpha=.8,
        staplewex=.2,cex=0,outwex=0,whisklty = 1,
        cex.main=2,
        main='Top 12 Predictors: Multi-species Analysis')
# box(col="black",lwd=2)
title(ylab="Relative Importance (%)", mgp=c(3,1,0),
      cex.lab=1.7)

##### Marginal Effect Curves #####
par(mfrow=c(3,4))
par(mgp=c(4,1,0))
par(mar=c(3,6,3,2))

# par(mfrow=c(4,4))
# par(mar=c(4,3,2,1))
for(k in 1:12){
  
  xvals<-matrix(nrow=100,ncol=pab.reps)
  for(i in 1:pab.reps){
    xvals[,i]<-plot(model_list[[i]]$model,which(ebola_relinf.ds$Var==ebola_ord.ri$Var[k]),return.grid=T)[,1]
  }
  yvals<-matrix(nrow=100,ncol=pab.reps)
  for(i in 1:pab.reps){
    yvals[,i]<-plot(model_list[[i]]$model,which(ebola_relinf.ds$Var==ebola_ord.ri$Var[k]),return.grid=T)[,2]
  }
  
  # plot all marginal effect curves for var with highest mean rel inf
  median.y<-rep(NA,length(yvals[,1]))
  mean.y<-rep(NA,length(yvals[,1]))
  upp.y<-rep(NA,length(yvals[,1]))
  lwr.y<-rep(NA,length(yvals[,1]))
  mean.x<-rep(NA,length(yvals[,1]))
  
  for(i in 1:length(yvals[,1])){
    median.y[i]<-median(yvals[i,])
    mean.y[i]<-mean(yvals[i,])
    mean.x[i]<-mean(xvals[i,])
    lwr.y[i]<-quantile(yvals[i,], probs = 0.025)
    upp.y[i]<-quantile(yvals[i,], probs = 0.975)
  }
  
  plot(mean.x,mean.y,lwd=3,col="white",type="l",ylim=c(min(lwr.y),max(upp.y)),
       las=1,cex.axis=1.8,cex.lab=2,xlab='',main='', ylab= 'Marg. Lik.',cex.main=1.6,tck=-.02,
       alpha=.8)
  mtext(varname[k],side=3,line=.3,cex=1.6)
  box(col="black",lwd=2)
  polygon(c(mean.x,rev(mean.x)),c(lwr.y,rev(upp.y)),col=col.vec[k],border=NA)
  lines(mean.x,mean.y,lwd=3)
}






