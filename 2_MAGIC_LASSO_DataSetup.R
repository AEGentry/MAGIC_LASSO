#######################################################
### MAGIC LASSO ###
### Script written by Amanda Elswick Gentry and B. Todd Webb
### For associated manuscript describing the method, see: https://www.biorxiv.org/content/10.1101/2021.04.29.442057v1
#######################################################

### Step 2 ###
### Further prepare UKB data for analysis

#######################################################

### Comparing missingness across top level variables as a function of missingingness in AUDIT ###

### UKB AUDIT working folder
workdir <- "your_directory_name_here"
plotdir <- paste(workdir, "Plots/", sep="")
preddir <- paste(workdir, "Scripts/Predictions/", sep="")
outdir <- paste(workdir, "Scripts/Output/", sep="")
######################################################################
### Load the table object with missingness counts
### Code below shows how it was computed

### The Comp.withAUD object has dimensions 631 x 9 (one row for every variable in the top-level set used for the Lasso application)
### The columns are:
###   Row.names     The UKB variable identifier
###   Both.obs      The count of subjects with both the variable AND AUDIT observed
###   Both.miss     The count of subjects with both the variable AND AUDIT missing
###   Obs.AUDmiss   The count of subjects with the variable observed and AUDIT missing
###   Miss.AUDobs   The count of subjects with the variable missing and AUDIT observed
###   Nmiss         The count of subjects with the variable missing
###   Nobs          The count of subjects with the variable observed
###   PctOverlap    Both.obs/Nobs = The proportion of subjects with the variable measured who also have AUDIT measured

### Explore some tables/plots to figure out how to filter variables with (a) lower missingness and (b) a good proportion of overlap with measured AUDIT

load(paste(outdir,"MissingnessComparison.RDS", sep=""))
library(tidyverse)
### Get the deciles
Comp.withAUD <- Comp.withAUD %>%
  mutate(quantile = ntile(Nobs, 10))
Comp.withAUD$quantile <- as.factor(Comp.withAUD$quantile)
### Look at the spread of percent overlap with measured AUDIT
Comp.sum.pct <- Comp.withAUD%>%
  group_by(quantile)%>% 
  summarise(Min=min(PctOverlap), Median=median(PctOverlap),Max=max(PctOverlap))
Comp.sum.nobs <- Comp.withAUD%>%
  group_by(quantile)%>% 
  summarise(Min=min(Nobs), Median=median(Nobs),Max=max(Nobs))         

Comp.sum <- cbind(Comp.sum.pct, Comp.sum.nobs[,2:4])
colnames(Comp.sum) <- c("Quantile", "Min.PropOverlap", "Med.PropOverlap", "Max.PropOverlap",
                        "Min.Nobs", "Med.Nobs", "Max.Nobs")
Comp.sum

######################################################################
### Starting here shows how the Comp.withAUD object was computed
### it's very slow, so it's not recommended to recreate it
### Just load the saved object, as shown above
load(paste(workdir,"23704_30302_DataDict.Rds", sep=""))
load(paste(workdir, "UKB_23704_30302Top_Labeled.RData", sep=""))
### Some dummy forgot to add the rownames to the bd.top dataframe
rownames(bd.top) <- bd.top$f.eid
load(paste(workdir, "UKB23704_30302_MissByCol_DF_Top.Rds", sep=""))
load(paste(workdir, "UKB23704_30302_Pairwise_Mat.Rds", sep=""))

load(paste(workdir, "UKB_AUDIT.Rds", sep=""))
AUDIT.sub <- AUDIT[, c("f.eid", "score")]
AUDIT.vars <- c("f.20414.0.0", "f.20403.0.0", "f.20416.0.0", "f.20413.0.0", "f.20407.0.0",
                "f.20412.0.0", "f.20409.0.0", "f.20408.0.0", "f.20411.0.0", "f.20405.0.0")

isID <- which(colnames(bd.top) == "f.eid")
bd.top.noID <- bd.top[, -isID]
AUDIT1 <- data.frame(AUDIT$AUDIT1)

colnames(AUDIT1) <- "AUDIT1"
rownames(AUDIT1) <- rownames(AUDIT)

### this takes forever
bd.topA <- merge(AUDIT1, bd.top.noID, by.x="row.names", by.y="row.names", all.y=TRUE)
rownames(bd.topA) <- bd.topA$Row.names
test <- bd.topA[,-1]
save(bd.topA, file=paste(workdir, "UKB23704_30302_Pairwise_Top_withAUDIT.Rds", sep=""))

bd.miss.top <- is.na(test)
chunkSz <- 100
nVars <- dim(bd.miss.top)[2]
nSets <- ceiling(nVars/chunkSz)
starting <- (0:(nSets - 1) * chunkSz + 1)
ending <- c(1:(nSets - 1) * chunkSz, nVars)
nChunks <- length(starting)
chunkSzA <- ending - starting + 1


for (ii in 1:nChunks){
  print(ii)
  aa <- bd.miss.top[ ,starting[ii]:ending[ii]]
  aa <- aa * 1
  if (ii == 1){
    bd.miss.ind.top <- aa
  } else {
    bd.miss.ind.top <- cbind(bd.miss.ind.top, aa)
  }
  print(dim(bd.miss.ind.top))
  rm(aa)
}
colnames(bd.miss.ind.top) <- colnames(test)

### Missing is 1
Both.miss <- Both.obs <- Obs.AUDmiss <- Miss.AUDobs <- numeric()
for (j in 1:(dim(bd.miss.ind.top)[2]-1)){
  Both.miss[j]    <- sum(bd.miss.ind.top[,1] * bd.miss.ind.top[,j +1])
  Both.obs[j]     <- sum( (bd.miss.ind.top[,1] + bd.miss.ind.top[,j +1])==0 )
  Obs.AUDmiss[j]  <- sum( (bd.miss.ind.top[,1] - bd.miss.ind.top[,j +1])==1 )
  Miss.AUDobs[j]  <- sum( (bd.miss.ind.top[,1] - bd.miss.ind.top[,j +1])==-1 )
}
Comp.withAUD <- cbind(Both.obs, Both.miss, Obs.AUDmiss, Miss.AUDobs)
rownames(Comp.withAUD) <- colnames(bd.miss.ind.top)[2:dim(bd.miss.ind.top)[2]]
Comp.withAUD <- data.frame(Comp.withAUD)
Comp.withAUD <- merge(Comp.withAUD, miss.total.topDF, by.x="row.names", by.y="Col")
Comp.withAUD$Nobs <- dim(bd.top)[1] - Comp.withAUD$Nmiss
Comp.withAUD$PctOverlap <- Comp.withAUD$Both.obs/Comp.withAUD$Nobs
save(Comp.withAUD, file=paste(outdir, "MissingnessComparison.RDS", sep=""))

### Once the objects are computer, you can reload them and use them to compute metrics for filtering

infile1<-paste(outdir, "MissingnessComparison.RDS", sep="") 
inName<-load(infile1)
missSum<-Comp.withAUD
rm(Comp.withAUD)
infile2<-paste(inputdir, "UKB23704_30302_Pairwise_Top_withAUDIT.Rds", sep="") ### object created in MissingnessCompwithAUD.R by Amanda 2020_05_22
### This object is just the data plus AUDIT Q1 added
load(infile2)
nrow(bd.topA)

##################
table(is.na(bd.topA$AUDIT.AUDIT1))

###Make summary of non missingness (completeness)
naSum<-apply(bd.topA,2,function(x) {sum(!is.na(x))}) ### count the non-missing (ie the PRESENT) observations in each column
naSum_obs<-apply(subset(bd.topA,!is.na(AUDIT.AUDIT1)),2,function(x) {sum(!is.na(x))}) ### count the non-missing in the subset of subjects WITH AUDIT
naSum_unobs<-apply(subset(bd.topA,is.na(AUDIT.AUDIT1)),2,function(x) {sum(!is.na(x))}) ### count the non-missing in the subset of subjects WITHOUT AUDIT

###Check names before cbinding### good to know Todd is equally paranoid about this as I am
table(names(naSum)==names(naSum_obs))
table(names(naSum)==names(naSum_unobs))
table(names(naSum_unobs)==names(naSum_obs))

###Combine NA counts into summary table
naSum<-as.data.frame(cbind(naSum,naSum_obs,naSum_unobs),stringsAsFactors=F)
names(naSum)<-c("nAll","nObs","nUnobs")

###Make completeness (1-missingness) variables
naSum$compAll<-naSum$nAll/max(naSum$nAll) ### proportion missing out of the whole sample
naSum$compObs<-naSum$nObs/max(naSum$nObs) ### proportion missing in the observed AUDIT set
naSum$compUnobs<-naSum$nUnobs/max(naSum$nUnobs) ### proporion missing in the unobserved AUDIT set

###What is minimum completeness in full, obs, and unobs sets?
fullCut<-0.01
obsCut<-0.01
unobsCut<-0.01

table(round(naSum$compAll[naSum$compAll>fullCut],2))
table(round(naSum$compObs[naSum$compAll>obsCut],2))
table(round(naSum$compUnobs[naSum$compAll>unobsCut],2))

###What are the minimum completeness observed across Observed(training), Unobserved (prediction) sets
min(naSum$compAll)
min(naSum$compObs)
min(naSum$compUnobs)

###Find ratio of completeness across available Training and Prediction subsets
naSum$Predict_V_Train<-naSum$compUnobs/naSum$compObs  ###0 indicates not informative in predicted set.
###Other ratios to explore
naSum$Train_V_Predict<-naSum$compObs/naSum$compUnobs
naSum$Predict_V_All<-naSum$compUnobs/naSum$compAll
naSum$Train_V_All<-naSum$compObs/naSum$compAll

###Check variables with low predcition to training completeness ratio
naSum[naSum$Predict_V_Train<0.1,]

###Summarize by ratio bin
table(round(naSum$Predict_V_Train,1))

###What other variables can be derived  
#missSum$predict_complet<-

### Save the table
save(naSum, file=paste(workdir, "MissingnessComparison_vBTW.Rds"))



