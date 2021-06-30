### MAGIC-LASSO Prediction for AUDIT in the UK Biobank
### Written by Amanda E. Gentry, Ph.D.
### Methods described in "Missingness Adapted Group Informed Clustered (MAGIC)-LASSO: A novel paradigm for prediction in data with widespread non-random missingness"
### by Amanda E. Gentry, Robert M. Kirkpatrick, Roseann E. Peterson, and Bradley Todd Webb

### Objects loaded here were mostly created using the "LoadUKB_23704_30302.Rmd" script

### Necessary Libraries ###
library(Matrix) ### v1.2.17
library(fastDummies) ### v1.5.0
library(grpreg) ### v3.2.1

### Set Directories ###############################################################################

workdir <- "/filepath/" 
plotdir <- "/filepath/"
inputdir <- "/filepath/"
preddir <- "/filepath/"
outdir <- "/filepath/"
PubOutdir <- "/filepath/"

### Load preconstructed objects ###################################################################

### Data Dictionary ###
### Object containing the UKB data dictionary, trimmed to include onnly the variables in our download baskets
load(paste(inputdir,"23704_30302_DataDict.Rds", sep="")) ### object "data.dict" (dim 9613 x  21)

### Data object with all subjects and all top-level variables ###
### This is a dataframe of all the UKB data for our filtered variables, with all the factor-class variables labeled and ordered according to ordering and labels provided by UKB
load(paste(inputdir, "UKB_23704_30302Top_Labeled.RData", sep="")) ### object "bd.top" (dim 502536  x  632)
print(dim(bd.top)) ### [1] 502536    632
rownames(bd.top) <- bd.top$f.eid ### Add rownames
Nobs <- dim(bd.top)[1] ### 502536

### Load the pre-computed missingness metrics ###
### First, a dataframe containing counts of missingness for each variable
load(paste(inputdir, "UKB23704_30302_MissByCol_DF_Top.Rds", sep="")) ### object "miss.total.topDF" (dim 631 x  2) 
### Second, a dataframe of counts of pairwise observations PRESENT for each pair of variables
load(paste(inputdir, "UKB23704_30302_Pairwise_Mat.Rds", sep="")) ### object bd.miss.pairwise.topS (dim 631 x 631) 

### Load the AUDIT data, total score, consumption, and problems computed separately and loaded here ###
load(paste(inputdir, "UKB_AUDIT.Rds", sep="")) ### object AUDIT (dim 157162  x   15)
AUDIT.sub <- AUDIT[, c("f.eid", "score")]
AUDIT.vars <- c("f.20414.0.0", "f.20403.0.0", "f.20416.0.0", "f.20413.0.0", "f.20407.0.0",
                "f.20412.0.0", "f.20409.0.0", "f.20408.0.0", "f.20411.0.0", "f.20405.0.0")

### Load balancing metrics ###
### Object precomputed to calculate tau for each variable in the set
load(paste(outdir, "MissingnessComparison_vBTW.Rds", sep="")) ### object "naSum" (dim 633 x 10)
### Drop Pred Vs Train set ratio (tau) where t<=0.7 (277 variables)
drop <- rownames(naSum[round(naSum$Predict_V_Train,1) <= 0.7, ])
drop.ID <- which(rownames(bd.miss.pairwise.topS) %in% drop)
length(drop.ID) # 277
bd.miss.pairwise.topS <- bd.miss.pairwise.topS[-drop.ID, -drop.ID]
dim(bd.miss.pairwise.topS) # 354 354

######################################################################## Initial Cluster/Lasso ####################################
########################### Cluster ###########################
dist.pair <- dist(bd.miss.pairwise.topS, method="euclidean") ### calculate the distance matrix
modSet1 <- hclust(dist.pair, method="average") ### cluster the data
modSet1$height <- round(modSet1$height, 6) # round the height, this is necessary in order to cut the tree
modSet1 <- cutree(modSet1,k=12) # cut the tree
nClus <- length(unique(modSet1)) # find the number of clusters
### Make an object to describe the variables in each cluster
Clus <- matrix(nrow=nClus, ncol=7)
colnames(Clus) <- c("Cluster", "Mean", "Min", "Median", "Max","N.items", "Comp.Cases")
for (i in 1:nClus){
  print(i)
  aa <- Nobs - miss.total.topDF$Nmiss[miss.total.topDF$Col %in% names(modSet1[modSet1 == i])]
  Clus[i,1] <- i
  Clus[i,2] <- mean(aa)
  Clus[i,3] <- min(aa)
  Clus[i,4] <- median(aa)
  Clus[i,5] <- max(aa)
  Clus[i,6] <- length(aa)
  bb <- bd.top[, colnames(bd.top) %in% names(modSet1[modSet1 == i])]
  # print(dim(bb))
  # print(dim(bb[complete.cases(bb), ]))
  if (length(modSet1[modSet1 == i]) == 1){
    Clus[i,7] <- 1} else {
      Clus[i,7] <- dim(bb[complete.cases(bb), ])[1]
    }
}
Clus <- Clus[order(Clus[,2]),]
Clus
#       Cluster      Mean    Min   Median    Max N.items Comp.Cases
# [1,]       9 116487.8 102975 120047.0 120981      11        720
# [2,]      11 119543.5 116105 119543.5 122982       2      84297
# [3,]       8 144905.0 130091 151979.0 152905       5          0
# [4,]      10 168629.6 142150 171563.0 172713      31     123961
# [5,]       6 203219.5 199347 203219.5 207092       2      83624
# [6,]       7 227301.6 226958 226958.0 228676       5     226958
# [7,]       4 270151.0 259833 272749.0 277009       6      43729
# [8,]       3 288447.6 284831 287170.0 298454       7     133756
# [9,]       2 344260.2 328840 346515.0 353330      12      66114
# [10,]       5 390229.1 375654 389625.5 415237      14      73959
# [11,]      12 428563.4 425232 429590.0 430098       9     383718
# [12,]       1 490068.1 448376 495155.0 502536     250     227415
### See any variables being dropped
# data.dict$Field[data.dict$UDI %in% names(modSet1[modSet1 == 8])]

########################### Dummy-Code ###########################
### Create the dummy-coded variable sets in each cluster
### Uses command "dummy_cols" to create the dummy-coded categorical variables in a format required by
### the Group-LASSO. Extra "data wrangling" steps included here to keep track of what "group" each
### variable belongs to so that the dummy-coded variables can be re-collapsed for clustering again
  for (i in 1:nClus){
    cat("Set ", i, "\n")
    mod1 <- modSet1
    bb <- bd.top[, colnames(bd.top) %in% names(mod1[mod1 == i])] ### subset the data for that cluster
    bb <- bb[complete.cases(bb), ] ### use complete cases
    factors <- sapply(bb, is.factor) ### identify which are factors (all will be ordered, by default)
    factor.names <- names(factors[factors])
    classes <- sapply(sapply(bb, class), function(x) x[1])
    cat("Complete case dimensions are:", dim(bb), "\n")
    if (length(grep("character", classes)) > 0){
      cat("Danger, Will Robinson! For i = ", i, "\n")
    } else if (length(grep("ordered", classes)) < 1) {
      cat("No factors in i = ", i, "\n")
      eval(parse(text = paste("bdSet", i, " <- ", "bb", sep="")))
      cc <- data.frame(group=seq(1:(dim(bb)[2])), item=sapply(strsplit(colnames(bb), "\\."), "[", 2), orig=colnames(bb))
      eval(parse(text = paste("rownames(bdSet", i, ") <- rownames(bb)", sep="")))
      eval(parse(text = paste("dd <- data.frame(catCol=colnames(bdSet", i, "), item=sapply(strsplit(colnames(bdSet",i, "), \"\\\\.\"), \"[\", 2))", sep="")))
      eval(parse(text = paste("group",i," <- merge(dd, cc, by=\"item\", all.x=TRUE)", sep="")))
    } else {
      bb[,factors] <- sapply(bb[,factors], factor, ordered=FALSE) ### remove the ordered from the factors
      ### Make the dummy-variable amplified matrix
      eval(parse(text = paste("bdSet",i, " <- ", "dummy_cols(bb, remove_first_dummy = TRUE, remove_most_frequent_dummy = FALSE, ignore_na = TRUE)", sep="")))
      ### drop the original factor variables
      eval(parse(text = paste("bdSet",i, " <- bdSet", i, "[, -which(colnames(bdSet", i,") %in% factor.names)]", sep="")))
      ### Assign the rownames (the identifiers)
      eval(parse(text = paste("rownames(bdSet", i, ") <- rownames(bb)", sep="")))
      cc <- data.frame(group=seq(1:(dim(bb)[2])), item=sapply(strsplit(colnames(bb), "\\."), "[", 2), orig=colnames(bb))
      eval(parse(text = paste("dd <- data.frame(catCol=colnames(bdSet",i, "), item=sapply(strsplit(colnames(bdSet",i, "), \"\\\\.\"), \"[\", 2))", sep="")))
      eval(parse(text = paste("group",i," <- merge(dd, cc, by=\"item\", all.x=TRUE)", sep="")))
      ### Print the dimensions of the groupi object, the first of which gives the new number of columns, after dummy expansion
      eval(parse(text = paste("print(dim(group", i, "))", sep="")))
    }
  }

# Set  1 
# Complete case dimensions are: 227415 250 
# [1] 550   4
# Set  2 
# Complete case dimensions are: 66114 12 
# [1] 12  4
# Set  3 
# Complete case dimensions are: 133756 7 
# [1] 19  4
# Set  4 
# Complete case dimensions are: 43729 6 
# [1] 16  4
# Set  5 
# Complete case dimensions are: 73959 14 
# [1] 26  4
# Set  6 
# Complete case dimensions are: 83624 2 
# [1] 7 4
# Set  7 
# Complete case dimensions are: 226958 5 
# [1] 17  4
# Set  8 
# Complete case dimensions are: 0 5 
# [1] 5 4
# Set  9 
# Complete case dimensions are: 720 11 
# [1] 20  4
# Set  10 
# Complete case dimensions are: 123961 31 
# [1] 120   4
# Set  11 
# Complete case dimensions are: 84297 2 
# [1] 6 4
# Set  12 
# Complete case dimensions are: 383718 9 
# No factors in i =  12 

########################### Group-LASSO ###########################
### This initial Group-LASSO run is the most computationally expensive, it's recommended not
### to run interactively, but to run via separate script
### Uses command "cv.grpreg" to fit the cross-validated Group-LASSO model
### Saves the output from each model so that they can be re-loaded later
rundate<-Sys.Date()
sets <- c("score", "AudC", "AudP")
setL <- c("AUDIT_T", "AUDIT_C", "AUDIT_P")

for (j in 1:length(setL)){
  for (i in 1:nClus){
    eval(parse(text = paste("bdSet <- bdSet", i, sep=""))) # set the data
    eval(parse(text = paste("group <- group", i, sep=""))) # set the group
    
    if (dim(bdSet)[1] ==0){
      cat("Set ", i, " has zero complete cases...", "\n")
    } else {
      cat("Dimensions of complete case matrix: ", dim(bdSet), "\n")
      aud <- AUDIT[AUDIT$f.eid %in% rownames(bdSet), c("f.eid",sets[j])]
      grpSet <- merge(aud, bdSet, by.x="f.eid", by.y="row.names") # merge with AUDIT
      cat("Dimensions after merging with AUDIT: ", dim(grpSet), "\n")
      if (dim(grpSet)[1] == 0){
        cat("No observations with AUDIT in this set. Abort.", "\n")
      } else
        groupdf <- group # keep the grouping dataframe
      group <- group[match(colnames(grpSet[3:dim(grpSet)[2]]), as.character(group$catCol)),] # order the group set
      X <- as.matrix(grpSet[, 3:(dim(grpSet)[2])]) # put data into a matrix
      y <- grpSet[, colnames(grpSet) %in% sets[j]] # set the AUDIT score outcome
      group <- group$group # set the group vector
      cv.testG <- cv.grpreg(X=X, y=y, group=group, penalty="grLasso", family="gaussian", nfolds=5, seed=51220) # cross-validated model
      retained <- unique(as.character(
        groupdf$orig[as.character(groupdf$catCol) %in% names(coef(cv.testG, s=lambda.min)[coef(cv.testG, s=lambda.min) != 0]) ])) # variables retained
      eval(parse(text = paste("cat(length(retained), \" variables retained in set ",i,"\", \"\n\")", sep="")))
      eval(parse(text = paste("retained", i, "<-retained", sep="")))
      #eval(parse(text = paste("cor", setL[[j]], "_", i, "<- cor(y, predict(cv.testG, X, lambda=cv.testG$lambda.min))", sep="")))
      eval(parse(text = paste("save(retained", i, ", file=paste(PubOutdir, \"GrpLas_retained_", setL[j], "_", i, "_Rd1_", rundate, ".Rds\", sep=\"\"))", sep="")))
      cat("Saved",setL[j], "Set", i, "\n")
    }  
  }
}

# Dimensions of complete case matrix:  227415 550 
# Dimensions after merging with AUDIT:  80839 552 
# 28  variables retained in set 1 
# Saved AUDIT_T Set 1 
# Dimensions of complete case matrix:  66114 12 
# Dimensions after merging with AUDIT:  19292 14 
# 9  variables retained in set 2 
# Saved AUDIT_T Set 2 
# Dimensions of complete case matrix:  133756 19 
# Dimensions after merging with AUDIT:  44941 21 
# 7  variables retained in set 3 
# Saved AUDIT_T Set 3 
# Dimensions of complete case matrix:  43729 16 
# Dimensions after merging with AUDIT:  16436 18 
# 1  variables retained in set 4 
# Saved AUDIT_T Set 4 
# Dimensions of complete case matrix:  73959 26 
# Dimensions after merging with AUDIT:  24888 28 
# 9  variables retained in set 5 
# Saved AUDIT_T Set 5 
# Dimensions of complete case matrix:  83624 7 
# Dimensions after merging with AUDIT:  27852 9 
# 2  variables retained in set 6 
# Saved AUDIT_T Set 6 
# Dimensions of complete case matrix:  226958 17 
# Dimensions after merging with AUDIT:  67766 19 
# 5  variables retained in set 7 
# Saved AUDIT_T Set 7 
# Set  8  has zero complete cases... 
# Dimensions of complete case matrix:  720 20 
# Dimensions after merging with AUDIT:  225 22 
# 5  variables retained in set 9 
# Saved AUDIT_T Set 9 
# Dimensions of complete case matrix:  123961 120 
# Dimensions after merging with AUDIT:  45638 122 
# 22  variables retained in set 10 
# Saved AUDIT_T Set 10 
# Dimensions of complete case matrix:  84297 6 
# Dimensions after merging with AUDIT:  31086 8 
# 2  variables retained in set 11 
# Saved AUDIT_T Set 11 
# Dimensions of complete case matrix:  383718 9 
# Dimensions after merging with AUDIT:  120601 11 
# 9  variables retained in set 12 
# Saved AUDIT_T Set 12 
# Dimensions of complete case matrix:  227415 550 
# Dimensions after merging with AUDIT:  80839 552 
# 24  variables retained in set 1 
# Saved AUDIT_C Set 1 
# Dimensions of complete case matrix:  66114 12 
# Dimensions after merging with AUDIT:  19292 14 
# 11  variables retained in set 2 
# Saved AUDIT_C Set 2 
# Dimensions of complete case matrix:  133756 19 
# Dimensions after merging with AUDIT:  44941 21 
# 7  variables retained in set 3 
# Saved AUDIT_C Set 3 
# Dimensions of complete case matrix:  43729 16 
# Dimensions after merging with AUDIT:  16436 18 
# 5  variables retained in set 4 
# Saved AUDIT_C Set 4 
# Dimensions of complete case matrix:  73959 26 
# Dimensions after merging with AUDIT:  24888 28 
# 8  variables retained in set 5 
# Saved AUDIT_C Set 5 
# Dimensions of complete case matrix:  83624 7 
# Dimensions after merging with AUDIT:  27852 9 
# 2  variables retained in set 6 
# Saved AUDIT_C Set 6 
# Dimensions of complete case matrix:  226958 17 
# Dimensions after merging with AUDIT:  67766 19 
# 5  variables retained in set 7 
# Saved AUDIT_C Set 7 
# Set  8  has zero complete cases... 
# Dimensions of complete case matrix:  720 20 
# Dimensions after merging with AUDIT:  225 22 
# 6  variables retained in set 9 
# Saved AUDIT_C Set 9 
# Dimensions of complete case matrix:  123961 120 
# Dimensions after merging with AUDIT:  45638 122 
# 28  variables retained in set 10 
# Saved AUDIT_C Set 10 
# Dimensions of complete case matrix:  84297 6 
# Dimensions after merging with AUDIT:  31086 8 
# 1  variables retained in set 11 
# Saved AUDIT_C Set 11 
# Dimensions of complete case matrix:  383718 9 
# Dimensions after merging with AUDIT:  120601 11 
# 9  variables retained in set 12 
# Saved AUDIT_C Set 12 
# Dimensions of complete case matrix:  227415 550 
# Dimensions after merging with AUDIT:  80839 552 
# 52  variables retained in set 1 
# Saved AUDIT_P Set 1 
# Dimensions of complete case matrix:  66114 12 
# Dimensions after merging with AUDIT:  19292 14 
# 12  variables retained in set 2 
# Saved AUDIT_P Set 2 
# Dimensions of complete case matrix:  133756 19 
# Dimensions after merging with AUDIT:  44941 21 
# 7  variables retained in set 3 
# Saved AUDIT_P Set 3 
# Dimensions of complete case matrix:  43729 16 
# Dimensions after merging with AUDIT:  16436 18 
# 0  variables retained in set 4 
# Saved AUDIT_P Set 4 
# Dimensions of complete case matrix:  73959 26 
# Dimensions after merging with AUDIT:  24888 28 
# 8  variables retained in set 5 
# Saved AUDIT_P Set 5 
# Dimensions of complete case matrix:  83624 7 
# Dimensions after merging with AUDIT:  27852 9 
# 2  variables retained in set 6 
# Saved AUDIT_P Set 6 
# Dimensions of complete case matrix:  226958 17 
# Dimensions after merging with AUDIT:  67766 19 
# 5  variables retained in set 7 
# Saved AUDIT_P Set 7 
# Set  8  has zero complete cases... 
# Dimensions of complete case matrix:  720 20 
# Dimensions after merging with AUDIT:  225 22 
# 5  variables retained in set 9 
# Saved AUDIT_P Set 9 
# Dimensions of complete case matrix:  123961 120 
# Dimensions after merging with AUDIT:  45638 122 
# 21  variables retained in set 10 
# Saved AUDIT_P Set 10 
# Dimensions of complete case matrix:  84297 6 
# Dimensions after merging with AUDIT:  31086 8 
# 2  variables retained in set 11 
# Saved AUDIT_P Set 11 
# Dimensions of complete case matrix:  383718 9 
# Dimensions after merging with AUDIT:  120601 11 
# 9  variables retained in set 12 
# Saved AUDIT_P Set 12 

#################################### Reload Objects ####################################
### Reload the computed Group-LASSO fit objects
suffixes <- c("T", "C", "P")
for (j in 1:length(suffixes)){
    ToRem <-ls(pattern="retain")
    rm(list=ToRem)
    last.rundate <- "2021-03-03"
    ### After this is run, reload the objects
    ToLoad <- list.files(paste(PubOutdir, sep=""), pattern="\\.Rds")
    ### Reload all 3 sets, total, consumption, and problems
    ToLoad <- ToLoad[grep(paste("GrpLas_retained_AUDIT_", suffixes[j], sep=""), ToLoad)]
    ToLoad <- ToLoad[grep("Rd1", ToLoad)]
    ToLoad <- ToLoad[grep(last.rundate, ToLoad)]
    ToLoad <- paste(PubOutdir, ToLoad, sep="")
    lapply(ToLoad,load,.GlobalEnv)
    
    ### Combine the retained variables into a single vector

    ret.list <- ls(pattern="retained")
    retained.all <-character()
    for (i in 1:length(ret.list)){
      if (i==1){
        retained.all <- retained1
        }else{
        eval(parse(text = paste("retained.all <- c(retained.all,",ret.list[i],")", sep="")))
      }
    }
    eval(parse(text = paste("ret", suffixes[j], "<-retained.all", sep="")))
    #eval(parse(text = paste("print(length(ret", suffixes[j], "))", sep="")))
    eval(parse(text = paste("cat(length(ret", suffixes[j], "), \" variables retained in set ",suffixes[j],"\n\")", sep="")))
}

### 99  variables retained in set T
### 106  variables retained in set C
### 123  variables retained in set P
    
### Subset the dataset and the observation count matrix for just the retained variables
  ### Do this once for each of the 3 sets (total, consump, prob)
  ### Total
  bd.grpT <- bd.top[, colnames(bd.top) %in% retT]
  grpT.index <- which(colnames(bd.miss.pairwise.topS) %in% retT)
  bd.miss.pairwise.topST <- bd.miss.pairwise.topS[grpT.index, grpT.index] 
  ### Consumption
  bd.grpC <- bd.top[, colnames(bd.top) %in% retC]
  grpC.index <- which(colnames(bd.miss.pairwise.topS) %in% retC)
  bd.miss.pairwise.topSC <- bd.miss.pairwise.topS[grpC.index, grpC.index]
  ### Problems
  bd.grpP <- bd.top[, colnames(bd.top) %in% retP]
  grpP.index <- which(colnames(bd.miss.pairwise.topS) %in% retP)
  bd.miss.pairwise.topSP <- bd.miss.pairwise.topS[grpP.index, grpP.index]

######################################################################## Cluster/Lasso Round 2 ####################################
### From this point, clustering and Group-LASSO application must be done individually for Total/Consumpt/Probs
### because these each retained different sets of variables from the initial step
  
  ################# Clustering ################################
  #############################################################
  ### Total
  dist.pairT <- dist(bd.miss.pairwise.topST, method="euclidean")
  modT <- hclust(dist.pairT, method="average")
    modT$height <- round(modT$height, 6) # this is necessary in order to cut the tree
    modT <- cutree(modT, k=5)
    nClus <- length(unique(modT))
    Clus <- matrix(nrow=nClus, ncol=7)
    colnames(Clus) <- c("Cluster", "Mean", "Min", "Median", "Max","N.items", "Comp.Cases")
    for (i in 1:nClus){
      aa <- dim(bd.top)[1] - miss.total.topDF$Nmiss[miss.total.topDF$Col %in% names(modT[modT == i])]
      Clus[i,1] <- i
      Clus[i,2] <- mean(aa)
      Clus[i,3] <- min(aa)
      Clus[i,4] <- median(aa)
      Clus[i,5] <- max(aa)
      Clus[i,6] <- length(aa)
      bb <- bd.top[, colnames(bd.top) %in% names(modT[modT == i])]
      Clus[i,7] <- dim(bb[complete.cases(bb), ])[1]
    }
    Clus <- Clus[order(Clus[,2]),]
    nClusT <- nClus
    ClusT <- Clus
    ClusT
    ### For total, use k=5 ###confirmed 3/5/2021
    # Cluster     Mean    Min   Median    Max N.items Comp.Cases
    # [1,]       4 119694.6 115483 120981.0 120981       5     114595
    # [2,]       5 165235.9 116105 172138.0 172713      24      71196
    # [3,]       3 220421.0 199347 226958.0 228676       7      41434
    # [4,]       2 317754.4 277009 329677.0 353330      17      13389
    # [5,]       1 455899.9 375654 467790.5 502536      46      46292
    
    #############################################################
    ### Consumption
    dist.pairC <- dist(bd.miss.pairwise.topSC, method="euclidean")
      modC <- hclust(dist.pairC, method="average")
      modC$height <- round(modC$height, 6) # this is necessary in order to cut the tree
      modC <- cutree(modC, k=6)
      nClus <- length(unique(modC))
      Clus <- matrix(nrow=nClus, ncol=7)
      colnames(Clus) <- c("Cluster", "Mean", "Min", "Median", "Max","N.items", "Comp.Cases")
      for (i in 1:nClus){
        aa <- dim(bd.top)[1] - miss.total.topDF$Nmiss[miss.total.topDF$Col %in% names(modC[modC == i])]
        Clus[i,1] <- i
        Clus[i,2] <- mean(aa)
        Clus[i,3] <- min(aa)
        Clus[i,4] <- median(aa)
        Clus[i,5] <- max(aa)
        Clus[i,6] <- length(aa)
        bb <- bd.top[, colnames(bd.top) %in% names(modC[modC == i])]
        Clus[i,7] <- dim(bb[complete.cases(bb), ])[1]
      }
      Clus <- Clus[order(Clus[,2]),]
      nClusC <- nClus
      ClusC <- Clus
      ClusC
      ### For consumption, use k=6 ###confirmed 3/5/2021
      # Cluster     Mean    Min   Median    Max N.items Comp.Cases
      # [1,]       5 116908.0 102975 120514.0 120981       6      21258
      # [2,]       6 166860.7 122982 171563.0 172713      29      89440
      # [3,]       4 220421.0 199347 226958.0 228676       7      41434
      # [4,]       3 281200.6 259833 285999.5 298454      12      17423
      # [5,]       2 364005.2 328840 353330.0 415237      19      16697
      # [6,]       1 468807.6 425232 469584.0 502536      33     246876

      #############################################################
      ### Problems
      dist.pairP <- dist(bd.miss.pairwise.topSP, method="euclidean")
      modP <- hclust(dist.pairP, method="average")
      modP$height <- round(modP$height, 4) # this is necessary in order to cut the tree
      modP <- cutree(modP, k=4)
      nClus <- length(unique(modP))
      Clus <- matrix(nrow=nClus, ncol=7)
      colnames(Clus) <- c("Cluster", "Mean", "Min", "Median", "Max","N.items", "Comp.Cases")
      for (i in 1:nClus){
        aa <- dim(bd.top)[1] - miss.total.topDF$Nmiss[miss.total.topDF$Col %in% names(modP[modP == i])]
        Clus[i,1] <- i
        Clus[i,2] <- mean(aa)
        Clus[i,3] <- min(aa)
        Clus[i,4] <- median(aa)
        Clus[i,5] <- max(aa)
        Clus[i,6] <- length(aa)
        bb <- bd.top[, colnames(bd.top) %in% names(modP[modP == i])]
        Clus[i,7] <- dim(bb[complete.cases(bb), ])[1]
      }
      Clus <- Clus[order(Clus[,2]),]
      nClusP <- nClus
      ClusP <- Clus
      ClusP
      ### For problems, use k=4 ###confirmed 3/5/2021   
      # Cluster     Mean    Min Median    Max N.items Comp.Cases
      # [1,]       4 156937.4 115483 170606 172713      28      17373
      # [2,]       3 220421.0 199347 226958 228676       7      41434
      # [3,]       2 323697.7 284831 336769 353330      19      29090
      # [4,]       1 469228.2 375654 478200 502536      69      47643
      
      
################# Dummy-Coding ################################
suffixes <- c("T", "C", "P")
Mods <- list(modT, modC, modP)
Clusters <- list(nClusT, nClusC, nClusP)

for (k in 1:length(suffixes)){
  for (i in 1:Clusters[[k]]){
    cat(suffixes[[k]], ", Set ", i, "\n")
    mod1 <- Mods[[k]]
    
    bb <- bd.top[, colnames(bd.top) %in% names(mod1[mod1 == i])] ### subset the data for that cluster
    bb <- bb[complete.cases(bb), ] ### use complete cases
    factors <- sapply(bb, is.factor) ### identify which are factors (all will be ordered, by default)
    factor.names <- names(factors[factors])
    classes <- sapply(sapply(bb, class), function(x) x[1])
    cat("Complete case dimensions are:", dim(bb), "\n")
    if (length(grep("character", classes)) > 0){
      cat("Danger, Will Robinson! For i = ", i, "\n")
    } else if (length(grep("ordered", classes)) < 1) {
      cat("No factors in i = ", i, "\n")
      eval(parse(text = paste("bdSet",suffixes[k], i, " <- ", "bb", sep="")))
      cc <- data.frame(group=seq(1:(dim(bb)[2])), item=sapply(strsplit(colnames(bb), "\\."), "[", 2), orig=colnames(bb))
      eval(parse(text = paste("rownames(bdSet",suffixes[k],  i, ") <- rownames(bb)", sep="")))
      eval(parse(text = paste("dd <- data.frame(catCol=colnames(bdSet",suffixes[k], i, "), item=sapply(strsplit(colnames(bdSet",suffixes[k], i, "), \"\\\\.\"), \"[\", 2))", sep="")))
      eval(parse(text = paste("group",suffixes[k], i," <- merge(dd, cc, by=\"item\", all.x=TRUE)", sep="")))
    } else {
      bb[,factors] <- sapply(bb[,factors], factor, ordered=FALSE) ### remove the ordered from the factors
      ### Make the dummy-variable amplified matrix
      eval(parse(text = paste("bdSet",suffixes[k], i, " <- ", "dummy_cols(bb, remove_first_dummy = TRUE, remove_most_frequent_dummy = FALSE, ignore_na = TRUE)", sep="")))
      ### drop the original factor variables
      eval(parse(text = paste("bdSet",suffixes[k], i, " <- bdSet",suffixes[k], i, "[, -which(colnames(bdSet",suffixes[k], i,") %in% factor.names)]", sep="")))
      ### Assign the rownames (the identifiers)
      eval(parse(text = paste("rownames(bdSet",suffixes[k],  i, ") <- rownames(bb)", sep="")))
      cc <- data.frame(group=seq(1:(dim(bb)[2])), item=sapply(strsplit(colnames(bb), "\\."), "[", 2), orig=colnames(bb))
      eval(parse(text = paste("dd <- data.frame(catCol=colnames(bdSet",suffixes[k], i, "), item=sapply(strsplit(colnames(bdSet",suffixes[k], i, "), \"\\\\.\"), \"[\", 2))", sep="")))
      eval(parse(text = paste("group",suffixes[k], i," <- merge(dd, cc, by=\"item\", all.x=TRUE)", sep="")))
      ### Print the dimensions of the groupi object, the first of which gives the new number of columns, after dummy expansion
      eval(parse(text = paste("print(dim(group",suffixes[k],  i, "))", sep="")))
    }
  }
}

# T , Set  1 
# Complete case dimensions are: 46292 46 
# [1] 69  4
# T , Set  2 
# Complete case dimensions are: 13389 17 
# [1] 29  4
# T , Set  3 
# Complete case dimensions are: 41434 7 
# [1] 24  4
# T , Set  4 
# Complete case dimensions are: 114595 5 
# [1] 9 4
# T , Set  5 
# Complete case dimensions are: 71196 24 
# [1] 102   4
# C , Set  1 
# Complete case dimensions are: 246876 33 
# [1] 43  4
# C , Set  2 
# Complete case dimensions are: 16697 19 
# [1] 31  4
# C , Set  3 
# Complete case dimensions are: 17423 12 
# [1] 34  4
# C , Set  4 
# Complete case dimensions are: 41434 7 
# [1] 24  4
# C , Set  5 
# Complete case dimensions are: 21258 6 
# [1] 12  4
# C , Set  6 
# Complete case dimensions are: 89440 29 
# [1] 119   4
# P , Set  1 
# Complete case dimensions are: 47643 69 
# [1] 119   4
# P , Set  2 
# Complete case dimensions are: 29090 19 
# [1] 31  4
# P , Set  3 
# Complete case dimensions are: 41434 7 
# [1] 24  4
# P , Set  4 
# Complete case dimensions are: 17373 28 
# [1] 102   4

################# Group-LASSO ################################
rundate<-Sys.Date()
sets <- c("score", "AudC", "AudP" )
setL <- c("AUDIT_T", "AUDIT_C", "AUDIT_P")
setI <- which(colnames(AUDIT) == sets)
suffixes <- c("T", "C", "P")
Clusters <- list(nClusT, nClusC, nClusP)
for (j in 1:length(sets)){
  for (i in 1:Clusters[[j]]){
    eval(parse(text = paste("bdSet <- bdSet",suffixes[j], i, sep=""))) # set the data
    eval(parse(text = paste("group <- group",suffixes[j], i, sep=""))) # set the group
    
    if (dim(bdSet)[1] ==0){
      cat("Set ", i, " has zero complete cases...", "\n")
    } else {
      cat("Dimensions of complete case matrix: ", dim(bdSet), "\n")
      aud <- AUDIT[AUDIT$f.eid %in% rownames(bdSet), c("f.eid",sets[j])]
      grpSet <- merge(aud, bdSet, by.x="f.eid", by.y="row.names") # merge with AUDIT
      cat("Dimensions after merging with AUDIT: ", dim(grpSet), "\n")
      if (dim(grpSet)[1] == 0){
        cat("No observations with AUDIT in this set. Abort.", "\n")
      } else
        groupdf <- group # keep the grouping dataframe
      group <- group[match(colnames(grpSet[3:dim(grpSet)[2]]), as.character(group$catCol)),] # order the group set
      X <- as.matrix(grpSet[, 3:(dim(grpSet)[2])]) # put data into a matrix
      y <- grpSet[, colnames(grpSet) %in% sets[j]] # set the AUDIT score outcome
      group <- group$group # set the group vector
      cv.testG <- cv.grpreg(X=X, y=y, group=group, penalty="grLasso", family="gaussian", nfolds=5, seed=51220) # cross-validated model
      retained <- unique(as.character(
        groupdf$orig[as.character(groupdf$catCol) %in% names(coef(cv.testG, s=lambda.min)[coef(cv.testG, s=lambda.min) != 0]) ])) # variables retained
      eval(parse(text = paste("cat(length(retained), \" variables retained in set ",i,"\", \"\n\")", sep="")))
      eval(parse(text = paste("retained", i, "<-retained", sep="")))
      #eval(parse(text = paste("cor", setL[[j]], "_", i, "<- cor(y, predict(cv.testG, X, lambda=cv.testG$lambda.min))", sep="")))
      eval(parse(text = paste("save(retained", i, ", file=paste(PubOutdir, \"GrpLas_retained_", setL[j], i, "_Rd2_", rundate, ".Rds\", sep=\"\"))", sep="")))
      cat("Saved",setL[j], "_", i, "\n")
    }  
  }
}

# Dimensions of complete case matrix:  46292 69 
# Dimensions after merging with AUDIT:  16363 71 
# 21  variables retained in set 1 
# Saved AUDIT_T _ 1 
# Dimensions of complete case matrix:  13389 29 
# Dimensions after merging with AUDIT:  4031 31 
# 13  variables retained in set 2 
# Saved AUDIT_T _ 2 
# Dimensions of complete case matrix:  41434 24 
# Dimensions after merging with AUDIT:  12680 26 
# 7  variables retained in set 3 
# Saved AUDIT_T _ 3 
# Dimensions of complete case matrix:  114595 9 
# Dimensions after merging with AUDIT:  36004 11 
# 5  variables retained in set 4 
# Saved AUDIT_T _ 4 
# Dimensions of complete case matrix:  71196 102 
# Dimensions after merging with AUDIT:  27737 104 
# 19  variables retained in set 5 
# Saved AUDIT_T _ 5 
# Dimensions of complete case matrix:  246876 43 
# Dimensions after merging with AUDIT:  85210 45 
# 17  variables retained in set 1 
# Saved AUDIT_C _ 1 
# Dimensions of complete case matrix:  16697 31 
# Dimensions after merging with AUDIT:  5169 33 
# 18  variables retained in set 2 
# Saved AUDIT_C _ 2 
# Dimensions of complete case matrix:  17423 34 
# Dimensions after merging with AUDIT:  6312 36 
# 8  variables retained in set 3 
# Saved AUDIT_C _ 3 
# Dimensions of complete case matrix:  41434 24 
# Dimensions after merging with AUDIT:  12680 26 
# 7  variables retained in set 4 
# Saved AUDIT_C _ 4 
# Dimensions of complete case matrix:  21258 12 
# Dimensions after merging with AUDIT:  6608 14 
# 5  variables retained in set 5 
# Saved AUDIT_C _ 5 
# Dimensions of complete case matrix:  89440 119 
# Dimensions after merging with AUDIT:  33174 121 
# 25  variables retained in set 6 
# Saved AUDIT_C _ 6 
# Dimensions of complete case matrix:  47643 119 
# Dimensions after merging with AUDIT:  16922 121 
# 14  variables retained in set 1 
# Saved AUDIT_P _ 1 
# Dimensions of complete case matrix:  29090 31 
# Dimensions after merging with AUDIT:  8244 33 
# 10  variables retained in set 2 
# Saved AUDIT_P _ 2 
# Dimensions of complete case matrix:  41434 24 
# Dimensions after merging with AUDIT:  12680 26 
# 7  variables retained in set 3 
# Saved AUDIT_P _ 3 
# Dimensions of complete case matrix:  17373 102 
# Dimensions after merging with AUDIT:  6612 104 
# 23  variables retained in set 4 
# Saved AUDIT_P _ 4 

#################################### Reload Objects from Round 2 ####################################
### After this is run, reload the objects if needed
suffixes <- c("T", "C", "P")
for (j in 1:length(suffixes)){
  ToRem <-ls(pattern="retain")
  rm(list=ToRem)
  last.rundate <- "2021-03-10"
  ### After this is run, reload the objects
  ToLoad <- list.files(paste(PubOutdir, sep=""), pattern="\\.Rds")
  ### Reload all 3 sets, total, consumption, and problems
  ToLoad <- ToLoad[grep(paste("GrpLas_retained_AUDIT_", suffixes[j], sep=""), ToLoad)]
  ToLoad <- ToLoad[grep("Rd2", ToLoad)]
  ToLoad <- ToLoad[grep(last.rundate, ToLoad)]
  ToLoad <- paste(PubOutdir, ToLoad, sep="")
  lapply(ToLoad,load,.GlobalEnv)
  
  ### Combine the retained variables into a single vector
  
  ret.list <- ls(pattern="retained")
  retained.all <-character()
  for (i in 1:length(ret.list)){
    if (i==1){
      retained.all <- retained1
    }else{
      eval(parse(text = paste("retained.all <- c(retained.all,",ret.list[i],")", sep="")))
    }
  }
  eval(parse(text = paste("ret", suffixes[j], "<-retained.all", sep="")))
  #eval(parse(text = paste("print(length(ret", suffixes[j], "))", sep="")))
  eval(parse(text = paste("cat(length(ret", suffixes[j], "), \" variables retained in set ",suffixes[j],"\n\")", sep="")))
}
### 65  variables retained in set T
### 80  variables retained in set C
### 54  variables retained in set P ### verified 3/1/2021
   
 ### Subset the dataset and the observation count matrix for just the retained variables
 ### Do this once for each of the 3 sets (total, consump, prob)

 ### Total
 bd.grpT <- bd.top[, colnames(bd.top) %in% retT]
 grpT.index <- which(colnames(bd.miss.pairwise.topS) %in% retT)
 bd.miss.pairwise.topST <- bd.miss.pairwise.topS[grpT.index, grpT.index] 
 ### Consumption
 bd.grpC <- bd.top[, colnames(bd.top) %in% retC]
 grpC.index <- which(colnames(bd.miss.pairwise.topS) %in% retC)
 bd.miss.pairwise.topSC <- bd.miss.pairwise.topS[grpC.index, grpC.index]
 ### Problems
 bd.grpP <- bd.top[, colnames(bd.top) %in% retP]
 grpP.index <- which(colnames(bd.miss.pairwise.topS) %in% retP)
 bd.miss.pairwise.topSP <- bd.miss.pairwise.topS[grpP.index, grpP.index] 
 
 ##########################################################################################################################################
 ### Now that's we're down to a reasonable number of variables in each set, we can assess remaining measures for balance between
 ### the "prediction" and the "training" sets to prioritize measures for the final model according to their balance between sets
 
 ########### TOTAL ################
 naSum.T <- naSum[rownames(naSum) %in% retT,]
 naSum.T <- naSum.T[order(naSum.T$nAll, decreasing=TRUE),]
 Rnames <- data.frame(rownames(naSum.T))
 Rnames <- data.frame(Rnames[seq(dim(Rnames)[1],1),])
 colnames(Rnames) <- "naSum.T.rownames"
 naSum.Ta <- merge(Rnames, naSum.T, by.x="naSum.T.rownames", by.y="row.names", sort=FALSE)
 rownames(naSum.Ta) <- naSum.Ta$naSum.T.rownames
 
 for (i in 1:dim(naSum.T)[1]){
   ### Decreasing order of completeness
   aa <- bd.grpT[, colnames(bd.grpT) %in% rownames(naSum.T)[1:i]] ### Starting with the var with most obs, add one var each iteration
   bb <- bd.grpT[, -which(colnames(bd.grpT) %in% rownames(naSum.T)[i])]
   ### Increasing order of completeness
   aa2 <- bd.grpT[, colnames(bd.grpT) %in% rownames(naSum.Ta)[1:i]] ### Starting with the var with least obs, add one var each iteration
   if (i == 1) {
     ### Decreasing order
     names(aa) <- rownames(bd.grpT)
     naSum.T$CompCase.Add[i] <- length(aa) - sum(is.na(aa))
     naSum.T$CompCase.Add.Obs[i] <- sum(names(aa[!is.na(aa)]) %in% rownames(AUDIT))
     naSum.T$CompCase.Add.UnObs[i] <- naSum.T$CompCase.Add[i] - naSum.T$CompCase.Add.Obs[i]
     naSum.T$CompCase.LeaveOut[i] <- dim(bb[complete.cases(bb),])[1]
     ### Increasing order
     names(aa2) <- rownames(bd.grpT)
     naSum.Ta$CompCase.Add.I[i] <- length(aa2) - sum(is.na(aa2))
     naSum.Ta$CompCase.Add.Obs.I[i] <- sum(names(aa2[!is.na(aa2)]) %in% rownames(AUDIT))
     naSum.Ta$CompCase.Add.UnObs.I[i] <- naSum.Ta$CompCase.Add.I[i] - naSum.Ta$CompCase.Add.Obs.I[i]
   } else{
     ### Decreasing Order
     rownames(aa) <- rownames(bd.grpT)
     naSum.T$CompCase.Add[i] <- dim(aa[complete.cases(aa),])[1]
     naSum.T$CompCase.Add.Obs[i] <- sum(rownames(aa[complete.cases(aa),]) %in% rownames(AUDIT))
     naSum.T$CompCase.Add.UnObs[i] <- naSum.T$CompCase.Add[i] - naSum.T$CompCase.Add.Obs[i]
     naSum.T$CompCase.LeaveOut[i] <- dim(bb[complete.cases(bb),])[1]
     ### Increasing Order
     rownames(aa2) <- rownames(bd.grpT)
     naSum.Ta$CompCase.Add.I[i] <- dim(aa2[complete.cases(aa2),])[1]
     naSum.Ta$CompCase.Add.Obs.I[i] <- sum(rownames(aa2[complete.cases(aa2),]) %in% rownames(AUDIT))
     naSum.Ta$CompCase.Add.UnObs.I[i] <- naSum.Ta$CompCase.Add.I[i] - naSum.Ta$CompCase.Add.Obs.I[i]
   }
 }
 
 naSum.T[,11:14]
 naSum.Ta[,12:14]
 naSum.Ta <- naSum.Ta[,12:14]
 naSum.T <- merge(naSum.T, naSum.Ta, by="row.names", sort=FALSE)
 save(naSum.T, file=paste(PubOutdir, "MissingnessComparison_TotalSet.Rds", sep=""))
 ########### CONSUMPTION ################
 naSum.C <- naSum[rownames(naSum) %in% retC,]
 naSum.C <- naSum.C[order(naSum.C$nAll, decreasing=TRUE),]
 Rnames <- data.frame(rownames(naSum.C))
 Rnames <- data.frame(Rnames[seq(dim(Rnames)[1],1),])
 colnames(Rnames) <- "naSum.C.rownames"
 naSum.Ca <- merge(Rnames, naSum.C, by.x="naSum.C.rownames", by.y="row.names", sort=FALSE)
 rownames(naSum.Ca) <- naSum.Ca$naSum.C.rownames
 
 for (i in 1:dim(naSum.C)[1]){
   ### Decreasing order of completeness
   aa <- bd.grpC[, colnames(bd.grpC) %in% rownames(naSum.C)[1:i]] ### Starting with the var with most obs, add one var each iteration
   bb <- bd.grpC[, -which(colnames(bd.grpC) %in% rownames(naSum.C)[i])]
   ### Increasing order of completeness
   aa2 <- bd.grpC[, colnames(bd.grpC) %in% rownames(naSum.Ca)[1:i]] ### Starting with the var with least obs, add one var each iteration
   if (i == 1) {
     ### Decreasing order
     names(aa) <- rownames(bd.grpC)
     naSum.C$CompCase.Add[i] <- length(aa) - sum(is.na(aa))
     naSum.C$CompCase.Add.Obs[i] <- sum(names(aa[!is.na(aa)]) %in% rownames(AUDIT))
     naSum.C$CompCase.Add.UnObs[i] <- naSum.C$CompCase.Add[i] - naSum.C$CompCase.Add.Obs[i]
     naSum.C$CompCase.LeaveOut[i] <- dim(bb[complete.cases(bb),])[1]
     ### Increasing order
     names(aa2) <- rownames(bd.grpC)
     naSum.Ca$CompCase.Add.I[i] <- length(aa2) - sum(is.na(aa2))
     naSum.Ca$CompCase.Add.Obs.I[i] <- sum(names(aa2[!is.na(aa2)]) %in% rownames(AUDIT))
     naSum.Ca$CompCase.Add.UnObs.I[i] <- naSum.Ca$CompCase.Add.I[i] - naSum.Ca$CompCase.Add.Obs.I[i]
   } else{
     ### Decreasing Order
     rownames(aa) <- rownames(bd.grpC)
     naSum.C$CompCase.Add[i] <- dim(aa[complete.cases(aa),])[1]
     naSum.C$CompCase.Add.Obs[i] <- sum(rownames(aa[complete.cases(aa),]) %in% rownames(AUDIT))
     naSum.C$CompCase.Add.UnObs[i] <- naSum.C$CompCase.Add[i] - naSum.C$CompCase.Add.Obs[i]
     naSum.C$CompCase.LeaveOut[i] <- dim(bb[complete.cases(bb),])[1]
     ### Increasing Order
     rownames(aa2) <- rownames(bd.grpC)
     naSum.Ca$CompCase.Add.I[i] <- dim(aa2[complete.cases(aa2),])[1]
     naSum.Ca$CompCase.Add.Obs.I[i] <- sum(rownames(aa2[complete.cases(aa2),]) %in% rownames(AUDIT))
     naSum.Ca$CompCase.Add.UnObs.I[i] <- naSum.Ca$CompCase.Add.I[i] - naSum.Ca$CompCase.Add.Obs.I[i]
   }
 }
 
 naSum.C[,11:14]
 naSum.Ca[,12:14]
 naSum.Ca <- naSum.Ca[,12:14]
 naSum.C <- merge(naSum.C, naSum.Ca, by="row.names", sort=FALSE)
 save(naSum.C, file=paste(PubOutdir, "MissingnessComparison_ConsumpSet.Rds", sep=""))
 
 ########### PROBLEMS ################
 naSum.P <- naSum[rownames(naSum) %in% retP,]
 naSum.P <- naSum.P[order(naSum.P$nAll, decreasing=TRUE),]
 Rnames <- data.frame(rownames(naSum.P))
 Rnames <- data.frame(Rnames[seq(dim(Rnames)[1],1),])
 colnames(Rnames) <- "naSum.P.rownames"
 naSum.Pa <- merge(Rnames, naSum.P, by.x="naSum.P.rownames", by.y="row.names", sort=FALSE)
 rownames(naSum.Pa) <- naSum.Pa$naSum.P.rownames
 
 for (i in 1:dim(naSum.P)[1]){
   ### Decreasing order of completeness
   aa <- bd.grpP[, colnames(bd.grpP) %in% rownames(naSum.P)[1:i]] ### Starting with the var with most obs, add one var each iteration
   bb <- bd.grpP[, -which(colnames(bd.grpP) %in% rownames(naSum.P)[i])]
   ### Increasing order of completeness
   aa2 <- bd.grpP[, colnames(bd.grpP) %in% rownames(naSum.Pa)[1:i]] ### Starting with the var with least obs, add one var each iteration
   if (i == 1) {
     ### Decreasing order
     names(aa) <- rownames(bd.grpP)
     naSum.P$CompCase.Add[i] <- length(aa) - sum(is.na(aa))
     naSum.P$CompCase.Add.Obs[i] <- sum(names(aa[!is.na(aa)]) %in% rownames(AUDIT))
     naSum.P$CompCase.Add.UnObs[i] <- naSum.P$CompCase.Add[i] - naSum.P$CompCase.Add.Obs[i]
     naSum.P$CompCase.LeaveOut[i] <- dim(bb[complete.cases(bb),])[1]
     ### Increasing order
     names(aa2) <- rownames(bd.grpP)
     naSum.Pa$CompCase.Add.I[i] <- length(aa2) - sum(is.na(aa2))
     naSum.Pa$CompCase.Add.Obs.I[i] <- sum(names(aa2[!is.na(aa2)]) %in% rownames(AUDIT))
     naSum.Pa$CompCase.Add.UnObs.I[i] <- naSum.Pa$CompCase.Add.I[i] - naSum.Pa$CompCase.Add.Obs.I[i]
   } else{
     ### Decreasing Order
     rownames(aa) <- rownames(bd.grpP)
     naSum.P$CompCase.Add[i] <- dim(aa[complete.cases(aa),])[1]
     naSum.P$CompCase.Add.Obs[i] <- sum(rownames(aa[complete.cases(aa),]) %in% rownames(AUDIT))
     naSum.P$CompCase.Add.UnObs[i] <- naSum.P$CompCase.Add[i] - naSum.P$CompCase.Add.Obs[i]
     naSum.P$CompCase.LeaveOut[i] <- dim(bb[complete.cases(bb),])[1]
     ### Increasing Order
     rownames(aa2) <- rownames(bd.grpP)
     naSum.Pa$CompCase.Add.I[i] <- dim(aa2[complete.cases(aa2),])[1]
     naSum.Pa$CompCase.Add.Obs.I[i] <- sum(rownames(aa2[complete.cases(aa2),]) %in% rownames(AUDIT))
     naSum.Pa$CompCase.Add.UnObs.I[i] <- naSum.Pa$CompCase.Add.I[i] - naSum.Pa$CompCase.Add.Obs.I[i]
   }
 }
 
 naSum.P[,11:14]
 naSum.Pa[,12:14]
 naSum.Pa <- naSum.Pa[,12:14]
 naSum.P <- merge(naSum.P, naSum.Pa, by="row.names", sort=FALSE)
 save(naSum.P, file=paste(PubOutdir, "MissingnessComparison_ProbSet.Rds", sep=""))
 ###################
 ### If necessary, reload the objects
 load(paste(PubOutdir, "MissingnessComparison_TotalSet.Rds", sep=""))
 load(paste(PubOutdir, "MissingnessComparison_ConsumpSet.Rds", sep=""))
 load(paste(PubOutdir, "MissingnessComparison_ProbSet.Rds", sep=""))
 
 naSum.T$Row <- seq(1:dim(naSum.T)[1])
 naSum.C$Row <- seq(1:dim(naSum.C)[1])
 naSum.P$Row <- seq(1:dim(naSum.P)[1])

###################################### Incremental Rerun ##################################
### Run the Lasso model for
### every subset of variables in a step-wise fashion
sets <- c("score", "AudC", "AudP")
setL <- c("AUDIT_T", "AUDIT_C", "AUDIT_P")
setI <- which(colnames(AUDIT) %in% sets)
suffixes <- c("T", "C", "P")
N <- c(which(naSum.T$CompCase.Add == 0)[1] - 1, which(naSum.C$CompCase.Add == 0)[1] - 1, which(naSum.P$CompCase.Add == 0)[1] - 1)
up.limit <- c(40, 12, 28)
low.limit <- c(0, 0, 0)

for (k in 1:3){
  eval(parse(text = paste("naSum.", suffixes[k], "$PredCor <- NA", sep="")))
  for (i in 2:N[k]){
    cat(setL[k],", Set ", i, "\n")
    eval(parse(text = paste("ret <- naSum.", suffixes[k], "$Row.names[1:", i, "]", sep="")))
    bb <- bd.top[, colnames(bd.top) %in% ret] ### subset the data for that cluster
    bb.cols <- rownames(bb[complete.cases(bb), ]) ### use complete cases
    factors <- sapply(bb, is.factor) ### identify which are factors (all will be ordered, by default)
    factor.names <- names(factors[factors])
    classes <- sapply(sapply(bb, class), function(x) x[1])
    if (length(grep("character", classes)) > 0){
      cat("Danger, Will Robinson! For i = ", i, "\n")
    } else if (length(grep("ordered", classes)) < 1) {
      cat("No factors in i = ", i, "\n")
      eval(parse(text = paste("bdSet", suffixes[k], "_", i, " <- ", "bb", sep="")))
      cc <- data.frame(group=seq(1:(dim(bb)[2])), item=sapply(strsplit(colnames(bb), "\\."), "[", 2), orig=colnames(bb))
      eval(parse(text = paste("rownames(bdSet", suffixes[k], "_", i, ") <- rownames(bb)", sep="")))
      eval(parse(text = paste("dd <- data.frame(catCol=colnames(bdSet",suffixes[k], "_", i, "), item=sapply(strsplit(colnames(bdSet",suffixes[k], "_", i, "), \"\\\\.\"), \"[\", 2))", sep="")))
      eval(parse(text = paste("group",suffixes[k], "_", i," <- merge(dd, cc, by=\"item\", all.x=TRUE)", sep="")))
      eval(parse(text = paste("PredSet",suffixes[k], "_",  i, "<- bdSet",suffixes[k], "_", i, sep="")))
      eval(parse(text = paste("bdSet",suffixes[k], "_", i, "<- bdSet",suffixes[k], "_", i, "[rownames(bdSet",suffixes[k], "_", i, ") %in% bb.cols, ]", sep="")))
    } else {
      bb[,factors] <- sapply(bb[,factors], factor, ordered=FALSE) ### remove the ordered from the factors
      ### Make the dummy-variable amplified matrix
      eval(parse(text = paste("bdSet",suffixes[k], "_", i, " <- ", "dummy_cols(bb, remove_first_dummy = TRUE, remove_most_frequent_dummy = FALSE, ignore_na = TRUE)", sep="")))
      ### drop the original factor variables
      eval(parse(text = paste("bdSet",suffixes[k], "_", i, " <- bdSet",suffixes[k], "_", i, "[, -which(colnames(bdSet",suffixes[k], "_", i,") %in% factor.names)]", sep="")))
      ### Assign the rownames (the identifiers)
      eval(parse(text = paste("rownames(bdSet",suffixes[k], "_", i, ") <- rownames(bb)", sep="")))
      cc <- data.frame(group=seq(1:(dim(bb)[2])), item=sapply(strsplit(colnames(bb), "\\."), "[", 2), orig=colnames(bb))
      eval(parse(text = paste("dd <- data.frame(catCol=colnames(bdSet",suffixes[k], "_", i, "), item=sapply(strsplit(colnames(bdSet",suffixes[k], "_", i, "), \"\\\\.\"), \"[\", 2))", sep="")))
      eval(parse(text = paste("group",suffixes[k], "_", i," <- merge(dd, cc, by=\"item\", all.x=TRUE)", sep="")))
      ### Print the dimensions of the groupi object, the first of which gives the new number of columns, after dummy expansion
      # eval(parse(text = paste("print(dim(group",suffixes[k], "_", i, "))", sep="")))
      eval(parse(text = paste("PredSet",suffixes[k], "_", i, "<- bdSet",suffixes[k], "_", i, sep="")))
      eval(parse(text = paste("bdSet",suffixes[k], "_", i, "<- bdSet",suffixes[k], "_", i, "[rownames(bdSet",suffixes[k], "_", i, ") %in% bb.cols, ]", sep="")))
    }
    
    eval(parse(text = paste("bdSet <- bdSet",suffixes[k], "_", i, sep="")))
    eval(parse(text = paste("rownames(bdSet) <- rownames(bdSet",suffixes[k], "_", i, ")", sep="")))
    eval(parse(text = paste("group <- group",suffixes[k], "_", i, sep="")))
    # cat("Set ", i, "\n")
    cat("Dimensions of complete case matrix: ", dim(bdSet), "\n")
    aud <- AUDIT[AUDIT$f.eid %in% rownames(bdSet), c("f.eid",sets[k])] 
    grpSet <- merge(aud, bdSet, by.x="f.eid", by.y="row.names") # merge with AUDIT
    cat("Dimensions after merging with AUDIT: ", dim(grpSet), "\n")
    groupdf <- group # keep the grouping dataframe
    group <- group[match(colnames(grpSet[3:dim(grpSet)[2]]), as.character(group$catCol)),] # order the group set
    
    X <- as.matrix(grpSet[, 3:(dim(grpSet)[2])]) # put data into a matrix
    y <- grpSet[, colnames(grpSet) %in% sets[k]]
    
    set.seed(93019)
    test <- sample(seq(1:length(y)), floor(0.25 * length(y)))
    X.test <- X[test, ]
    y.test <- y[test]
    X <- X[-test, ]
    y <- y[-test]
    group <- group$group # set the group vector
    cv.testG <- cv.grpreg(X=X, y=y, group=group, penalty="grLasso", family="gaussian", nfolds=5, seed=51320) # cross-validated model
    pred.y <- predict(cv.testG, X.test, lambda=cv.testG$lambda.min)
    pred.y[pred.y > up.limit[k]] <- up.limit[k]
    pred.y[pred.y < low.limit[k]] <- low.limit[k]
    cat("Test set correlation btw predicted and actual: ", cor(pred.y, y.test), "\n")
    eval(parse(text = paste("naSum.", suffixes[k], "$PredCor[", i, "] <- cor(pred.y, y.test)", sep="")))
    retained <- unique(as.character(
      groupdf$orig[as.character(groupdf$catCol) %in% names(coef(cv.testG, s=lambda.min)[coef(cv.testG, s=lambda.min) != 0]) ])) # variables retained
    cat(length(retained), " variables retained", "\n")
    X <- as.matrix(grpSet[, 3:(dim(grpSet)[2])]) # put data into a matrix
    y <- grpSet[, colnames(grpSet) %in% sets[k]]
    pred.y.full <- predict(cv.testG, X, lambda=cv.testG$lambda.min)
    pred.y.full[pred.y.full > up.limit[k]] <- up.limit[k]
    pred.y.full[pred.y.full < low.limit[k]] <- low.limit[k]
    cat("Full set correlation btw predicted and actual: ",cor(pred.y.full, y), "\n")
    eval(parse(text = paste("pred.y.",suffixes[k], "_", i, "<- pred.y", sep="")))
    eval(parse(text = paste("pred.y.full.",suffixes[k], "_", i,"<- pred.y.full", sep="")))
    eval(parse(text = paste("retained",suffixes[k], "_", i, "<- retained", sep="")))
    eval(parse(text = paste("grpSet",suffixes[k], "_", i, "<- grpSet", sep="")))
    eval(parse(text = paste("Coefs", suffixes[k], "_",i, "<- coef(cv.testG, s=lambda.min)", sep="")))
  }
} 

### Save the tables with the correlations added as Rds and csv files (for tables)
# naSum.T$PredCor[(N[1]+1):dim(naSum.T)[1]] <- NA
# naSum.C$PredCor[(N[2]+1):dim(naSum.T)[1]] <- NA
save(naSum.T, file=paste(PubOutdir, "MissingnessComparison_TotalSet_wCorr.Rds", sep=""))
save(naSum.C, file=paste(PubOutdir, "MissingnessComparison_ConsumpSet_wCorr.Rds", sep=""))
save(naSum.P, file=paste(PubOutdir, "MissingnessComparison_ProbSet_wCorr.Rds", sep=""))

naSum.T.pub <- naSum.T[,c("Row", "CompCase.Add", "CompCase.Add.Obs", "CompCase.Add.UnObs","PropObs", "PropUnobs", "PropRatio", "PredCor")]
naSum.T.pub[, c(5:8)] <- apply(naSum.T.pub[,c(5:8)], 2, round, 4)
write.csv(naSum.T.pub, file=paste(PubOutdir, "MissingnessComparison_TotalSet_wCorrPub.csv", sep=""))
naSum.C.pub <- naSum.C[,c("Row", "CompCase.Add", "CompCase.Add.Obs", "CompCase.Add.UnObs","PropObs", "PropUnobs", "PropRatio", "PredCor")]
naSum.C.pub[, c(5:8)] <- apply(naSum.C.pub[,c(5:8)], 2, round, 4)
write.csv(naSum.C.pub, file=paste(PubOutdir, "MissingnessComparison_ConsumpSet_wCorrPub.csv", sep=""))
naSum.P.pub <- naSum.P[,c("Row", "CompCase.Add", "CompCase.Add.Obs", "CompCase.Add.UnObs","PropObs", "PropUnobs", "PropRatio", "PredCor")]
naSum.P.pub[, c(5:8)] <- apply(naSum.P.pub[,c(5:8)], 2, round, 4)
write.csv(naSum.P.pub, file=paste(PubOutdir, "MissingnessComparison_ProbSet_wCorrPub.csv", sep=""))

### Reload if needed
load(paste(PubOutdir, "MissingnessComparison_TotalSet_wCorr.Rds", sep=""))
load(paste(PubOutdir, "MissingnessComparison_ConsumpSet_wCorr.Rds", sep=""))
load(paste(PubOutdir, "MissingnessComparison_ProbSet_wCorr.Rds", sep=""))

###################### Plotting #####################################
### Making plots of the stepwise performance is helpful

### Total
naSum.T$PropObs <- naSum.T$CompCase.Add.Obs / 157162
naSum.T$PropUnobs <- naSum.T$CompCase.Add.UnObs / (502536 - 157162)
naSum.T$PropRatio <- naSum.T$PropUnobs/naSum.T$PropObs

png(paste(plotdir, "Pub_AUDIT_T_BalanceVsPrediction.png", sep=""))
plot(naSum.T$Row, naSum.T$PropRatio, xlim=c(0,40), ylim=c(0.2,1.5),
     ylab="", xlab="Number of Phenotypes", main="AUDIT-Total")
lines(naSum.T$Row, naSum.T$PropRatio, xlim=c(0,30))
points(naSum.T$Row, naSum.T$PredCor, xlim=c(0,30), col="red")
lines(naSum.T$Row, naSum.T$PredCor, xlim=c(0,30), col="red")
points(naSum.T$Row[which(naSum.T$PredCor == max(naSum.T$PredCor,na.rm=TRUE))],
       max(naSum.T$PredCor,na.rm=TRUE), col="blue", pch=16)
legend("topleft", 
       c("Ratio of Prop Unmeas/Meas", "Correlation between Obs/Pred"),
       fill=c("black", "red"))
dev.off()
T.point <- which(naSum.T$PredCor == max(naSum.T$PredCor,na.rm=TRUE))
T.Coefs <- naSum.T$Row.names[1:T.point]
save(T.Coefs, file=paste(preddir, "Round3_Pub/GrLassoAUDIT_FinalVars_T.Rds", sep=""))
### ^ This above will save the names of the variables that go INTO the final model, but not the retained variables

### Consumption
naSum.C$PropObs <- naSum.C$CompCase.Add.Obs / 157162
naSum.C$PropUnobs <- naSum.C$CompCase.Add.UnObs / (502536 - 157162)
naSum.C$PropRatio <- naSum.C$PropUnobs/naSum.C$PropObs

png(paste(plotdir, "Pub_AUDIT_C_BalanceVsPrediction.png", sep=""))
plot(naSum.C$Row, naSum.C$PropRatio, xlim=c(0,45), ylim=c(0.2,1.5),
     ylab="", xlab="Number of Phenotypes", main="AUDIT-C")
lines(naSum.C$Row, naSum.C$PropRatio, xlim=c(0,30))
points(naSum.C$Row, naSum.C$PredCor, xlim=c(0,30), col="red")
lines(naSum.C$Row, naSum.C$PredCor, xlim=c(0,30), col="red")
points(naSum.C$Row[which(naSum.C$PredCor == max(naSum.C$PredCor,na.rm=TRUE))],
       max(naSum.C$PredCor,na.rm=TRUE), col="blue", pch=16)
legend("topleft", 
       c("Ratio of Prop Unmeas/Meas", "Correlation between Obs/Pred"),
       fill=c("black", "red"))
dev.off()
C.point <- which(naSum.C$PredCor == max(naSum.C$PredCor,na.rm=TRUE))
C.Coefs <- naSum.C$Row.names[1:C.point]
save(C.Coefs, file=paste(preddir, "Round3_Pub/GrLassoAUDIT_FinalVars_C.Rds", sep=""))

### Problems
naSum.P$PropObs <- naSum.P$CompCase.Add.Obs / 157162
naSum.P$PropUnobs <- naSum.P$CompCase.Add.UnObs / (502536 - 157162)
naSum.P$PropRatio <- naSum.P$PropUnobs/naSum.P$PropObs

png(paste(plotdir, "Pub_AUDIT_P_BalanceVsPrediction.png", sep=""))
plot(naSum.P$Row, naSum.P$PropRatio, xlim=c(0,45), ylim=c(0.2,1.5),
     ylab="", xlab="Number of Phenotypes", main="AUDIT-P")
lines(naSum.P$Row, naSum.P$PropRatio, xlim=c(0,30))
points(naSum.P$Row, naSum.P$PredCor, xlim=c(0,30), col="red")
lines(naSum.P$Row, naSum.P$PredCor, xlim=c(0,30), col="red")
points(naSum.P$Row[which(naSum.P$PredCor == max(naSum.P$PredCor,na.rm=TRUE))],
       max(naSum.P$PredCor,na.rm=TRUE), col="blue", pch=16)
legend("topleft", 
       c("Ratio of Prop Unmeas/Meas", "Correlation between Obs/Pred"),
       fill=c("black", "red"))
dev.off()
P.point <- which(naSum.P$PredCor == max(naSum.P$PredCor,na.rm=TRUE))
P.Coefs <- naSum.P$Row.names[1:P.point]
save(P.Coefs, file=paste(preddir, "Round3_Pub/GrLassoAUDIT_FinalVars_P.Rds", sep=""))

########################################################################
################ Make Prediction Sets ##################################
### Calculate the index for the model with best performance
index<-c(which(naSum.T$PredCor == max(naSum.T$PredCor,na.rm=TRUE)),
         which(naSum.C$PredCor == max(naSum.C$PredCor,na.rm=TRUE)),
         which(naSum.P$PredCor == max(naSum.P$PredCor,na.rm=TRUE)))

for (i in 1:3){
eval(parse(text = paste("aa", "<- Coefs",suffixes[i], "_", index[i], sep="")))
bb <- names(aa[aa != 0])
eval(parse(text = paste(suffixes[i], ".NonZero <- unique(group",suffixes[i], "_", index[i],"$orig[group",suffixes[i], "_", index[i],"$catCol %in% bb])", sep="")))
eval(parse(text = paste("save(", suffixes[i], ".NonZero, file=paste(preddir, \"/Round3_Pub/GrLassoAUDIT_FinalNonZeroVars_", suffixes[i], ".Rds\", sep=\"\"))", sep="")))
}

for (i in 1:3){
### Create the predictions for the full UKB set
eval(parse(text = paste("PredSet", "<- PredSet",suffixes[i], "_", index[i], sep="")))
PredSet <- data.matrix(PredSet)
PredSet <- cbind(rep(1,dim(PredSet)[1]), PredSet)
PredSet[is.na(PredSet)] <- 0
eval(parse(text = paste("full.pred <- PredSet %*% Coefs", suffixes[i], "_", index[i], sep="")))
full.pred[full.pred > up.limit[i]] <- up.limit[i]
full.pred[full.pred < low.limit[i]] <- low.limit[i]
eval(parse(text = paste("full.pred", suffixes[i], "_", index[i], "<- full.pred", sep="")))

### Create the predictions for the set with AUDIT measured
eval(parse(text = paste("PredSet", "<- PredSet",  suffixes[i], "_", index[i], sep="")))
PredSet <- PredSet[rownames(PredSet) %in% AUDIT$f.eid,]
PredSet <- data.matrix(PredSet)
PredSet <- cbind(rep(1,dim(PredSet)[1]), PredSet)
PredSet[is.na(PredSet)] <- 0
eval(parse(text = paste("full.pred <- PredSet %*% Coefs", suffixes[i], "_", index[i], sep="")))
full.pred[full.pred > up.limit[i]] <- up.limit[i]
full.pred[full.pred < low.limit[i]] <- low.limit[i]
eval(parse(text = paste("full.pred.measAUDIT.", suffixes[i], "_", index[i]," <- full.pred", sep="")))

### Merge this all into a dataframe to save and output for GWAS
### Make a set with the predicted values
eval(parse(text = paste("GrLassoAUDIT <- data.frame(cbind(grpSet", suffixes[i], "_", index[i], "$f.eid, pred.y.full.", suffixes[i], "_", index[i], "))", sep="")))
colnames(GrLassoAUDIT) <- c("f.eid", "AUDIThat_subset")
GrLassoAUDIT$f.eid <- as.character(GrLassoAUDIT$f.eid)
aud <- AUDIT[, c("f.eid",sets[i])]
### Merge the AUDIT observations
GrLassoAUDIT <- merge(aud, GrLassoAUDIT, by="f.eid", all.x=TRUE)
### Merge the predictions where AUDIT is measured
eval(parse(text = paste("GrLassoAUDIT <- merge(GrLassoAUDIT, full.pred.measAUDIT.", suffixes[i], "_", index[i], ",by.x=\"f.eid\", by.y=\"row.names\", all.y=TRUE)", sep="")))
colnames(GrLassoAUDIT) <- c("f.eid", "AUDIT", "AUDIThat_subset", "AUDIThat_full")
eval(parse(text = paste("GrLassoAUDIT <- merge(GrLassoAUDIT, full.pred", suffixes[i], "_", index[i], ",by.x=\"f.eid\", by.y=\"row.names\", all.y=TRUE)", sep="")))
colnames(GrLassoAUDIT) <- c("f.eid", "AUDIT", "AUDIThat_subset", "AUDIThat_full", "AUDIThat_total")
eval(parse(text = paste("colnames(GrLassoAUDIT) <- c(\"f.eid\", \"AUDIT_", suffixes[i],"\", \"AUDIThat_subset\", \"AUDIThat_full\", \"AUDIThat_total\")", sep="")))
eval(parse(text = paste("GrLassoAUDIT_", suffixes[i], " <- GrLassoAUDIT", sep="")))
rundate <- Sys.Date()
eval(parse(text = paste("save(GrLassoAUDIT_", suffixes[i], ", file=paste(preddir, \"/Round3_Pub/GrLassoAUDITpred_Rerun_", suffixes[i], "_", rundate, ".Rds\", sep=\"\"))", sep="")))
}

######## FINAL variables entering each model ##############
data.dict$Field[data.dict$UDI %in% T.Coefs]
### Total Score Variables:
# [1] Episodes containing "Dates of operations" data           
# [2] Apolipoprotein A                                         
# [3] Cystatin C                                               
# [4] Gamma glutamyltransferase                                
# [5] Urate                                                    
# [6] Arm fat-free mass (right)                                
# [7] Leg predicted mass (left)                                
# [8] Waist circumference                                      
# [9] Forced expiratory volume in 1-second (FEV1), Best measure
# [10] Mother's age at death                                    
# [11] Alcohol usually taken with meals                         
# [12] Average weekly spirits intake                            
# [13] Average weekly beer plus cider intake                    
# [14] Average weekly fortified wine intake                     
# [15] Average weekly red wine intake                           
# [16] Alcohol intake frequency.                                
# [17] Average weekly champagne plus white wine intake          
# [18] Cereal intake                                            
# [19] Age first had sexual intercourse                         
# [20] Past tobacco smoking                                     
# [21] Smoking status                                           
# [22] Ever smoked                                              
# [23] Neuroticism score                                        
# [24] Irritability                                             
# [25] Length of working week for main job                      
# [26] Mean corpuscular volume                                  
# [27] Mean corpuscular haemoglobin                             
# [28] Potassium in urine                                       
# [29] Year of birth                                            
# [30] Sex                                  

data.dict$Field[data.dict$UDI %in% C.Coefs]
### Consumption Variables:
# [1] Apolipoprotein A            
# [2] Cystatin C                  
# [3] Gamma glutamyltransferase   
# [4] Urate                       
# [5] Arm fat-free mass (right)   
# [6] Leg predicted mass (left)   
# [7] Hand grip strength (left)   
# [8] Alcohol intake frequency.   
# [9] Cereal intake               
# [10] Cereal type                 
# [11] Smoking status              
# [12] Ever smoked                 
# [13] Mean corpuscular volume     
# [14] Mean corpuscular haemoglobin
# [15] Year of birth               
# [16] Sex                         
# [17] Genetic sex                 
# [18] Heterozygosity 

data.dict$Field[data.dict$UDI %in% P.Coefs]
# [1] Aspartate aminotransferase                     
# [2] Gamma glutamyltransferase                      
# [3] Urate                                          
# [4] Leg predicted mass (left)                      
# [5] Waist circumference                            
# [6] Mother's age at death                          
#  [7] Alcohol usually taken with meals               
#  [8] Average weekly spirits intake                  
#  [9] Average weekly beer plus cider intake          
# [10] Average weekly fortified wine intake           
# [11] Average weekly red wine intake                 
# [12] Alcohol intake frequency.                      
# [13] Average weekly champagne plus white wine intake
# [14] Cereal intake                                  
# [15] Smoking status                                 
# [16] Ever smoked                                    
# [17] Neuroticism score                              
# [18] Guilty feelings                                
# [19] Mean corpuscular volume                        
# [20] Year of birth      

######## FINAL NON-ZERO variables in each model ##############
data.dict$Field[data.dict$UDI %in% T.NonZero]
# [1] Episodes containing "Dates of operations" data           
# [2] Apolipoprotein A                                         
# [3] Cystatin C                                               
# [4] Gamma glutamyltransferase                                
# [5] Urate                                                    
# [6] Waist circumference                                      
# [7] Forced expiratory volume in 1-second (FEV1), Best measure
# [8] Alcohol usually taken with meals                         
# [9] Average weekly spirits intake                            
# [10] Average weekly beer plus cider intake                    
# [11] Average weekly fortified wine intake                     
# [12] Average weekly red wine intake                           
# [13] Alcohol intake frequency.                                
# [14] Average weekly champagne plus white wine intake          
# [15] Cereal intake                                            
# [16] Age first had sexual intercourse                         
# [17] Past tobacco smoking                                     
# [18] Smoking status                                           
# [19] Ever smoked                                              
# [20] Neuroticism score                                        
# [21] Irritability                                             
# [22] Length of working week for main job                      
# [23] Mean corpuscular volume                                  
# [24] Mean corpuscular haemoglobin                             
# [25] Potassium in urine                                       
# [26] Year of birth                                            
# [27] Sex                   

data.dict$Field[data.dict$UDI %in% C.NonZero]
# [1] Apolipoprotein A             
# [2] Cystatin C                  
# [3] Gamma glutamyltransferase    
# [4] Urate                       
# [5] Arm fat-free mass (right)    
# [6] Leg predicted mass (left)   
# [7] Hand grip strength (left)    
# [8] Alcohol intake frequency.   
# [9] Cereal intake                
# [10] Cereal type                 
# [11] Smoking status               
# [12] Ever smoked                 
# [13] Mean corpuscular volume      
# [14] Mean corpuscular haemoglobin
# [15] Year of birth                
# [16] Sex                         
# [17] Genetic sex                  
# [18] Heterozygosity      

data.dict$Field[data.dict$UDI %in% P.NonZero]
# [1] Gamma glutamyltransferase                      
# [2] Alcohol usually taken with meals               
# [3] Average weekly spirits intake                  
# [4] Average weekly beer plus cider intake          
# [5] Average weekly fortified wine intake           
# [6] Average weekly red wine intake                 
# [7] Average weekly champagne plus white wine intake
# [8] Cereal intake                                  
# [9] Smoking status                                 
# [10] Ever smoked                                    
# [11] Neuroticism score                              
# [12] Guilty feelings                                
# [13] Mean corpuscular volume                        
# [14] Year of birth   

####################################################################
### Output from the stepwise model runs

# AUDIT_T , Set  2 
# Dimensions of complete case matrix:  502536 2 
# Dimensions after merging with AUDIT:  157162 4 
# Test set correlation btw predicted and actual:  0.2683643 
# 2  variables retained 
# Full set correlation btw predicted and actual:  0.2751769 
# AUDIT_T , Set  3 
# Dimensions of complete case matrix:  501645 5 
# Dimensions after merging with AUDIT:  157089 7 
# Test set correlation btw predicted and actual:  0.3418764 
# 3  variables retained 
# Full set correlation btw predicted and actual:  0.3414631 
# AUDIT_T , Set  4 
# Dimensions of complete case matrix:  501639 11 
# Dimensions after merging with AUDIT:  157088 13 
# Test set correlation btw predicted and actual:  0.60442 
# 4  variables retained 
# Full set correlation btw predicted and actual:  0.6031113 
# AUDIT_T , Set  5 
# Dimensions of complete case matrix:  501632 14 
# Dimensions after merging with AUDIT:  157088 16 
# Test set correlation btw predicted and actual:  0.6057935 
# 5  variables retained 
# Full set correlation btw predicted and actual:  0.6051382 
# AUDIT_T , Set  6 
# Dimensions of complete case matrix:  501515 15 
# Dimensions after merging with AUDIT:  157067 17 
# Test set correlation btw predicted and actual:  0.6110423 
# 6  variables retained 
# Full set correlation btw predicted and actual:  0.6089528 
# AUDIT_T , Set  7 
# Dimensions of complete case matrix:  499725 16 
# Dimensions after merging with AUDIT:  156837 18 
# Test set correlation btw predicted and actual:  0.6160497 
# 7  variables retained 
# Full set correlation btw predicted and actual:  0.6117077 
# AUDIT_T , Set  8 
# Dimensions of complete case matrix:  497760 17 
# Dimensions after merging with AUDIT:  156548 19 
# Test set correlation btw predicted and actual:  0.608098 
# 8  variables retained 
# Full set correlation btw predicted and actual:  0.612161 
# AUDIT_T , Set  9 
# Dimensions of complete case matrix:  489704 18 
# Dimensions after merging with AUDIT:  154505 20 
# 
# 
# Test set correlation btw predicted and actual:  0.6089419 
# 9  variables retained 
# Full set correlation btw predicted and actual:  0.6123988 
# AUDIT_T , Set  10 
# Dimensions of complete case matrix:  489654 19 
# Dimensions after merging with AUDIT:  154491 21 
# Test set correlation btw predicted and actual:  0.6128231 
# 10  variables retained 
# Full set correlation btw predicted and actual:  0.6123957 
# AUDIT_T , Set  11 
# Dimensions of complete case matrix:  472443 20 
# Dimensions after merging with AUDIT:  149925 22 
# Test set correlation btw predicted and actual:  0.6129868 
# 10  variables retained 
# Full set correlation btw predicted and actual:  0.6130932 
# AUDIT_T , Set  12 
# Dimensions of complete case matrix:  453342 21 
# Dimensions after merging with AUDIT:  144398 23 
# Test set correlation btw predicted and actual:  0.6173864 
# 11  variables retained 
# Full set correlation btw predicted and actual:  0.6173119 
# AUDIT_T , Set  13 
# Dimensions of complete case matrix:  453339 22 
# Dimensions after merging with AUDIT:  144397 24 
# Test set correlation btw predicted and actual:  0.6186652 
# 12  variables retained 
# Full set correlation btw predicted and actual:  0.617555 
# AUDIT_T , Set  14 
# Dimensions of complete case matrix:  433878 23 
# Dimensions after merging with AUDIT:  138101 25 
# Test set correlation btw predicted and actual:  0.6188697 
# 14  variables retained 
# Full set correlation btw predicted and actual:  0.6196578 
# AUDIT_T , Set  15 
# Dimensions of complete case matrix:  433378 24 
# Dimensions after merging with AUDIT:  137942 26 
# Test set correlation btw predicted and actual:  0.6196756 
# 15  variables retained 
# Full set correlation btw predicted and actual:  0.6222682 
# AUDIT_T , Set  16 
# Dimensions of complete case matrix:  432764 25 
# Dimensions after merging with AUDIT:  137767 27 
# Test set correlation btw predicted and actual:  0.6280607 
# 16  variables retained 
# Full set correlation btw predicted and actual:  0.6248224 
# AUDIT_T , Set  17 
# Dimensions of complete case matrix:  399391 29 
# Dimensions after merging with AUDIT:  131141 31 
# Test set correlation btw predicted and actual:  0.6229661 
# 16  variables retained 
# Full set correlation btw predicted and actual:  0.6253834 
# AUDIT_T , Set  18 
# Dimensions of complete case matrix:  364427 30 
# Dimensions after merging with AUDIT:  125509 32 
# Test set correlation btw predicted and actual:  0.6035287 
# 12  variables retained 
# Full set correlation btw predicted and actual:  0.6053793 
# AUDIT_T , Set  19 
# Dimensions of complete case matrix:  331556 31 
# Dimensions after merging with AUDIT:  113971 33 
# Test set correlation btw predicted and actual:  0.627982 
# 19  variables retained 
# Full set correlation btw predicted and actual:  0.6270715 
# AUDIT_T , Set  20 
# Dimensions of complete case matrix:  271230 32 
# Dimensions after merging with AUDIT:  95947 34 
# Test set correlation btw predicted and actual:  0.6155192 
# 15  variables retained 
# Full set correlation btw predicted and actual:  0.6144059 
# AUDIT_T , Set  21 
# Dimensions of complete case matrix:  215930 36 
# Dimensions after merging with AUDIT:  79510 38 
# Test set correlation btw predicted and actual:  0.5544056 
# 17  variables retained 
# Full set correlation btw predicted and actual:  0.5518448 
# AUDIT_T , Set  22 
# Dimensions of complete case matrix:  163047 37 
# Dimensions after merging with AUDIT:  62004 39 
# Test set correlation btw predicted and actual:  0.5469444 
# 19  variables retained 
# Full set correlation btw predicted and actual:  0.5572768 
# AUDIT_T , Set  23 
# Dimensions of complete case matrix:  149450 38 
# Dimensions after merging with AUDIT:  56980 40 
# Test set correlation btw predicted and actual:  0.4953301 
# 10  variables retained 
# Full set correlation btw predicted and actual:  0.4959403 
# AUDIT_T , Set  24 
# Dimensions of complete case matrix:  149450 39 
# Dimensions after merging with AUDIT:  56980 41 
# Test set correlation btw predicted and actual:  0.5284406 
# 8  variables retained 
# Full set correlation btw predicted and actual:  0.5317964 
# AUDIT_T , Set  25 
# Dimensions of complete case matrix:  149450 40 
# Dimensions after merging with AUDIT:  56980 42 
# Test set correlation btw predicted and actual:  0.6178606 
# 12  variables retained 
# Full set correlation btw predicted and actual:  0.6224632 
# AUDIT_T , Set  26 
# Dimensions of complete case matrix:  149450 41 
# Dimensions after merging with AUDIT:  56980 43 
# Test set correlation btw predicted and actual:  0.6375945 
# 14  variables retained 
# Full set correlation btw predicted and actual:  0.6411493 
# AUDIT_T , Set  27 
# Dimensions of complete case matrix:  149450 42 
# Dimensions after merging with AUDIT:  56980 44 
# Test set correlation btw predicted and actual:  0.6375945 
# 14  variables retained 
# Full set correlation btw predicted and actual:  0.6411493 
# AUDIT_T , Set  28 
# Dimensions of complete case matrix:  91677 43 
# Dimensions after merging with AUDIT:  33330 45 
# Test set correlation btw predicted and actual:  0.6264875 
# 14  variables retained 
# Full set correlation btw predicted and actual:  0.6368451 
# AUDIT_T , Set  29 
# Dimensions of complete case matrix:  55671 44 
# Dimensions after merging with AUDIT:  19386 46 
# Test set correlation btw predicted and actual:  0.6365635 
# 26  variables retained 
# Full set correlation btw predicted and actual:  0.6600155 
                                          # AUDIT_T , Set  30 
                                          # Dimensions of complete case matrix:  26307 45 
                                          # Dimensions after merging with AUDIT:  9862 47 
                                          # Test set correlation btw predicted and actual:  0.6390275 
                                          # 27  variables retained 
                                          # Full set correlation btw predicted and actual:  0.6486998 
# AUDIT_T , Set  31 
# Dimensions of complete case matrix:  26307 50 
# Dimensions after merging with AUDIT:  9862 52 
# Test set correlation btw predicted and actual:  0.6370106 
# 26  variables retained 
# Full set correlation btw predicted and actual:  0.6469072 
# AUDIT_T , Set  32 
# Dimensions of complete case matrix:  26307 55 
# Dimensions after merging with AUDIT:  9862 57 
# Test set correlation btw predicted and actual:  0.6370106 
# 26  variables retained 
# Full set correlation btw predicted and actual:  0.6469072 
# AUDIT_T , Set  33 
# Dimensions of complete case matrix:  26307 56 
# Dimensions after merging with AUDIT:  9862 58 
# Test set correlation btw predicted and actual:  0.6370106 
# 26  variables retained 
# Full set correlation btw predicted and actual:  0.6469072 
# AUDIT_T , Set  34 
# Dimensions of complete case matrix:  12636 57 
# Dimensions after merging with AUDIT:  4904 59 
# Test set correlation btw predicted and actual:  0.6241357 
# 25  variables retained 
# Full set correlation btw predicted and actual:  0.6464378 
# AUDIT_T , Set  35 
# Dimensions of complete case matrix:  5844 60 
# Dimensions after merging with AUDIT:  2062 62 
# Test set correlation btw predicted and actual:  0.6271389 
# 20  variables retained 
# Full set correlation btw predicted and actual:  0.662173 
# AUDIT_T , Set  36 
# Dimensions of complete case matrix:  5844 64 
# Dimensions after merging with AUDIT:  2062 66 
# Test set correlation btw predicted and actual:  0.6271389 
# 20  variables retained 
# Full set correlation btw predicted and actual:  0.662173 
# AUDIT_T , Set  37 
# Dimensions of complete case matrix:  5844 68 
# Dimensions after merging with AUDIT:  2062 70 
# Test set correlation btw predicted and actual:  0.6271389 
# 20  variables retained 
# Full set correlation btw predicted and actual:  0.662173 
# AUDIT_T , Set  38 
# Dimensions of complete case matrix:  5844 73 
# Dimensions after merging with AUDIT:  2062 75 
# Test set correlation btw predicted and actual:  0.6271389 
# 20  variables retained 
# Full set correlation btw predicted and actual:  0.662173 
# AUDIT_T , Set  39 
# Dimensions of complete case matrix:  5844 74 
# Dimensions after merging with AUDIT:  2062 76 
# Test set correlation btw predicted and actual:  0.6271389 
# 20  variables retained 
# Full set correlation btw predicted and actual:  0.662173 
# AUDIT_T , Set  40 
# Dimensions of complete case matrix:  2612 80 
# Dimensions after merging with AUDIT:  852 82 
# Test set correlation btw predicted and actual:  0.5727765 
# 26  variables retained 
# Full set correlation btw predicted and actual:  0.6260669 
#############################################################
# AUDIT_C , Set  2 
# Dimensions of complete case matrix:  502536 2 
# Dimensions after merging with AUDIT:  157162 4 
# Test set correlation btw predicted and actual:  0.2883792 
# 2  variables retained 
# Full set correlation btw predicted and actual:  0.294723 
# AUDIT_C , Set  3 
# Dimensions of complete case matrix:  501645 5 
# Dimensions after merging with AUDIT:  157089 7 
# Test set correlation btw predicted and actual:  0.3470688 
# 3  variables retained 
# Full set correlation btw predicted and actual:  0.3480297 
# AUDIT_C , Set  4 
# Dimensions of complete case matrix:  501639 11 
# Dimensions after merging with AUDIT:  157088 13 
# Test set correlation btw predicted and actual:  0.6880834 
# 4  variables retained 
# Full set correlation btw predicted and actual:  0.6862633 
# AUDIT_C , Set  5 
# Dimensions of complete case matrix:  501522 12 
# Dimensions after merging with AUDIT:  157067 14 
# Test set correlation btw predicted and actual:  0.6870141 
# 5  variables retained 
# Full set correlation btw predicted and actual:  0.6886524 
# AUDIT_C , Set  6 
# Dimensions of complete case matrix:  499529 13 
# Dimensions after merging with AUDIT:  156778 15 
# Test set correlation btw predicted and actual:  0.6881191 
# 6  variables retained 
# Full set correlation btw predicted and actual:  0.6889924 
# AUDIT_C , Set  7 
# Dimensions of complete case matrix:  496519 14 
# Dimensions after merging with AUDIT:  156224 16 
# Test set correlation btw predicted and actual:  0.6921933 
# 7  variables retained 
# Full set correlation btw predicted and actual:  0.689085 
# AUDIT_C , Set  8 
# Dimensions of complete case matrix:  488461 15 
# Dimensions after merging with AUDIT:  154161 17 
# Test set correlation btw predicted and actual:  0.68984 
# 8  variables retained 
# Full set correlation btw predicted and actual:  0.6902993 
# AUDIT_C , Set  9 
# Dimensions of complete case matrix:  488411 16 
# Dimensions after merging with AUDIT:  154147 18 
# Test set correlation btw predicted and actual:  0.6901305 
# 9  variables retained 
# Full set correlation btw predicted and actual:  0.6903247 
# AUDIT_C , Set  10 
# Dimensions of complete case matrix:  475775 17 
# Dimensions after merging with AUDIT:  150955 19 
# Test set correlation btw predicted and actual:  0.6927375 
# 10  variables retained 
# Full set correlation btw predicted and actual:  0.6903496 
# AUDIT_C , Set  11 
# Dimensions of complete case matrix:  475775 18 
# Dimensions after merging with AUDIT:  150955 20 
# Test set correlation btw predicted and actual:  0.6934636 
# 11  variables retained 
# Full set correlation btw predicted and actual:  0.6910805 
# AUDIT_C , Set  12 
# Dimensions of complete case matrix:  461337 19 
# Dimensions after merging with AUDIT:  146357 21 
# Test set correlation btw predicted and actual:  0.6896062 
# 12  variables retained 
# Full set correlation btw predicted and actual:  0.6939171 
# AUDIT_C , Set  13 
# Dimensions of complete case matrix:  461334 20 
# Dimensions after merging with AUDIT:  146356 22 
# Test set correlation btw predicted and actual:  0.6932464 
# 13  variables retained 
# Full set correlation btw predicted and actual:  0.6940466 
# AUDIT_C , Set  14 
# Dimensions of complete case matrix:  442202 21 
# Dimensions after merging with AUDIT:  140174 23 
# Test set correlation btw predicted and actual:  0.6920688 
# 14  variables retained 
# Full set correlation btw predicted and actual:  0.6953755 
# AUDIT_C , Set  15 
# Dimensions of complete case matrix:  441698 22 
# Dimensions after merging with AUDIT:  140017 24 
# Test set correlation btw predicted and actual:  0.6978933 
# 14  variables retained 
# Full set correlation btw predicted and actual:  0.6968771 
# AUDIT_C , Set  16 
# Dimensions of complete case matrix:  441068 23 
# Dimensions after merging with AUDIT:  139835 25 
# Test set correlation btw predicted and actual:  0.6990753 
# 16  variables retained 
# Full set correlation btw predicted and actual:  0.6997134 
# AUDIT_C , Set  17 
# Dimensions of complete case matrix:  401599 24 
# Dimensions after merging with AUDIT:  127090 26 
# Test set correlation btw predicted and actual:  0.7017683 
# 17  variables retained 
# Full set correlation btw predicted and actual:  0.7025164 
                                  # AUDIT_C , Set  18 
                                  # Dimensions of complete case matrix:  332697 30 
                                  # Dimensions after merging with AUDIT:  107054 32 
                                  # Test set correlation btw predicted and actual:  0.7052558 
                                  # 18  variables retained 
                                  # Full set correlation btw predicted and actual:  0.7001291 
# AUDIT_C , Set  19 
# Dimensions of complete case matrix:  267747 31 
# Dimensions after merging with AUDIT:  89546 33 
# Test set correlation btw predicted and actual:  0.6976437 
# 19  variables retained 
# Full set correlation btw predicted and actual:  0.6983257 
# AUDIT_C , Set  20 
# Dimensions of complete case matrix:  229101 32 
# Dimensions after merging with AUDIT:  77294 34 
# Test set correlation btw predicted and actual:  0.693831 
# 19  variables retained 
# Full set correlation btw predicted and actual:  0.6939683 
# AUDIT_C , Set  21 
# Dimensions of complete case matrix:  181950 36 
# Dimensions after merging with AUDIT:  64355 38 
# Test set correlation btw predicted and actual:  0.6073834 
# 21  variables retained 
# Full set correlation btw predicted and actual:  0.6117308 
# AUDIT_C , Set  22 
# Dimensions of complete case matrix:  147982 41 
# Dimensions after merging with AUDIT:  53641 43 
# Test set correlation btw predicted and actual:  0.6106492 
# 22  variables retained 
# Full set correlation btw predicted and actual:  0.6083669 
# AUDIT_C , Set  23 
# Dimensions of complete case matrix:  109901 42 
# Dimensions after merging with AUDIT:  38096 44 
# Test set correlation btw predicted and actual:  0.5997893 
# 22  variables retained 
# Full set correlation btw predicted and actual:  0.6062214 
# AUDIT_C , Set  24 
# Dimensions of complete case matrix:  84491 43 
# Dimensions after merging with AUDIT:  28471 45 
# Test set correlation btw predicted and actual:  0.6134048 
# 22  variables retained 
# Full set correlation btw predicted and actual:  0.6127586 
# AUDIT_C , Set  25 
# Dimensions of complete case matrix:  67376 44 
# Dimensions after merging with AUDIT:  22807 46 
# Test set correlation btw predicted and actual:  0.6026196 
# 14  variables retained 
# Full set correlation btw predicted and actual:  0.6001543 
# AUDIT_C , Set  26 
# Dimensions of complete case matrix:  50148 45 
# Dimensions after merging with AUDIT:  17732 47 
# Test set correlation btw predicted and actual:  0.6022882 
# 22  variables retained 
# Full set correlation btw predicted and actual:  0.6141533 
# AUDIT_C , Set  27 
# Dimensions of complete case matrix:  50148 46 
# Dimensions after merging with AUDIT:  17732 48 
# Test set correlation btw predicted and actual:  0.6020622 
# 23  variables retained 
# Full set correlation btw predicted and actual:  0.6142951 
# AUDIT_C , Set  28 
# Dimensions of complete case matrix:  46032 47 
# Dimensions after merging with AUDIT:  16357 49 
# Test set correlation btw predicted and actual:  0.5731675 
# 25  variables retained 
# Full set correlation btw predicted and actual:  0.5916907 
# AUDIT_C , Set  29 
# Dimensions of complete case matrix:  46032 48 
# Dimensions after merging with AUDIT:  16357 50 
# Test set correlation btw predicted and actual:  0.5975308 
# 25  variables retained 
# Full set correlation btw predicted and actual:  0.6141345 
# AUDIT_C , Set  30 
# Dimensions of complete case matrix:  46032 49 
# Dimensions after merging with AUDIT:  16357 51 
# Test set correlation btw predicted and actual:  0.6476782 
# 26  variables retained 
# Full set correlation btw predicted and actual:  0.6580099 
# AUDIT_C , Set  31 
# Dimensions of complete case matrix:  46032 50 
# Dimensions after merging with AUDIT:  16357 52 
# Test set correlation btw predicted and actual:  0.6609421 
# 26  variables retained 
# Full set correlation btw predicted and actual:  0.6716674 
# AUDIT_C , Set  32 
# Dimensions of complete case matrix:  46032 51 
# Dimensions after merging with AUDIT:  16357 53 
# Test set correlation btw predicted and actual:  0.6613014 
# 27  variables retained 
# Full set correlation btw predicted and actual:  0.6724205 
# AUDIT_C , Set  33 
# Dimensions of complete case matrix:  30054 52 
# Dimensions after merging with AUDIT:  8731 54 
# Test set correlation btw predicted and actual:  0.6677409 
# 27  variables retained 
# Full set correlation btw predicted and actual:  0.6695221 
# AUDIT_C , Set  34 
# Dimensions of complete case matrix:  25683 53 
# Dimensions after merging with AUDIT:  7353 55 
# Test set correlation btw predicted and actual:  0.6628387 
# 28  variables retained 
# Full set correlation btw predicted and actual:  0.6687983 
# AUDIT_C , Set  35 
# Dimensions of complete case matrix:  14602 54 
# Dimensions after merging with AUDIT:  4503 56 
# Test set correlation btw predicted and actual:  0.6663018 
# 27  variables retained 
# Full set correlation btw predicted and actual:  0.6609701 
# AUDIT_C , Set  36 
# Dimensions of complete case matrix:  9876 55 
# Dimensions after merging with AUDIT:  2979 57 
# Test set correlation btw predicted and actual:  0.6642751 
# 28  variables retained 
# Full set correlation btw predicted and actual:  0.6596457 
# AUDIT_C , Set  37 
# Dimensions of complete case matrix:  7380 56 
# Dimensions after merging with AUDIT:  2215 58 
# Test set correlation btw predicted and actual:  0.6626375 
# 20  variables retained 
# Full set correlation btw predicted and actual:  0.6422458 
# AUDIT_C , Set  38 
# Dimensions of complete case matrix:  7380 61 
# Dimensions after merging with AUDIT:  2215 63 
# Test set correlation btw predicted and actual:  0.6632329 
# 22  variables retained 
# Full set correlation btw predicted and actual:  0.6435876 
# AUDIT_C , Set  39 
# Dimensions of complete case matrix:  7380 66 
# Dimensions after merging with AUDIT:  2215 68 
# Test set correlation btw predicted and actual:  0.6629536 
# 23  variables retained 
# Full set correlation btw predicted and actual:  0.6437024 
# AUDIT_C , Set  40 
# Dimensions of complete case matrix:  7380 71 
# Dimensions after merging with AUDIT:  2215 73 
# Test set correlation btw predicted and actual:  0.6629536 
# 23  variables retained 
# Full set correlation btw predicted and actual:  0.6437024 
# AUDIT_C , Set  41 
# Dimensions of complete case matrix:  3311 72 
# Dimensions after merging with AUDIT:  1066 74 
# Test set correlation btw predicted and actual:  0.5482012 
# 23  variables retained 
# Full set correlation btw predicted and actual:  0.6497146 
# AUDIT_C , Set  42 
# Dimensions of complete case matrix:  1586 73 
# Dimensions after merging with AUDIT:  574 75 
# Test set correlation btw predicted and actual:  0.5792567 
# 10  variables retained 
# Full set correlation btw predicted and actual:  0.5636586 
# AUDIT_C , Set  43 
# Dimensions of complete case matrix:  760 81 
# Dimensions after merging with AUDIT:  256 83 
# Test set correlation btw predicted and actual:  0.6511354 
# 11  variables retained 
# Full set correlation btw predicted and actual:  0.7186735 
# AUDIT_P , Set  2 
# Dimensions of complete case matrix:  501645 4 
# Dimensions after merging with AUDIT:  157089 6 
# Test set correlation btw predicted and actual:  0.2075831 
# 2  variables retained 
# Full set correlation btw predicted and actual:  0.210108 
# AUDIT_P , Set  3 
# Dimensions of complete case matrix:  501639 10 
# Dimensions after merging with AUDIT:  157088 12 
# Test set correlation btw predicted and actual:  0.3199296 
# 3  variables retained 
# Full set correlation btw predicted and actual:  0.321301 
# AUDIT_P , Set  4 
# Dimensions of complete case matrix:  501629 13 
# Dimensions after merging with AUDIT:  157088 15 
# Test set correlation btw predicted and actual:  0.3301084 
# 4  variables retained 
# Full set correlation btw predicted and actual:  0.3305902 
# AUDIT_P , Set  5 
# Dimensions of complete case matrix:  501512 14 
# Dimensions after merging with AUDIT:  157067 16 
# Test set correlation btw predicted and actual:  0.3386491 
# 5  variables retained 
# Full set correlation btw predicted and actual:  0.3363997 
# AUDIT_P , Set  6 
# Dimensions of complete case matrix:  499724 15 
# Dimensions after merging with AUDIT:  156837 17 
# Test set correlation btw predicted and actual:  0.3547507 
# 6  variables retained 
# Full set correlation btw predicted and actual:  0.3487327 
# AUDIT_P , Set  7 
# Dimensions of complete case matrix:  497759 16 
# Dimensions after merging with AUDIT:  156548 18 
# Test set correlation btw predicted and actual:  0.3438425 
# 7  variables retained 
# Full set correlation btw predicted and actual:  0.349338 
# AUDIT_P , Set  8 
# Dimensions of complete case matrix:  489703 17 
# Dimensions after merging with AUDIT:  154505 19 
# Test set correlation btw predicted and actual:  0.3494441 
# 8  variables retained 
# Full set correlation btw predicted and actual:  0.3512235 
# AUDIT_P , Set  9 
# Dimensions of complete case matrix:  467267 18 
# Dimensions after merging with AUDIT:  148214 20 
# Test set correlation btw predicted and actual:  0.3589755 
# 9  variables retained 
# Full set correlation btw predicted and actual:  0.3596084 
# AUDIT_P , Set  10 
# Dimensions of complete case matrix:  446935 19 
# Dimensions after merging with AUDIT:  141689 21 
# Test set correlation btw predicted and actual:  0.3706367 
# 10  variables retained 
# Full set correlation btw predicted and actual:  0.3669333 
# AUDIT_P , Set  11 
# Dimensions of complete case matrix:  446291 20 
# Dimensions after merging with AUDIT:  141504 22 
# Test set correlation btw predicted and actual:  0.3655089 
# 11  variables retained 
# Full set correlation btw predicted and actual:  0.3690363 
# AUDIT_P , Set  12 
# Dimensions of complete case matrix:  444751 21 
# Dimensions after merging with AUDIT:  141046 23 
# Test set correlation btw predicted and actual:  0.3718494 
# 12  variables retained 
# Full set correlation btw predicted and actual:  0.3690657 
# AUDIT_P , Set  13 
# Dimensions of complete case matrix:  357287 22 
# Dimensions after merging with AUDIT:  117944 24 
# Test set correlation btw predicted and actual:  0.3797169 
# 13  variables retained 
# Full set correlation btw predicted and actual:  0.3790589 
# AUDIT_P , Set  14 
# Dimensions of complete case matrix:  280039 26 
# Dimensions after merging with AUDIT:  97110 28 
# Test set correlation btw predicted and actual:  0.3865632 
# 14  variables retained 
# Full set correlation btw predicted and actual:  0.3875286 
# AUDIT_P , Set  15 
# Dimensions of complete case matrix:  252356 27 
# Dimensions after merging with AUDIT:  88187 29 
# Test set correlation btw predicted and actual:  0.4084047 
# 15  variables retained 
# Full set correlation btw predicted and actual:  0.4103257 
# AUDIT_P , Set  16 
# Dimensions of complete case matrix:  252356 28 
# Dimensions after merging with AUDIT:  88187 30 
# Test set correlation btw predicted and actual:  0.4376156 
# 16  variables retained 
# Full set correlation btw predicted and actual:  0.4409822 
# AUDIT_P , Set  17 
# Dimensions of complete case matrix:  252356 29 
# Dimensions after merging with AUDIT:  88187 31 
# Test set correlation btw predicted and actual:  0.456881 
# 17  variables retained 
# Full set correlation btw predicted and actual:  0.4667126 
# AUDIT_P , Set  18 
# Dimensions of complete case matrix:  252356 30 
# Dimensions after merging with AUDIT:  88187 32 
# Test set correlation btw predicted and actual:  0.4674318 
# 18  variables retained 
# Full set correlation btw predicted and actual:  0.4770713 
# AUDIT_P , Set  19 
# Dimensions of complete case matrix:  252356 31 
# Dimensions after merging with AUDIT:  88187 33 
# Test set correlation btw predicted and actual:  0.4684872 
# 19  variables retained 
# Full set correlation btw predicted and actual:  0.4778064 
                                # AUDIT_P , Set  20 
                                # Dimensions of complete case matrix:  147063 32 
                                # Dimensions after merging with AUDIT:  48890 34 
                                # Test set correlation btw predicted and actual:  0.4784857 
                                # 14  variables retained 
                                # Full set correlation btw predicted and actual:  0.4647833 
# AUDIT_P , Set  21 
# Dimensions of complete case matrix:  71542 33 
# Dimensions after merging with AUDIT:  26002 35 
# Test set correlation btw predicted and actual:  0.4623286 
# 14  variables retained 
# Full set correlation btw predicted and actual:  0.4523023 
# AUDIT_P , Set  22 
# Dimensions of complete case matrix:  71542 38 
# Dimensions after merging with AUDIT:  26002 40 
# Test set correlation btw predicted and actual:  0.4623286 
# 14  variables retained 
# Full set correlation btw predicted and actual:  0.4523023 
# AUDIT_P , Set  23 
# Dimensions of complete case matrix:  71542 39 
# Dimensions after merging with AUDIT:  26002 41 
# Test set correlation btw predicted and actual:  0.4623286 
# 14  variables retained 
# Full set correlation btw predicted and actual:  0.4523023 
# AUDIT_P , Set  24 
# Dimensions of complete case matrix:  71542 40 
# Dimensions after merging with AUDIT:  26002 42 
# Test set correlation btw predicted and actual:  0.4623286 
# 14  variables retained 
# Full set correlation btw predicted and actual:  0.4523023 
# AUDIT_P , Set  25 
# Dimensions of complete case matrix:  39805 43 
# Dimensions after merging with AUDIT:  13484 45 
# Test set correlation btw predicted and actual:  0.4662938 
# 14  variables retained 
# Full set correlation btw predicted and actual:  0.4598631 
# AUDIT_P , Set  26 
# Dimensions of complete case matrix:  39805 47 
# Dimensions after merging with AUDIT:  13484 49 
# Test set correlation btw predicted and actual:  0.4662938 
# 14  variables retained 
# Full set correlation btw predicted and actual:  0.4598631 
# AUDIT_P , Set  27 
# Dimensions of complete case matrix:  39805 51 
# Dimensions after merging with AUDIT:  13484 53 
# Test set correlation btw predicted and actual:  0.4662938 
# 14  variables retained 
# Full set correlation btw predicted and actual:  0.4598631 
# AUDIT_P , Set  28 
# Dimensions of complete case matrix:  39805 56 
# Dimensions after merging with AUDIT:  13484 58 
# Test set correlation btw predicted and actual:  0.4662938 
# 14  variables retained 
# Full set correlation btw predicted and actual:  0.4598631 
# AUDIT_P , Set  29 
# Dimensions of complete case matrix:  39805 57 
# Dimensions after merging with AUDIT:  13484 59 
# Test set correlation btw predicted and actual:  0.4662938 
# 14  variables retained 
# Full set correlation btw predicted and actual:  0.4598631 
# AUDIT_P , Set  30 
# Dimensions of complete case matrix:  17562 63 
# Dimensions after merging with AUDIT:  5429 65 
# Test set correlation btw predicted and actual:  0.3982243 
# 26  variables retained 
# Full set correlation btw predicted and actual:  0.4607969 
