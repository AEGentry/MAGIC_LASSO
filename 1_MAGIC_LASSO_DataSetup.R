#######################################################
### MAGIC LASSO ###
### Script written by Amanda Elswick Gentry and B. Todd Webb
### For associated manuscript describing the method, see: https://www.biorxiv.org/content/10.1101/2021.04.29.442057v1
#######################################################

### Step 1 ###
### Prepare UKB data for analysis

#######################################################
### Point to your directories
workdir <- "your_directory_name_here"
UKBdir23704 <- paste(workdir, "ukb23704/", sep="")
UKBdir30302 <- paste(workdir, "ukb30302/", sep="")
plotdir <- paste(workdir, "Plots/", sep="")

### Read in the UKB data showcase, as downloaded from the showcase webpage
data.show <- read.csv(paste(workdir, "Data_Dictionary_Showcase_downloaded072419.csv", sep=""))

### Read the first line from the R.tab files containing the phenotypes. 
### For our purposes here, we're interested in reading in the "top level" variables only, ie, 
### the variables from the first measurement instance of each phenotype. We do some 
### transformation on the variable names from these files in order to match what's in the data 
### showcase. This is done multiple times here, once for each batch of phenotypes we received from UKB, 
### and the variable names are combined into a data dictionary.

### Here, we are reading in data from two separate "sends" from UKB

phenoFile1<-paste(UKBdir30302,"ukb30302_R.tab",sep="")
bdKey<-names(fread(phenoFile1, nrows=1))
bdNames <- gsub("f.", "", bdKey)
bdNames <- gsub("\\..*", "", bdNames)
data.dict30302 <- data.frame(UDI=bdKey, FieldID=bdNames)
data.dict30302 <- merge(data.dict30302, data.show, by="FieldID", all.x=TRUE)
data.dict30302$release <- 30302
data.dict30302$top.lev <- 0 # initialize a column
data.dict30302$top.lev[grep("\\.0\\.", data.dict30302$UDI)] <- 1

phenoFile2<-paste(UKBdir23704,"ukb23704_R.tab",sep="")
bdKey<-names(fread(phenoFile2, nrows=1))
bdNames <- gsub("f.", "", bdKey)
bdNames <- gsub("\\..*", "", bdNames)
data.dict23704 <- data.frame(UDI=bdKey, FieldID=bdNames)
data.dict23704 <- merge(data.dict23704, data.show, by="FieldID", all.x=TRUE)
data.dict23704$release <- 23704
data.dict23704$top.lev <- 0 # initialize a column
data.dict23704$top.lev[grep("\\.0\\.", data.dict23704$UDI)] <- 1

### Combine the two data dictionaries
data.dict <- rbind(data.dict23704, data.dict30302)
### Identify any duplicate columns (should only be the eid column)
which(duplicated(data.dict$UDI))
data.dict <- data.dict[-which(duplicated(data.dict$UDI)), ]
data.dict$top.lev[grep("eid", data.dict$UDI)] <- 1

### Load the data.
### The fread() function reads in the data from the original files most efficiently. Then we convert 
### those objects to dataframes and merge them together and save that dataframe as an R object. Note 
### that this object contains all the data in its original, numeric, unlabeled format.

bd23704<-fread(phenoFile2)
setDF(bd23704)

bd30302<-fread(phenoFile1)
setDF(bd30302)

bd <- merge(bd23704, bd30302, by="f.eid")
save(bd, file=paste(workdir, "UKB_23704_30302noLabels.RData", sep=""))

load(paste(workdir, "UKB_23704_30302noLabels.RData", sep=""))
dim(bd)
### 502536   9613

### Quantify missingness.
### For analysis purposes, we want to quantify the missingness for each variable.
### Here, simply count the missing values for each measure.

### Count Missing Variables 
nVars <- dim(bd)[2]
### Divide the variables into chunks, for efficiency
chunkSz <- 100
nSets <- ceiling(nVars/chunkSz)
starting <- (0:(nSets - 1) * chunkSz + 1)
ending <- c(1:(nSets - 1) * chunkSz, nVars)
nChunks <- length(starting)
chunkSzA <- ending - starting + 1

miss.by.col <- numeric()
kk <- 0
for (ii in 1:nChunks){
  for (jj in starting[ii]:ending[ii]){
    print(jj)
    print(ii)
    aa <- bd[,jj]
    miss.by.col[1 + kk] <- sum(is.na(aa))
    kk <- kk+1
  }
  eval(parse(text = paste("miss.by.col", ii, " <- ", "miss.by.col", sep="")))
  kk<-0
  if (ii==1){
    miss.total <- miss.by.col
  } else {
    miss.total <- c(miss.total, miss.by.col)
  }
  rm(miss.by.col)
  miss.by.col <- numeric()
}

save(miss.total, file=paste(workdir, "UKB23704_30302_MissByCol.Rds", sep=""))
### Make a dataframe so that the column names may be saved
miss.total.DF <- data.frame(Col=colnames(bd), Nmiss=miss.total, stringsAsFactors=FALSE)
save(miss.total.DF, file=paste(workdir, "UKB23704_30302_MissByCol_DF.Rds", sep=""))

### For the AUDIT Lasso application, we want to segregate the top-level variables for the initial analysis. 
### This includes first-instance measures (where longitudinal measurements exist) and variables which will 
### be appropriate to include in a Lasso model. Most variables in UKB are categorical as factors or ordered 
### factors (ie, any multiple choice question with a single response allowed will be coded as a factor, or 
### an ordered factor for scale items.) Where measures are R-class "categorical," these are free-answer type 
### responses and not appropriate for inclusion in our model. Other variables from UKB appear as UKB-class 
### "Array" measures. These are questions which allow multiple ("check all that apply") responses, among other 
### things. These, in their current coding, are not appropriate for model inclusion, they must first be 
### transformed into a series of Y/N variables, which we have not attempted yet.

### The ICD codes and related variables are in this UKB "Array" format, so we separate and save these separately 
### for later use. Additionally, we separate the AUDIT variables and save these separately.

### Separate ICD codes
bd.ICD <- cbind(bd$f.eid, bd[, colnames(bd) %in% data.dict$UDI[grep("ICD", data.dict$Field)]])
save(bd.ICD, file=paste(workdir, "UKB_ICD.Rds", sep=""))

AUDIT.vars <- c("f.20414.0.0", "f.20403.0.0", "f.20416.0.0", "f.20413.0.0", "f.20407.0.0",
                "f.20412.0.0", "f.20409.0.0", "f.20408.0.0", "f.20411.0.0", "f.20405.0.0")
AUDIT <- bd[, colnames(bd) %in% AUDIT.vars]
rownames(AUDIT) <- bd$f.eid

###dim(bd)
### [1] 502536   9613
bd.top <- bd[, miss.total < (dim(bd)[1] - 100000)]
dim(bd.top)
### [1] 502536   1122
### Rm ICDs
isICD <- which(colnames(bd.top) %in% data.dict$UDI[grep("ICD", data.dict$Field)])
bd.top <- bd.top[, -isICD]
dim(bd.top)
### [1] 502536   1074
### Char/Date classes
isChar <- which(sapply(sapply(bd.top, class), function(x) x[1]) == "character")
bd.top <- bd.top[, -isChar]
dim(bd.top)
### [1] 502536   1010
# isDate <- which(sapply(sapply(bd.top, class), function(x) x[1]) == "Date")
# bd.top <- bd.top[, -isDate]
### AUDIT vars
AUDIT.vars <- c("f.20414.0.0", "f.20403.0.0", "f.20416.0.0", "f.20413.0.0", "f.20407.0.0",
                "f.20412.0.0", "f.20409.0.0", "f.20408.0.0", "f.20411.0.0", "f.20405.0.0")
isAUDIT <-which(colnames(bd.top) %in% AUDIT.vars)
bd.top <- bd.top[, -isAUDIT]
dim(bd.top)
### [1] 502536   1005
### Top level, as define above
isNOTtop <- which(colnames(bd.top) %in% data.dict$UDI[data.dict$top.lev == 0])
bd.top <-bd.top[, -isNOTtop]
dim(bd.top)
### [1] 502536    977
### Array variables
isArray <- which(colnames(bd.top) %in% data.dict$UDI[data.dict$Array > 1])
# data.dict$Field[data.dict$UDI %in% colnames(bd.top)[colnames(bd.top) %in% data.dict$UDI[data.dict$Array > 1]]]
bd.top <- bd.top[, -isArray]
dim(bd.top)
### [1] 502536    638
### Zero variance
isZero <- apply(bd.top, 2, var, na.rm=TRUE)
bd.top <- bd.top[, isZero > 0]
dim(bd.top)
# [1] 502536    632

### Label UKB measures.
### The raw UKB data is largely numeric. Most of the variables, however, are categorical and the category 
### labels and (where appropriate) category ordering is stored in a separate file. Here, we extract the 
### labels/levels of these variables and apply the coding to the variables in our dataset. This is 
### computationaly burdensome and takes some time. Isolate and run this portion in a separate job, if necessary.

###Where is master annotation file
masterFile1<-paste(UKBdir23704,"ukb23704_R.r",sep="")
masterFile2<-paste(UKBdir30302,"ukb30302_R.r",sep="")

###Pull annotate commands
varKey23704 <- colnames(bd.top)[colnames(bd.top) %in% data.dict$UDI[data.dict$release == 23704]]
varKey30302 <- colnames(bd.top)[colnames(bd.top) %in% data.dict$UDI[data.dict$release == 30302]]
nVars1<-length(varKey23704)
nVars2<-length(varKey30302)

for (ii in 1:nVars1){
  print(ii)
  #comT<-system(paste("grep -F -B1 ",varKeyWW[ii]," ", masterFile1,sep=""), intern=T)
  comT<-system(paste("grep -F ",varKey23704[ii]," ", masterFile1,sep=""), intern=T)
  if (ii==1) comW<-comT else comW<-rbind(comW,comT)
}
### No factor variables in this set
# for (ii in 1:nVars2){
# print(ii)
# #comT<-system(paste("grep -F -B1 ",varKeyWW[ii]," ", masterFile1,sep=""), intern=T)
# comT2<-system(paste("grep -F ",varKey30302[ii]," ", masterFile2,sep=""), intern=T)
# if (ii==1) comW2<-comT2 else comW2<-rbind(comW2,comT2)
# }
###Find Levels and Lables
crapT<-do.call(rbind,strsplit(comW," ",fixed=TRUE))
lvlT<-gsub("levels=","",crapT[,4][grep("levels=", crapT[,4])])
lvlT<-gsub(",","",lvlT)
lblT<-gsub("labels=","",crapT[,5][grep("labels=", crapT[,5])])
lblT<-gsub(")","",lblT)
lblW<-c(lvlT,lblT)
lblW<-names(table(lblW))

###Pull Level and Label commands
nlbl<-length(lblW)
for (jj in 1:nlbl){
  print(jj)
  lblComT<-system(paste("awk \'$1== \"",lblW[jj],"\"{print $0}\' ", masterFile1,sep=""), intern=T)
  if (jj==1) lblComW<-lblComT else lblComW<-rbind(lblComW,lblComT)
}

###Update variables with labels and levels
comAll<-c(lblComW,comW)
nComs<-length(comAll)

for (kk in 1:nComs){
  print(kk)
  eval(parse(text=comAll[kk]))
}

###The code above executes to label the appropriate factors (a subset of variables we're interested in) from the full bd dataset
###To achieve bd.top (subsetted above), with labels, just subset the actual dataframe again

bd.top <- bd[, colnames(bd) %in% colnames(bd.top)]
bd.top.classes <- sapply(sapply(bd.top, class), function(x) x[1])

###Save the file
save(bd.top, file=paste(workdir, "UKB_23704_30302Top_Labeled.RData", sep=""))

### Create the AUDIT dataset
### We isolated the AUDIT variables from the subsetted data above because we do not intend to use 
### the AUDIT variables to predict themselves (obviously.) We take the AUDIT measures here and 
### label the responses using the same process as above. We also score total AUDIT, AUDIT-C, 
### and AUDIT-P. For our analysis, we consider AUDIT to be "present" if AUDIT Question 1 was 
### answered (this means excluding subjects who were not administered the AUDIT questionairre 
### at all [N=345,179])

### Repeat the same process for the AUDIT data, which was part of the 23704 send
nVarsA<-length(AUDIT.vars)
for (ii in 1:nVarsA){
  print(ii)
  comT<-system(paste("grep -F ",AUDIT.vars[ii]," ", masterFile1,sep=""), intern=T)
  if (ii==1) comW<-comT else comW<-rbind(comW,comT)
}
###Find Levels and Lables
crapT<-do.call(rbind,strsplit(comW," ",fixed=TRUE))
lvlT<-gsub("levels=","",crapT[,4][grep("levels=", crapT[,4])])
lvlT<-gsub(",","",lvlT)
lblT<-gsub("labels=","",crapT[,5][grep("labels=", crapT[,5])])
lblT<-gsub(")","",lblT)
lblW<-c(lvlT,lblT)
lblW<-names(table(lblW))
###Pull Level and Label commands
nlbl<-length(lblW)
for (jj in 1:nlbl){
  print(jj)
  lblComT<-system(paste("awk \'$1== \"",lblW[jj],"\"{print $0}\' ", masterFile1,sep=""), intern=T)
  if (jj==1) lblComW<-lblComT else lblComW<-rbind(lblComW,lblComT)
}
comAll<-c(lblComW,comW)
nComs<-length(comAll)
comAll <- gsub("bd", "AUDIT", comAll)

for (kk in 1:nComs){
  print(kk)
  eval(parse(text=comAll[kk]))
}
### Save AUDIT before scoring or removing missing variables
RawAUDIT <- AUDIT
RawAUDIT <- data.frame(AUDIT1 = RawAUDIT$f.20414.0.0,
                       AUDIT2 = RawAUDIT$f.20403.0.0,
                       AUDIT3 = RawAUDIT$f.20416.0.0,
                       AUDIT4 = RawAUDIT$f.20413.0.0,
                       AUDIT5 = RawAUDIT$f.20407.0.0,
                       AUDIT6 = RawAUDIT$f.20412.0.0,
                       AUDIT7 = RawAUDIT$f.20409.0.0,
                       AUDIT8 = RawAUDIT$f.20408.0.0,
                       AUDIT9 = RawAUDIT$f.20411.0.0,
                       AUDIT10 = RawAUDIT$f.20405.0.0,
                       row.names = rownames(AUDIT))
save(RawAUDIT, file=paste(workdir, "RawAUDIT.Rds",sep=""))
load(paste(workdir, "RawAUDIT.Rds", sep=""))

### Look at some variables
### AUDIT1 and AUDIT9
table(as.numeric(AUDIT$f.20414.0.0), as.numeric(AUDIT$f.20411.0.0), useNA="ifany")
### AUDIT1 and AUDIT10
table(as.numeric(AUDIT$f.20414.0.0), as.numeric(AUDIT$f.20405.0.0), useNA="ifany")
AUDIT1 <- ifelse(as.numeric(AUDIT$f.20414.0.0) == 1,
                 -1, as.numeric(AUDIT$f.20414.0.0) - 2)
AUDIT2 <- ifelse(as.numeric(AUDIT$f.20403.0.0) == 1, 
                 -1, as.numeric(AUDIT$f.20403.0.0) - 2)
AUDIT3 <- ifelse(as.numeric(AUDIT$f.20416.0.0) == 1,
                 -1, as.numeric(AUDIT$f.20416.0.0) - 2)
AUDIT4 <- ifelse(as.numeric(AUDIT$f.20413.0.0) == 1,
                 -1, as.numeric(AUDIT$f.20413.0.0) - 2)
AUDIT5 <- ifelse(as.numeric(AUDIT$f.20407.0.0) == 1,
                 -1, as.numeric(AUDIT$f.20407.0.0) - 2)
AUDIT6 <- ifelse(as.numeric(AUDIT$f.20412.0.0) == 1,
                 -1, as.numeric(AUDIT$f.20412.0.0) - 2)
AUDIT7 <- ifelse(as.numeric(AUDIT$f.20409.0.0) == 1,
                 -1, as.numeric(AUDIT$f.20409.0.0) - 2)
AUDIT8 <- ifelse(as.numeric(AUDIT$f.20408.0.0) == 1,
                 -1, as.numeric(AUDIT$f.20408.0.0) - 2)
AUDIT9 <- ifelse(as.numeric(AUDIT$f.20411.0.0) == 1, -1,
                 ifelse(as.numeric(AUDIT$f.20411.0.0) == 2, 0,
                        ifelse(as.numeric(AUDIT$f.20411.0.0) == 3, 2,
                               as.numeric(AUDIT$f.20411.0.0))))
AUDIT10 <- ifelse(as.numeric(AUDIT$f.20405.0.0) == 1, -1,
                  ifelse(as.numeric(AUDIT$f.20405.0.0) == 2, 0,
                         ifelse(as.numeric(AUDIT$f.20405.0.0) == 3, 2,
                                as.numeric(AUDIT$f.20405.0.0))))
AUDIT <- data.frame(AUDIT1=AUDIT1,
                    AUDIT2=AUDIT2,
                    AUDIT3=AUDIT3,
                    AUDIT4=AUDIT4,
                    AUDIT5=AUDIT5,
                    AUDIT6=AUDIT6,
                    AUDIT7=AUDIT7,
                    AUDIT8=AUDIT8,
                    AUDIT9=AUDIT9,
                    AUDIT10=AUDIT10,
                    row.names=rownames(AUDIT))
### Add the ID's to the AUDIT dataframe
AUDIT$f.eid <- rownames(AUDIT)

### Drop the Q1 skips
# CountSkips <- AUDIT == -1
# CountSkips[is.na(CountSkips)] <- FALSE
# summary(CountSkips[,1])
#     #    Mode   FALSE    TRUE 
#     # logical  502341     195 
# AUDIT$NoSkips <- 0
# AUDIT$NoSkips <- apply(CountSkips, 1, sum)
# sum(AUDIT$NoSkips > 0) ### 1320

AUDIT$measured <- 0
AUDIT$measured <- ifelse(!is.na(AUDIT$AUDIT1), 1,0) ### AUDIT Q1 present, at least
keep <- which(AUDIT$measured == 1)

# keep.a <- which(AUDIT$measured == 1 & AUDIT$NoSkips < 1)
test<-AUDIT
AUDIT <- AUDIT[keep,]
keep2 <- which(AUDIT$AUDIT1 == -1)
AUDIT <- AUDIT[-keep2,]
AUDIT.old[! (as.character(AUDIT.old$f.eid) %in% AUDIT$f.eid),]

AUDIT[,2:10][is.na(AUDIT[,2:10])] <- 0
AUDIT[AUDIT == -1] <- 0
### Calculate the AUDIT score
for (i in 1:dim(AUDIT)[1]){
  if (AUDIT$AUDIT1[i] == 0) {
    AUDIT$score[i] <- AUDIT$AUDIT9[i] + AUDIT$AUDIT10[i]
  } else if (AUDIT$AUDIT2[i] + AUDIT$AUDIT3[i] < 1) {
    AUDIT$score[i] <- AUDIT$AUDIT1[i] + AUDIT$AUDIT9[i] + AUDIT$AUDIT10[i]
  } else {
    AUDIT$score[i] <- sum(AUDIT[i,1:10])
  }
}

AUDIT$AudC <- AUDIT$AUDIT1 + AUDIT$AUDIT2 + AUDIT$AUDIT3
AUDIT$AudP <- AUDIT$AUDIT4 + AUDIT$AUDIT5 + AUDIT$AUDIT6 +
  AUDIT$AUDIT7 + AUDIT$AUDIT8 + AUDIT$AUDIT9 +
  AUDIT$AUDIT10

save(AUDIT, file=paste(workdir, "UKB_AUDIT.Rds", sep=""))

### Quantify Missingness in the top level variables
### For computation of missingness, remove the ID variable
### EID variable
isID <- which(colnames(bd.top) == "f.eid")
bd.top.noID <- bd.top[, -isID]

nVars <- dim(bd.top.noID)[2]
### Divide the variables into chunks
chunkSz <- 100
nSets <- ceiling(nVars/chunkSz)
starting <- (0:(nSets - 1) * chunkSz + 1)
ending <- c(1:(nSets - 1) * chunkSz, nVars)
nChunks <- length(starting)
chunkSzA <- ending - starting + 1

miss.by.col <- numeric()
kk <- 0
for (ii in 1:nChunks){
  for (jj in starting[ii]:ending[ii]){
    print(jj)
    print(ii)
    aa <- bd.top.noID[,jj]
    miss.by.col[1 + kk] <- sum(is.na(aa))
    kk <- kk+1
  }
  #eval(parse(text = paste("miss.by.col", ii, " <- ", "miss.by.col", sep="")))
  kk<-0
  if (ii==1){
    miss.total <- miss.by.col
  } else {
    miss.total <- c(miss.total, miss.by.col)
  }
  rm(miss.by.col)
  miss.by.col <- numeric()
}
miss.total.top<-miss.total
save(miss.total.top, file=paste(workdir, "UKB23704_30302_MissByCol_Top.Rds", sep=""))
### Make a dataframe so that the column names may be saved
miss.total.topDF <- data.frame(Col=colnames(bd.top.noID), Nmiss=miss.total.top, stringsAsFactors=FALSE)
save(miss.total.topDF, file=paste(workdir, "UKB23704_30302_MissByCol_DF_Top.Rds", sep=""))

### Create Missingness Indicator Matrix 
bd.miss.top <- is.na(bd.top.noID)
chunkSz <- 100
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
colnames(bd.miss.ind.top) <- colnames(bd.top.noID)
save(bd.miss.ind.top, file=paste(workdir, "UKB23704_30302_MissIndMatrix_Top.Rds", sep=""))

### Create Missingness Correlation Matrix 
### Top-Level Vars 
### Create the missingness correlation matrix here, but highly recommend running this portion as
### a separate job because it's computationally burdensome

workdir <- "/home/projects/UKB_23704/20181126/Phenos/working/AUDITproject/"
UKBdir23704 <- "/home/projects/UKB_23704/20181126/Phenos/ukb23704/"
UKBdir30302 <- "/home/projects/UKB_23704/20181126/Phenos/ukb30302/"
plotdir <- "/home/projects/UKB_23704/20181126/Phenos/working/AUDITproject/Plots/"
load(paste(workdir, "UKB23704_30302_MissIndMatrix_Top.Rds", sep=""))
### To avoid NA values in the correlation calculations, add two dummy subjects
dummy.all <- rep(0, dim(bd.miss.ind.top)[2])
dummy.none <- rep(1, dim(bd.miss.ind.top)[2])
bd.miss.indA <- rbind(dummy.all, dummy.none, bd.miss.ind.top)
library(propagate)
bd.miss.cor <- bigcor(bd.miss.indA, fun="cor", size=100)
save(bd.miss.cor, file=paste(workdir, "UKB23704_30302_MissCor.Rds", sep=""))

############################################################
### In the indicator matrix, 1 indicates missing and 0 indicates present
### Here, count how many complete cases exist for each pair of variables by
### adding two indicator columns together and counting how many zeros there are
### Code below split and run in two separate jobs
### Part 1
workdir <- "/home/projects/UKB_23704/20181126/Phenos/working/AUDITproject/"
UKBdir23704 <- "/home/projects/UKB_23704/20181126/Phenos/ukb23704/"
UKBdir30302 <- "/home/projects/UKB_23704/20181126/Phenos/ukb30302/"
plotdir <- "/home/projects/UKB_23704/20181126/Phenos/working/AUDITproject/Plots/"
load(paste(workdir, "UKB23704_30302_MissIndMatrix_Top.Rds", sep=""))
### Make a matrix for pairwise missingness counts
nVars <- dim(bd.miss.ind.top)[2]
cut.point <- 300
pairwise.missA <- matrix(nrow=cut.point, ncol=nVars)
colnames(pairwise.missA) <-  colnames(bd.miss.ind.top)
for (ii in 1:cut.point){
  cat("Calculating row ", ii, "\n", sep="")
  pairwise.missA[ii, 1:ii] <- 0
  for (jj in 1:(nVars - ii)){
    pairwise.missA[ii, ii + jj] <- sum(bd.miss.ind.top[,ii] + bd.miss.ind.top[, ii + jj] == 0)
  }
}
### Save the output
save(pairwise.missA, file=paste(workdir, "UKB23704_30302_Pairwise_Top_Pt1.Rds", sep=""))

### Part 2
workdir <- "/home/projects/UKB_23704/20181126/Phenos/working/AUDITproject/"
UKBdir23704 <- "/home/projects/UKB_23704/20181126/Phenos/ukb23704/"
UKBdir30302 <- "/home/projects/UKB_23704/20181126/Phenos/ukb30302/"
plotdir <- "/home/projects/UKB_23704/20181126/Phenos/working/AUDITproject/Plots/"
load(paste(workdir, "UKB23704_30302_MissIndMatrix_Top.Rds", sep=""))
### Make a matrix for pairwise missingness counts
nVars <- dim(bd.miss.ind.top)[2]
cut.point <- 300
pairwise.missB <- matrix(0, nrow=(nVars-cut.point), ncol=nVars)
colnames(pairwise.missB) <-  colnames(bd.miss.ind.top)
for (ii in (cut.point + 1):(nVars-1)){
  cat("Calculating row ", ii, "\n", sep="")
  #pairwise.missB[(ii - cut.point), 1:ii] <- 0
  for (jj in 1:(nVars - ii)){
    # print(jj)
    pairwise.missB[(ii - cut.point), ii + jj] <- sum(bd.miss.ind.top[,ii] + bd.miss.ind.top[, ii + jj] == 0)
  }
}
### Save the output
save(pairwise.missB, file=paste(workdir, "UKB23704_30302_Pairwise_Top_Pt2.Rds", sep=""))

### Make a new matrix with correlations on the lower triangle and
### pairwise present (not-missing) counts on the upper triangle
### These are calculated in two parts
### "pairwise.missA" and "pairwise.missB"
load(paste(workdir, "UKB23704_30302_Pairwise_Top_Pt1.Rds", sep=""))
load(paste(workdir, "UKB23704_30302_Pairwise_Top_Pt2.Rds", sep=""))
### Combine these into a single object
bd.miss.pairwise.top <- rbind(pairwise.missA, pairwise.missB)

### Reorder the missing values DF to match the matrix
miss.total.DF <- miss.total.topDF[match(colnames(bd.miss.pairwise.top), miss.total.topDF$Col), ]
### check
all.equal(colnames(bd.miss.pairwise.top), miss.total.topDF$Col)
### Set the diagonal of the pairwise-present count matrix to be the number of present values for that variable
### why this wasn't calculated in the first place, I cannot say
diag(bd.miss.pairwise.top) <- max(miss.total.topDF$Nmiss) - miss.total.topDF$Nmiss
library(Matrix)
bd.miss.pairwise.topS <- forceSymmetric(bd.miss.pairwise.top, uplo="U") ### make symmetric from the upper
save(bd.miss.pairwise.topS, file=paste(workdir, "UKB23704_30302_Pairwise_Mat.Rds", sep=""))


