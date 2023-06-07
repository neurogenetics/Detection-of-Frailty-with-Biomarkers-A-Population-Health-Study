#ARIC Data Analysis 6-6-23


getwd()
library(tidyverse)
library(haven)
library(rio)
library(foreign)
library(readxl)
library(magrittr)
library(aod)

#number of NA in each column remove data >15% and data check
getwd()
full_data <- read.csv("Full_data_frame_ARIC.csv", header = T)
head(full_data)
full_data %<>% mutate_if(is.factor,as.numeric)
str(full_data)
nrow(full_data)
sum(is.na(full_data))
sum(is.na(full_data$V2PTH))
sum(is.na(full_data$CESD51))
sum(is.na(full_data$parkinson))
biomarker.missing.data <- map_dbl(full_data, ~sum(is.na(.)))
head(biomarker.missing.data)
write.csv(biomarker.missing.data, file="Missing_Data_Full_data_frame_ARIC.csv")

#Frailty number of NA in each column remove data >15% and data check
getwd()
Frailty_full_data <- read.csv("Frailty_Full_data_frame_ARIC.csv", header = T)
head(Frailty_full_data)
str(Frailty_full_data)
Frailty_full_data %<>% mutate_if(is.factor,as.numeric)
nrow(Frailty_full_data)
sum(is.na(Frailty_full_data))
Frailty.missing.data <- map_dbl(Frailty_full_data, ~sum(is.na(.)))
head(Frailty.missing.data)
write.csv(Frailty.missing.data, file="Missing_Data_Frailty_data_frame_ARIC.csv")

########################write dataframe for genoML######################################
getwd()
Frailty_full_data <- read.csv("Frailty_Full_data_frame_ARIC.csv", header = T)
head(Frailty_full_data)
Frailty_full_data %<>% mutate_if(is.factor,as.numeric)
str(Frailty_full_data)
#######Count for Frail dataframe with the total n with exact match variables
Frail_exact_match <- Frailty_full_data[,c("ID","AGE","V3_VITD","FOL_WO","FOLA","B6","CHMA08","HDLSIU51",
                                          "FRAIL","CHM45","LIP23","THYR8","ALT_V4","HMEB2","IL6","V2PTH","tnf","ACB_total_clean","CESD_Code")]
head(Frail_exact_match)
Frail.drop.na <- Frail_exact_match %>% drop_na(Frail_exact_match)

#make frail dataframe write out as csv and table
Frail_full_data1 <- Frailty_full_data[,c("ID","AGE","EDUCATIO","HMTB13","V1_VITD","FOLA","V3_VITD","B6","E_IU","CHMA08","HDLSIU51",
                                          "FRAIL","CHM45","LIP23","THYR8","ALT_V4","TESTOSTE","MCP1","IL6","IL1","V2PTH","tnf","ACB_total_clean","CESD_Code")]
head(Frail_full_data1)
write.csv(Frail_full_data1, file="Frail_Model4.csv")
Frail <- read.csv("Frail_Model4.csv", header = TRUE)
head(Frail)
Frail.drop.na <- Frail %>% drop_na(FRAIL)
head(Frail.drop.na)
dim(Frail.drop.na)
write.csv(Frail.drop.na, file="Frail.drop.na.model4.csv")
write.table(Frail.drop.na, file = "Frail.drop.na")
frail_drop <- read.table("Frail.drop.na", header=TRUE)
head(frail_drop)

#make Prefrail dataframe write out as csv and table
PreFrail_full_data1 <- Frailty_full_data[,c("ID","AGE","EDUCATIO","HMTB13","V1_VITD","FOLA","V3_VITD","B6","E_IU","CHMA08","HDLSIU51",
                                         "PRE.FRAIL","CHM45","LIP23","THYR8","ALT_V4","TESTOSTE","MCP1","IL6","IL1","V2PTH","tnf","ACB_total_clean","CESD_Code")]
head(PreFrail_full_data1)
write.csv(PreFrail_full_data1, file="PreFrail_Model4.csv")
PreFrail <- read.csv("PreFrail_Model4.csv", header = TRUE)
head(PreFrail)
PreFrail.drop.na <- PreFrail %>% drop_na(PRE.FRAIL)
head(PreFrail.drop.na)
dim(PreFrail.drop.na)
write.csv(PreFrail.drop.na, file="PreFrail.drop.na.model4.csv")
write.table(PreFrail.drop.na, file = "PreFrail.drop.na")
Prefrail_drop <- read.table("PreFrail.drop.na", header=TRUE)
head(Prefrail_drop)

#one hot encoding gender, race
Gender <- within(full_data,GENDER2 <- match(full_data$GENDER,unique(full_data$GENDER)))
write.csv(Gender, file="Gender_one_hot_code.csv")
Race <- within(full_data,RACEGRP <- match(full_data$RACEGRP,unique(full_data$RACEGRP)))
write.csv(Race, file="Race_one_hot_code.csv")

str(full_data)
#Convert data from int to factor for phenotype Frail and pre-frail
full_data$PRE.FRAIL <- factor(full_data$PRE.FRAIL)
full_data$FRAIL <- factor(full_data$FRAIL)

#not sure if needed can delete later 
full_data$CHM45 <- as.numeric(full_data$CHM45)
full_data$CHMA08 <- as.numeric(full_data$CHMA08)
full_data$FOL_WO <- as.numeric(full_data$FOL_WO)
full_data$ACB_total_clean <- as.numeric(full_data$ACB_total_clean)
full_data$ACB_Total_adj <- as.numeric(full_data$ACB_Total_adj)
full_data$HMTB13 <- as.numeric(full_data$HMTB13)
full_data$LIP23 <- as.numeric(full_data$LIP23)

############Build Frail model then build Pre-frail model off of the Frail parameters#############

##Full_data_Boosting Model##
Full_Final_Dataframe <- read.csv("Full_data_frame_ARIC.csv", header = T)
head(Full_Final_Dataframe)
Final_Dataframe <- Full_Final_Dataframe[,c(1,3,5,8,10:25,29,30,31:34,38:56,59,62:64)]
head(Final_Dataframe)
ncol(Final_Dataframe)
nrow(Final_Dataframe)
Final_Dataframe$seed <- rep(1:10,length.out = length(Final_Dataframe$ID))
full_train <- subset(Final_Dataframe, seed <= 7)
full_test <- subset(Final_Dataframe, seed >7)
full_train.ids <- paste(full_train$ID, sep = " ")
full_test.ids <- paste(full_test$ID, sep = " ")
#write.table(full_train.ids, file = "full_train.ids", quote = F, row.names = F, col.names = F)
#write.table(full_test.ids, file = "full_test.ids", quote = F, row.names = F, col.names = F)
#write out the data.frames as test and train
#write.table (full_train, file ="full_train.data", quote = F, sep = "\t", row.names = F)
#write.table (full_test, file ="full_test.data", quote = F, sep = "\t", row.names = F)

##START HERE###FINAL MODELS################################################################################################################################################################
####now quit R, restart and get modeling####
require(parallel)
require(data.table)
require(caret)
require(corrplot)
require(Rtsne)
require(xgboost)
require(stats)
require(knitr)
require(ggplot2)
require(pROC)
library(tidyr)
library(xgboostExplainer)

#ncores <- detectCores() - 4
ncores <- 12

##############################################################Model 1 Frail Run model with only variables predictive for specific outcome (pre-frail/frail) Variables#################################################################

#####Frail model 1
#training set
train.load7 <- read.table("full_train.data", header = T)
head(train.load7)
dim(train.load7)
#TOOK OUT FOLATE AND VIT D V1
train.load1 <- train.load7[,c("ID","AGE","FRAIL","V3_VITD","FOLA","B6","CHMA08","HDLSIU51","CHM45","LIP23","THYR8","ALT_V4","V2PTH","ACB_total_clean","CESD_Code")]
head(train.load1)
dim(train.load1)
sum(is.na(train.load1$FRAIL))
train.load <- train.load1 %>% drop_na(FRAIL)
sum(is.na(train.load$FRAIL))
####Rename columns for final model graphs and tables
names(train.load) <- c("ID","Age","FRAIL","Vitamin D","Folate","Vitamin B6","Blood Urea Nitrogen","HDL Cholesterol","Urine Creatinine","Fasting Glucose","Free T4","Alanine Aminotransferase ALT","Parathyoid Hormone","Anticholinergic Drug Burden","Depression")
head(train.load)
dim(train.load)
train <- as.data.table(train.load)
train[, c("ID"):=NULL]
#train.outcome <- as.matrix(as.integer(train$FRAIL))
train.outcome <- as.integer(train$FRAIL)
train[, c("FRAIL"):=NULL]
train.matrix = as.matrix(train)
mode(train.matrix) = "numeric"
head(train.matrix)
#test set
test.load3 <- read.table("full_test.data", header = T)
test.load1 <- test.load3[,c("ID","AGE","FRAIL","V3_VITD","FOLA","B6","CHMA08","HDLSIU51","CHM45","LIP23","THYR8","ALT_V4","V2PTH","ACB_total_clean","CESD_Code")]
head(test.load1)
dim(test.load1)
sum(is.na(test.load1$FRAIL))
test.load <- test.load1 %>% drop_na(FRAIL)
sum(is.na(test.load$FRAIL))
#use these to create classes for confusion matrix 
FRAIL.full_test <- test.load$FRAIL
####Rename columns for final model graphs and tables
names(test.load) <- c("ID","Age","FRAIL","Vitamin D","Folate","Vitamin B6","Blood Urea Nitrogen","HDL Cholesterol","Urine Creatinine","Fasting Glucose","Free T4","Alanine Aminotransferase ALT","Parathyoid Hormone","Anticholinergic Drug Burden","Depression")
head(test.load)
dim(test.load)
test <- as.data.table(test.load)
test[, c("ID"):=NULL]
#test.outcome <- as.matrix(as.integer(test$FRAIL))
test.outcome <- as.integer(test$FRAIL)
test[, c("FRAIL"):=NULL]
test.matrix = as.matrix(test)
mode(test.matrix) = "numeric"
head(test.matrix)

####now run testing a variety of rounds and depths####
##testing parameters for the modeling###
Frail7 <- matrix(NA,90,3)
count <- 0
for(i in 1:10)
{
  for(j in c(5,10,20,30,40,50,75,100,200))
  {
    count <- count + 1    
    boosted <- xgboost(data = train.matrix, label = train.outcome, max.depth = i, nround = j, nthread = ncores, objective = "binary:logistic", missing = NA, eval_metric = "auc",scale_pos_weight = 2)
    predicted <- predict(boosted, test.matrix, missing = NA)
    testAuc <- auc(test.outcome, predicted)
    #print(paste("tested max.depth =",i,"and nround =",j,"resulting in a test AUC =",testAuc, sep = " "))
    Frail7[count,1] <- i
    Frail7[count,2] <- j
    Frail7[count,3] <- testAuc  
  }
}
#q("no")
##95% CIs for AUC 
ci.auc(testAuc)
getwd()
max(Frail7[,3])
save(Frail7,file = "Rounds_depths_Final.AUC.Rdata")
load("Rounds_depths_Final.AUC.Rdata")
write.csv(Frail7,file = "Rounds_depths_Final.AUC.csv")
#saveimage

#################################################################################
Frail7 <- matrix(NA,1,3)
count <- 0
for(i in 3)
{
  for(j in c(10))
  {
    count <- count + 1    
    boosted <- xgboost(data = train.matrix, label = train.outcome, max.depth = i, nround = j, nthread = ncores, objective = "binary:logistic", missing = NA, eval_metric = "auc", scale_pos_weight = 6)
    predicted <- predict(boosted, test.matrix, missing = NA)
    testAuc <- auc(test.outcome, predicted)
    #print(paste("tested max.depth =",i,"and nround =",j,"resulting in a test AUC =",testAuc, sep = " "))
    Frail7[count,1] <- i
    Frail7[count,2] <- j
    Frail7[count,3] <- testAuc  
  }
}
#q("no")
ci.auc(testAuc)
getwd()
max(Frail7[,3])
save(boosted, file = "Final_Frail_boosted_Model")
head(Frail7)
save(Frail7,file = "Frail_Final.AUC.Rdata")
load("Frail_Final.AUC.Rdata")
write.csv(Frail7,file = "Frail_Final_Model.AUC.csv")
##get the trained model and save##
xgb.dump(boosted, 'xgb.model.dump', with_stats = TRUE)
# get the feature real names
names = dimnames(train.matrix)[[2]]
# compute feature importance matrix
Adjusted_Frail.model.importance = xgb.importance(names, model=boosted)
print(Adjusted_Frail.model.importance) 
# plot
library(Ckmeans.1d.dp)
library(caret)
xgb.ggplot.importance(Adjusted_Frail.model.importance, top_n = 15, measure = "Gain")
#plot trees
library("DiagrammeR")
#plots 2 trees
xgb.plot.tree(model = boosted, trees = 1, feature_names = colnames(names))
#plots multitrees
xgb.plot.multi.trees(model = boosted, n_first_tree = 1, feature_names = colnames(names))
write.csv(Adjusted_Frail.model.importance, file = "Adjusted_Frail.model_AGE.importance.BEST.11.1.20.csv")

#test out the ranges of the predicted classes
(predicted[FRAIL.full_test == 1])
range(predicted[FRAIL.full_test == 1])
(predicted[FRAIL.full_test == 0])
range(predicted[FRAIL.full_test == 0])
range(predicted)
#set the cutoff threshold based on the mean of the classes
mean(predicted[FRAIL.full_test == 1])
mean(predicted[FRAIL.full_test == 0])

# Set our cutoff threshold
pred.resp <- factor(ifelse(predicted >= mean(predicted[FRAIL.full_test == 1]), 1, 0)) # makes this set: takes exactly the mean and uses the mean as the cut off or class 1 

# Create the confusion matrix
confusionMatrix(pred.resp, factor(FRAIL.full_test), positive="1")
#F1, precision, recall
confusionMatrix(pred.resp, factor(FRAIL.full_test), mode = "prec_recall", positive="1")
levels(pred.resp)
levels(FRAIL.full_test)


########################################################################## Model Pre-frail Run model###############################################

Prefrail8_train.load3 <- read.table("full_train.data", header = T)
head(Prefrail8_train.load3)
dim(Prefrail8_train.load3)
#TOOK OUT VITD_V1 AND FOLATE 
PreFrail8_train.load1 <- Prefrail8_train.load3[,c("ID","AGE","PRE.FRAIL","HMTB13","FOLA","V3_VITD","B6","E_IU","TESTOSTE","CHM45","ACB_total_clean","CESD_Code")]
head(PreFrail8_train.load1)
dim(PreFrail8_train.load1)
sum(is.na(PreFrail8_train.load1$PRE.FRAIL))
PreFrail8_train.load <- PreFrail8_train.load1 %>% drop_na(PRE.FRAIL)
str(PreFrail8_train.load)
names(PreFrail8_train.load) <- c("ID","Age","PRE.FRAIL","Mean Corpuscular Volume","Folate","Vitamin D","Vitamin B6","Vitamin E","Testosterone","Urine Creatinine","Anticholinergic Drug Burden","Depression")
head(PreFrail8_train.load)
dim(PreFrail8_train.load)
train <- as.data.table(PreFrail8_train.load)
train[, c("ID"):=NULL]
#train.outcome <- as.matrix(as.integer(train$PRE.FRAIL))
train.outcome <- as.integer(train$PRE.FRAIL)
train[, c("PRE.FRAIL"):=NULL]
train.matrix = as.matrix(train)
mode(train.matrix) = "numeric"
head(train.matrix)
#test set
PreFrail8_test.load3 <- read.table("full_test.data", header = T)
PreFrail8_test.load1 <- PreFrail8_test.load3[,c("ID","AGE","PRE.FRAIL","HMTB13","FOLA","V3_VITD","B6","E_IU","TESTOSTE","CHM45","ACB_total_clean","CESD_Code")]
head(PreFrail8_test.load1)
dim(PreFrail8_test.load1)
sum(is.na(PreFrail8_test.load1$PRE.FRAIL))
test.load <- PreFrail8_test.load1 %>% drop_na(PRE.FRAIL)
#use these to create classes for confusion matrix 
PRE.FRAIL.full_test <- test.load$PRE.FRAIL
names(test.load) <- c("ID","Age","PRE.FRAIL","Mean Corpuscular Volume","Folate","Vitamin D","Vitamin B6","Vitamin E","Testosterone","Urine Creatinine","Anticholinergic Drug Burden","Depression")
head(test.load)
test <- as.data.table(test.load)
test[, c("ID"):=NULL]
#test.outcome <- as.matrix(as.integer(test$PRE.FRAIL))
test.outcome <- as.integer(test$PRE.FRAIL)
test[, c("PRE.FRAIL"):=NULL]
test.matrix = as.matrix(test)
mode(test.matrix) = "numeric"

############################################Prefrail rounds and depths only used for stepwise#######################################
Prefrail.stepwise <- matrix(NA,90,3)
count <- 0
for(i in 1:10)
{
  for(j in c(5,10,20,30,40,50,75,100,200))
  {
    count <- count + 1    
    boosted <- xgboost(data = train.matrix, label = train.outcome, max.depth = i, nround = j, nthread = ncores, objective = "binary:logistic", missing = NA, eval_metric = "auc")
    predicted <- predict(boosted, test.matrix, missing = NA)
    testAuc <- auc(test.outcome, predicted)
    #print(paste("tested max.depth =",i,"and nround =",j,"resulting in a test AUC =",testAuc, sep = " "))
    Prefrail.stepwise[count,1] <- i
    Prefrail.stepwise[count,2] <- j
    Prefrail.stepwise[count,3] <- testAuc  
  }
}
#q("no")
##95% CIs for AUC = See model stepwise xgboost document
ci.auc(testAuc)
getwd()
max(Prefrail.stepwise[,3])

##testing parameters for the modeling###
PreFrail8 <- matrix(NA,1,3)
count <- 0
for(i in 2)
{
  for(j in c(20))
  {
    count <- count + 1    
    boosted <- xgboost(data = train.matrix, label = train.outcome, max.depth = i, nround = j, nthread = ncores, objective = "binary:logistic", missing = NA, eval_metric = "auc")
    predicted <- predict(boosted, test.matrix, missing = NA)
    testAuc <- auc(test.outcome, predicted)
    #print(paste("tested max.depth =",i,"and nround =",j,"resulting in a test AUC =",testAuc, sep = " "))
    PreFrail8[count,1] <- i
    PreFrail8[count,2] <- j
    PreFrail8[count,3] <- testAuc  
  }
}
#q("no")
ci.auc(testAuc)
getwd()
max(PreFrail8[,3])
save(boosted, file = "Final_PreFrail_boosted_Model")
head(PreFrail8)
save(PreFrail8,file = "PreFrail_Final_Model.AUC.Rdata")
load("PreFrail_Final_Model.AUC.Rdata")
write.csv(PreFrail8,file = "PreFrail_Final_Model.AUC.csv")
##get the trained model and save##
xgb.dump(boosted, 'xgb.model.dump', with_stats = TRUE)
# get the feature real names
names = dimnames(train.matrix)[[2]]
# compute feature importance matrix
importance_matrix = xgb.importance(names, model=boosted)
# plot
Adjusted_PreFrail.model.importance = xgb.plot.importance(importance_matrix)
print(Adjusted_PreFrail.model.importance) 
# plot
library(Ckmeans.1d.dp)
library(caret)
xgb.ggplot.importance(Adjusted_PreFrail.model.importance, top_n = 15, measure = "Gain")
#plot trees
write.csv(importance_matrix, file = "Adjusted_PreFrail.model.importance.BEST.9.10.19.csv")


#test out the ranges of the predicted classes
(predicted[PRE.FRAIL.full_test == 1])
range(predicted[PRE.FRAIL.full_test == 1])
(predicted[PRE.FRAIL.full_test == 0])
range(predicted[PRE.FRAIL.full_test == 0])
range(predicted)
#set the cutoff threshold based on the mean of the classes
mean(predicted[PRE.FRAIL.full_test == 1])
mean(predicted[PRE.FRAIL.full_test == 0])

# Set our cutoff threshold
pred.resp <- factor(ifelse(predicted >= mean(predicted[PRE.FRAIL.full_test == 1]), 1, 0)) # makes this set: takes exactly the mean and uses the mean as the cut off or class 1 

# Create the confusion matrix
confusionMatrix(pred.resp, factor(PRE.FRAIL.full_test), positive="1")
#F1, precision, recall
confusionMatrix(pred.resp, factor(PRE.FRAIL.full_test), mode = "prec_recall", positive="1")
levels(pred.resp)
levels(PRE.FRAIL.full_test)

#####Frail model ##########################################################Build Frail model for race = black/AA############################################################################## 

#training set
train.load <- read.table("B.full_train.data", header = T)
head(train.load)
dim(train.load)
train.load1 <- train.load[,c("ID","AGE","V3_VITD","FOLA","B6","CHMA08","HDLSIU51",
                             "FRAIL","CHM45","LIP23","THYR8","ALT_V4","V2PTH","ACB_total_clean","CESD_Code")]

head(train.load1)
dim(train.load1)
sum(is.na(train.load1$FRAIL))
train.load <- train.load1 %>% drop_na(FRAIL)
sum(is.na(train.load$FRAIL))
head(train.load)
dim(train.load)
train <- as.data.table(train.load)
train[, c("ID"):=NULL]
#train.outcome <- as.matrix(as.integer(train$FRAIL))
train.outcome <- as.integer(train$FRAIL)
train[, c("FRAIL"):=NULL]
train.matrix = as.matrix(train)
mode(train.matrix) = "numeric"
head(train.matrix)
#test set
test.load <- read.table("B.full_test.data", header = T)
head(test.load)
dim(test.load)
test.load1 <- test.load[,c("ID","AGE","V3_VITD","FOLA","B6","CHMA08","HDLSIU51",
                           "FRAIL","CHM45","LIP23","THYR8","ALT_V4","V2PTH","ACB_total_clean","CESD_Code")]
head(test.load1)
dim(test.load1)
sum(is.na(test.load1$FRAIL))
test.load <- test.load1 %>% drop_na(FRAIL)
sum(is.na(test.load$FRAIL))
#use these to create classes for confusion matrix 
B.FRAIL.full_test <- test.load$FRAIL
head(test.load)
dim(test.load)
test <- as.data.table(test.load)
test[, c("ID"):=NULL]
#test.outcome <- as.matrix(as.integer(test$FRAIL))
test.outcome <- as.integer(test$FRAIL)
test[, c("FRAIL"):=NULL]
test.matrix = as.matrix(test)
mode(test.matrix) = "numeric"
head(test.matrix)

####now run testing a variety of rounds and depths####
##testing parameters for the modeling###
Frail.B <- matrix(NA,90,3)
count <- 0
for(i in 1:10)
{
  for(j in c(5,10,20,30,40,50,75,100,200))
  {
    count <- count + 1    
    boosted <- xgboost(data = train.matrix, label = train.outcome, max.depth = i, nround = j, nthread = ncores, objective = "binary:logistic", missing = NA, eval_metric = "auc")
    predicted <- predict(boosted, test.matrix, missing = NA)
    testAuc <- auc(test.outcome, predicted)
    #print(paste("tested max.depth =",i,"and nround =",j,"resulting in a test AUC =",testAuc, sep = " "))
    Frail.B[count,1] <- i
    Frail.B[count,2] <- j
    Frail.B[count,3] <- testAuc  
  }
}
#q("no")
##95% CIs for AUC 
ci.auc(testAuc)
getwd()
max(Frail.B[,3])
save(Frail.B,file = "Black.Rounds_depths_Final.AUC.Rdata")
load("Black.Rounds_depths_Final.AUC.Rdata")
write.csv(Frail.B,file = "Frail.B.Rounds_depths_Final.AUC.csv")
#saveimage

##get the trained model and save##
xgb.dump(boosted, 'xgb.model.dump', with_stats = TRUE)
# get the feature real names
names = dimnames(train.matrix)[[2]]
# compute feature importance matrix
importance_matrix = xgb.importance(names, model=boosted)
# plot
Frail.modelB.importance = xgb.plot.importance(importance_matrix)
print(Frail.modelB.importance) 
write.csv(Frail.modelB.importance, file = "ModelB.Rounds_depth_Frail.importance.csv")
#make histogram of Frail model and save data###

##testing parameters for the modeling###

Frail.B <- matrix(NA,1,3)
count <- 0
for(i in 1)
{
  for(j in c(200))
  {
    count <- count + 1    
    boosted <- xgboost(data = train.matrix, label = train.outcome, max.depth = i, nround = j, nthread = ncores, objective = "binary:logistic", missing = NA, eval_metric = "auc")
    predicted <- predict(boosted, test.matrix, missing = NA)
    testAuc <- auc(test.outcome, predicted)
    #print(paste("tested max.depth =",i,"and nround =",j,"resulting in a test AUC =",testAuc, sep = " "))
    Frail.B[count,1] <- i
    Frail.B[count,2] <- j
    Frail.B[count,3] <- testAuc  
  }
}
#q("no")
##95% CIs for Max AUC 
ci.auc(testAuc)
getwd()
max(Frail.B[,3])
save(boosted, file = "Final_Frail.Black_boosted")
head(Frail.B)
save(Frail.B,file = "Frail_Black.AUC.Rdata")
load("Frail_Black.AUC.Rdata")
write.csv(Frail.B,file = "Frail_Black_Final.AUC.csv")
##get the trained model and save##
xgb.dump(boosted, 'xgb.model.dump', with_stats = TRUE)
# get the feature real names
names = dimnames(train.matrix)[[2]]
# compute feature importance matrix
importance_matrix = xgb.importance(names, model=boosted)
# plot
Adjusted_Frail.Black.importance = xgb.plot.importance(importance_matrix)
print(Adjusted_Frail.Black.importance) 
write.csv(importance_matrix, file = "Adjusted_Frail.Black.importance.csv")

#test out the ranges of the predicted classes
(predicted[B.FRAIL.full_test == 1])
range(predicted[B.FRAIL.full_test == 1])
(predicted[B.FRAIL.full_test == 0])
range(predicted[B.FRAIL.full_test == 0])
range(predicted)
#set the cutoff threshold based on the mean of the classes
mean(predicted[B.FRAIL.full_test == 1])
mean(predicted[B.FRAIL.full_test == 0])

# Set our cutoff threshold
pred.resp <- factor(ifelse(predicted >= mean(predicted[B.FRAIL.full_test == 1]), 1, 0)) # makes this set: takes exactly the mean and uses the mean as the cut off or class 1 

# Create the confusion matrix
confusionMatrix(pred.resp, factor(B.FRAIL.full_test), positive="1")
#F1, precision, recall
confusionMatrix(pred.resp, factor(B.FRAIL.full_test), mode = "prec_recall", positive="1")
levels(pred.resp)
levels(B.FRAIL.full_test)

############################################################################## PreFrail Black/AA#############################

Prefrail.B_train.load <- read.table("B.full_train.data", header = T)
head(Prefrail.B_train.load)
dim(Prefrail.B_train.load)
Prefrail.B_train.load1 <- Prefrail.B_train.load[,c("ID","AGE","HMTB13","FOLA","V3_VITD","B6","E_IU",
                                                  "PRE.FRAIL","CHM45","TESTOSTE","ACB_total_clean","CESD_Code")]
head(Prefrail.B_train.load1)
dim(Prefrail.B_train.load1)
sum(is.na(Prefrail.B_train.load1$PRE.FRAIL))
Prefrail.B_train.load <- Prefrail.B_train.load1 %>% drop_na(PRE.FRAIL)
str(Prefrail.B_train.load)
head(Prefrail.B_train.load)
dim(Prefrail.B_train.load)
train <- as.data.table(Prefrail.B_train.load)
train[, c("ID"):=NULL]
#train.outcome <- as.matrix(as.integer(train$PRE.FRAIL))
train.outcome <- as.integer(train$PRE.FRAIL)
train[, c("PRE.FRAIL"):=NULL]
train.matrix = as.matrix(train)
mode(train.matrix) = "numeric"
head(train.matrix)
#test set
PreFrail.B_test.load <- read.table("B.full_test.data", header = T)
PreFrail.B_test.load1 <- PreFrail.B_test.load[,c("ID","AGE","HMTB13","FOLA","V3_VITD","B6","E_IU",
                                                "PRE.FRAIL","CHM45","TESTOSTE","ACB_total_clean","CESD_Code")]
head(PreFrail.B_test.load1)
dim(PreFrail.B_test.load1)
sum(is.na(PreFrail.B_test.load1$PRE.FRAIL))
test.load <- PreFrail.B_test.load1 %>% drop_na(PRE.FRAIL)
head(test.load)
#use these to create classes for confusion matrix 
B.PRE.FRAIL.full_test <- test.load$PRE.FRAIL
test <- as.data.table(test.load)
test[, c("ID"):=NULL]
#test.outcome <- as.matrix(as.integer(test$PRE.FRAIL))
test.outcome <- as.integer(test$PRE.FRAIL)
test[, c("PRE.FRAIL"):=NULL]
test.matrix = as.matrix(test)
mode(test.matrix) = "numeric"

##################################################
Prefrail.B <- matrix(NA,90,3)
count <- 0
for(i in 1:10)
{
  for(j in c(5,10,20,30,40,50,75,100,200))
  {
    count <- count + 1    
    boosted <- xgboost(data = train.matrix, label = train.outcome, max.depth = i, nround = j, nthread = ncores, objective = "binary:logistic", missing = NA, eval_metric = "auc")
    predicted <- predict(boosted, test.matrix, missing = NA)
    testAuc <- auc(test.outcome, predicted)
    #print(paste("tested max.depth =",i,"and nround =",j,"resulting in a test AUC =",testAuc, sep = " "))
    Prefrail.B[count,1] <- i
    Prefrail.B[count,2] <- j
    Prefrail.B[count,3] <- testAuc  
  }
}
#q("no")
ci.auc(testAuc)
getwd()
max(Prefrail.B[,3])
save(Prefrail.B,file = "Black.Rounds_depths_Final.AUC.Rdata")
load("Black.Rounds_depths_Final.AUC.Rdata")
write.csv(Prefrail.B,file = "Prefrail.B.Rounds_depths_Final.AUC.csv")
#saveimage


##testing parameters for the modeling################################################################Using parameters from rounds.it ###################################
PreFrail.B <- matrix(NA,1,3)
count <- 0
for(i in 2)
{
  for(j in c(75))
  {
    count <- count + 1    
    boosted <- xgboost(data = train.matrix, label = train.outcome, max.depth = i, nround = j, nthread = ncores, objective = "binary:logistic", missing = NA, eval_metric = "auc")
    predicted <- predict(boosted, test.matrix, missing = NA)
    testAuc <- auc(test.outcome, predicted)
    #print(paste("tested max.depth =",i,"and nround =",j,"resulting in a test AUC =",testAuc, sep = " "))
    PreFrail.B[count,1] <- i
    PreFrail.B[count,2] <- j
    PreFrail.B[count,3] <- testAuc  
  }
}
#q("no")
##95% CIs for Max AUC
ci.auc(testAuc)
getwd()
max(PreFrail.B[,3])
save(boosted, file = "Final_PreFrail_Black.boosted")
head(PreFrail.B)
save(PreFrail.B,file = "PreFrail_Final_Black.AUC.Rdata")
load("PreFrail_Final_Black.AUC.Rdata")
write.csv(PreFrail.B,file = "PreFrail_Final_Black.AUC.csv")
##get the trained model and save##
xgb.dump(boosted, 'xgb.model.dump', with_stats = TRUE)
# get the feature real names
names = dimnames(train.matrix)[[2]]
# compute feature importance matrix
importance_matrix = xgb.importance(names, model=boosted)
# plot
Adjusted_PreFrail.Black.importance = xgb.plot.importance(importance_matrix)
print(Adjusted_PreFrail.Black.importance) 
write.csv(importance_matrix, file = "Adjusted_PreFrail.black.importance.csv")

#test out the ranges of the predicted classes
(predicted[B.PRE.FRAIL.full_test == 1])
range(predicted[B.PRE.FRAIL.full_test == 1])
(predicted[B.PRE.FRAIL.full_test == 0])
range(predicted[B.PRE.FRAIL.full_test == 0])
range(predicted)
#set the cutoff threshold based on the mean of the classes
mean(predicted[B.PRE.FRAIL.full_test == 1])
mean(predicted[B.PRE.FRAIL.full_test == 0])

# Set our cutoff threshold
pred.resp <- factor(ifelse(predicted >= mean(predicted[B.PRE.FRAIL.full_test == 1]), 1, 0)) # makes this set: takes exactly the mean and uses the mean as the cut off or class 1 

# Create the confusion matrix
confusionMatrix(pred.resp, factor(B.PRE.FRAIL.full_test), positive="1")
#F1, precision, recall
confusionMatrix(pred.resp, factor(B.PRE.FRAIL.full_test), mode = "prec_recall", positive="1")
levels(pred.resp)
levels(B.PRE.FRAIL.full_test)
###############################################################################Build Frail model for race = White###################################################################

##Train,test data for race = white 
train.load <- read.table("W.full_train.data", header = T)
head(train.load)
dim(train.load)
train.load1 <- train.load[,c("ID","AGE","V3_VITD","FOLA","B6","CHMA08","HDLSIU51",
                             "FRAIL","CHM45","LIP23","THYR8","ALT_V4","V2PTH","ACB_total_clean","CESD_Code")]
head(train.load1)
dim(train.load1)
sum(is.na(train.load1$FRAIL))
train.load <- train.load1 %>% drop_na(FRAIL)
sum(is.na(train.load$FRAIL))
head(train.load)
dim(train.load)
train <- as.data.table(train.load)
train[, c("ID"):=NULL]
#train.outcome <- as.matrix(as.integer(train$FRAIL))
train.outcome <- as.integer(train$FRAIL)
train[, c("FRAIL"):=NULL]
train.matrix = as.matrix(train)
mode(train.matrix) = "numeric"
head(train.matrix)
#test set
test.load <- read.table("W.full_test.data", header = T)
head(test.load)
dim(test.load)
test.load1 <- test.load[,c("ID","AGE","V3_VITD","FOLA","B6","CHMA08","HDLSIU51",
                           "FRAIL","CHM45","LIP23","THYR8","ALT_V4","V2PTH","ACB_total_clean","CESD_Code")]
head(test.load1)
dim(test.load1)
sum(is.na(test.load1$FRAIL))
test.load <- test.load1 %>% drop_na(FRAIL)
sum(is.na(test.load$FRAIL))
#use these to create classes for confusion matrix 
W.FRAIL.full_test <- test.load$FRAIL
head(test.load)
dim(test.load)
test <- as.data.table(test.load)
test[, c("ID"):=NULL]
#test.outcome <- as.matrix(as.integer(test$FRAIL))
test.outcome <- as.integer(test$FRAIL)
test[, c("FRAIL"):=NULL]
test.matrix = as.matrix(test)
mode(test.matrix) = "numeric"
head(test.matrix)

####now run testing a variety of rounds and depths####
##testing parameters for the modeling###
Frail.W <- matrix(NA,90,3)
count <- 0
for(i in 1:10)
{
  for(j in c(5,10,20,30,40,50,75,100,200))
  {
    count <- count + 1    
    boosted <- xgboost(data = train.matrix, label = train.outcome, max.depth = i, nround = j, nthread = ncores, objective = "binary:logistic", missing = NA, eval_metric = "auc")
    predicted <- predict(boosted, test.matrix, missing = NA)
    testAuc <- auc(test.outcome, predicted)
    #print(paste("tested max.depth =",i,"and nround =",j,"resulting in a test AUC =",testAuc, sep = " "))
    Frail.W[count,1] <- i
    Frail.W[count,2] <- j
    Frail.W[count,3] <- testAuc  
  }
}
#q("no")
##95% CIs for AUC 
ci.auc(testAuc)
getwd()
max(Frail.W[,3])
save(Frail.W,file = "White.Rounds_depths_Final.AUC.Rdata")
load("White.Rounds_depths_Final.AUC.Rdata")
write.csv(Frail.W,file = "Frail.W.Rounds_depths_Final.AUC.csv")
#saveimage

##get the trained model and save##
xgb.dump(boosted, 'xgb.model.dump', with_stats = TRUE)
# get the feature real names
names = dimnames(train.matrix)[[2]]
# compute feature importance matrix
importance_matrix = xgb.importance(names, model=boosted)
# plot
Frail.modelW.importance = xgb.plot.importance(importance_matrix)
print(Frail.modelW.importance) 
write.csv(Frail.modelW.importance, file = "ModelW.Rounds_depth_Frail.importance.csv")
#make histogram of Frail model and save data###

##testing parameters for the modeling###

Frail.W <- matrix(NA,1,3)
count <- 0
for(i in 1)
{
  for(j in c(30))
  {
    count <- count + 1    
    boosted <- xgboost(data = train.matrix, label = train.outcome, max.depth = i, nround = j, nthread = ncores, objective = "binary:logistic", missing = NA, eval_metric = "auc")
    predicted <- predict(boosted, test.matrix, missing = NA)
    testAuc <- auc(test.outcome, predicted)
    #print(paste("tested max.depth =",i,"and nround =",j,"resulting in a test AUC =",testAuc, sep = " "))
    Frail.W[count,1] <- i
    Frail.W[count,2] <- j
    Frail.W[count,3] <- testAuc  
  }
}
#q("no")
##95% CIs for Max AUC 
ci.auc(testAuc)
getwd()
max(Frail.W[,3])
save(boosted, file = "Final_Frail.White_boosted")
head(Frail.W)
save(Frail.W,file = "Frail_White.AUC.Rdata")
load("Frail_White.AUC.Rdata")
write.csv(Frail.W,file = "Frail_White_Final.AUC.csv")
##get the trained model and save##
xgb.dump(boosted, 'xgb.model.dump', with_stats = TRUE)
# get the feature real names
names = dimnames(train.matrix)[[2]]
# compute feature importance matrix
importance_matrix = xgb.importance(names, model=boosted)
# plot
Adjusted_Frail.White.importance = xgb.plot.importance(importance_matrix)
print(Adjusted_Frail.White.importance) 
write.csv(importance_matrix, file = "Adjusted_Frail.White.importance.csv")

#test out the ranges of the predicted classes
(predicted[W.FRAIL.full_test == 1])
range(predicted[W.FRAIL.full_test == 1])
(predicted[W.FRAIL.full_test == 0])
range(predicted[W.FRAIL.full_test == 0])
range(predicted)
#set the cutoff threshold based on the mean of the classes
mean(predicted[W.FRAIL.full_test == 1])
mean(predicted[W.FRAIL.full_test == 0])

# Set our cutoff threshold
pred.resp <- factor(ifelse(predicted >= mean(predicted[W.FRAIL.full_test == 1]), 1, 0)) # makes this set: takes exactly the mean and uses the mean as the cut off or class 1 

# Create the confusion matrix
confusionMatrix(pred.resp, factor(W.FRAIL.full_test), positive="1")
#F1, precision, recall
confusionMatrix(pred.resp, factor(W.FRAIL.full_test), mode = "prec_recall", positive="1")
levels(pred.resp)
levels(W.FRAIL.full_test)
########################################################################################PreFrail Race = White############################

Prefrail4_train.load3 <- read.table("W.full_train.data", header = T)
head(Prefrail4_train.load3)
dim(Prefrail4_train.load3)
PreFrail4_train.load1 <- Prefrail4_train.load3[,c("ID","AGE","HMTB13","FOLA","V3_VITD","B6","E_IU",
                                                  "PRE.FRAIL","CHM45","TESTOSTE","ACB_total_clean","CESD_Code")]
head(PreFrail4_train.load1)
dim(PreFrail4_train.load1)
sum(is.na(PreFrail4_train.load1$PRE.FRAIL))
PreFrail4_train.load <- PreFrail4_train.load1 %>% drop_na(PRE.FRAIL)
str(PreFrail4_train.load)
head(PreFrail4_train.load)
dim(PreFrail4_train.load)
train <- as.data.table(PreFrail4_train.load)
train[, c("ID"):=NULL]
#train.outcome <- as.matrix(as.integer(train$PRE.FRAIL))
train.outcome <- as.integer(train$PRE.FRAIL)
train[, c("PRE.FRAIL"):=NULL]
train.matrix = as.matrix(train)
mode(train.matrix) = "numeric"
head(train.matrix)
#test set
PreFrail4_test.load3 <- read.table("W.full_test.data", header = T)
PreFrail4_test.load1 <- PreFrail4_test.load3[,c("ID","AGE","HMTB13","FOLA","V3_VITD","B6","E_IU",
                                                "PRE.FRAIL","CHM45","TESTOSTE","ACB_total_clean","CESD_Code")]
head(PreFrail4_test.load1)
dim(PreFrail4_test.load1)
sum(is.na(PreFrail4_test.load1$PRE.FRAIL))
test.load <- PreFrail4_test.load1 %>% drop_na(PRE.FRAIL)
head(test.load)
#use these to create classes for confusion matrix 
W.PRE.FRAIL.full_test <- test.load$PRE.FRAIL
test <- as.data.table(test.load)
test[, c("ID"):=NULL]
#test.outcome <- as.matrix(as.integer(test$PRE.FRAIL))
test.outcome <- as.integer(test$PRE.FRAIL)
test[, c("PRE.FRAIL"):=NULL]
test.matrix = as.matrix(test)
mode(test.matrix) = "numeric"

##testing parameters for the modeling################################################Based off of white pre-frail rounds.it. #############################################
PreFrail.W <- matrix(NA,1,3)
count <- 0
for(i in 1)
{
  for(j in c(20))
  {
    count <- count + 1    
    boosted <- xgboost(data = train.matrix, label = train.outcome, max.depth = i, nround = j, nthread = ncores, objective = "binary:logistic", missing = NA, eval_metric = "auc")
    predicted <- predict(boosted, test.matrix, missing = NA)
    testAuc <- auc(test.outcome, predicted)
    #print(paste("tested max.depth =",i,"and nround =",j,"resulting in a test AUC =",testAuc, sep = " "))
    PreFrail.W[count,1] <- i
    PreFrail.W[count,2] <- j
    PreFrail.W[count,3] <- testAuc  
  }
}
#q("no")
##95% CIs for Max AUC
ci.auc(testAuc)
getwd()
max(PreFrail.W[,3])
save(boosted, file = "Final_PreFrail_White.boosted")
head(PreFrail.W)
save(PreFrail.W,file = "PreFrail_Final_White.AUC.Rdata")
load("PreFrail_Final_White.AUC.Rdata")
write.csv(PreFrail.W,file = "PreFrail_Final_White.AUC.csv")
##get the trained model and save##
xgb.dump(boosted, 'xgb.model.dump', with_stats = TRUE)
# get the feature real names
names = dimnames(train.matrix)[[2]]
# compute feature importance matrix
importance_matrix = xgb.importance(names, model=boosted)
# plot
Adjusted_PreFrail.White.importance = xgb.plot.importance(importance_matrix)
print(Adjusted_PreFrail.White.importance) 
write.csv(importance_matrix, file = "Adjusted_PreFrail.White.importance.csv")

#test out the ranges of the predicted classes
(predicted[W.PRE.FRAIL.full_test == 1])
range(predicted[W.PRE.FRAIL.full_test == 1])
(predicted[W.PRE.FRAIL.full_test == 0])
range(predicted[W.PRE.FRAIL.full_test == 0])
range(predicted)
#set the cutoff threshold based on the mean of the classes
mean(predicted[W.PRE.FRAIL.full_test == 1])
mean(predicted[W.PRE.FRAIL.full_test == 0])

# Set our cutoff threshold
pred.resp <- factor(ifelse(predicted >= mean(predicted[W.PRE.FRAIL.full_test == 1]), 1, 0)) # makes this set: takes exactly the mean and uses the mean as the cut off or class 1 

# Create the confusion matrix
confusionMatrix(pred.resp, factor(W.PRE.FRAIL.full_test), positive="1")
#F1, precision, recall
confusionMatrix(pred.resp, factor(W.PRE.FRAIL.full_test), mode = "prec_recall", positive="1")
levels(pred.resp)
levels(W.PRE.FRAIL.full_test)

############################################################################PRE-FRAIL TO FRAIL MODEL###################################################################################

getwd()
Frailty_full_data <- read.csv("Frailty_Full_data_frame_ARIC.csv", header = T)
head(Frailty_full_data)
Frailty_full_data %<>% mutate_if(is.factor,as.numeric)
str(Frailty_full_data)

getwd()
wbc <- read.csv("WBC-ARIC.csv", header = T)
head(wbc)
Pre.frail_Frail.wbc <- merge(Frailty_full_data,wbc, by= "ID", all.x=TRUE)
head(Pre.frail_Frail.wbc)
write.csv(Pre.frail_Frail.wbc, file="WBC.Pre.frail_Frail.Full_data_frame_ARIC.csv")

####first build the training and test dataset####
ordinal.dataframe <- read.csv("WBC.Pre.frail_Frail.Full_data_frame_ARIC.csv", header = T)
head(ordinal.dataframe)
ordinal.dataframe1 <- ordinal.dataframe[,c(1,3,5,8,10:25,29,30,31:35,39:57,59:66)]
head(ordinal.dataframe1)
ncol(ordinal.dataframe1)
nrow(ordinal.dataframe1)
ordinal.dataframe1$seed <- rep(1:10,length.out = length(ordinal.dataframe1$ID))
full_train <- subset(ordinal.dataframe1, seed <= 7)
full_test <- subset(ordinal.dataframe1, seed >7)
full_train.ids <- paste(full_train$ID, sep = " ")
full_test.ids <- paste(full_test$ID, sep = " ")
write.table(full_train.ids, file = "full_train.ids", quote = F, row.names = F, col.names = F)
write.table(full_test.ids, file = "full_test.ids", quote = F, row.names = F, col.names = F)
#write out the data.frames as test and train
write.table (full_train, file ="full_train.data", quote = F, sep = "\t", row.names = F)
write.table (full_test, file ="full_test.data", quote = F, sep = "\t", row.names = F)
#use these to create classes for confusion matrix 
sum(is.na(ordinal.dataframe$AXNPre.Frail))

#Prefrail =0 Frail =1 data.frame paremeters
pre.and.frail <- ordinal.dataframe[,c("AXNPre.Frail")]
table(pre.and.frail)

############################################################START BUILDING MODEL####################################
###load up the data####
#training set
train.load2 <- read.table("full_train.data", header = TRUE)
head(train.load2)
train.load1 <- train.load2[,c("ID","AGE","AXNPre.Frail","GLUSIU21","B6","THYR8","WBC")]
head(train.load1)
sum(is.na(train.load1$AXNPre.Frail))
train.load <- train.load1 %>% drop_na(AXNPre.Frail)
#use these to create classes for confusion matrix (only reporting out test but this can be used to run the confusion matrix on the train data if needed)
full_train.AXNPre.Frail <- train.load$AXNPre.Frail
str(train.load)
head(train.load)
dim(train.load)
train <- as.data.table(train.load)
train[, c("ID"):=NULL]
#train.outcome <- as.matrix(as.integer(train$AXNPre.Frail))
train.outcome <- as.integer(train$AXNPre.Frail)
train[, c("AXNPre.Frail"):=NULL]
train.matrix = as.matrix(train)
mode(train.matrix) = "numeric"
head(train.matrix)
#test set
test.load2 <- read.table("full_test.data", header = TRUE)
test.load1 <- test.load2[,c("ID","AGE","AXNPre.Frail","GLUSIU21","B6","THYR8","WBC")]
sum(is.na(test.load1$AXNPre.Frail))
test.load <- test.load1 %>% drop_na(AXNPre.Frail)
#use these to create classes for confusion matrix (only reporting out test but this can be used to run the confusion matrix on the train data if needed)
full_test.AXNPre.Frail <- test.load$AXNPre.Frail
head(test.load)
dim(test.load)
test <- as.data.table(test.load)
test[, c("ID"):=NULL]
#test.outcome <- as.matrix(as.integer(test$AXNPre.Frail))
test.outcome <- as.integer(test$AXNPre.Frail)
test[, c("AXNPre.Frail"):=NULL]
test.matrix = as.matrix(test)
mode(test.matrix) = "numeric"

####now run testing a variety of rounds and depths####
##testing parameters for the modeling###
pre.frail <- matrix(NA,90,3)
count <- 0
for(i in 1:10)
{
  for(j in c(5,10,20,30,40,50,75,100,200))
  {
    count <- count + 1    
    boosted <- xgboost(data = train.matrix, label = train.outcome, max.depth = i, nround = j, nthread = ncores, objective = "binary:logistic", missing = NA, eval_metric = "auc")
    predicted <- predict(boosted, test.matrix, missing = NA)
    testAuc <- auc(test.outcome, predicted)
    #print(paste("tested max.depth =",i,"and nround =",j,"resulting in a test AUC =",testAuc, sep = " "))
    pre.frail[count,1] <- i
    pre.frail[count,2] <- j
    pre.frail[count,3] <- testAuc  
  }
}
#q("no")
##95% CIs for AUC
ci.auc(testAuc)
max(pre.frail[,3])
save(boosted, file = "Adjusted_pre.Frail_boosted_Model")
#make histogram of PreFrail model and save data### 
getwd()
head(pre.frail)
save(pre.frail,file = "pre.frail_Final.AUC.Rdata")
load("pre.frail_Final.AUC.Rdata")
write.csv(pre.frail,file = "pre.frail_Final.AUC.csv")


################################################################################################################################
pre.frail <- matrix(NA,1,3)
count <- 0
for(i in 2)
{
  for(j in c(30))
  {
    count <- count + 1    
    boosted <- xgboost(data = train.matrix, label = train.outcome, max.depth = i, nround = j, nthread = ncores, objective = "binary:logistic", missing = NA, eval_metric = "auc")
    predicted <- predict(boosted, test.matrix, missing = NA)
    testAuc <- auc(test.outcome, predicted)
    #print(paste("tested max.depth =",i,"and nround =",j,"resulting in a test AUC =",testAuc, sep = " "))
    pre.frail[count,1] <- i
    pre.frail[count,2] <- j
    pre.frail[count,3] <- testAuc  
  }
}
#q("no")
##95% CIs for AUC
ci.auc(testAuc)
max(pre.frail[,3])
save(boosted, file = "Final_pre.frail_boosted_Model")
head(pre.frail)
save(pre.frail,file = "pre.frail_Final.AUC.Rdata")
write.csv(pre.frail,file = "pre.frail_Final.AUC.csv")
#saveimage
#read.table("pre.frail.COG1.AUC.txt")
##get the trained model and save##
xgb.dump(boosted, 'xgb.model.dump', with_stats = TRUE)
# get the feature real names
names = dimnames(train.matrix)[[2]]
# compute feature importance matrix
importance_matrix = xgb.importance(names, model=boosted)
# plot
Adjusted_pre.frail.model.importance = xgb.plot.importance(importance_matrix)
print(Adjusted_pre.frail.model.importance) 
write.csv(Adjusted_pre.frail.model.importance, file = "Final_AUC_Adjusted_pre.frail.model.importance.csv")

#test out the ranges of the predicted classes
(predicted[full_test.AXNPre.Frail == 1])
range(predicted[full_test.AXNPre.Frail == 1])
(predicted[full_test.AXNPre.Frail == 0])
range(predicted[full_test.AXNPre.Frail == 0])
range(predicted)
#set the cutoff threshold based on the mean of the classes
mean(predicted[full_test.AXNPre.Frail == 1])
mean(predicted[full_test.AXNPre.Frail == 0])

# Set our cutoff threshold
pred.resp <- factor(ifelse(predicted >= mean(predicted[full_test.AXNPre.Frail == 1]), 1, 0)) # makes this set: takes exactly the mean and uses the mean as the cut off or class 1 

# Create the confusion matrix
confusionMatrix(pred.resp, factor(full_test.AXNPre.Frail), positive="1")
#F1, precision, recall
confusionMatrix(pred.resp, factor(full_test.AXNPre.Frail), mode = "prec_recall", positive="1")
levels(pred.resp)
levels(full_test.AXNPre.Frail)

#####################################################################################RACE MODELS FOR PRE.FRAIL TO FRAIL #############################################################

ordinal.dataframe <- read.csv("B.WBC.Pre.frail_Frail.Full_data_frame_ARIC.csv", header = T)
head(ordinal.dataframe)
ordinal.dataframe1 <- ordinal.dataframe[,c(1,3,5,8,10:25,29,30,31:35,39:57,59:66)]
head(ordinal.dataframe1)
ncol(ordinal.dataframe1)
nrow(ordinal.dataframe1)
ordinal.dataframe1$seed <- rep(1:10,length.out = length(ordinal.dataframe1$ID))
full_train <- subset(ordinal.dataframe1, seed <= 7)
full_test <- subset(ordinal.dataframe1, seed >7)
full_train.ids <- paste(full_train$ID, sep = " ")
full_test.ids <- paste(full_test$ID, sep = " ")
write.table(full_train.ids, file = "B.full_train.ids", quote = F, row.names = F, col.names = F)
write.table(full_test.ids, file = "B.full_test.ids", quote = F, row.names = F, col.names = F)
#write out the data.frames as test and train
write.table (full_train, file ="B.full_train.data", quote = F, sep = "\t", row.names = F)
write.table (full_test, file ="B.full_test.data", quote = F, sep = "\t", row.names = F)
#use these to create classes for confusion matrix 
sum(is.na(ordinal.dataframe$AXNPre.Frail))

#Prefrail =0 Frail =1 data.frame paremeters
pre.and.frail <- ordinal.dataframe[,c("AXNPre.Frail")]
table(pre.and.frail)
#pre.and.frail

############################################################START BUILDING MODEL####################################
###load up the data####
#training set
train.load2 <- read.table("B.full_train.data", header = TRUE)
head(train.load2)
train.load1 <- train.load2[,c("ID","AGE","AXNPre.Frail","GLUSIU21","B6","THYR8","WBC")]
head(train.load1)
sum(is.na(train.load1$AXNPre.Frail))
train.load <- train.load1 %>% drop_na(AXNPre.Frail)
#use these to create classes for confusion matrix (only reporting out test but this can be used to run the confusion matrix on the train data if needed)
full_train.AXNPre.Frail <- train.load$AXNPre.Frail
str(train.load)
head(train.load)
dim(train.load)
train <- as.data.table(train.load)
train[, c("ID"):=NULL]
#train.outcome <- as.matrix(as.integer(train$AXNPre.Frail))
train.outcome <- as.integer(train$AXNPre.Frail)
train[, c("AXNPre.Frail"):=NULL]
train.matrix = as.matrix(train)
mode(train.matrix) = "numeric"
head(train.matrix)
#test set
test.load2 <- read.table("B.full_test.data", header = TRUE)
test.load1 <- test.load2[,c("ID","AGE","AXNPre.Frail","GLUSIU21","B6","THYR8","WBC")]
sum(is.na(test.load1$AXNPre.Frail))
test.load <- test.load1 %>% drop_na(AXNPre.Frail)
#use these to create classes for confusion matrix (only reporting out test but this can be used to run the confusion matrix on the train data if needed)
full_test.AXNPre.Frail <- test.load$AXNPre.Frail
head(test.load)
dim(test.load)
test <- as.data.table(test.load)
test[, c("ID"):=NULL]
#test.outcome <- as.matrix(as.integer(test$AXNPre.Frail))
test.outcome <- as.integer(test$AXNPre.Frail)
test[, c("AXNPre.Frail"):=NULL]
test.matrix = as.matrix(test)
mode(test.matrix) = "numeric"

####now run testing a variety of rounds and depths####
##testing parameters for the modeling###
pre.frail <- matrix(NA,90,3)
count <- 0
for(i in 1:10)
{
  for(j in c(5,10,20,30,40,50,75,100,200))
  {
    count <- count + 1    
    boosted <- xgboost(data = train.matrix, label = train.outcome, max.depth = i, nround = j, nthread = ncores, objective = "binary:logistic", missing = NA, eval_metric = "auc")
    predicted <- predict(boosted, test.matrix, missing = NA)
    testAuc <- auc(test.outcome, predicted)
    #print(paste("tested max.depth =",i,"and nround =",j,"resulting in a test AUC =",testAuc, sep = " "))
    pre.frail[count,1] <- i
    pre.frail[count,2] <- j
    pre.frail[count,3] <- testAuc  
  }
}
#q("no")
##95% CIs for AUC
ci.auc(testAuc)
max(pre.frail[,3])
save(boosted, file = "Adjusted_pre.Frail_boosted_Model")
#make histogram of PreFrail model and save data### 
getwd()
head(pre.frail)
save(pre.frail,file = "B.pre.frail_Final.AUC.Rdata")
load("B.pre.frail_Final.AUC.Rdata")
write.csv(pre.frail,file = "B.pre.frail_Final.AUC.csv")

################################################################################################################    
pre.frail <- matrix(NA,1,3)
count <- 0
for(i in 4)
{
  for(j in c(75))
  {
    count <- count + 1    
    boosted <- xgboost(data = train.matrix, label = train.outcome, max.depth = i, nround = j, nthread = ncores, objective = "binary:logistic", missing = NA, eval_metric = "auc")
    predicted <- predict(boosted, test.matrix, missing = NA)
    testAuc <- auc(test.outcome, predicted)
    #print(paste("tested max.depth =",i,"and nround =",j,"resulting in a test AUC =",testAuc, sep = " "))
    pre.frail[count,1] <- i
    pre.frail[count,2] <- j
    pre.frail[count,3] <- testAuc  
  }
}
#q("no")
##95% CIs for AUC
ci.auc(testAuc)
max(pre.frail[,3])
save(boosted, file = "Final_pre.frail_boosted_Model")
head(pre.frail)
save(pre.frail,file = "pre.frail_Final.AUC.Rdata")
write.csv(pre.frail,file = "pre.frail_Final.AUC.csv")
#saveimage
#read.table("pre.frail.COG1.AUC.txt")
##get the trained model and save##
xgb.dump(boosted, 'xgb.model.dump', with_stats = TRUE)
# get the feature real names
names = dimnames(train.matrix)[[2]]
# compute feature importance matrix
importance_matrix = xgb.importance(names, model=boosted)
# plot
Adjusted_pre.frail.model.importance = xgb.plot.importance(importance_matrix)
print(Adjusted_pre.frail.model.importance) 
write.csv(Adjusted_pre.frail.model.importance, file = "Final_AUC_Adjusted_pre.frail.model.importance.csv")

#test out the ranges of the predicted classes
(predicted[full_test.AXNPre.Frail == 1])
range(predicted[full_test.AXNPre.Frail == 1])
(predicted[full_test.AXNPre.Frail == 0])
range(predicted[full_test.AXNPre.Frail == 0])
range(predicted)
#set the cutoff threshold based on the mean of the classes
mean(predicted[full_test.AXNPre.Frail == 1])
mean(predicted[full_test.AXNPre.Frail == 0])

# Set our cutoff threshold
pred.resp <- factor(ifelse(predicted >= mean(predicted[full_test.AXNPre.Frail == 1]), 1, 0)) # makes this set: takes exactly the mean and uses the mean as the cut off or class 1 

# Create the confusion matrix
confusionMatrix(pred.resp, factor(full_test.AXNPre.Frail), positive="1")
#F1, precision, recall
confusionMatrix(pred.resp, factor(full_test.AXNPre.Frail), mode = "prec_recall", positive="1")
levels(pred.resp)
levels(full_test.AXNPre.Frail)


#################################################################RACE WHITE PRE.FRAIL TO FRAIL MODEL###########################################################


####first build the training and test datasets####
ordinal.dataframe <- read.csv("W.WBC.Pre.frail_Frail.Full_data_frame_ARIC.csv", header = T)
head(ordinal.dataframe)
ordinal.dataframe1 <- ordinal.dataframe[,c(1,3,5,8,10:25,29,30,31:35,39:57,59:66)]
head(ordinal.dataframe1)
ncol(ordinal.dataframe1)
nrow(ordinal.dataframe1)
ordinal.dataframe1$seed <- rep(1:10,length.out = length(ordinal.dataframe1$ID))
full_train <- subset(ordinal.dataframe1, seed <= 7)
full_test <- subset(ordinal.dataframe1, seed >7)
full_train.ids <- paste(full_train$ID, sep = " ")
full_test.ids <- paste(full_test$ID, sep = " ")
write.table(full_train.ids, file = "W.full_train.ids", quote = F, row.names = F, col.names = F)
write.table(full_test.ids, file = "W.full_test.ids", quote = F, row.names = F, col.names = F)
#write out the data.frames as test and train
write.table (full_train, file ="W.full_train.data", quote = F, sep = "\t", row.names = F)
write.table (full_test, file ="W.full_test.data", quote = F, sep = "\t", row.names = F)
#use these to create classes for confusion matrix 
sum(is.na(ordinal.dataframe$AXNPre.Frail))

#Prefrail =0 Frail =1 data.frame paremeters
pre.and.frail <- ordinal.dataframe[,c("AXNPre.Frail")]
table(pre.and.frail)

############################################################START BUILDING MODEL####################################
###load up the data####
#training set
train.load2 <- read.table("W.full_train.data", header = TRUE)
head(train.load2)
train.load1 <- train.load2[,c("ID","AGE","AXNPre.Frail","GLUSIU21","B6","THYR8","WBC")]
head(train.load1)
table(train.load1$AXNPre.Frail)
sum(is.na(train.load1$AXNPre.Frail))
train.load <- train.load1 %>% drop_na(AXNPre.Frail)
#use these to create classes for confusion matrix (only reporting out test but this can be used to run the confusion matrix on the train data if needed)
full_train.AXNPre.Frail <- train.load3$AXNPre.Frail
str(train.load)
head(train.load)
dim(train.load)
train <- as.data.table(train.load)
train[, c("ID"):=NULL]
#train.outcome <- as.matrix(as.integer(train$AXNPre.Frail))
train.outcome <- as.integer(train$AXNPre.Frail)
train[, c("AXNPre.Frail"):=NULL]
train.matrix = as.matrix(train)
mode(train.matrix) = "numeric"
head(train.matrix)
#test set
test.load2 <- read.table("W.full_test.data", header = TRUE)
test.load1 <- test.load2[,c("ID","AGE","AXNPre.Frail","GLUSIU21","B6","THYR8","WBC")]
sum(is.na(test.load1$AXNPre.Frail))
test.load <- test.load1 %>% drop_na(AXNPre.Frail)
#use these to create classes for confusion matrix (only reporting out test but this can be used to run the confusion matrix on the train data if needed)
full_test.AXNPre.Frail <- test.load$AXNPre.Frail
head(test.load)
dim(test.load)
test <- as.data.table(test.load)
test[, c("ID"):=NULL]
#test.outcome <- as.matrix(as.integer(test$AXNPre.Frail))
test.outcome <- as.integer(test$AXNPre.Frail)
test[, c("AXNPre.Frail"):=NULL]
test.matrix = as.matrix(test)
mode(test.matrix) = "numeric"

####now run testing a variety of rounds and depths####
##testing parameters for the modeling###
pre.frail <- matrix(NA,90,3)
count <- 0
for(i in 1:10)
{
  for(j in c(5,10,20,30,40,50,75,100,200))
  {
    count <- count + 1    
    boosted <- xgboost(data = train.matrix, label = train.outcome, max.depth = i, nround = j, nthread = ncores, objective = "binary:logistic", missing = NA, eval_metric = "auc")
    predicted <- predict(boosted, test.matrix, missing = NA)
    testAuc <- auc(test.outcome, predicted)
    #print(paste("tested max.depth =",i,"and nround =",j,"resulting in a test AUC =",testAuc, sep = " "))
    pre.frail[count,1] <- i
    pre.frail[count,2] <- j
    pre.frail[count,3] <- testAuc  
  }
}
#q("no")
##95% CIs for AUC
ci.auc(testAuc)
max(pre.frail[,3])
save(boosted, file = "Adjusted_pre.Frail_boosted_Model")
#make histogram of PreFrail model and save data### 
getwd()
head(pre.frail)
save(pre.frail,file = "W.pre.frail_Final.AUC.Rdata")
load("W.pre.frail_Final.AUC.Rdata")
write.csv(pre.frail,file = "W.pre.frail_Final.AUC.csv")

#################################################################################################################    
pre.frail <- matrix(NA,1,3)
count <- 0
for(i in 6)
{
  for(j in c(5))
  {
    count <- count + 1    
    boosted <- xgboost(data = train.matrix, label = train.outcome, max.depth = i, nround = j, nthread = ncores, objective = "binary:logistic", missing = NA, eval_metric = "auc")
    predicted <- predict(boosted, test.matrix, missing = NA)
    testAuc <- auc(test.outcome, predicted)
    #print(paste("tested max.depth =",i,"and nround =",j,"resulting in a test AUC =",testAuc, sep = " "))
    pre.frail[count,1] <- i
    pre.frail[count,2] <- j
    pre.frail[count,3] <- testAuc  
  }
}
#q("no")
##95% CIs for AUC
ci.auc(testAuc)
max(pre.frail[,3])
save(boosted, file = "Final_pre.frail_boosted_Model")
head(pre.frail)
save(pre.frail,file = "pre.frail_Final.AUC.Rdata")
write.csv(pre.frail,file = "pre.frail_Final.AUC.csv")
#saveimage
#read.table("pre.frail.COG1.AUC.txt")
##get the trained model and save##
xgb.dump(boosted, 'xgb.model.dump', with_stats = TRUE)
# get the feature real names
names = dimnames(train.matrix)[[2]]
# compute feature importance matrix
importance_matrix = xgb.importance(names, model=boosted)
# plot
Adjusted_pre.frail.model.importance = xgb.plot.importance(importance_matrix)
print(Adjusted_pre.frail.model.importance) 
write.csv(Adjusted_pre.frail.model.importance, file = "W_pre.frail.to.frailmodel.importance.csv")

#test out the ranges of the predicted classes
(predicted[full_test.AXNPre.Frail == 1])
range(predicted[full_test.AXNPre.Frail == 1])
(predicted[full_test.AXNPre.Frail == 0])
range(predicted[full_test.AXNPre.Frail == 0])
range(predicted)
#set the cutoff threshold based on the mean of the classes
mean(predicted[full_test.AXNPre.Frail == 1])
mean(predicted[full_test.AXNPre.Frail == 0])

# Set our cutoff threshold
pred.resp <- factor(ifelse(predicted >= mean(predicted[full_test.AXNPre.Frail == 1]), 1, 0)) # makes this set: takes exactly the mean and uses the mean as the cut off or class 1 

# Create the confusion matrix
confusionMatrix(pred.resp, factor(full_test.AXNPre.Frail), positive="1")
#F1, precision, recall
confusionMatrix(pred.resp, factor(full_test.AXNPre.Frail), mode = "prec_recall", positive="1")
levels(pred.resp)
levels(full_test.AXNPre.Frail)



