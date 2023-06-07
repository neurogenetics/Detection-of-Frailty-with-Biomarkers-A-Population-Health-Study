#####################Ordinal model build in xgboost
#Robust and prefrail
#Robust and frail
#prefrail and robust

ordinal.dataframe <- read.csv("V1_All_Phenotypes_Recode.Final.Data.Frame.csv", header = T)
head(ordinal.dataframe)
ordinal.dataframe$seed <- rep(1:10,length.out = length(ordinal.dataframe$IID))
train <- subset(ordinal.dataframe, seed <= 7)
test <- subset(ordinal.dataframe, seed >7)
train.ids <- paste(train$IID, sep = " ")
test.ids <- paste(test$IID, sep = " ")
write.table(train.ids, file = "train.ids", quote = F, row.names = F, col.names = F)
write.table(test.ids, file = "test.ids", quote = F, row.names = F, col.names = F)
head(test)
#write out the data.frames as test and train
write.table (train, file ="train.data", quote = F, sep = "\t", row.names = F)
write.table (test, file ="test.data", quote = F, sep = "\t", row.names = F)
#use these to create classes for confusion matrix 
sum(is.na(ordinal.dataframe$AXNPre.Frail))

#PreFrail data.frame parameters#
PreFrail <- ordinal.dataframe[,c("AXNPREFRAIL")]
table(PreFrail)

#Frail data.frame parameters#
Frail <- ordinal.dataframe[,c("AXNFRAIL")]
table(Frail)

#Prefrail =0 Frail =1 data.frame paremeters
pre.and.frail <- ordinal.dataframe[,c("AXNPre.Frail")]
table(pre.and.frail)
pre.and.frail

##START HERE###
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

#ncores <- detectCores() - 4
ncores <- 12



######################################################################################Pre-frail and frail Model Boosting##########################################################################################
####load up the data####
#training set
train.load2 <- read.table("train.data", header = TRUE)
head(train.load2)
train.load1 <- train.load2[,c(1,3:298)]
head(train.load1)
sum(is.na(train.load1$AXNPre.Frail))
train.load <- train.load1 %>% drop_na(AXNPre.Frail)
#use these to create classes for confusion matrix (only reporting out test but this can be used to run the confusion matrix on the train data if needed)
full_train.AXNPre.Frail <- train.load$AXNPre.Frail
str(train.load)
head(train.load)
dim(train.load)
train <- as.data.table(train.load)
train[, c("IID"):=NULL]
#train.outcome <- as.matrix(as.integer(train$AXNPre.Frail))
train.outcome <- as.integer(train$AXNPre.Frail)
train[, c("AXNPre.Frail"):=NULL]
train.matrix = as.matrix(train)
mode(train.matrix) = "numeric"
head(train.matrix)
#test set
test.load2 <- read.table("test.data", header = TRUE)
test.load1 <- test.load2[,c(1,3:298)]
sum(is.na(test.load1$AXNPre.Frail))
test.load <- test.load1 %>% drop_na(AXNPre.Frail)
#use these to create classes for confusion matrix (only reporting out test but this can be used to run the confusion matrix on the train data if needed)
full_test.AXNPre.Frail <- test.load$AXNPre.Frail
head(test.load)
dim(test.load)
test <- as.data.table(test.load)
test[, c("IID"):=NULL]
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
for(i in 1)
{
  for(j in c(10))
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
write.csv(Adjusted_pre.frail.model.importance, file = "Final_pre.frail.model.importance.csv")

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


######################################################Feature selection based on significance###########################################################

#training set
train.load2 <- read.table("train.data", header = TRUE)
head(train.load2)
train.load1 <- train.load2[,c("IID","AXNPre.Frail","X_GB","IXAGE","X_GLU","X_UCRE24","X_FT4","X_IL1RA")]
head(train.load1)
sum(is.na(train.load1$AXNPre.Frail))
train.load <- train.load1 %>% drop_na(AXNPre.Frail)
#use these to create classes for confusion matrix (only reporting out test but this can be used to run the confusion matrix on the train data if needed)
full_train.AXNPre.Frail <- train.load$AXNPre.Frail
####Rename columns for final model graphs and tables
names(train.load) <- c("IID","Pre.Frail","WBC","AGE","Glucose","24-hour Urine Creatinine","Free T4","IL1RA")
str(train.load)
head(train.load)
dim(train.load)
train <- as.data.table(train.load)
train[, c("IID"):=NULL]
train.outcome <- as.integer(train$Pre.Frail)
train[, c("Pre.Frail"):=NULL]
train.matrix = as.matrix(train)
mode(train.matrix) = "numeric"
head(train.matrix)
#test set
test.load2 <- read.table("test.data", header = TRUE)
test.load1 <- test.load2[,c("IID","AXNPre.Frail","X_GB","IXAGE","X_GLU","X_UCRE24","X_FT4","X_IL1RA")]
sum(is.na(test.load1$AXNPre.Frail))
test.load <- test.load1 %>% drop_na(AXNPre.Frail)
#use these to create classes for confusion matrix (only reporting out test but this can be used to run the confusion matrix on the train data if needed)
full_test.AXNPre.Frail <- test.load$AXNPre.Frail
####Rename columns for final model graphs and tables
names(test.load) <- c("IID","Pre.Frail","WBC","AGE","Glucose","24-hour Urine Creatinine","Free T4","IL1RA")
head(test.load)
dim(test.load)
test <- as.data.table(test.load)
test[, c("IID"):=NULL]
test.outcome <- as.integer(test$Pre.Frail)
test[, c("Pre.Frail"):=NULL]
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
    boosted <- xgboost(data = train.matrix, label = train.outcome, max.depth = i, nround = j, nthread = ncores, objective = "binary:logistic", missing = NA, eval_metric = "auc", scale_pos_weight=)
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

pre.frail <- matrix(NA,1,3)
count <- 0
for(i in 1)
{
  for(j in c(10))
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
write.csv(Adjusted_pre.frail.model.importance, file = "Final_pre.frail.model.importance.csv")

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
