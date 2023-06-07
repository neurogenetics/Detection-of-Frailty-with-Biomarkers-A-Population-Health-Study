##Multinomial Classification ARIC
###6.4.23

##https://medium.com/@imoisharma/gentle-introduction-of-xgboost-library-2b1ac2669680
####even better https://rstudio-pubs-static.s3.amazonaws.com/233321_06bcdf2c8bc445dbb635740bb44f978b.html
###code used from https://rpubs.com/mharris/multiclass_xgboost
###balance data package = ROSE https://www.rdocumentation.org/packages/ROSE/versions/0.0-3/topics/ovun.sample

require(parallel)
require(data.table)
require(corrplot)
require(Rtsne)
require(stats)
require(knitr)
require(ggplot2)
require(pROC)
library(tidyr)


#multiclass packages
library("xgboost")  # the main algorithm
library("caret")    # for the confusionmatrix() function (also needs e1071 package)
library("dplyr")    # for some data preperation
library("Ckmeans.1d.dp") # for xgb.ggplot.importance
library(ROSE)
library(nnet)

Multinominal1 <- read.csv("Multinomial.WBC.Full_data_frame_ARIC.csv", header = T)
head(Multinominal1)
str(Multinominal1)
dim(Multinominal1)
Multinominal <- Multinominal1[,c(1,3,5,8,10:25,27:36,40:58,60:67)]
head(Multinominal)
dim(Multinominal)
table(Multinominal$AXNMULTIFRAIL)
class(Multinominal)
Multi2 <- Multinominal %>% drop_na(AXNMULTIFRAIL)
dim(Multi2)

set.seed(12345)

inTrain = createDataPartition(Multi4$AXNMULTIFRAIL, p=2/3, list = FALSE)

#Train
train.data=Multi4[inTrain,]
dim(train.data)
table(train.data$AXNMULTIFRAIL)
#split the train data with class 0 and 2 = dataframe1
#split data 0 and 2
dataframe1 <- train.data[train.data$AXNMULTIFRAIL %in% c(0,2),]
table(dataframe1$AXNMULTIFRAIL)
#write.csv(dataframe1, file="viewdataframe1.csv")
# then dataframe with class 1 = dataframe2

#now balance the data

# take dataframe1 and balance using the overmethod = train.balanced
train.balanced1 <- ovun.sample(AXNMULTIFRAIL ~ ., data=dataframe1, N=500,seed=1,method="both")$data
table(train.balanced1$AXNMULTIFRAIL)
#now balance the data 
dataframe2 <- train.data[train.data$AXNMULTIFRAIL %in% c(0,1),]
table(dataframe2$AXNMULTIFRAIL)

train.balanced2 <- ovun.sample(AXNMULTIFRAIL ~ ., data=dataframe2, N=500,seed=2,method="both")$data
table(train.balanced2$AXNMULTIFRAIL)

#take 1 out of train.balance2
dataframe3 <- train.balanced2[train.balanced2$AXNMULTIFRAIL %in% c(1),]
# now combine together
final.train.balanced <- merge(dataframe3,train.balanced1, all=TRUE)
head(final.train.balanced)
dim(final.train.balanced)
table(final.train.balanced$AXNMULTIFRAIL)
#classes after balance train data

# then do the rest of the train and test 
train.label=final.train.balanced[,"AXNMULTIFRAIL"]
#remove class col
final.train.balanced = subset(final.train.balanced,select = -c(ID,AXNMULTIFRAIL))
head(final.train.balanced)
dim(as.matrix(final.train.balanced))
#create xgboost matrix
train.matrix=xgb.DMatrix(data=as.matrix(final.train.balanced), label=train.label)

#Test
test.data=Multi4[-inTrain,]
test.label=test.data[,"AXNMULTIFRAIL"]
#remove class col
test.data = subset(test.data,select = -c(ID,AXNMULTIFRAIL))
head(test.data)
dim(as.matrix(test.data))
#create xgboost matrix
test.matrix=xgb.DMatrix(data=as.matrix(test.data), label=test.label)


#run multiclass model
numberOfClasses <- length(unique(Multi4$AXNMULTIFRAIL))
xgb_params <- list("objective" = "multi:softprob",
                   "eval_metric" = "mlogloss",
                   "num_class" = numberOfClasses)
nround    <- 100 # number of XGBoost rounds
cv.nfold  <- 5

# Fit cv.nfold * cv.nround XGB models and save OOF predictions
multiclass.frail <- xgb.cv(params = xgb_params,
                           data = train.matrix, 
                           nrounds = nround,
                           nfold = cv.nfold,
                           verbose = FALSE,
                           prediction = TRUE)

OOF_prediction <- data.frame(multiclass.frail$pred) %>%
  mutate(max_prob = max.col(., ties.method = "last"),
         label = train.label + 1)

head(OOF_prediction)


conf_frail <- confusionMatrix(factor(OOF_prediction$max_prob),
                              factor(OOF_prediction$label),
                              mode = "everything")
print(conf_frail)
##https://stackoverflow.com/questions/34842837/saving-output-of-confusionmatrix-as-a-csv-table
#write out confussionMatrix results
t <- as.table(conf_frail)
o <- as.matrix(conf_frail,what="overall")
c <- as.matrix(conf_frail, what = "classes")
write.csv(c, file="Multiclass.classes.stats.csv")

library(pROC)
library(plyr)
#do this process with train and test.data 
dim(Multi4)

#Pheno <- revalue(Multi4$AXNMULTIFRAIL,c(0="X1",1="X2",2="X3"))
Pheno <- train.label
Pheno[Pheno==0] <- "X1"
Pheno[Pheno==1] <- "X2"
Pheno[Pheno==2] <- "X3"
Pheno=as.factor(Pheno)
length(Pheno)
dim(OOF_prediction)
auc <- multiclass.roc(Pheno,OOF_prediction[,c("X1","X2","X3")], levels = c("X1","X2","X3")) 
print(auc)
#plot(auc, ylim=c(0,1), print.thres=TRUE, main=paste('AUC:',round(auc$auc[[1]],2)))
#abline(h=1,col='blue',lwd=2)
#abline(h=0,col='red',lwd=2)

#Train Full Model and Assess Test Set Error on Test Data 
bst_model <- xgb.train(params = xgb_params,
                       data = train.matrix,
                       nrounds = nround)

# Predict hold-out test set
test_pred <- predict(bst_model, newdata = test.matrix)
test_prediction <- matrix(test_pred, nrow = numberOfClasses,
                          ncol=length(test_pred)/numberOfClasses) %>%
  t() %>%
  data.frame() %>%
  mutate(label = test.label + 1,
         max_prob = max.col(., "last"))
# confusion matrix of test set
confusionMatrix(factor(test_prediction$max_prob),
                factor(test_prediction$label),
                mode = "everything")



# get the feature real names
names = dimnames(train.matrix)[[1]]
# compute feature importance matrix
importance_matrix = xgb.importance(feature_names = names, model = bst_model)
head(importance_matrix)
write.csv(importance_matrix, file = "final.importance.multi.allvars.csv")

# plot
gp = xgb.ggplot.importance(importance_matrix)
print(gp) 

#################################################################################Race Multinomial models start here#################################################################

##############Race White
Multinominal1 <- read.csv("W.Multinomial.WBC.Full_data_frame_ARIC.csv", header = T)
head(Multinominal1)
str(Multinominal1)
dim(Multinominal1)
Multinominal <- Multinominal1[,c(1,3,5,8,10:25,27:36,40:58,60:67)]
head(Multinominal)
dim(Multinominal)
table(Multinominal$AXNMULTIFRAIL)
class(Multinominal)
Multi2 <- Multinominal %>% drop_na(AXNMULTIFRAIL)
dim(Multi2)

set.seed(12345)

inTrain = createDataPartition(Multi4$AXNMULTIFRAIL, p=2/3, list = FALSE)

#Train
train.data=Multi4[inTrain,]
dim(train.data)
table(train.data$AXNMULTIFRAIL)
#split the train data with class 0 and 2 = dataframe1
#split data 0 and 2
dataframe1 <- train.data[train.data$AXNMULTIFRAIL %in% c(0,2),]
table(dataframe1$AXNMULTIFRAIL)
#write.csv(dataframe1, file="viewdataframe1.csv")
# then dataframe with class 1 = dataframe2

#now balance the data

# take dataframe1 and balance using the overmethod = train.balanced
train.balanced1 <- ovun.sample(AXNMULTIFRAIL ~ ., data=dataframe1, N=500,seed=1,method="both")$data
table(train.balanced1$AXNMULTIFRAIL)
#now balance the data 
dataframe2 <- train.data[train.data$AXNMULTIFRAIL %in% c(0,1),]
table(dataframe2$AXNMULTIFRAIL)

train.balanced2 <- ovun.sample(AXNMULTIFRAIL ~ ., data=dataframe2, N=500,seed=2,method="both")$data
table(train.balanced2$AXNMULTIFRAIL)

#take 1 out of train.balance2
dataframe3 <- train.balanced2[train.balanced2$AXNMULTIFRAIL %in% c(1),]
# now combine together
final.train.balanced <- merge(dataframe3,train.balanced1, all=TRUE)
head(final.train.balanced)
dim(final.train.balanced)
table(final.train.balanced$AXNMULTIFRAIL)
#classes after balance train data

# then do the rest of the train and test 
train.label=final.train.balanced[,"AXNMULTIFRAIL"]
#remove class col
final.train.balanced = subset(final.train.balanced,select = -c(ID,AXNMULTIFRAIL))
head(final.train.balanced)
dim(as.matrix(final.train.balanced))
#create xgboost matrix
train.matrix=xgb.DMatrix(data=as.matrix(final.train.balanced), label=train.label)

#Test
test.data=Multi4[-inTrain,]
test.label=test.data[,"AXNMULTIFRAIL"]
#remove class col
test.data = subset(test.data,select = -c(ID,AXNMULTIFRAIL))
head(test.data)
dim(as.matrix(test.data))
#create xgboost matrix
test.matrix=xgb.DMatrix(data=as.matrix(test.data), label=test.label)

#run multiclass model
numberOfClasses <- length(unique(Multi4$AXNMULTIFRAIL))
xgb_params <- list("objective" = "multi:softprob",
                   "eval_metric" = "mlogloss",
                   "num_class" = numberOfClasses)
nround    <- 100 # number of XGBoost rounds
cv.nfold  <- 5

# Fit cv.nfold * cv.nround XGB models and save OOF predictions
multiclass.frail <- xgb.cv(params = xgb_params,
                           data = train.matrix, 
                           nrounds = nround,
                           nfold = cv.nfold,
                           verbose = FALSE,
                           prediction = TRUE)

OOF_prediction <- data.frame(multiclass.frail$pred) %>%
  mutate(max_prob = max.col(., ties.method = "last"),
         label = train.label + 1)

head(OOF_prediction)


conf_frail <- confusionMatrix(factor(OOF_prediction$max_prob),
                              factor(OOF_prediction$label),
                              mode = "everything")
print(conf_frail)
##https://stackoverflow.com/questions/34842837/saving-output-of-confusionmatrix-as-a-csv-table
#write out confussionMatrix results
t <- as.table(conf_frail)
o <- as.matrix(conf_frail,what="overall")
c <- as.matrix(conf_frail, what = "classes")
write.csv(c, file="Multiclass.classes.stats.csv")

library(pROC)
library(plyr)
#do this process with train and test.data 
dim(Multi4)

#Pheno <- revalue(Multi4$AXNMULTIFRAIL,c(0="X1",1="X2",2="X3"))
Pheno <- train.label
Pheno[Pheno==0] <- "X1"
Pheno[Pheno==1] <- "X2"
Pheno[Pheno==2] <- "X3"
Pheno=as.factor(Pheno)
length(Pheno)
dim(OOF_prediction)
auc <- multiclass.roc(Pheno,OOF_prediction[,c("X1","X2","X3")], levels = c("X1","X2","X3")) 
print(auc)
#plot(auc, ylim=c(0,1), print.thres=TRUE, main=paste('AUC:',round(auc$auc[[1]],2)))
#abline(h=1,col='blue',lwd=2)
#abline(h=0,col='red',lwd=2)

#Train Full Model and Assess Test Set Error on Test Data 
bst_model <- xgb.train(params = xgb_params,
                       data = train.matrix,
                       nrounds = nround)

# Predict hold-out test set
test_pred <- predict(bst_model, newdata = test.matrix)
test_prediction <- matrix(test_pred, nrow = numberOfClasses,
                          ncol=length(test_pred)/numberOfClasses) %>%
  t() %>%
  data.frame() %>%
  mutate(label = test.label + 1,
         max_prob = max.col(., "last"))
# confusion matrix of test set
confusionMatrix(factor(test_prediction$max_prob),
                factor(test_prediction$label),
                mode = "everything")

# get the feature real names
names = dimnames(train.matrix)[[1]]
# compute feature importance matrix
importance_matrix = xgb.importance(feature_names = names, model = bst_model)
head(importance_matrix)
write.csv(importance_matrix, file = "W.final.importance.multi.allvars.csv")

# plot
gp = xgb.ggplot.importance(importance_matrix)
print(gp) 


############################################################################RACE Multinominal Model Black#######################################################################33

##############Race Black
Multinominal1 <- read.csv("B.Multinomial.WBC.Full_data_frame_ARIC.csv", header = T)
head(Multinominal1)
str(Multinominal1)
dim(Multinominal1)
Multinominal <- Multinominal1[,c(1,3,5,8,10:25,27:36,40:58,60:67)]
head(Multinominal)
dim(Multinominal)
table(Multinominal$AXNMULTIFRAIL)
class(Multinominal)
Multi2 <- Multinominal %>% drop_na(AXNMULTIFRAIL)
dim(Multi2)

set.seed(12345)

inTrain = createDataPartition(Multi4$AXNMULTIFRAIL, p=2/3, list = FALSE)

#Train
train.data=Multi4[inTrain,]
dim(train.data)
table(train.data$AXNMULTIFRAIL)
#split the train data with class 0 and 2 = dataframe1
#split data 0 and 2
dataframe1 <- train.data[train.data$AXNMULTIFRAIL %in% c(0,2),]
table(dataframe1$AXNMULTIFRAIL)
#write.csv(dataframe1, file="viewdataframe1.csv")
# then dataframe with class 1 = dataframe2

#now balance the data

# take dataframe1 and balance using the overmethod = train.balanced
train.balanced1 <- ovun.sample(AXNMULTIFRAIL ~ ., data=dataframe1, N=160,seed=1,method="both")$data
table(train.balanced1$AXNMULTIFRAIL)
#now balance the data 
dataframe2 <- train.data[train.data$AXNMULTIFRAIL %in% c(0,1),]
table(dataframe2$AXNMULTIFRAIL)

train.balanced2 <- ovun.sample(AXNMULTIFRAIL ~ ., data=dataframe2, N=160,seed=2,method="both")$data
table(train.balanced2$AXNMULTIFRAIL)

#take 1 out of train.balance2
dataframe3 <- train.balanced2[train.balanced2$AXNMULTIFRAIL %in% c(1),]
# now combine together
final.train.balanced <- merge(dataframe3,train.balanced1, all=TRUE)
head(final.train.balanced)
dim(final.train.balanced)
table(final.train.balanced$AXNMULTIFRAIL)
#classes after balance train data

# then do the rest of the train and test 
train.label=final.train.balanced[,"AXNMULTIFRAIL"]
#remove class col
final.train.balanced = subset(final.train.balanced,select = -c(ID,AXNMULTIFRAIL))
head(final.train.balanced)
dim(as.matrix(final.train.balanced))
#create xgboost matrix
train.matrix=xgb.DMatrix(data=as.matrix(final.train.balanced), label=train.label)

#Test
test.data=Multi4[-inTrain,]
test.label=test.data[,"AXNMULTIFRAIL"]
#remove class col
test.data = subset(test.data,select = -c(ID,AXNMULTIFRAIL))
head(test.data)
dim(as.matrix(test.data))
#create xgboost matrix
test.matrix=xgb.DMatrix(data=as.matrix(test.data), label=test.label)


#run multiclass model
numberOfClasses <- length(unique(Multi4$AXNMULTIFRAIL))
xgb_params <- list("objective" = "multi:softprob",
                   "eval_metric" = "mlogloss",
                   "num_class" = numberOfClasses)
nround    <- 100 # number of XGBoost rounds
cv.nfold  <- 5

# Fit cv.nfold * cv.nround XGB models and save OOF predictions
multiclass.frail <- xgb.cv(params = xgb_params,
                           data = train.matrix, 
                           nrounds = nround,
                           nfold = cv.nfold,
                           verbose = FALSE,
                           prediction = TRUE)

OOF_prediction <- data.frame(multiclass.frail$pred) %>%
  mutate(max_prob = max.col(., ties.method = "last"),
         label = train.label + 1)

head(OOF_prediction)


conf_frail <- confusionMatrix(factor(OOF_prediction$max_prob),
                              factor(OOF_prediction$label),
                              mode = "everything")
print(conf_frail)
##https://stackoverflow.com/questions/34842837/saving-output-of-confusionmatrix-as-a-csv-table
#write out confussionMatrix results
t <- as.table(conf_frail)
o <- as.matrix(conf_frail,what="overall")
c <- as.matrix(conf_frail, what = "classes")
write.csv(c, file="Multiclass.classes.stats.csv")

library(pROC)
library(plyr)
#do this process with train and test.data 

#Pheno <- revalue(Multi4$AXNMULTIFRAIL,c(0="X1",1="X2",2="X3"))
Pheno <- train.label
Pheno[Pheno==0] <- "X1"
Pheno[Pheno==1] <- "X2"
Pheno[Pheno==2] <- "X3"
Pheno=as.factor(Pheno)
length(Pheno)
dim(OOF_prediction)
auc <- multiclass.roc(Pheno,OOF_prediction[,c("X1","X2","X3")], levels = c("X1","X2","X3")) 
print(auc)
#plot(auc, ylim=c(0,1), print.thres=TRUE, main=paste('AUC:',round(auc$auc[[1]],2)))
#abline(h=1,col='blue',lwd=2)
#abline(h=0,col='red',lwd=2)

#Train Full Model and Assess Test Set Error on Test Data 
bst_model <- xgb.train(params = xgb_params,
                       data = train.matrix,
                       nrounds = nround)

# Predict hold-out test set
test_pred <- predict(bst_model, newdata = test.matrix)
test_prediction <- matrix(test_pred, nrow = numberOfClasses,
                          ncol=length(test_pred)/numberOfClasses) %>%
  t() %>%
  data.frame() %>%
  mutate(label = test.label + 1,
         max_prob = max.col(., "last"))
# confusion matrix of test set
confusionMatrix(factor(test_prediction$max_prob),
                factor(test_prediction$label),
                mode = "everything")

# get the feature real names
names = dimnames(train.matrix)[[1]]
# compute feature importance matrix
importance_matrix = xgb.importance(feature_names = names, model = bst_model)
head(importance_matrix)
write.csv(importance_matrix, file = "B.final.importance.multi.allvars.csv")

# plot
gp = xgb.ggplot.importance(importance_matrix)
print(gp) 





