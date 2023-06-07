##Multinomial Classification InCHIANTI Multiclass random forest classifier to classify characteristics of frailty
###06.4.23

#precision = how few false positives and recall = sensitivity 
#F1 harmonic mean 
#high precision high recall (but can handle lower recall)
#xgb.max_precision(pred, dtrain)

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
#install.packages("SHAPforxgboost")
library(SHAPforxgboost)

Multinominal1 <- read.csv("V1_multi_With_Performance.csv", header = T)
head(Multinominal1)
str(Multinominal1)
dim(Multinominal1)
Multinominal <- Multinominal1[,c(1,4,6:8,10:301)]
head(Multinominal)
dim(Multinominal)
table(Multinominal$AXNMULTIFRAIL)
class(Multinominal)
Multi2 <- Multinominal %>% drop_na(AXNMULTIFRAIL)
table(Multinominal$AXNMULTIFRAIL)
dim(Multi2)

#IXCESD_T_Depression.1 = feature = CESD continuous 


#############################################################split and create train test data sets########################################################################

inTrain = createDataPartition(Multi2$AXNMULTIFRAIL, p=2/3, list = FALSE)

#Train
train.data=Multi2[inTrain,]
train.label=train.data[,"AXNMULTIFRAIL"]
#remove class col
train.data=train.data[,-2]
dim(as.matrix(train.data))
length(train.label)
#create xgboost matrix
train.matrix=xgb.DMatrix(data=as.matrix(train.data), label=train.label)

#Test
test.data=Multi2[-inTrain,]
test.label=test.data[,"AXNMULTIFRAIL"]
#remove class col
test.data=test.data[,-2]
test.matrix=xgb.DMatrix(data=as.matrix(test.data), label=test.label)

table(train.data$AXNMULTIFRAIL)
 
table(test.data$AXNMULTIFRAIL)
#transfer to xgboost matrix

############################################################################Sort data down to features from original binary model ############################################################################

#check all features for NAs Table on all variables NAs 
table(Multi2$X_25OH_D)
sum(is.na(Multi2$X_25OH_D))
myfeatures <- c("IID","AXNMULTIFRAIL","IXAGE","X_UCOR24","X_25OH_D","X_CPK","X_CL24","X_VES","X_C24_0B","X_FOLICG","X_FREETS","X_HOMCYS","X_IL1B","X_IL6","X_VGM","X_MCP1_B","X_RETINL","X_TNFAR1","X_TNFAR2",
                "X_U_PRO","X_VITB6G","X_GAMTOC","IXCESD_T_Depression.1","Total_ACB","X_UCRE24","X_GLU","X_BUN","X_FT4","X_GPT","X_COLHDL","X_LYCOPN","X_GB","X_SCD14","X_IL1RA")

for (x in names(Multi2)){
  print(x)
  print(sum(is.na(Multi2[,x])))
  }
  
for (x in myfeatures){
  print(x)
  print(sum(is.na(Multi2[,x])))
}

###############################################################################################################################################################

##########################################################################FINAL Model Data run model ##################################################

Multi4 <- Multi2[,c("AXNMULTIFRAIL","X_FREETS","X_CPK","X_TNFAR1","X_GLU","X_PTH","X_GPT","X_BUN","X_COLHDL",
                    "X_LYCOPN","X_FT4","IXCESD_T_Depression.1","X_IL6","X_HOMCYS","X_VITB6G","X_VES","X_FOLICG","X_TNFAR2","X_UCRE24","X_25OH_D","IXAGE")]

###Age alone
Multi4 <- Multi2[,c("AXNMULTIFRAIL","IXAGE")]

#Age and CESD continuous
Multi4 <- Multi2[,c("AXNMULTIFRAIL","IXAGE","IXCESD_T_Depression.1")]

#Age and CESD continuous and ACB---not completed
#Multi4 <- Multi2[,c("AXNMULTIFRAIL","IXAGE","Total_ACB")]

head(Multi4)

set.seed(123)

#sample of 2/3rds of each level training then 1/3 for remaining testing 
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
train.balanced1 <- ovun.sample(AXNMULTIFRAIL ~ ., data=dataframe1, N=180,seed=1,method="both")$data
table(train.balanced1$AXNMULTIFRAIL)
#now balance the data 
dataframe2 <- train.data[train.data$AXNMULTIFRAIL %in% c(0,1),]
table(dataframe2$AXNMULTIFRAIL)

train.balanced2 <- ovun.sample(AXNMULTIFRAIL ~ ., data=dataframe2, N=180,seed=2,method="both")$data
table(train.balanced2$AXNMULTIFRAIL)

#take 1 out of train.balance2
dataframe3 <- train.balanced2[train.balanced2$AXNMULTIFRAIL %in% c(1),]
# now combine together
final.train.balanced <- merge(dataframe3,train.balanced1, all=TRUE)
dim(final.train.balanced)
table(final.train.balanced$AXNMULTIFRAIL)

# then do the rest of the train and test 
train.label=final.train.balanced[,"AXNMULTIFRAIL"]
#remove class col
final.train.balanced=final.train.balanced[,-1]
dim(as.matrix(final.train.balanced))
#create xgboost matrix
train.matrix=xgb.DMatrix(data=as.matrix(final.train.balanced), label=train.label)

summary(final.train.balanced)

#Test
test.data=Multi4[-inTrain,]
test.label=test.data[,"AXNMULTIFRAIL"]
#remove class col
test.data=test.data[,-1]
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
#length(y_pred)
#y_pred<-as.ordered(OOF_prediction)
table(Multi4$AXNMULTIFRAIL)
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

#AUC test
Pheno <- test.label
Pheno[Pheno==0] <- "X1"
Pheno[Pheno==1] <- "X2"
Pheno[Pheno==2] <- "X3"
Pheno=as.factor(Pheno)
length(Pheno)
dim(test_prediction)
auc <- multiclass.roc(Pheno,test_prediction[,c("X1","X2","X3")], levels = c("X1","X2","X3")) 
print(auc)


# get the feature real names
names <-  colnames(Multi4[,-1])
#names = dimnames(train.matrix)[[1]]
# compute feature importance matrix
importance_matrix = xgb.importance(feature_names = names, model = bst_model)
head(importance_matrix)
write.csv(importance_matrix, file = "final.importance.multi.allvars.csv")

# plot
gp = xgb.ggplot.importance(importance_matrix)
print(gp) 




