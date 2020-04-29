# Clear workspace
rm(list=ls())

# Close any open graphics devices
graphics.off()

library("readxl")
library(caret)
library(Rmisc)

source("CreateObj.R")
source("accuracy.R")

# load xlsx file
df <- read_excel("chickenThighFillet_FTIR_raw_v3.xlsx")

#######################################################################################
# set the row names and remove the first column with the row IDs after this treatment
sampleNames <- c()
for(i in 1:nrow(df[,1])){
  sampleNames[i] <- df[i,1]
}
row.names(df) <- sampleNames
df <- df[,-1]

############## method to round wavelengths ###########################
###### df
#round colnames (wavelengths)
roundedCol_df <- round(as.numeric(colnames(df[,4:length(colnames(df))])))
colnames(df) <- c(colnames(df[,1:3]), roundedCol_df)
df_orgin <- df
data <- df
d <- c()
for(i in 4:(length(colnames(df_orgin))-1)){
  if(as.numeric(colnames(df_orgin[,i])) == as.numeric(colnames(df_orgin[,i+1]))){
    for(j in 1:nrow(df_orgin)){
      d[j] <- mean(as.numeric(df_orgin[j,i]), as.numeric(df_orgin[j,i+1]))
    }
    d <- matrix(d, ncol = 1)
    d <- as.data.frame(d)
    rownames(d) <- rownames(df_orgin)
    colnames(d) <- colnames(df_orgin[,i])
    for(z in 4:(length(colnames(data))-1)){
      if(as.numeric(colnames(data[z])) == as.numeric(colnames(df_orgin[,i]))){
        df <- CreateObj(data[,1:(z-1)], d)
        df <- CreateObj(df, data[,(z+2):length(data)])
        data <- df
      }
    }
    d <- c()
    i = i+1
  }
}

######################################################################
#feature extraction (for wavelength between 0 and 1000 nm)
removeID <- c()
features <- colnames(df)
for(i in 6:ncol(df)){
  if(as.numeric(features[i]) <= 1000){
    removeID[length(removeID)+1] <- i
  }
}
df <- df[,-removeID]

#####################################################################################################
############################# feature selection step 1 ##############################################
############################# remove redundant features #############################################

#ensure the results are repeatable
set.seed(123)
# create specra dataframe
spectra <- df[,4:ncol(df)]
#calculate correlation matrix
correlationMatrix <- cor(spectra)
#summarize the correlation matrix
print(correlationMatrix)
# find attributes which are highly correlated (ideally > 0.75)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff = 0.999)
# print index of highly correlated attributes
print(highlyCorrelated)
#remove highly correlated attributes
spectra2 <- spectra[, - highlyCorrelated]

######################################################
# that data is highly correlated 
# create a dataframe with removed features
df2 <- CreateObj(df[2:3], spectra2)

#######################################################################################################
############################# feature selection step 2 ################################################
############################# rank features by importance #############################################
#prerate training scheme
control <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
# create df for diffrent bacterial counts separetly
df_TSA_wf <- df2[,-2]
df_CFC_wf <- df2[,-c(1)]

#######################################################################################################
################################# Total Viable Counts (TVC) ###########################################
#######################################################################################################

######################### select features for the whole dataset #######################################
#ensure results are repeatable
set.seed(123)
#train the model
model_varImp_tvc <- train(TVC~., data = df_TSA_wf, method = "rf", trControl = control, importance = TRUE)
#estimate viable importance
importance <- varImp(model_varImp_tvc, scale = FALSE)
#sumarize importance
print(importance)
#plot importance
plot(importance)
# get an index of important variables
feature_ID_TSA <- c(1)
for(i in 1:nrow(importance$importance)){
  if(importance$importance[i,1] >= 0.5){
    feature_ID_TSA[length(feature_ID_TSA)+1] <- i+1
  }
}
# create a dataframe eith selected variables
df_TSA <- df_TSA_wf[,feature_ID_TSA]

###############################################################################################################
################################### Prediction for TVC ########################################################
###############################################################################################################

#split data
set.seed(123) #caret package

#data split into 70% for training and 30% for testing, based on df (whole data set)
trainIndex_tsa1 <- createDataPartition(df_TSA$TVC, p = .7, 
                                       list = FALSE, 
                                       times = 1, groups =3) #caret package
#split the whole dataset
train_tsa1 <- df_TSA[trainIndex_tsa1,]
test_tsa1 <- df_TSA[-trainIndex_tsa1,]

# resembling (method to avoid overfitting)
fitControl <- trainControl(## 10-fold CV
  method = "cv",
  number = 10)

#############################################################################################
############################# feature selection step 3 ################################################
############################# feature selection #############################################

# rfe control (define the control using a random forest selection function)
rctrl1 <- rfeControl(method = "cv",
                     number = 10, #10-fold cv
                     functions = rfFuncs)

######################### select features for the whole dataset #######################################
# run the RFE algorithm
results_TSA <- rfe(TVC ~ ., data = train_tsa1,
                   sizes = c(2:ncol(train_tsa1)), rfeControl = rctrl1)
#sumarize the results
print(results_TSA)
#list the chosen features
predictors_tsa <- predictors(results_TSA)
#plot the results
plot(results_TSA, type = c("g", "o"))
#create a dataframe with selected features (predictors)
predictors_tsa_id <- c(1)
for(i in 2:ncol(train_tsa1)){
  for(j in 1:length(predictors_tsa)){
    if( as.numeric(colnames(train_tsa1[i]))  == as.numeric(gsub("`", "", (as.name(predictors_tsa[j]))))){
      predictors_tsa_id[length(predictors_tsa_id)+1] <- i
    }
  }
}
train_tsa2 <- train_tsa1[,predictors_tsa_id] #training set with just selected features
test_tsa2 <- test_tsa1[,predictors_tsa_id] #testing set with just selected features
df2_TSA <- df_TSA[,predictors_tsa_id] #not-splitted data set with just selected features

#save selected features in one vector
predictors_tsa_toSave <- c() # vector with selected feature names (rounded wavelengths)
len_predictors_tsa <- 0 # number of selected features
for(i in predictors_tsa){ # to get number of selected features
  if(as.numeric(lengths(predictors_tsa[i]))== 1){
    len_predictors_tsa <- len_predictors_tsa +1
  }
}
for(i in 1:len_predictors_tsa){ # to get vector with selected feature names (rounded wavelengths)
  predictors_tsa_toSave[i] <- as.numeric(gsub("`", "", (as.name(predictors_tsa[i]))))
}
saveRDS(predictors_tsa_toSave, "predictors_TVC.rds") #save model as an rds object

##################################################################################################
##### Random Forest for Total Viable Counts (TVC) by chickenThighFillet_FTIR_raw_v3.xlsx file ####

model_Rf_TSA1 <- train(TVC ~ ., data = train_tsa2,
                       method = "rf", trControl = fitControl) #train a model
pred_Rf_TSA1<- predict(model_Rf_TSA1,test_tsa2) # predict bacterial counts
RMSE.Rf_TSA1 <- RMSE(test_tsa2$TVC,pred_Rf_TSA1) #compute RMSE  #RMSE =1.32
# Random Forest model Accuracy for TVC
accuracy_Rf_TSA <- accuracy(pred_Rf_TSA1, test_tsa2$TVC)

# Random Forest model plot predicted vs Observed
plot(pred_Rf_TSA1,test_tsa2$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Random Forest model for Total Viable Counts (TVC) \nRMSE:",round(RMSE.Rf_TSA1,digits = 2), 
                "\nAccuracy",accuracy_Rf_TSA,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(model_Rf_TSA1, "model_RF_TVC.rds") #save model as an rds object

######################################################################################################################
### k-nearest neighbors (k-NN) algorithm for Total Viable Counts (TVC) by chickenThighFillet_FTIR_raw_v3.xlsx file ###

model_knn_TSA1 <- train(TVC ~ ., data = train_tsa2,
                        method = "knn", trControl = fitControl) #train a model
pred_knn_TSA1<- predict(model_knn_TSA1,test_tsa2) # predict bacterial counts
RMSE.knn_TSA1 <- RMSE(test_tsa2$TVC,pred_knn_TSA1) # compute RMSE #RMSE =1.37
# k-nearest neighbors (k-NN) model Accuracy for TVC
accuracy_knn_TSA <- accuracy(pred_knn_TSA1, test_tsa2$TVC)

# k-nearest neighbors (k-NN) model plot predicted vs Observed
plot(pred_knn_TSA1,test_tsa2$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("k-nearest neighbors (k-NN) model for Total Viable Counts (TVC) \nRMSE:",round(RMSE.knn_TSA1,digits = 2), 
                "\nAccuracy",accuracy_knn_TSA,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(model_knn_TSA1, "model_knn_TVC.rds") #save model as an rds object

#####################################################################################################################
### Linear support-vector machine (SVM) for Total Viable Counts (TVC) by chickenThighFillet_FTIR_raw_v3.xlsx file ###

model_SVMLM_TSA1 <- train(TVC ~ ., data = train_tsa2,
                          method = "svmLinear", trControl = fitControl) # model training
pred_SVMLM_TSA1<- predict(model_SVMLM_TSA1,test_tsa2) # predict bacterial counts
RMSE.SVMLM_TSA1 <- RMSE(test_tsa2$TVC,pred_SVMLM_TSA1) #calculate RMSE #RMSE =1.25
# linear support-vector machine (SVM) model Accuracy for TVC
accuracy_SVMLM_TSA <- accuracy(pred_SVMLM_TSA1, test_tsa2$TVC)

# Linear support-vector machine (SVM) model plot predicted vs Observed
plot(pred_SVMLM_TSA1,test_tsa2$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Linear support-vector machine (SVM) model for Total Viable Counts (TVC) \nRMSE:",round(RMSE.SVMLM_TSA1,digits = 2), 
                "\nAccuracy",accuracy_SVMLM_TSA,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(model_SVMLM_TSA1, "model_SVMLM_TVC.rds") #save model as an rds object

#####################################################################################################################
### Radial support-vector machine (SVM) for Total Viable Counts (TVC) by chickenThighFillet_FTIR_raw_v3.xlsx file ###

##################
# the first method
model_SVMR_TSA1 <- train(TVC ~ ., data = train_tsa2,
                         method = "svmRadial", trControl = fitControl) # model training
pred_SVMR_TSA1<- predict(model_SVMR_TSA1,test_tsa2) # predict bacterial counts
RMSE.SVMR_TSA1 <- RMSE(test_tsa2$TVC,pred_SVMR_TSA1) #compute RMSE #RMSE =1.26
# Radial support-vector machine (SVM) model Accuracy for TVC
accuracy_SVMR_TSA <- accuracy(pred_SVMR_TSA1, test_tsa2$TVC)

# Radial support-vector machine (SVM) model plot predicted vs Observed
plot(pred_SVMR_TSA1,test_tsa2$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Radial support-vector machine (SVM) model for Total Viable Counts (TVC) \nRMSE:",round(RMSE.SVMR_TSA1,digits = 2), 
                "\nAccuracy",accuracy_SVMR_TSA,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(model_SVMR_TSA1, "model_SVMR_TVC.rds") #save model as an rds object

#########################################################################################################################
### Polynomial support-vector machine (SVM) for Total Viable Counts (TVC) by chickenThighFillet_FTIR_raw_v3.xlsx file ###

##################
# the first method
model_SVMP_TSA1 <- train(TVC ~ ., data = train_tsa2,
                         method = "svmPoly", trControl = fitControl) # model training
pred_SVMP_TSA1<- predict(model_SVMP_TSA1,test_tsa2) # predict bacterial counts
RMSE.SVMP_TSA1 <- RMSE(test_tsa2$TVC,pred_SVMP_TSA1) #compute RMSE #RMSE =1.26
# Polynomial support-vector machine (SVM) model Accuracy for TVC
accuracy_SVMP_TSA <- accuracy(pred_SVMP_TSA1, test_tsa2$TVC)

# Polynomial support-vector machine (SVM) model plot predicted vs Observed
plot(pred_SVMP_TSA1,test_tsa2$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Polynomial support-vector machine (SVM) model for Total Viable Counts (TVC) \nRMSE:",round(RMSE.SVMP_TSA1,digits = 2), 
                "\nAccuracy",accuracy_SVMP_TSA,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(model_SVMP_TSA1, "model_SVMP_TVC.rds") #save model as an rds object

##############################################################################################
### Linear Model for Total Viable Counts (TVC) by chickenThighFillet_FTIR_raw_v3.xlsx file ###

model_lm_TSA1 <- train(TVC ~ ., data = train_tsa2,
                       method = "lm", trControl = fitControl) # model training
pred_lm_TSA1<- predict(model_lm_TSA1,test_tsa2) # predict bacterial counts
RMSE.lm_TSA1 <- RMSE(test_tsa2$TVC,pred_lm_TSA1) #compute RMSE #RMSE =1.4
# Linear model Accuracy for TVC
accuracy_lm_TSA <- accuracy(pred_lm_TSA1, test_tsa2$TVC)

# Linear model plot predicted vs Observed
plot(pred_lm_TSA1,test_tsa2$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Linear Model for Total Viable Counts (TVC) \nRMSE:",round(RMSE.lm_TSA1,digits = 2), 
                "\nAccuracy",accuracy_lm_TSA,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(model_lm_TSA1, "model_LM_TVC.rds") #save model as an rds object

#######################################################################################################
################### Pseudomonas spp. colony forming units (CFU) counts ################################
#######################################################################################################

######################### select features for the whole dataset #######################################
#ensure results are repeatable
set.seed(123)
#train the model
model_varImp_cfc <- train(Ps~., data = df_CFC_wf, method = "rf", trControl = control, importance = TRUE)
#estimate viable importance
importance_CFC <- varImp(model_varImp_cfc, scale = FALSE)
#sumarize importance
print(importance_CFC)
#plot importance
plot(importance_CFC)
# get an index of important variables
feature_ID_CFC <- c(1)
for(i in 1:nrow(importance_CFC$importance)){
  if(importance_CFC$importance[i,1] >= 0.5){
    feature_ID_CFC[length(feature_ID_CFC)+1] <- i+1
  }
}
# create a dataframe eith selected variables
df_CFC <- df_CFC_wf[,feature_ID_CFC]

###############################################################################################################
############################## Prediction for Pseudomonas #####################################################
###############################################################################################################

#split data
set.seed(123) #caret package

#data split into 70% for training and 30% for testing, based on df (whole data set)
trainIndex_cfc1 <- createDataPartition(df_CFC$Ps, p = .7, 
                                       list = FALSE, 
                                       times = 1, groups =3) #caret package
#split the whole dataset
train_cfc1 <- df_CFC[trainIndex_cfc1,]
test_cfc1 <- df_CFC[-trainIndex_cfc1,]

# resembling (method to avoid overfitting)
fitControl <- trainControl(## 10-fold CV
  method = "cv",
  number = 10)
## repeated ten times

#############################################################################################
############################# feature selection step 3 ################################################
############################# feature selection #############################################

# rfe control (define the control using a random forest selection function)
rctrl1 <- rfeControl(method = "cv",
                     number = 10, #10-fold cv
                     functions = rfFuncs)

######################### select features for the whole dataset #######################################
# run the RFE algorithm
results_CFC <- rfe(Ps ~ ., data = train_cfc1,
                   sizes = c(2:ncol(train_cfc1)), rfeControl = rctrl1)
#sumarize the results
print(results_CFC)
#list the chosen features
predictors_cfc <- predictors(results_CFC)
#plot the results
plot(results_CFC, type = c("g", "o"))
#create a dataframe with selected features
predictors_cfc_id <- c(1)
for(i in 2:ncol(train_cfc1)){
  for(j in 1:length(predictors_cfc)){
    if( as.numeric(colnames(train_cfc1[i]))  == as.numeric(gsub("`", "", (as.name(predictors_cfc[j]))))){
      predictors_cfc_id[length(predictors_cfc_id)+1] <- i
    }
  }
}
train_cfc2 <- train_cfc1[,predictors_cfc_id]
test_cfc2 <- test_cfc1[,predictors_cfc_id]
df2_CFC <- df_CFC[,predictors_cfc_id]

#save selected features in one vector
predictors_cfc_toSave <- c()
len_predictors_cfc <- 0
for(i in predictors_cfc){
  if(as.numeric(lengths(predictors_cfc[i]))== 1){
    len_predictors_cfc <- len_predictors_cfc +1
  }
}
for(i in 1:len_predictors_cfc){
  predictors_cfc_toSave[i] <- as.numeric(gsub("`", "", (as.name(predictors_cfc[i]))))
}
saveRDS(predictors_cfc_toSave, "predictors_PS.rds") #save model as an rds object

#################################################################################################
### Random Forest for Pseudomonas spp. CFU counts by chickenThighFillet_FTIR_raw_v3.xlsx file ###

model_Rf_CFC1 <- train(Ps ~ ., data = train_cfc2,
                       method = "rf", trControl = fitControl) # model training
pred_Rf_CFC1<- predict(model_Rf_CFC1,test_cfc2) # predict bacterial counts
RMSE.Rf_CFC1 <- RMSE(test_cfc2$Ps,pred_Rf_CFC1) #compute RMSE #RMSE =1.32
# Random Forest model Accuracy for Pseudomonas spp. colony forming units (CFU) counts
accuracy_Rf_CFC <- accuracy(pred_Rf_CFC1, test_cfc2$Ps)

# Random Forest model plot predicted vs Observed
plot(pred_Rf_CFC1,test_cfc2$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 Pseudomonas spp. CFU counts/g",
     ylab="Actual log10 Pseudomonas spp. CFU counts/g",col = "blue", 
     main=paste("Random Forest model for Pseudomonas spp. colony forming units (CFU) counts \nRMSE:",
                round(RMSE.Rf_CFC1,digits = 2), 
                "\nAccuracy",accuracy_Rf_CFC,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(model_Rf_CFC1, "model_RF_PS.rds") #save model as an rds object

########################################################################################################################
### k-nearest neighbors (k-NN) algorithm for Pseudomonas spp. CFU counts by chickenThighFillet_FTIR_raw_v3.xlsx file ###

model_knn_CFC1 <- train(Ps ~ ., data = train_cfc2,
                        method = "knn", trControl = fitControl) # model training
pred_knn_CFC1<- predict(model_knn_CFC1,test_cfc2) # predict bacterial counts
RMSE.knn_CFC1 <- RMSE(test_cfc2$Ps,pred_knn_CFC1) #compute RMSE #RMSE =1.33
# k-nearest neighbors (k-NN) model Accuracy for Pseudomonas spp. colony forming units (CFU) counts
accuracy_knn_CFC <- accuracy(pred_knn_CFC1, test_cfc2$Ps)

# k-nearest neighbors (k-NN) model plot predicted vs Observed
plot(pred_knn_CFC1,test_cfc2$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 Pseudomonas spp. CFU counts/g",
     ylab="Actual log10 Pseudomonas spp. CFU counts/g",col = "blue", 
     main=paste("k-nearest neighbors (k-NN) model for Pseudomonas spp. colony forming units (CFU) counts \nRMSE:",
                round(RMSE.knn_CFC1,digits = 2), 
                "\nAccuracy",accuracy_knn_CFC,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(model_knn_CFC1, "model_knn_PS.rds") #save model as an rds object


#######################################################################################################################
### Linear support-vector machine (SVM) for Pseudomonas spp. CFU counts by chickenThighFillet_FTIR_raw_v3.xlsx file ###

model_SVMLM_CFC1 <- train(Ps ~ ., data = train_cfc2,
                          method = "svmLinear", trControl = fitControl) # model training
pred_SVMLM_CFC1<- predict(model_SVMLM_CFC1,test_cfc2) # predict bacterial counts
RMSE.SVMLM_CFC1 <- RMSE(test_cfc2$Ps,pred_SVMLM_CFC1) #compute RMSE #RMSE =1.24
# Linear support-vector machine (SVM) model Accuracy for Pseudomonas spp. colony forming units (CFU) counts
accuracy_SVMLM_CFC <- accuracy(pred_SVMLM_CFC1, test_cfc2$Ps)

# Linear support-vector machine (SVM) model plot predicted vs Observed
plot(pred_SVMLM_CFC1,test_cfc2$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 Pseudomonas spp. CFU counts/g",
     ylab="Actual log10 Pseudomonas spp. CFU counts/g",col = "blue", 
     main=paste("Linear support-vector machine (SVM) model for Pseudomonas spp. colony forming units (CFU) counts \nRMSE:",
                round(RMSE.SVMLM_CFC1,digits = 2), 
                "\nAccuracy",accuracy_SVMLM_CFC,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(model_SVMLM_CFC1, "model_SVMLM_PS.rds") #save model as an rds object


#######################################################################################################################
### Radial support-vector machine (SVM) for Pseudomonas spp. CFU counts by chickenThighFillet_FTIR_raw_v3.xlsx file ###

model_SVMR_CFC1 <- train(Ps ~ ., data = train_cfc2,
                         method = "svmRadial", trControl = fitControl) # model training
pred_SVMR_CFC1<- predict(model_SVMR_CFC1,test_cfc2) # predict bacterial counts
RMSE.SVMR_CFC1 <- RMSE(test_cfc2$Ps,pred_SVMR_CFC1) #compute RMSE #RMSE =1.29
# Radial support-vector machine (SVM) model Accuracy for Pseudomonas spp. colony forming units (CFU) counts
accuracy_SVMR_CFC <- accuracy(pred_SVMR_CFC1, test_cfc2$Ps)

# Radial support-vector machine (SVM) model plot predicted vs Observed
plot(pred_SVMR_CFC1,test_cfc2$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 Pseudomonas spp. CFU counts/g",
     ylab="Actual log10 Pseudomonas spp. CFU counts/g",col = "blue", 
     main=paste("Radial support-vector machine (SVM) model for Pseudomonas spp. colony forming units (CFU) counts \nRMSE:",
                round(RMSE.SVMR_CFC1,digits = 2), 
                "\nAccuracy",accuracy_SVMR_CFC,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(model_SVMR_CFC1, "model_SVMR_PS.rds") #save model as an rds object

###########################################################################################################################
### Polynomial support-vector machine (SVM) for Pseudomonas spp. CFU counts by chickenThighFillet_FTIR_raw_v3.xlsx file ###

model_SVMP_CFC1 <- train(Ps ~ ., data = train_cfc2,
                         method = "svmPoly", trControl = fitControl) # model training
pred_SVMP_CFC1<- predict(model_SVMP_CFC1,test_cfc2)# predict bacterial counts
RMSE.SVMP_CFC1 <- RMSE(test_cfc2$Ps,pred_SVMP_CFC1) #compute RMSE #RMSE =1.26
# Polynomial support-vector machine (SVM) model Accuracy for Pseudomonas spp. colony forming units (CFU) counts
accuracy_SVMP_CFC <- accuracy(pred_SVMP_CFC1, test_cfc2$Ps)

# Polynomial support-vector machine (SVM) model plot predicted vs Observed
plot(pred_SVMP_CFC1,test_cfc2$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 Pseudomonas spp. CFU counts/g",
     ylab="Actual log10 Pseudomonas spp. CFU counts/g",col = "blue", 
     main=paste("Polynomial support-vector machine (SVM) model for Pseudomonas spp. colony forming units (CFU) counts \nRMSE:",
                round(RMSE.SVMP_CFC1,digits = 2), 
                "\nAccuracy",accuracy_SVMP_CFC,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(model_SVMP_CFC1, "model_SVMP_PS.rds") #save model as an rds object

################################################################################################
### Linear model for Pseudomonas spp. CFU counts by chickenThighFillet_FTIR_raw_v3.xlsx file ###

model_lm_CFC1 <- train(Ps ~ ., data = train_cfc2,
                       method = "lm", trControl = fitControl) # model training
pred_lm_CFC1<- predict(model_lm_CFC1,test_cfc2) # predict bacterial counts
RMSE.lm_CFC1 <- RMSE(test_cfc2$Ps,pred_lm_CFC1) #compute RMSE #RMSE =1.4
# Linear model Accuracy for Pseudomonas spp. colony forming units (CFU) counts
accuracy_lm_CFC <- accuracy(pred_lm_CFC1, test_cfc2$Ps)

# Linear model plot predicted vs Observed
plot(pred_lm_CFC1,test_cfc2$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 Pseudomonas spp. CFU counts/g",
     ylab="Actual log10 Pseudomonas spp. CFU counts/g",col = "blue", 
     main=paste("Linear Model for Pseudomonas spp. colony forming units (CFU) counts \nRMSE:",
                round(RMSE.lm_CFC1,digits = 2), 
                "\nAccuracy",accuracy_lm_CFC,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(model_lm_CFC1, "model_LM_PS.rds") #save model as an rds object
