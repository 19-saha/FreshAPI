predict.bacCounts <- function(bacteria, method, entry, info, product, platform){
  #
  # Function to predict selected bacteria counts  
  #
  # bacteria = information which bactera type was chosen to predict (TVC, Ps, Bth or LAB)
  # entry = a single FTIR entry row
  # info = data frame with informations about models (specification in a technical document)
  # method = information which machine learning method was chosen for prediction (rf, knn, svmLinear, tsvmRadial, svmPoly, lm)
  # platform = analytical platform which was used for measurement: FTIR or MSI
  # product = orgin of the sample: Chicken Thigh Fillet (CTF) or Chicken Burger (CB)
  #
  # Hint - function requires "FTIR_models.xlsx" file loaded as model.info variable 
  #        (model.info <- read_excel("FTIR_models.xlsx"))
  #      - prospectr package is required
  #      - readxl package is required
  #      - jsonlite package is required
  # 
  data <- entry
  model.info <- info
  prediction <- 0
  for(i in 1:nrow(model.info)){
    if(model.info$product[i] == product && model.info$bacterialCounts[i] == bacteria && model.info$ML[i] == method && model.info$platform[i] == platform){
      if(model.info$preprocessing[i] == "SG"){ #Savitzky-Golay filter
        data <- preprocess.sg(data) #appropriate pre-processing
        model <- readRDS(model.info$model[i])
        prediction <- stats::predict(model, data) #predict
      }
      else if(model.info$preprocessing[i] == "nSG"){ #normalisation + Savitzky-Golay filter
        spectra <- readRDS(model.info$rawSpectra[i])
        data <- preprocess.nor_sg(spectra, data) #appropriate pre-processing
        model <- readRDS(model.info$model[i])
        prediction <- stats::predict(model, data) #predict
      }
      else{ #raw data
        model <- readRDS(model.info$model[i])
        prediction <- stats::predict(model, data) #predict
      }
    }
  }
  prediction <- prediction
}


###################################################################################################################
###################################################################################################################
predict.bestTVC <- function(entry, info, product, platform){
  #
  # Function to predict selected bacteria counts  
  #
  # bacteria = information which bactera type was chosen to predict (TVC, Ps, Bth or LAB)
  # entry = a single FTIR entry row
  # info = data frame with informations about models (specification in a technical document)
  # method = information which machine learning method was chosen for prediction (rf, knn, svmLinear, tsvmRadial, svmPoly, lm)
  # platform = analytical platform which was used for measurement: FTIR or MSI
  # product = orgin of the sample: Chicken Thigh Fillet (CTF) or Chicken Burger (CB)
  #
  # Hint - function requires "models.xlsx" file loaded as model.info variable 
  #        (model.info <- read_excel("models.xlsx"))
  #      - prospectr package is required
  #      - readxl package is required
  #      - jsonlite package is required
  # 
  data <- entry
  model.info <- info
  prediction <- 0
  model_name <- 0
  for(i in 1:nrow(model.info)){
    if(model.info$product[i] == product && model.info$bacterialCounts[i]=="TVC" && model.info$place[i]=="1" && model.info$platform[i] == platform){
      if(model.info$preprocessing[i] == "SG" && model.info$platform == "FTIR"){ #Savitzky-Golay filter
        data <- preprocess.sg(data) #appropriate pre-processing
        model <- readRDS(model.info$model[i])
        prediction <- stats::predict(model, data) #predict
        if(model.info$ML[i] == "knn"){
          model_name <- "k-Nearest Neighbors (k-NN) algorithm"
        } else if(model.info$ML[i] == "rf"){
          model_name <- "Random Forest algorithm"
        } else if(model.info$ML[i] == "svmLinear"){
          model_name <- "Support-Vector Machine (SVM) with linear kernel"
        } else if(model.info$ML[i] == "svmRadial"){
          model_name <- "Support-Vector Machine (SVM) with radial kernel"
        } else if(model.info$ML[i] == "svmPoly"){
          model_name <- "Support-Vector Machine (SVM) with polynomial kernel"
        } else {
          model_name <- "Linear Model"
        }
      }
      else if(model.info$preprocessing[i] == "nSG" && model.info$platform == "FTIR"){ #normalisation + Savitzky-Golay filter
        spectra <- readRDS(model.info$rawSpectra[i])
        data <- preprocess.nor_sg(spectra, data) #appropriate pre-processing
        model <- readRDS(model.info$model[i])
        prediction <- stats::predict(model, data) #predict
        if(model.info$ML[i] == "knn"){
          model_name <- "k-Nearest Neighbors (k-NN) algorithm"
        } else if(model.info$ML[i] == "rf"){
          model_name <- "Random Forest algorithm"
        } else if(model.info$ML[i] == "svmLinear"){
          model_name <- "Support-Vector Machine (SVM) with linear kernel"
        } else if(model.info$ML[i] == "svmRadial"){
          model_name <- "Support-Vector Machine (SVM) with radial kernel"
        } else if(model.info$ML[i] == "svmPoly"){
          model_name <- "Support-Vector Machine (SVM) with polynomial kernel"
        } else {
          model_name <- "Linear Model"
        }
      }
      else{ #raw data
        model <- readRDS(model.info$model[i])
        prediction <- stats::predict(model, data) #predict
        if(model.info$ML[i] == "knn"){
          model_name <- "k-Nearest Neighbors (k-NN) algorithm"
        } else if(model.info$ML[i] == "rf"){
          model_name <- "Random Forest algorithm"
        } else if(model.info$ML[i] == "svmLinear"){
          model_name <- "Support-Vector Machine (SVM) with linear kernel"
        } else if(model.info$ML[i] == "svmRadial"){
          model_name <- "Support-Vector Machine (SVM) with radial kernel"
        } else if(model.info$ML[i] == "svmPoly"){
          model_name <- "Support-Vector Machine (SVM) with polynomial kernel"
        } else {
          model_name <- "Linear Model"
        }
      }
    }
  }
  return(paste0("The most accurate prediction using ", model_name, ": ", round(prediction,digits = 3)))
}

###################################################################################################################
###################################################################################################################
preprocess.sg <- function(singleRow){
  #
  # Function to apply Savitzky Golay filter on a single row spectra
  #
  # spectra = data frame with FTIR results based on which machine learning model was built
  # singleRow = a single row FTIR entry 
  #
  # Hint - prospectr package is required
  #
  singleRow <- singleRow[,order(as.numeric(colnames(singleRow)))] #order by wavelengths from least to greatest
  columns <- colnames(singleRow[,-c((ncol(singleRow)-9):ncol(singleRow))]) #write columnnames which will be used after Savitzky-Golay smoothing
  singleRow <- prospectr::savitzkyGolay(singleRow, w=11, p=3, m=1) #apply SG
  singleRow <- as.data.frame(singleRow) #write received single row as a dataframe
  colnames(singleRow) <- columns #set data frame columns 
  singleRow <- singleRow
}


###################################################################################################################
###################################################################################################################
preprocess.nor_sg <- function(spectra, singleRow){
  #
  # Function to apply Standard Normal Variate (SNV) Normalization and Savitzky Golay filter on a single row spectra
  #
  # spectra = data frame with FTIR results based on which machine learning model was built
  # singleRow = a single row FTIR entry 
  #
  # Hint - prospectr package is required
  #
  features <- c(1001:4000) #create a vector with wavelengths which are suitable for results which can be received from different machines  
  
  #feature extraction in spectra dataset
  selected_spectra <- spectra[, c(as.character(features))] #select features
  selected_spectra <- selected_spectra[, order(as.numeric(colnames(selected_spectra)))] #order by wavelengths from least to greatest
  
  #remove features in a single row
  selected_singleRow <- singleRow[,c(as.character(features))] #select features
  selected_singleRow <- selected_singleRow[,order(as.numeric(colnames(selected_singleRow)))] #order by wavelengths from least to greatest
  
  #merge spectra
  merged <- rbind.data.frame(selected_spectra, selected_singleRow) #merge a single row and spectra together
  #SNV
  snv <- scale(merged, center=TRUE, scale=TRUE) #normalise the data
  output <- snv[nrow(snv),] #extract the single row entry, which is the last row
  #apply SG
  columns <- colnames(snv[,-c((ncol(snv)-9):ncol(snv))]) #write columnnames which will be used after Savitzky-Golay smoothing
  output <- prospectr::savitzkyGolay(output, w=11, p=3, m=1) #apply SG
  output <- matrix(output, nrow=1) # #write received single row as a dataframe
  output <- as.data.frame(output) #write matrix as a dataframe
  colnames(output) <- columns #set column names
  singleRow <- output
}


####################################################################################################################
####################################################################################################################
roundWavelengths <- function(df){
  #
  # function to round column names/wavelengs 
  #
  # df = data frame where column names are numeric
  #
  # Hint - function works just for numeric column names
  #      - CreateObj function is required
  #
  roundedCol_df <- round(as.numeric(colnames(df))) #round column names/wavelengths
  colnames(df) <- roundedCol_df #set rounded wavelengths as column names
  df_orgin <- df 
  data <- df
  d <- c()
  for (i in 1:(ncol(df_orgin)-1)){
    if(as.numeric(colnames(df_orgin[i])) == as.numeric(colnames(df_orgin[i+1]))){
      for(j in 1:nrow(df_orgin)){ #calculate means based on columns of which wavelengths are the same after rounding
        d[j] <- mean(as.numeric(df_orgin[j,i]), as.numeric(df_orgin[j,i+1]))
      }
      d <- matrix(d, ncol = 1)
      d <- as.data.frame(d)
      rownames(d) <- rownames(df_orgin)
      colnames(d) <- colnames(df_orgin[i])
      for(z in 1:(length(colnames(data))-1)){ #create a dataframe where wavelengths are not repeated
        if(colnames(data[z]) == colnames(df_orgin[i])){
          df <- CreateObj(data[,1:(z-1)], d)
          df <- CreateObj(df, data[,(z+2):length(data)])
          data <- df
        }
      }
      d <- c()
      i = i+1
    }
  }
  dataFrame <- df
  
}


###################################################################################################################
###################################################################################################################
CreateObj <- function(data1, data2){
  #
  # Data pretreatment function to create combined data sets
  #
  # data1 = the first data set to combining
  # data2 = the second data set to combining
  #
  # Hint - each data set has to have the same number of rows
  #
  #combine all rows from the both dataset 
  merged <- merge(data1, data2, by = 'row.names')
  rownames(merged) = merged[,1]
  # remove the row names column, which was added th the merged dataset during mergining (additional one)
  as.data.frame(merged[,-1])
}

