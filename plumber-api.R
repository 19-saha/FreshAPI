#* @apiTitle Food Spoilage API
#* @apiDescription Predicting the Bacterial counts to determine level of spoilage in meat samples

model.info <- readxl::read_excel("models.xlsx")
source("functions.R")

#* parse JSON
#* @post /predict_MSI
function(req, bacteria, model){
  #
  # Function to give an user repond for REST API request 
  #
  # req = user request (a single row MSI entry)
  # bacteria = information which bactera type was chosen to predict (TVC, Ps, Bth or LAB)
  # method = information which machine learning method was chosen for prediction (rf, knn, svmLinear, tsvmRadial, svmPoly, lm)
  #
  # Hint - function requires "models.xlsx" file loaded as model.info variable 
  #        (model.info <- read_excel("models.xlsx"))
  #      - readxl package is required
  #      - jsonlite package is required
  #      - function.R script is required
  #        (source("functions.R"))
  #
  
  platform <- "MSI" #specify an analytical platform
  product <- "CTF" #specify a product (platform shall adopt entries just for Chicken Thigh Fillet)
  #load single row entry
  data <- jsonlite::fromJSON(req$postBody)
  data <- data[,1:18]
  colnames(data) <- c("Mean.01", "Mean.02", "Mean.03", "Mean.04", "Mean.05", "Mean.06",
                      "Mean.07", "Mean.08", "Mean.09", "Mean.10", "Mean.11", "Mean.12",
                      "Mean.13", "Mean.14", "Mean.15", "Mean.16", "Mean.17", "Mean.18")
  
  prediction <- 0
  # predict bacterial counts
  if(bacteria == "TVC"){
    if(model == "rf"){
      prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
      return(paste0("Prediction of Total Viable Counts using randomForest: ", round(prediction,digits = 3)))
    }
    else if (model=="knn"){
      prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
      return(paste0("Prediction of Total Viable Counts using KNN: ", round(prediction,digits = 3)))
    }
    else if (model=="svmLinear"){
      prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
      return(paste0("Prediction of Total Viable Counts using svmLinear: ", round(prediction,digits = 3)))
    }
    else if (model=="svmRadial"){
      prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
      return(paste0("Prediction of Total Viable Counts using svmRadial: ", round(prediction,digits = 3)))
    }
    else if (model=="svmPoly"){
      prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
      return(paste0("Prediction of Total Viable Counts using svmPoly: ", round(prediction,digits = 3)))
    }
    else if (model=="lm"){
      prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
      return(paste0("Prediction of Total Viable Counts using Linear Regression: ", round(prediction,digits = 3)))
    }
    else{
      print("Error occured: Please choose the correct model: rf | knn | svmLinear | svmRadial | svmPoly | lm")
    }
  } 
  ##########################
  else if(bacteria == "Ps"){
    if(model == "rf"){
      prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
      return(paste0("Prediction of Pseudomonas count using randomForest: ", round(prediction,digits = 3)))
    }
    else if (model=="knn"){
      prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
      return(paste0("Prediction of Pseudomonas count using KNN: ", round(prediction,digits = 3)))
    }
    else if (model=="svmLinear"){
      prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
      return(paste0("Prediction of Pseudomonas count using svmLinear: ", round(prediction,digits = 3)))
    }
    else if (model=="svmRadial"){
      prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
      return(paste0("Prediction of Pseudomonas count using svmRadial: ", round(prediction,digits = 3)))
    }
    else if (model=="svmPoly"){
      prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
      return(paste0("Prediction of Pseudomonas count using svmPoly: ", round(prediction,digits = 3)))
    }
    else if (model=="lm"){
      prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
      return(paste0("Prediction of Pseudomonas count using Linear Regression: ", round(prediction,digits = 3)))
    }
    else{
      print("Error occured: Please choose the correct model: rf | knn | svmLinear | svmRadial | svmPoly | lm")
    }
  }
  else{
    print("Error occured: Please choose the correct bacteria: TVC | Ps")
  }
}

###################################################################################################################
###################################################################################################################
#* parse JSON
#* @post /predict_FTIR
function(req, product, bacteria, model){
  #
  # Function to give an user repond for REST API request
  #
  # req = user request (a single row FTIR entry)
  # product = chicken burger (paramiter: CB) or chicken thigh filet (parameter: CTF)
  # bacteria = information which bactera type was chosen to predict (TVC, Ps, Bth or LAB)
  # method = information which machine learning method was chosen for prediction (rf, knn, svmLinear, tsvmRadial, svmPoly, lm)
  #
  # Hint - function requires "models.xlsx" file loaded as model.info variable 
  #        (model.info <- read_excel("models.xlsx"))
  #      - prospectr package is required
  #      - readxl package is required
  #      - jsonlite package is required
  #      - function.R script is required
  #        (source("functions.R"))
  #
  
  #load single row entry
  data <- jsonlite::fromJSON(req$postBody)
  platform <- "FTIR" #specify an analytical platform
  
  # round if this is required
  row_colnames<- colnames(data)
  row_col_integer <- 0
  for(i in 1:length(row_colnames)){
    if(as.numeric(row_colnames[i])==round(as.numeric(row_colnames[i]), 0)){
      row_col_integer <- row_col_integer+1
    }
    else{
      row_col_integer <- row_col_integer
    }
  }
  if(length(row_colnames) == row_col_integer){
    data<-data  
  } else{
    data <- data
    data <- roundWavelengths(data)
  }
  
  # predict bacterial counts
  if(product=="CTF"){
    if(bacteria=="TVC"){
      if (model=="rf"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Total Viable Counts using randomForest: ", round(prediction,digits = 3)))
      }
      else if (model=="knn"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Total Viable Counts using KNN: ", round(prediction,digits = 3)))
      }
      else if (model=="svmLinear"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Total Viable Counts using svmLinear: ", round(prediction,digits = 3)))
      }
      else if (model=="svmRadial"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Total Viable Counts using svmRadial: ", round(prediction,digits = 3)))
      }
      else if (model=="svmPoly"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Total Viable Counts using svmPoly: ", round(prediction,digits = 3)))
      }
      else if (model=="lm"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Total Viable Counts using Linear Regression: ", round(prediction,digits = 3)))
      }
      else{
        print("Error occured: Please choose the correct model: rf | knn | svmLinear | svmRadial | svmPoly | lm")
      }
    } 
    
    #########################
    else if(bacteria=="Ps"){
      if (model=="rf"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Pseudomonas count using randomForest: ", round(prediction,digits = 3)))
      }
      else if (model=="knn"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Pseudomonas count using KNN: ", round(prediction,digits = 3)))
      }
      else if (model=="svmLinear"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Pseudomonas count using svmLinear: ", round(prediction,digits = 3)))
      }
      else if (model=="svmRadial"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Pseudomonas count using svmRadial: ", round(prediction,digits = 3)))
      }
      else if (model=="svmPoly"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Pseudomonas count using svmPoly: ", round(prediction,digits = 3)))
      }
      else if (model=="lm"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Pseudomonas count using Linear Regression: ", round(prediction,digits = 3)))
      }
      else{
        print("Error occured: Please choose the correct model: rf | knn | svmLinear | svmRadial | svmPoly | lm")
      }
    }
    else{
      print("Error occured: Please choose the correct bacteria: TVC | Ps")
    }
  } 
  
  ########################
  else if(product=="CB"){
    if(bacteria=="TVC"){
      if (model=="rf"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Total Viable Counts using randomForest: ", round(prediction,digits = 3)))
      }
      else if (model=="knn"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Total Viable Counts using KNN: ", round(prediction,digits = 3)))
      }
      else if (model=="svmLinear"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Total Viable Counts using svmLinear: ", round(prediction,digits = 3)))
      }
      else if (model=="svmRadial"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Total Viable Counts using svmRadial: ", round(prediction,digits = 3)))
      }
      else if (model=="svmPoly"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Total Viable Counts using svmPoly: ", round(prediction,digits = 3)))
      }
      else if (model=="lm"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Total Viable Counts using Linear Regression: ", round(prediction,digits = 3)))
      }
      else{
        print("Error occured: Please choose the correct model: rf | knn | svmLinear | svmRadial | svmPoly | lm")
      }
    }
    
    ###########################
    else if(bacteria == "Ps"){
      if (model=="rf"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Pseudomonas count using randomForest: ", round(prediction,digits = 3)))
      }
      else if (model=="knn"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Pseudomonas count using KNN: ", round(prediction,digits = 3)))
      }
      else if (model=="svmLinear"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Pseudomonas count using svmLinear: ", round(prediction,digits = 3)))
      }
      else if (model=="svmRadial"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Pseudomonas count using svmRadial: ", round(prediction,digits = 3)))
      }
      else if (model=="svmPoly"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Pseudomonas count using svmPoly: ", round(prediction,digits = 3)))
      }
      else if (model=="lm"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Pseudomonas count using Linear Regression: ", round(prediction,digits = 3)))
      }
      else{
        print("Error occured: Please choose the correct model: rf | knn | svmLinear | svmRadial | svmPoly | lm")
      }
    }
    
    ###########################
    else if(bacteria == "Bth"){
      if (model=="rf"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Brochothrix thermosphacta count using randomForest: ", round(prediction,digits = 3)))
      }
      else if (model=="knn"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Brochothrix thermosphacta count using KNN: ", round(prediction,digits = 3)))
      }
      else if (model=="svmLinear"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Brochothrix thermosphacta count using svmLinear: ", round(prediction,digits = 3)))
      }
      else if (model=="svmRadial"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Brochothrix thermosphacta count using svmRadial: ", round(prediction,digits = 3)))
      }
      else if (model=="svmPoly"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Brochothrix thermosphacta count using svmPoly: ", round(prediction,digits = 3)))
      }
      else if (model=="lm"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Brochothrix thermosphacta count using Linear Regression: ", round(prediction,digits = 3)))
      }
      else{
        print("Error occured: Please choose the correct model: rf | knn | svmLinear | svmRadial | svmPoly | lm")
      }
    }
    
    ###########################
    else if(bacteria == "LAB"){
      if (model=="rf"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Lactic acid bacteria count using randomForest: ", round(prediction,digits = 3)))
      }
      else if (model=="knn"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Lactic acid bacteria count using KNN: ", round(prediction,digits = 3)))
      }
      else if (model=="svmLinear"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Lactic acid bacteria count using svmLinear: ", round(prediction,digits = 3)))
      }
      else if (model=="svmRadial"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Lactic acid bacteria count using svmRadial: ", round(prediction,digits = 3)))
      }
      else if (model=="svmPoly"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Lactic acid bacteria count using svmPoly: ", round(prediction,digits = 3)))
      }
      else if (model=="lm"){
        prediction <- predict.bacCounts(bacteria, model, data, model.info, product, platform)
        return(paste0("Prediction of Lactic acid bacteria count using Linear Regression: ", round(prediction,digits = 3)))
      }
      else{
        print("Error occured: Please choose the correct model: | rf | knn | svmLinear | svmRadial | svmPoly | lm")
      }
    }
    else{
      print("Error occured: Please choose the correct bacteria: TVC | Ps | Bth | LAB")
    }
  } 
  else {
    print("Error occured: Please choose the correct product: CB | CTF")
  }
}


##################################################################################################################
##################################################################################################################
#* parse JSON
#* @post /predict
function(req, platform, product){
  #
  # Function to give an user repond for REST API request
  #
  # req = user request (a single row FTIR or MSI entry)
  # product = chicken burger (parameter: CB) or chicken thigh filet (parameter: CTF)
  # platform = analytical platform which was used for measurement: FTIR or MSI
  #
  # Hint - function requires "models.xlsx" file loaded as model.info variable 
  #        (model.info <- read_excel("models.xlsx"))
  #      - prospectr package is required
  #      - readxl package is required
  #      - jsonlite package is required
  #      - function.R script is required
  #        (source("functions.R"))
  #
  data <- jsonlite::fromJSON(req$postBody) #load single row entry
  
  if(platform == "MSI"){
    if(product == "CTF"){
      data <- data[,1:18]
      colnames(data) <- c("Mean.01", "Mean.02", "Mean.03", "Mean.04", "Mean.05", "Mean.06",
                          "Mean.07", "Mean.08", "Mean.09", "Mean.10", "Mean.11", "Mean.12",
                          "Mean.13", "Mean.14", "Mean.15", "Mean.16", "Mean.17", "Mean.18")
      predict.bestTVC(data, model.info, product, platform)
    }
    else{
      print("Error occured: For MSI platform shall adopt entries just for Chicken Thigh Fillet. Please choose parameter: product=CTF")
    }
  }
  #############################
  else if(platform == "FTIR"){
    if(product == "CB" | product == "CTF"){
      # round if this is required
      row_colnames<- colnames(data)
      row_col_integer <- 0
      for(i in 1:length(row_colnames)){
        if(as.numeric(row_colnames[i])==round(as.numeric(row_colnames[i]), 0)){
          row_col_integer <- row_col_integer+1
        }
        else{
          row_col_integer <- row_col_integer
        }
      }
      if(length(row_colnames) == row_col_integer){
        data<-data  
      } else{
        data <- data
        data <- roundWavelengths(data)
      }
      predict.bestTVC(data, model.info, product, platform)
    }
    else{
      print("Error occured: Please choose the correct product: CB or CTF")
    }
  }
  else{
    print("Error occured: Please choose the correct platform: MSI or FTIR")
  }
}

###################################################################################################################
###################################################################################################################

#* @serializer contentType list(type="application/pdf")
#* @get /report

function(platform, product) {
  
  if (platform=="MSI"){
    rmarkdown::render("MSI_Report.Rmd", output_format = "pdf_document")
  } 
  else if (platform=="FTIR"){
    if(product == "CB"){
      rmarkdown::render("FTIR_CB_Report.Rmd", output_format = "pdf_document")
    }
    else if(product == "CTF"){
      rmarkdown::render("FTIR_CTF_Report.Rmd", output_format = "pdf_document")
    }
    else {
      print("Error occured: Please type the correct product: CB or CTF")
    }
  }
  else{
    print("Error occured: Please type the correct platform: MSI or FTIR")
  }
  
  
  
}
