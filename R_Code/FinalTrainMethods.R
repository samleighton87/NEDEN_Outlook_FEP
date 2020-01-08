#native libraries
library(haven)
library(foreign)
library(doParallel)
library(dplyr)
library(readr)
library(labelled)
library(caret)
library(pROC)
library(clinfun)
library(gtools)
library(devtools)
library(mice)
library(corrplot)
library(car)
library(multcompView)
library(emmeans)
library(multcomp)
library(stats)
library(png)
library(PostcodesioR)
library(AppliedPredictiveModeling)
library(rmda)
library(pmsampsize)
#github libraries
library(CalibrationCurves)
library(dca)

#enable multicore (windows) which roughly halfs time for analysis runs
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)

options(max.print=1000000)

#don't use scientific notation (revert back with options(scipen=0)
options(scipen=999)
options(digits = 4)

setwd("/Users/sam_l/Desktop/Outlook/FinalPostReview/")

makeVlist <- function(dta) { 
  labels <- sapply(dta, function(x) attr(x, "label"))
  tibble(name = names(labels),
         label = labels)
}

#haven to factor
all_to_factor <- function(haven_x)
{
  for (i in seq(1:length(haven_x)))
  {
    if (!is.null(attr(haven_x[[i]], "labels")))
      haven_x[[i]] = to_factor(haven_x[[i]])
  }
  return(format.data.frame(haven_x))
}

#change char cols to factor
csv_to_factor <- function(imported_csv)
{
  cols_char_csv = colnames(imported_csv[, sapply(imported_csv, class) == 'character'])
  for (i in seq(1:length(cols_char_csv)))
  {
    imported_csv[[cols_char_csv[i]]] = as.factor(imported_csv[[cols_char_csv[i]]])
  }
  return(imported_csv)
}

#don't truncate list output of str
str_all <- function(dataframe)
{
  str(dataframe,list.len = ncol(dataframe))
}

#https://github.com/nogueirs/JMLR2018
getStability <- function(X,alpha=0.05) {
  ## the input X is a binary matrix of size M*d where:
  ## M is the number of bootstrap replicates
  ## d is the total number of features
  ## alpha is the level of significance (e.g. if alpha=0.05, we will get 95% confidence intervals)
  ## it's an optional argument and is set to 5% by default
  ### first we compute the stability
  
  M<-nrow(X)
  d<-ncol(X)
  hatPF<-colMeans(X)
  kbar<-sum(hatPF)
  v_rand=(kbar/d)*(1-kbar/d)
  stability<-1-(M/(M-1))*mean(hatPF*(1-hatPF))/v_rand ## this is the stability estimate
  
  ## then we compute the variance of the estimate
  ki<-rowSums(X)
  phi_i<-rep(0,M)
  for(i in 1:M){ 
    phi_i[i]<-(1/v_rand)*((1/d)*sum(X[i,]*hatPF)-(ki[i]*kbar)/d^2-(stability/2)*((2*kbar*ki[i])/d^2-ki[i]/d-kbar/d+1))
  }
  phi_bar=mean(phi_i)
  var_stab=(4/M^2)*sum((phi_i-phi_bar)^2) ## this is the variance of the stability estimate
  
  ## then we calculate lower and upper limits of the confidence intervals
  z<-qnorm(1-alpha/2) # this is the standard normal cumulative inverse at a level 1-alpha/2
  upper<-stability+z*sqrt(var_stab) ## the upper bound of the (1-alpha) confidence interval
  lower<-stability-z*sqrt(var_stab) ## the lower bound of the (1-alpha) confidence interval
  
  return(list("stability"=stability,"variance"=var_stab,"lower"=lower,"upper"=upper))
  
}

processData = function(dataset, outcomeVariable, siteOutcomeColumnNames, maxMissing = .25, site = "Site", removeMissing = T)
{
  datasetOutcome = dataset[ ,!(colnames(dataset) %in% siteOutcomeColumnNames)]
  #Dummy code (not outcome or site)
  dummies = dummyVars(~ ., data = datasetOutcome, fullRank = T)
  datasetOutcome <- data.frame(predict(dummies, newdata = datasetOutcome))
  #Add factor site and outcome back in (removing unused levels)
  datasetOutcome[, site] = droplevels(get(site, dataset))
  datasetOutcome[, outcomeVariable] = get(outcomeVariable, dataset)
  #remove na from outcome
  datasetOutcome = datasetOutcome[which(!is.na(datasetOutcome[, outcomeVariable])),]
  if (removeMissing)
  {
    #remove columns with more than 25% missing data default
    datasetOutcome = datasetOutcome[, colMeans(is.na(datasetOutcome)) <= maxMissing]
  }
  #remove zero and near zero variance columns
  nzv_cols <- nearZeroVar(datasetOutcome)
  if(length(nzv_cols) > 0) datasetOutcome <- datasetOutcome[, -nzv_cols]
}

#match up columns in yourData to otherData, return yourData
matchData = function(otherData, yourData)
{
  print(setdiff(colnames(yourData),colnames(otherData)))
  cols_not_in_otherData = setdiff(colnames(yourData),colnames(otherData))
  #remove
  for (i in seq(1:length(cols_not_in_otherData)))
  {
    yourData[[cols_not_in_otherData[i]]] = NULL
  }
  print(setdiff(colnames(yourData),colnames(otherData)))
  
  return(yourData)
}

#dataset must have site column - default to "Site"
#better to use multicore
#requires caret
nestedSiteCV = function(datasetWithSites, outcomeVariable, grid = T, tuneGrid, seed = 987, site = "Site", method = "glmnet", tuneLength = 10,
                        control = trainControl(method="cv", number=10, classProbs=TRUE, summaryFunction=twoClassSummary, 
                                               selectionFunction="best"), shrinkage = 1, colsToSave = NULL, multipleGLM = F)
{
  results = list()
  mods = list()
  siteColumn = which( colnames(datasetWithSites)==site )
  outcomeColumn = which( colnames(datasetWithSites)==outcomeVariable )
  
  #nested over number of Sites
  for(i in seq(1:length(levels(datasetWithSites[[siteColumn]]))))
  {
    set.seed(seed)
    
    #set up leave one site out cv
    test_set = datasetWithSites[ which(datasetWithSites[[siteColumn]] == levels(datasetWithSites[[siteColumn]])[i]), ]
    train_set = datasetWithSites[ -which(datasetWithSites[[siteColumn]] == levels(datasetWithSites[[siteColumn]])[i]), ]
    
    #remove site column
    test_set = test_set[ ,-siteColumn]
    train_set = train_set[ ,-siteColumn]
    
    #Get the observed outcome classes for this test set
    result = data.frame(obs=get(outcomeVariable, test_set)) 
    
    #train model default over grid of lambda alpha, knn impute, standardise
    if(grid)
    {
      mod <- train(as.formula(paste(outcomeVariable, "~ .")), data=train_set, method=method, metric="ROC", 
                   tuneGrid = tuneGrid, preProc = c("center", "scale","knnImpute"), trControl=control, 
                   na.action = na.pass)
    }else
    {
      mod <- train(as.formula(paste(outcomeVariable, "~ .")), data=train_set, method=method, metric="ROC", 
                   tuneLength = tuneLength, preProc = c("center", "scale","knnImpute"), trControl=control, 
                   na.action = na.pass)
      
      #Multiple imputation derived model
      if(multipleGLM)
      {
        train_set_imp = train_set[,-outcomeColumn]
        tempData <- mice(train_set_imp,m=10,seed=987)
        models = list()
        for (j in seq(1:tempData$m))
        {
          #Get imputed data
          train_imp = complete(tempData,j)
          #standardise the columns before building model
          preProcValues = preProcess(train_imp, method = c("center", "scale"))
          train_imp_stand = predict(preProcValues, train_imp)
          #Add factor outcome back in
          train_imp_stand[,outcomeVariable] = train_set[,outcomeVariable]
          models[[j]] = glm(as.formula(paste(outcomeVariable, "~ .")), data = train_imp_stand, family = "binomial")
        }
        #Pool results to get predictor estimates based on Rubin's rule and apply shrinkage from internal validation
        mod$finalModel$coefficients = summary(pool(models))[[1]]*shrinkage
      }else
      {
        mod$finalModel$coefficients = mod$finalModel$coefficients*shrinkage
      }
    }
    
    result$pred = predict(mod, test_set, type = "prob", na.action = na.pass)
    if (!is.null(colsToSave))
    {
      result[,colsToSave] = test_set[,colsToSave]
    }
    results[[i]] = result
    mods[[i]] = mod
    print(paste0("Site Round: ", i))
  }
  results_mods = list(results = results, mods = mods)
}

combineResultsInSequence = function(results_mods, predictYes = T, colsToSave = NULL)
{
  results_seq = NULL
  for (i in seq(1:length(results_mods[[1]])))
  {
    if(predictYes)
    {
      results_seq$pred = c(results_seq$pred, results_mods$results[[i]]$pred$Yes)
    }else
    {
      results_seq$pred = c(results_seq$pred, results_mods$results[[i]]$pred$No)
    }
    results_seq$obs = c(results_seq$obs, as.character(results_mods$results[[i]]$obs))
    
    if (!is.null(colsToSave))
    {
      results_seq[[colsToSave]] = c(results_seq[[colsToSave]], results_mods$results[[i]][,colsToSave])
    }
  }
  results_seq$obs = as.factor(results_seq$obs)
  return(results_seq)
}

permutationPValue = function(results_seq, auc_seq, seed = 987, predictYes=T, seq=T, dataSet = NULL, outcomeVariable = NULL)
{
  set.seed(seed)
  auc_null = NULL
  #significance level <0.0001
  if (!seq)
  {
    if (is.null(dataSet)|is.null(outcomeVariable))
    {
      print("You need to define a dataset and an outcome variable!")
      pPerm = NULL
    }else
    {
      for(i in seq (1:10001))
      {
        perm = permute(dataSet[,outcomeVariable])
        #set direction explicitly so not biased towards higher roc values (just makes it less likely for things to be significant, however)
        #https://www.rdocumentation.org/packages/pROC/versions/1.15.3/topics/roc
        if(predictYes)
        {
          auc_null = c(auc_null, roc(predictor = results_seq$Yes, response = perm, levels=c("No", "Yes"), direction="<")$auc)
        }else
        {
          auc_null = c(auc_null, roc(predictor = results_seq$No, response = perm, levels=c("No", "Yes"), direction=">")$auc)
        }
      }
      pPerm = (1+sum(auc_null >= auc_seq))/10001
    }
  }else
  {
    for(i in seq (1:10001))
    {
      perm = permute(results_seq$obs)
      #set direction explicitly so not biased towards higher roc values (just makes it less likely for things to be significant, however)
      #https://www.rdocumentation.org/packages/pROC/versions/1.15.3/topics/roc
      if(predictYes)
      {
        auc_null = c(auc_null, roc(predictor = results_seq$pred, response = perm, levels=c("No", "Yes"), direction="<")$auc)
      }else
      {
        auc_null = c(auc_null, roc(predictor = results_seq$pred, response = perm, levels=c("No", "Yes"), direction=">")$auc)
      }
    }
    pPerm = (1+sum(auc_null >= auc_seq))/10001    
  }
  #get p value by taking proportion of permutated values greater or equal to the actual value
  return(pPerm)
}

coefEvaluation = function(datasetWithSites, results_mods, site = "Site", isGLM = F)
{
  siteColumn = which( colnames(datasetWithSites)==site )
  coefs = NULL
  for (i in seq(1:length(levels(datasetWithSites[[siteColumn]]))))
  {
    if(isGLM)
    {
      coefs = rbind(coefs, coef(results_mods$mods[[i]]$finalModel))
    }else
    {
      coefs = c(coefs, coef(results_mods$mods[[i]]$finalModel, results_mods$mods[[i]]$bestTune$lambda))
    }
  }
  lengthC = NULL
  if (isGLM)
  {
    lengthC = length(coefs[1,])
  }
  else
  {
    lengthC = length(coefs[[1]])
  }
  
  #just get numbers
  coefs_extract = NULL
  for(i in seq(1:length(levels(datasetWithSites[[siteColumn]]))))
  {
    if(isGLM)
    {
      coefs_extract = coefs
    }else
    {
      coefs_extract = rbind(coefs_extract, coefs[[i]][1:lengthC])
    }
  }
  
  #get matrix of coefficients presence (1) or absence (0)
  #Presence or absence of predictors across all 14 LOSOCV models
  coefs_presence = NULL
  coefs_presence = coefs_extract[1:length(levels(datasetWithSites[[siteColumn]])),1:lengthC]
  coefs_presenceint = coefs_extract[1:length(levels(datasetWithSites[[siteColumn]])),1:lengthC]
  coefs_presence[coefs_presence != 0] <- 1
  coefs_presenceint[coefs_presenceint != 0] <- 1
  
  #stability of feature selection http://jmlr.org/papers/volume18/17-514/17-514.pdf
  print(getStability(coefs_presence))
  
  #get rank of coef by importance as in sports ranking
  coefs_rank = NULL
  
  for(i in seq(1:length(levels(datasetWithSites[[siteColumn]]))))
  {
    #rank absolute value excluding the intercept for each model
    coefs_rank = rbind(coefs_rank, rank(abs(coefs_extract[i,1:lengthC]), ties.method = "min"))
  }
  
  # rank the mean ranks of each column across all models
  coefs_rank_mean = colMeans(coefs_rank)
  
  #Invert order of rank to identify top models
  coefs_order = rank(-coefs_rank_mean)
  
  #Get the column names (not the intercept)
  if(isGLM)
  {
    coef_names = colnames(coefs[,1:lengthC])
  }else
  {
    coef_names = dimnames(coefs[[1]])[[1]][1:lengthC]
  }
  
  coefs_means = colMeans(coefs_extract)[1:lengthC]
  
  df = data.frame(coef_names, coefs_order, coefs_means, colMeans(coefs_presenceint))
}

smileyDiagram = function(probPercent)
{
  if(probPercent > 100 | probPercent < 0)
  {
    print("Enter a percentage probability between 0 and 100")
  }else
  {
    #Smiley face diagram
    #data with fixed ratio and random distribution
    probPercent = floor(probPercent + 0.5)
    x=sample(as.factor(rep(c("A","B"),c(as.integer(probPercent), as.integer(100-probPercent)))))
    # input parameters - nr * nc should equal length(x)
    
    #colours are actually only used for mapping - they just have to be called a valid colour
    #here blue is happy and red is sad but the colours only change if you change the pictures
    cols = NULL
    if(probPercent == 0)
    {
      cols <- c("red", "blue")
    }else
    {
      cols <- c("blue", "red")
    }
    
    nr <- 10
    nc <- 10
    # create data.frame of positions and colors
    m <- matrix(cols[x], nr, nc)
    DF <- data.frame(row = c(row(m)), col = c(col(m)[, nc:1]), value = c(m), 
                     stringsAsFactors = FALSE)
    
    # blank graph to insert man icons into
    xp <- 1.25
    plot(col ~ row, DF, col = DF$value, asp = 1,
         xlim = c(0, xp * nr), ylim = c(0, xp * nc),
         axes = FALSE, xlab = "", ylab = "", type = "n")
    #images must be in working directory
    smile_face<- readPNG("smile_face.png")
    sad_face<- readPNG("sad_face.png")
    
    G <- subset(transform(DF, row = xp * row, col = xp * col), value == "blue")
    with(G, rasterImage(sad_face,
                        row - .5, col - .5, row + .5, col + .5, 
                        xlim = c(0, xp * nr), ylim = c(0, xp * nc),
                        xlab = "", ylab = ""))
    
    R <- subset(transform(DF, row = xp * row, col = xp * col), value == "red")
    with(R, rasterImage(smile_face,
                        row - .5, col - .5, row + .5, col + .5, 
                        xlim = c(0, xp * nr), ylim = c(0, xp * nc),
                        xlab = "", ylab = ""))
  }
}

#takes postcode and 4 column table of ccg code, ccg name, raw imd score and standardised score
getCCGFromPostcode = function(postcode, ccgTable)
{
  letters_only <- function(x) !grepl("[^A-Za-z]", x)
  numbers_only <- function(x) !grepl("\\D", x)
  outCodeToPostCode = random_postcode(postcode)
  completePostcodeBool = postcode_validation(postcode)
  completePostCode = NULL
  CCGName = NULL
  DepZScore = NULL
  if(!is.null(outCodeToPostCode))
  {
    #is a valid outcode
    print("Valid partial postcode outcode")
    if(outCodeToPostCode$country != "England")
    {
      print("Sorry only valid for NHS England")
    }else
    {
      ccg = ccgTable[which(ccgTable[[1]] == outCodeToPostCode$codes$ccg),]
      print(paste0("Your postocde outcode corresponds to ",ccg[[2]]))
      print(paste0("The Deprivation Z-Score for ",ccg[[2]]," is ", ccg[[4]], "."))
      completePostCode = outCodeToPostCode$postcode
      CCGName = ccg[[2]]
      DepZScore = ccg[[4]]
    }
  }else if (completePostcodeBool)
  {
    print("Valid complete postcode")
    #is a valid complete postcode
    completePostCode = postcode_lookup(postcode)
    if(completePostCode$country != "England")
    {
      print("Sorry only valid for NHS England")
      completePostCode = NULL
    }else
    {
      ccg = ccgTable[which(ccgTable[[1]] ==  completePostCode$ccg_code),]
      print(paste0("Your postocde corresponds to ",ccg[[2]]))
      print(paste0("The Deprivation Z-Score for ",ccg[[2]]," is ", ccg[[4]], "."))
      CCGName = ccg[[2]]
      DepZScore = ccg[[4]]
    }
  }else if (!letters_only(postcode) & #not all letters
            !numbers_only(postcode) & #not all numbers
            #don't autocomplete partial outcodes 
            #so at least 3 characters if second character is a number (if it is a complete outcode it is handled above)
            ((nchar(postcode)>2 & numbers_only(strsplit(as.character(postcode), "")[[1]][2]))|
             #or at least 4 characters if second character is a letter (if it is a complete outcode it is handled above)
             (nchar(postcode)>3 & letters_only(strsplit(as.character(postcode), "")[[1]][2])))
  )
  {
    tryCatch({#not an outcode but partial postcode greater than outcode
      autoCompletePostCode = postcode_autocomplete(postcode, 1)[[1]]
      completePostCode = postcode_lookup(autoCompletePostCode)
      print("Valid partial postcode longer than outcode")
      if(completePostCode$country != "England")
      {
        print("Sorry only valid for NHS England")
        completePostCode = NULL
      }else
      {
        ccg = ccgTable[which(ccgTable[[1]] ==  completePostCode$ccg_code),]
        print(paste0("Your partial postocde corresponds to ",ccg[[2]]))
        print(paste0("The Deprivation Z-Score for ",ccg[[2]]," is ", ccg[[4]], "."))
        CCGName = ccg[[2]]
        DepZScore = ccg[[4]]
      }
    },
    #not a postcode
    error=function(error_message) {
      print("Neither a valid complete postcode nor a valid partial postcode")
    })
  }else
  {
    #too short and meaningless
    print("Please enter at least a valid postcode outcode")
  }
  if(is.null(completePostCode))
  {
    print("Nothing to return")
  }else
  {
    x = list("completePostCode" = completePostCode, "CCGName" = CCGName, "DepZScore" = DepZScore)
    return(x)
  }
}

#Custom method for internal validation glm only
customMultipleImputeCV = function(dataset, #dataset with outcome
                            outcomeVariable, #outcome variable name String
                            repeats = 10, #defaults to 10 repeats as per steyerberg
                            cv = 10 , #defaults to 10 folds
                            control = trainControl(method="none", classProbs=TRUE, summaryFunction=twoClassSummary), # no tuning
                            preProcess = c("center", "scale","knnImpute"),#what preprocessing for model testing
                            multipleGLM = T,#multiple imputation for model building
                            shrinkage = 1)#any shrinkage - generally no this is to find amount shrinkage 
{
  resultsOuter = list()
  modsOuter = list()
  outcomeColumn = which( colnames(dataset)==outcomeVariable )
  
  for(j in seq(1:repeats)) 
  {
    #make replicable
    #
    set.seed(j)
    #createFolds splits the data into k groups (defaults to 10 groups, & as list) 
    #when returnTrain = TRUE, the values returned are the sample positions corresponding to the data used during training, 
    #returns a list or matrix of row position integers corresponding to the training data
    #
    splits <- createFolds(dataset[,outcomeVariable], returnTrain = TRUE, k = cv)
    
    #lapply returns a list of the same length as splits to results, each element of which is the result of applying function to the corresponding element of splits, 
    #the result of the function is the data.frame created for each split. This sorts out the ordering problem when collating the results I assume
    #split number (row position integer) becomes x in function
    #holdout is a vector of numbers from 1 to the number of rows in dat except those with where the row numbers are in x (no duplicates or attribute names due to unique() ). 
    #N.B. the rows in x are those used in the training data with the held out removed (returnTrain is True) so this just recreates held out
    #the data.frame is the combination of the index vector (the row numbers for the holdout data) and the obs vector (values corresponding to the index)
    #
    results = lapply(splits, function(x, dat) 
    {
      holdout <- (1:nrow(dat))[-unique(x)]
      data.frame(index = holdout, obs = dat[,outcomeVariable][holdout])
    },
    dat = dataset)
    
    #a vector to hold the models as 10 different ones will be created from the 10 splits
    #
    mods = vector(mode = "list", length = length(splits))
    
    #seq generates the sequence 1, 2, ..., length(splits) for each of the 10 separate 10ths of the data
    #having two square braces only gives the dataframe column
    #unclear why need to reset seed but Max says so
    #create your model with the 9/10 of data with the 1/10 removed (this will happen 100 times in the manual cv loop. 
    # a column of predictions (pred) created by the testing the model on the held out 1/10 is added to the corresponding vector of the same indices and the actual observations (obs)
    #model added to mods vector for use later
    #
    print(paste0("Repeat: ",j,"/",repeats))
    for(i in seq(along = splits)) 
    {
      in_train <- unique(splits[[i]])
      set.seed(j+1)
      mod <- train(as.formula(paste(outcomeVariable, "~ .")), data=dataset[in_train, ], method="glm", metric="ROC", preProc=preProcess, 
                   trControl=control, na.action=na.pass)
      
      #Multiple imputation derived model
      if(multipleGLM)
      {
        train_set = dataset[in_train, ]
        train_set_imp = train_set[,-outcomeColumn]
        tempData <- mice(train_set_imp,m=10,seed=987)
        models = list()
        for (j in seq(1:tempData$m))
        {
          #Get imputed data
          train_imp = complete(tempData,j)
          #standardise the columns before building model
          preProcValues = preProcess(train_imp, method = c("center", "scale"))
          train_imp_stand = predict(preProcValues, train_imp)
          #Add factor outcome back in
          train_imp_stand[,outcomeVariable] = train_set[,outcomeVariable]
          models[[j]] = glm(as.formula(paste(outcomeVariable, "~ .")), data = train_imp_stand, family = "binomial")
        }
        #Pool results to get predictor estimates based on Rubin's rule and apply shrinkage from internal validation
        mod$finalModel$coefficients = summary(pool(models))[[1]]*shrinkage
      }else
      {
        mod$finalModel$coefficients = mod$finalModel$coefficients*shrinkage
      }
      
      results[[i]]$pred = predict(mod, dataset[-in_train, ], type = "prob", na.action = na.pass)
      
      mods[[i]] = mod
      print(paste0("Cross Validation: ", i,"/", cv))
    }
    resultsOuter[[j]] = results
    modsOuter[[j]] = mods
  }
  results_mods = list(results = resultsOuter, mods = modsOuter)  
}

#Custom method to combine the results in a sequence in order to make a ROC curve for the whole lot
combineResultsInSequenceNested = function(results, predictYes = F)
{
  results_seq = NULL
  for (j in seq(1:length(results)))
  {
    for (i in seq(1:length(results[[j]])))
    {
      if(predictYes)
      {
        results_seq$pred = c(results_seq$pred, results[[j]][[i]]$pred$Yes)
      }else
      {
        results_seq$pred = c(results_seq$pred, results[[j]][[i]]$pred$No)
      }
      results_seq$obs = c(results_seq$obs, as.character(results[[j]][[i]]$obs))
    }
  }
  results_seq$obs = as.factor(results_seq$obs)
  return(results_seq)
}

#sample size
pmsampsize(type = "b",rsquared = 0.25, parameters = 20, shrinkage = 0.9, prevalence = 0.5)

#Impute 12 month outcome for EDEN and Outlook datasets separately, using only PANSS at 6 months and 12 months.
#EDEN
eden_6_12 = read_csv("EDEN_6_12_PANSS.csv")
eden_6_12 = csv_to_factor(eden_6_12)
eden_6_12$row = rownames(eden_6_12)
#remove empty rows (not including coloumn 64 which is the row number and never empty)
eden_6_12 = eden_6_12[rowSums(is.na(eden_6_12[,-64])) != ncol(eden_6_12[,-64]),]
#Impute the missing data (not using row number column 64)
tempData_eden_6_12 = mice(eden_6_12[,-64], m=1, seed=987)
#Take the imputed outcome
complete_eden_6_12 = complete(tempData_eden_6_12,1)
#preserve original rownames
rownames(complete_eden_6_12) = eden_6_12$row

#Outlook
outlook_6_12 = read_csv("Outlook_6_12_PANSS.csv")
outlook_6_12 = csv_to_factor(outlook_6_12)
outlook_6_12$row = rownames(outlook_6_12)
#remove empty rows (not including coloumn 64 which is the row number and never empty)
outlook_6_12 = outlook_6_12[rowSums(is.na(outlook_6_12[,-64])) != ncol(outlook_6_12[,-64]),]
#Impute the missing data (not using row number column 64)
tempData_outlook_6_12 = mice(outlook_6_12[,-64], m=1, seed=987)
#Take the imputed outcome
complete_outlook_6_12 = complete(tempData_outlook_6_12,1)
#preserve original rownames
rownames(complete_outlook_6_12) = outlook_6_12$row

#Outlookall = read_spss("OutlookData/PSYGRID 1E.sav")
#Outlookall_labels = makeVlist(Outlookall)
#write_csv(all_to_factor(Outlookall), path = "Outlookall.csv") #Have to manually remove "NA"s from csv file
#write_csv(format.data.frame(Outlookall_labels), path = "Outlookall_labels.csv")

eden_outlook_preproc = read_csv("outlook_eden_preproc_rem.csv")
eden_outlook_preproc_factor = csv_to_factor(eden_outlook_preproc)

eden_all = eden_outlook_preproc_factor[ which(eden_outlook_preproc_factor$Study == "EDEN"),!(colnames(eden_outlook_preproc_factor) %in% c("Study")) ]
outlook_all = eden_outlook_preproc_factor[ which(eden_outlook_preproc_factor$Study == "Outlook"),!(colnames(eden_outlook_preproc_factor) %in% c("Study")) ]

#Choose Final predictors based on expert knowledge
eden_all_final = eden_all[rownames(complete_eden_6_12),-c(3:9,17,19:24,113:119)]
outlook_all_final = outlook_all[rownames(complete_outlook_6_12),-c(3:9,17,19:24,113:119)]
eden_all_final_old = eden_all_final
outlook_all_final_old = outlook_all_final

#Replace missing outcome column with imputed outcome column from above
eden_all_final$M12_PANSS_Period_Rem = complete_eden_6_12$M12_PANSS_Period_Rem
outlook_all_final$M12_PANSS_Period_Rem = complete_outlook_6_12$M12_PANSS_Period_Rem

#Baseline stats
#summary(aov(y~group, data = data.frame(group=factor(rep(1:4, c(1027,901,399,278))),y=c(eden_all$ADJ_DUP,
#                                                                                           eden_all_final$ADJ_DUP,
#                                                                                           outlook_all$ADJ_DUP,
#                                                                                           outlook_all_final$ADJ_DUP))))

#chisq.test(as.table(rbind(c(245,399,262,98),c(216,347,228,91),c(89,130,92,69),c(58,96,68,47))), correct = F)

#tune over three separate grids of static alpha and varied lambda to allow us to vary selection criteria by tolerance
#which is based on lambda only
tuneGrid.1 = expand.grid(alpha=0.1, lambda = 10 ^ seq(-0.3, -5, length = 100))
tuneGrid.5 = expand.grid(alpha=0.5, lambda = 10 ^ seq(-0.3, -5, length = 100))
tuneGrid.9 = expand.grid(alpha=0.9, lambda = 10 ^ seq(-0.3, -5, length = 100))

#Remission Data set up
eden_all_final_Rem = processData(dataset = eden_all_final, outcomeVariable = "M12_PANSS_Period_Rem",
                           siteOutcomeColumnNames = c("Site","M12_PANSS_Period_Rem"))
#Don't have 25% cut off for test dataset predictors
outlook_all_final_Rem = processData(dataset = outlook_all_final, outcomeVariable = "M12_PANSS_Period_Rem", removeMissing = F,
                              siteOutcomeColumnNames = c("Site","M12_PANSS_Period_Rem"))
#same predictor columns in both datasets
eden_all_final_Rem = matchData(otherData = outlook_all_final_Rem, yourData = eden_all_final_Rem)
outlook_all_final_Rem = matchData(otherData = eden_all_final_Rem, yourData = outlook_all_final_Rem)

#Tune over three ranges an alpha (0.1, 0.5 & 0.9) and 100 lambdas - native caret doesn't work correctly for tuning here
results_mods_final_Rem.1 = nestedSiteCV(datasetWithSites = eden_all_final_Rem, outcomeVariable = "M12_PANSS_Period_Rem", tuneGrid = tuneGrid.1)
results_mods_final_Rem.5 = nestedSiteCV(datasetWithSites = eden_all_final_Rem, outcomeVariable = "M12_PANSS_Period_Rem", tuneGrid = tuneGrid.5)
results_mods_final_Rem.9 = nestedSiteCV(datasetWithSites = eden_all_final_Rem, outcomeVariable = "M12_PANSS_Period_Rem", tuneGrid = tuneGrid.9)

results_seq_final_Rem.1 = combineResultsInSequence(results_mods_final_Rem.1, predictYes = F)
results_seq_final_Rem.5 = combineResultsInSequence(results_mods_final_Rem.5, predictYes = F)
results_seq_final_Rem.9 = combineResultsInSequence(results_mods_final_Rem.9, predictYes = F)

#which has best combination of discrimation, predictor stability & sparsity?
roc_seq_final_Rem.1 = roc(predictor = results_seq_final_Rem.9$pred, response = results_seq_final_Rem.9$obs, ci = T, levels=c("No", "Yes"), direction=">")
auc_seq_final_Rem.1 = roc_seq_final_Rem.1$auc
roc_seq_final_Rem.1 #AUC 0.672
df_coefs_final_Rem.1 = coefEvaluation(datasetWithSites = eden_all_final_Rem, results_mods = results_mods_final_Rem.1)
View(df_coefs_final_Rem.1) #stability 0.5775
roc_seq_final_Rem.5 = roc(predictor = results_seq_final_Rem.9$pred, response = results_seq_final_Rem.9$obs, ci = T, levels=c("No", "Yes"), direction=">")
auc_seq_final_Rem.5 = roc_seq_final_Rem.5$auc
roc_seq_final_Rem.5 #AUC 0.672
df_coefs_final_Rem.5 = coefEvaluation(datasetWithSites = eden_all_final_Rem, results_mods = results_mods_final_Rem.5)
View(df_coefs_final_Rem.5) #stability 0.6435
roc_seq_final_Rem.9 = roc(predictor = results_seq_final_Rem.9$pred, response = results_seq_final_Rem.9$obs, ci = T, levels=c("No", "Yes"), direction=">")
auc_seq_final_Rem.9 = roc_seq_final_Rem.9$auc
roc_seq_final_Rem.9 #AUC 0.672
df_coefs_final_Rem.9 = coefEvaluation(datasetWithSites = eden_all_final_Rem, results_mods = results_mods_final_Rem.9)
View(df_coefs_final_Rem.9) #stability 0.6137
#results_mods_final_Rem.5 has best combination of discrimination and stability
val.prob.ci.2(p=results_seq_final_Rem.5$pred, y=results_seq_final_Rem.5$obs=="No", g=10, logistic.cal = T, lty.log=9,
              col.log="red", lwd.log=1.5, col.ideal="blue", lwd.ideal=0.5, dostats = T)

#Build GLM with predictors present across all 14 LOSOCV
#Tune coefficient inclusion by percentage presence - 100%, 90%, 80%, 70%, 50% - whichever performs best at internal validation
#(Steyerberg uses 50% in https://ajp.psychiatryonline.org/doi/10.1176/appi.ajp.2018.18050566)
#this further reduces the degrees of freedom
View(df_coefs_final_Rem.5)
#select coefficients at different perentage cut offs here
#80% is best at internal validation
eden_all_final_Rem_outlook = eden_all_final_Rem[,c(8,35,25,32,29,27,19,2,88,43,34,36,13,93,3,55,86,78,1,70)] #not outcome column as imputing first
#Multiple Imputation to Construct Final GLM model - not based on outcome column
#10 datasets
tempData_eden_outlook <- mice(eden_all_final_Rem_outlook,m=10,seed=987)
finalModels = list()
for (i in seq(1:tempData_eden_outlook$m))
{
 #Get imputed data
 eden_outlook_imp = complete(tempData_eden_outlook,i)
 #standardise the columns before building model
 preProcValues = preProcess(eden_outlook_imp, method = c("center", "scale"))
 eden_outlook_imp_stand = predict(preProcValues, eden_outlook_imp)
 #Add factor outcome back in
 eden_outlook_imp_stand$M12_PANSS_Period_Rem = eden_all_final_Rem$M12_PANSS_Period_Rem
 finalModels[[i]] = glm(M12_PANSS_Period_Rem ~ ., data = eden_outlook_imp_stand, family = "binomial")
}
#Pool results to get predictor estimates based on Rubin's rule (coefficients reversed as using to predict non-remission)
View(-summary(pool(finalModels), conf.int = T, exponentiate = F, conf.level = 0.95)[,c(1,6,7)])
View(1/summary(pool(finalModels), conf.int = T, exponentiate = T, conf.level = 0.95)[,c(1,6,7)])

#Internal Validation to determine best percentage presence of predictors
#Now 10x repeated 10x cross-validation for internal validation with mulitple imputation
#but not bootstrap as overly optimistic as per https://www.r-bloggers.com/part-2-optimism-corrected-bootstrapping-is-definitely-bias-further-evidence/
#https://stats.stackexchange.com/a/46344/114271
#https://intobioinformatics.wordpress.com/2018/12/25/optimism-corrected-bootstrapping-a-problematic-method/
#Add outcome back in
eden_all_final_Rem_outlook$M12_PANSS_Period_Rem = eden_all_final_Rem$M12_PANSS_Period_Rem
#Takes over 10 mins
eden_all_final_Rem_outlook_glm_internal_mods_results = customMultipleImputeCV(dataset = eden_all_final_Rem_outlook, outcomeVariable = "M12_PANSS_Period_Rem")
eden_all_final_Rem_outlook_glm_internal_results_seq = combineResultsInSequenceNested(eden_all_final_Rem_outlook_glm_internal_mods_results$results, predictYes = F)
#Check Internal Validation Calibration and Discrimination
pdf("figure_2.pdf", width = 7, height = 7)
val.prob.ci.2(p=eden_all_final_Rem_outlook_glm_internal_results_seq$pred, y=eden_all_final_Rem_outlook_glm_internal_results_seq$obs=="No", g=10, logistic.cal = T, lty.log=9,
              col.log="red", lwd.log=1.5, col.ideal="blue", lwd.ideal=0.5, dostats = T)
#AUC 0.7176573 for 50%,  0.7210798 for 70%, 0.7237540 for 80%, 0.721707 for 90%, 0.720897 for 100%
dev.off()
eden_all_final_Rem_outlook_glm_internal_results_seq_roc = roc(predictor = eden_all_final_Rem_outlook_glm_internal_results_seq$pred, 
                                                              response = eden_all_final_Rem_outlook_glm_internal_results_seq$obs, ci = T, levels=c("No", "Yes"), direction=">")
eden_all_final_Rem_outlook_glm_internal_results_seq_auc = eden_all_final_Rem_outlook_glm_internal_results_seq_roc$auc
eden_all_final_Rem_outlook_glm_internal_results_seq_roc
eden_all_final_Rem_outlook_glm_internal_results_seq_auc_p = permutationPValue(results_seq = eden_all_final_Rem_outlook_glm_internal_results_seq, 
                                                                              auc_seq = eden_all_final_Rem_outlook_glm_internal_results_seq_auc, predictYes = F)
#multiply by shrinkage factor
View(-summary(pool(finalModels), conf.int = T, exponentiate = F, conf.level = 0.95)[,c(1,6,7)]* 0.8727967)
View(exp(-summary(pool(finalModels), conf.int = T, exponentiate = F, conf.level = 0.95)[,c(1,6,7)]* 0.8727967))

#Internal-External Validation
#get a true nested leave one-site out performance metric (adding outcome back in) - previous just does random nested cv
eden_all_final_Rem_outlook$Site = eden_all_final_Rem$Site
#with shrinkage from internal and multiple imputation for model development
eden_all_final_Rem_outlook_glm_LOSOCV_mods_results = nestedSiteCV(datasetWithSites = eden_all_final_Rem_outlook, outcomeVariable = "M12_PANSS_Period_Rem", grid = F, method = "glm", 
                                                                  control = trainControl(method="none", classProbs=TRUE, summaryFunction=twoClassSummary), shrinkage =  0.8727967,
                                                                  colsToSave = "ADJ_DUP", multipleGLM = T)
eden_all_final_Rem_outlook_glm_LOSOCV_results_seq = combineResultsInSequence(eden_all_final_Rem_outlook_glm_LOSOCV_mods_results, predictYes = F,colsToSave = "ADJ_DUP")
#Calibration Curve Internal-External
pdf("figure_3.pdf", width = 7, height = 7)
val.prob.ci.2(p=eden_all_final_Rem_outlook_glm_LOSOCV_results_seq$pred, y=eden_all_final_Rem_outlook_glm_LOSOCV_results_seq$obs=="No", g=10, logistic.cal = T, lty.log=9,
              col.log="red", lwd.log=1.5, col.ideal="blue", lwd.ideal=0.5, dostats = T)
#0.71 vs 0.71 vs 0.72 vs 0.71 vs 0.71
dev.off()
eden_all_final_Rem_outlook_glm_LOSOCV_results_seq_roc = roc(predictor = eden_all_final_Rem_outlook_glm_LOSOCV_results_seq$pred, 
                                                              response = eden_all_final_Rem_outlook_glm_LOSOCV_results_seq$obs, ci = T, levels=c("No", "Yes"), direction=">")
eden_all_final_Rem_outlook_glm_LOSOCV_results_seq_auc = eden_all_final_Rem_outlook_glm_LOSOCV_results_seq_roc$auc
eden_all_final_Rem_outlook_glm_LOSOCV_results_seq_roc
eden_all_final_Rem_outlook_glm_LOSOCV_results_seq_auc_p = permutationPValue(results_seq = eden_all_final_Rem_outlook_glm_LOSOCV_results_seq, 
                                                                              auc_seq = eden_all_final_Rem_outlook_glm_LOSOCV_results_seq_auc, predictYes = F)
#Decision Curve Analysis Internal-External
dca_glm_Rem_int_ext = NULL
dca_glm_Rem_int_ext$M12_PANSS_Period_Rem = as.integer(as.character(eden_all_final_Rem_outlook_glm_LOSOCV_results_seq$obs)=="No")
dca_glm_Rem_int_ext$Model = eden_all_final_Rem_outlook_glm_LOSOCV_results_seq$pred
dca_glm_Rem_int_ext$ADJ_DUP = eden_all_final_Rem_outlook_glm_LOSOCV_results_seq$ADJ_DUP
#get data across all thresholds
pdf("figure_4.pdf", width = 7, height = 7)
dca(data = as.data.frame(dca_glm_Rem_int_ext), outcome = "M12_PANSS_Period_Rem", predictors = c("Model","ADJ_DUP"), smooth = "TRUE", loess.span = 0.25, probability = c(TRUE, FALSE), graph = T)
dev.off()
#Remove site as not needed for external validation
eden_all_final_Rem_outlook$Site = NULL

#External Validation
outlook_all_final_Rem_eden = matchData(otherData = eden_all_final_Rem_outlook, yourData = outlook_all_final_Rem)
#change coefficients glm from caret to match those from pooled multiple imputation model * shrinkage factor from internal validation
#But continue to use knn imputation rules from caret model object from training so as not to be overly optimistic (as would result if imputed test set separately)
set.seed(987)
mod_eden_all_final_Rem_outlook_glm = train(M12_PANSS_Period_Rem ~ ., data=eden_all_final_Rem_outlook, method="glm", metric="ROC", preProc = c("center", "scale","knnImpute"), 
                                            trControl = trainControl(method="none", classProbs=TRUE, summaryFunction=twoClassSummary), na.action = na.pass)
mod_eden_all_final_Rem_outlook_glm_shrink = mod_eden_all_final_Rem_outlook_glm
mod_eden_all_final_Rem_outlook_glm_shrink$finalModel$coefficients = summary(pool(finalModels))[[1]]* 0.8727967
outlook_all_final_Rem_glm_result_shrink = predict(mod_eden_all_final_Rem_outlook_glm_shrink, outlook_all_final_Rem_eden, type = "prob", na.action = na.pass)
pdf("figure_5.pdf", width = 7, height = 7)
val.prob.ci.2(p=outlook_all_final_Rem_glm_result_shrink$No, y=outlook_all_final_Rem_eden$M12_PANSS_Period_Rem=="No", g=10, logistic.cal = T, lty.log=9,
              col.log="red", lwd.log=1.5, col.ideal="blue", lwd.ideal=0.5)
dev.off()
outlook_all_final_Rem_glm_result_shrink_roc = roc(predictor = outlook_all_final_Rem_glm_result_shrink$No, 
                                                            response = outlook_all_final_Rem_eden$M12_PANSS_Period_Rem, ci = T, levels=c("No", "Yes"), direction=">")
outlook_all_final_Rem_glm_result_shrink_auc = outlook_all_final_Rem_glm_result_shrink_roc$auc
outlook_all_final_Rem_glm_result_shrink_roc
outlook_all_final_Rem_glm_result_shrink_auc_p = permutationPValue(results_seq = outlook_all_final_Rem_glm_result_shrink, seq = F, dataSet = outlook_all_final_Rem_eden,
                                                                            auc_seq = outlook_all_final_Rem_glm_result_shrink_auc, predictYes = F, 
                                                                            outcomeVariable = "M12_PANSS_Period_Rem")
#Decision Curve Analysis External
dca_glm_Rem_ext = NULL
dca_glm_Rem_ext$M12_PANSS_Period_Rem = as.integer(as.character(outlook_all_final_Rem_eden$M12_PANSS_Period_Rem)=="No")
dca_glm_Rem_ext$Model = outlook_all_final_Rem_glm_result_shrink$No
dca_glm_Rem_ext$ADJ_DUP = outlook_all_final_Rem_eden$ADJ_DUP
#get data across all thresholds
pdf("figure_6.pdf", width = 7, height = 7)
dca(data = as.data.frame(dca_glm_Rem_ext), outcome = "M12_PANSS_Period_Rem", predictors = c("Model","ADJ_DUP"), smooth = "TRUE", loess.span = 0.35, probability = c(TRUE, FALSE), graph = T, xstart = 0.32, xstop = 0.725)
dev.off()

##

#for online model, refit with prestandardised local concentration on based on whole range of values from IMD 2007 just to get correct preprocess info
#gives identical coefficients
write_csv(eden_all_final_Rem_outlook,"unstandardised_dep.csv", na = "")
#replace local concentration with population standardised
eden_all_final_Rem_outlook_new = read_csv("standardised_dep.csv")
eden_all_final_Rem_outlook_new$M12_PANSS_Period_Rem = as.factor(eden_all_final_Rem_outlook_new$M12_PANSS_Period_Rem)
set.seed(987)
mod_eden_all_final_Rem_outlook_glm_new = train(M12_PANSS_Period_Rem ~ ., data=eden_all_final_Rem_outlook_new, method="glm", metric="ROC", preProc = c("center", "scale","knnImpute"), 
                                               trControl = trainControl(method="none", classProbs=TRUE, summaryFunction=twoClassSummary), na.action = na.pass)
mod_eden_all_final_Rem_outlook_glm_new_shrink = mod_eden_all_final_Rem_outlook_glm_new
mod_eden_all_final_Rem_outlook_glm_new_shrink$finalModel$coefficients = summary(pool(finalModels))[[1]]*0.8727967
#remove training data for sharing model
mod_eden_all_final_Rem_outlook_glm_new_shrink$trainingData = NULL
save(mod_eden_all_final_Rem_outlook_glm_new_shrink, file = "mod_eden_all_final_Rem_outlook_glm_new_shrink.rda")
