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
library(pmsampsize)
#github libraries
library(CalibrationCurves)
library(dca)

#enable multicore (windows) which roughly halfs time for analysis runs
cl <- makeCluster(detectCores(), type = 'PSOCK')
registerDoParallel(cl)

options(max.print = 1000000)

#don't use scientific notation (revert back with options(scipen=0)
options(scipen = 999)
options(digits = 4)

makeVlist <- function(dta) {
  labels <- sapply(dta, function(x)
    attr(x, "label"))
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
  str(dataframe, list.len = ncol(dataframe))
}

#https://github.com/nogueirs/JMLR2018
getStability <- function(X, alpha = 0.05) {
  ## the input X is a binary matrix of size M*d where:
  ## M is the number of bootstrap replicates
  ## d is the total number of features
  ## alpha is the level of significance (e.g. if alpha=0.05, we will get 95% confidence intervals)
  ## it's an optional argument and is set to 5% by default
  ### first we compute the stability
  
  M <- nrow(X)
  d <- ncol(X)
  hatPF <- colMeans(X)
  kbar <- sum(hatPF)
  v_rand = (kbar / d) * (1 - kbar / d)
  stability <-
    1 - (M / (M - 1)) * mean(hatPF * (1 - hatPF)) / v_rand ## this is the stability estimate
  
  ## then we compute the variance of the estimate
  ki <- rowSums(X)
  phi_i <- rep(0, M)
  for (i in 1:M) {
    phi_i[i] <-
      (1 / v_rand) * ((1 / d) * sum(X[i,] * hatPF) - (ki[i] * kbar) / d ^ 2 -
                        (stability / 2) * ((2 * kbar * ki[i]) / d ^ 2 - ki[i] / d - kbar / d + 1))
  }
  phi_bar = mean(phi_i)
  var_stab = (4 / M ^ 2) * sum((phi_i - phi_bar) ^ 2) ## this is the variance of the stability estimate
  
  ## then we calculate lower and upper limits of the confidence intervals
  z <-
    qnorm(1 - alpha / 2) # this is the standard normal cumulative inverse at a level 1-alpha/2
  upper <-
    stability + z * sqrt(var_stab) ## the upper bound of the (1-alpha) confidence interval
  lower <-
    stability - z * sqrt(var_stab) ## the lower bound of the (1-alpha) confidence interval
  
  return(list(
    "stability" = stability,
    "variance" = var_stab,
    "lower" = lower,
    "upper" = upper
  ))
  
}

processData = function(dataset,
                       outcomeVariable,
                       siteOutcomeColumnNames,
                       maxMissing = .25,
                       site = "Site",
                       removeMissing = T)
{
  datasetOutcome = dataset[,!(colnames(dataset) %in% siteOutcomeColumnNames)]
  #Dummy code (not outcome or site)
  dummies = dummyVars(~ ., data = datasetOutcome, fullRank = T)
  datasetOutcome <-
    data.frame(predict(dummies, newdata = datasetOutcome))
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
  if (length(nzv_cols) > 0)
    datasetOutcome <- datasetOutcome[, -nzv_cols]
  return(datasetOutcome)
}

#match up columns in yourData to otherData, return yourData
matchData = function(otherData, yourData)
{
  print(setdiff(colnames(yourData), colnames(otherData)))
  cols_not_in_otherData = setdiff(colnames(yourData), colnames(otherData))
  #remove
  for (i in seq(1:length(cols_not_in_otherData)))
  {
    yourData[[cols_not_in_otherData[i]]] = NULL
  }
  print(setdiff(colnames(yourData), colnames(otherData)))
  
  return(yourData)
}

#dataset must have site column - default to "Site"
#better to use multicore
#requires caret
nestedSiteCVUpdate = function(datasetWithSites,
                              outcomeVariable,
                              tuneGrid,
                              seed = 987,
                              site = "Site",
                              method = "glmnet",
                              nestedSite = TRUE,
                              control = trainControl(
                                method = "cv",
                                number = 10,
                                classProbs = TRUE,
                                summaryFunction = twoClassSummary,
                                selectionFunction =
                                  "best"
                              ),
                              preProc = c("center", "scale", "knnImpute"),
                              colsToSave = NULL,
                              repeats = 10,
                              cv = 10)
{
  resultsOuter = list()
  modsOuter = list()
  siteColumn = which(colnames(datasetWithSites) == site)
  outcomeColumn = which(colnames(datasetWithSites) == outcomeVariable)
  
  if (nestedSite)
  {
    #nested over number of Sites
    for (i in seq(1:length(levels(datasetWithSites[[siteColumn]]))))
    {
      set.seed(seed)
      
      #set up leave one site out cv
      test_set = datasetWithSites[which(datasetWithSites[[siteColumn]] == levels(datasetWithSites[[siteColumn]])[i]), ]
      train_set = datasetWithSites[-which(datasetWithSites[[siteColumn]] == levels(datasetWithSites[[siteColumn]])[i]), ]
      
      #remove site column
      test_set = test_set[,-siteColumn]
      train_set = train_set[,-siteColumn]
      
      #Get the observed outcome classes for this test set
      result = data.frame(obs = get(outcomeVariable, test_set))
      
      #train model default over grid of lambda alpha, knn impute, standardise
      mod <-
        train(
          as.formula(paste(outcomeVariable, "~ .")),
          data = train_set,
          method = method,
          metric = "ROC",
          tuneGrid = tuneGrid,
          preProc = preProc,
          trControl = control,
          na.action = na.pass
        )
      
      
      result$pred = predict(mod, test_set, type = "prob", na.action = na.pass)
      if (!is.null(colsToSave))
      {
        result[, colsToSave] = test_set[, colsToSave]
      }
      resultsOuter[[i]] = result
      modsOuter[[i]] = mod
      print(paste0("Site Round: ", i))
    }
  } else
    #nested CV
  {
    #fit on whole dataset
    dataset = datasetWithSites[,-siteColumn]
    for (j in seq(1:repeats))
    {
      #make replicable
      #
      set.seed(j)
      #createFolds splits the data into k groups (defaults to 10 groups, & as list)
      #when returnTrain = TRUE, the values returned are the sample positions corresponding to the data used during training,
      #returns a list or matrix of row position integers corresponding to the training data
      #
      splits <-
        createFolds(dataset[, outcomeVariable], returnTrain = TRUE, k = cv)

      #lapply returns a list of the same length as splits to results, each element of which is the result of applying function to the corresponding element of splits,
      #the result of the function is the data.frame created for each split. This sorts out the ordering problem when collating the results I assume
      #split number (row position integer) becomes x in function
      #holdout is a vector of numbers from 1 to the number of rows in dat except those with where the row numbers are in x (no duplicates or attribute names due to unique() ).
      #N.B. the rows in x are those used in the training data with the held out removed (returnTrain is True) so this just recreates held out
      #the data.frame is the combination of the index vector (the row numbers for the holdout data) and the obs vector (values corresponding to the index)
      #
      results = lapply(splits, function(x, dat)
      {
        #holdout not in training
        holdout <- (1:nrow(dat))[-unique(x)]
        #test data obs
        data.frame(index = holdout, obs = dat[, outcomeVariable][holdout])
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
      print(paste0("Repeat: ", j, "/", repeats))
      for (i in seq(along = splits))
      {
        in_train <- unique(splits[[i]])
        set.seed(j + 1)
        mod <-
          train(
            as.formula(paste(outcomeVariable, "~ .")),
            data = dataset[in_train, ],
            method = "glmnet",
            metric = "ROC",
            tuneGrid = tuneGrid,
            preProc = preProc,
            trControl = control,
            na.action = na.pass
          )
        
        results[[i]]$pred = predict(mod, dataset[-in_train, ], type = "prob", na.action = na.pass)
        
        mods[[i]] = mod
        print(paste0("Cross Validation: ", i, "/", cv))
      }
      resultsOuter[[j]] = results
      modsOuter[[j]] = mods
    }
  }
  set.seed(987)
  #fit on whole dataset 
  train_set = datasetWithSites[,-siteColumn]
  finalMod = train(
    as.formula(paste(outcomeVariable, "~ .")),
    data = train_set,
    method = method,
    metric = "ROC",
    tuneGrid = tuneGrid,
    preProc = c("center", "scale", "knnImpute"),
    trControl = trainControl(
      method = "repeatedcv",
      number = 10,
      repeats = 10,
      classProbs = TRUE,
      summaryFunction = twoClassSummary,
      selectionFunction =
        "best"
    ),
    na.action = na.pass
  )
  
  results_mods = list(results = resultsOuter,
                      mods = modsOuter,
                      finalMod = finalMod)
}

combineResultsInSequence = function(results_mods,
                                    predictYes = T,
                                    colsToSave = NULL,
                                    #only works nested site
                                    nestedSite = T)
{
  if (nestedSite)
  {
    results_seq = NULL
    for (i in seq(1:length(results_mods[[1]])))
    {
      if (predictYes)
      {
        results_seq$pred = c(results_seq$pred, results_mods$results[[i]]$pred$Yes)
      } else
      {
        results_seq$pred = c(results_seq$pred, results_mods$results[[i]]$pred$No)
      }
      results_seq$obs = c(results_seq$obs,
                          as.character(results_mods$results[[i]]$obs))
      
      if (!is.null(colsToSave))
      {
        results_seq[[colsToSave]] = c(results_seq[[colsToSave]], results_mods$results[[i]][, colsToSave])
      }
    }
    results_seq$obs = as.factor(results_seq$obs)
    return(results_seq)
  } else
    #nestedCV
  {
    results = results_mods$results
    results_seq = NULL
    for (j in seq(1:length(results)))
    {
      for (i in seq(1:length(results[[j]])))
      {
        if (predictYes)
        {
          results_seq$pred = c(results_seq$pred, results[[j]][[i]]$pred$Yes)
        } else
        {
          results_seq$pred = c(results_seq$pred, results[[j]][[i]]$pred$No)
        }
        results_seq$obs = c(results_seq$obs, as.character(results[[j]][[i]]$obs))
      }
    }
    results_seq$obs = as.factor(results_seq$obs)
    return(results_seq)
  }
  
}

permutationPValue = function(results_seq,
                             auc_seq,
                             seed = 987,
                             predictYes = T,
                             seq = T,
                             dataSet = NULL,
                             outcomeVariable = NULL)
{
  set.seed(seed)
  auc_null = NULL
  #significance level <0.0001
  if (!seq)
  {
    if (is.null(dataSet) | is.null(outcomeVariable))
    {
      print("You need to define a dataset and an outcome variable!")
      pPerm = NULL
    } else
    {
      for (i in seq (1:10001))
      {
        perm = permute(dataSet[, outcomeVariable])
        #set direction explicitly so not biased towards higher roc values (just makes it less likely for things to be significant, however)
        #https://www.rdocumentation.org/packages/pROC/versions/1.15.3/topics/roc
        if (predictYes)
        {
          auc_null = c(
            auc_null,
            roc(
              predictor = results_seq$Yes,
              response = perm,
              levels = c("No", "Yes"),
              direction = "<"
            )$auc
          )
        } else
        {
          auc_null = c(
            auc_null,
            roc(
              predictor = results_seq$No,
              response = perm,
              levels = c("No", "Yes"),
              direction = ">"
            )$auc
          )
        }
      }
      pPerm = (1 + sum(auc_null >= auc_seq)) / 10001
    }
  } else
  {
    for (i in seq (1:10001))
    {
      perm = permute(results_seq$obs)
      #set direction explicitly so not biased towards higher roc values (just makes it less likely for things to be significant, however)
      #https://www.rdocumentation.org/packages/pROC/versions/1.15.3/topics/roc
      if (predictYes)
      {
        auc_null = c(
          auc_null,
          roc(
            predictor = results_seq$pred,
            response = perm,
            levels = c("No", "Yes"),
            direction = "<"
          )$auc
        )
      } else
      {
        auc_null = c(
          auc_null,
          roc(
            predictor = results_seq$pred,
            response = perm,
            levels = c("No", "Yes"),
            direction = ">"
          )$auc
        )
      }
    }
    pPerm = (1 + sum(auc_null >= auc_seq)) / 10001
  }
  #get p value by taking proportion of permutated values greater or equal to the actual value
  return(pPerm)
}

coefEvaluation = function(datasetWithSites,
                          results_mods,
                          site = "Site",
                          isNestedCV = F)
{
  siteColumn = which(colnames(datasetWithSites) == site)
  coefs = NULL
  if (isNestedCV)
  {
    #nested cv
    for (i in seq(1:length(results_mods$mods)))
    {
      for (j in seq(1:length(results_mods$mods[[i]])))
      {
        coefs = c(
          coefs,
          coef(
            results_mods$mods[[i]][[j]]$finalModel,
            results_mods$mods[[i]][[j]]$bestTune$lambda
          )
        )
      }
    }
  } else
    #nested site
  {
    for (i in seq(1:length(levels(datasetWithSites[[siteColumn]]))))
    {
      coefs = c(
        coefs,
        coef(
          results_mods$mods[[i]]$finalModel,
          results_mods$mods[[i]]$bestTune$lambda
        )
      )
      
    }
  }
  
  
  lengthC = NULL
  
  lengthC = length(coefs[[1]])
  
  #just get numbers
  coefs_extract = NULL
  for (i in seq(1:length(levels(datasetWithSites[[siteColumn]]))))
  {
    coefs_extract = rbind(coefs_extract, coefs[[i]][1:lengthC])
    
  }
  
  #get matrix of coefficients presence (1) or absence (0)
  #Presence or absence of predictors across all 14 LOSOCV models
  coefs_presence = NULL
  coefs_presence = coefs_extract[1:length(levels(datasetWithSites[[siteColumn]])), 1:lengthC]
  coefs_presenceint = coefs_extract[1:length(levels(datasetWithSites[[siteColumn]])), 1:lengthC]
  coefs_presence[coefs_presence != 0] <- 1
  coefs_presenceint[coefs_presenceint != 0] <- 1
  
  #stability of feature selection http://jmlr.org/papers/volume18/17-514/17-514.pdf
  print(getStability(coefs_presence))
  
  #get rank of coef by importance as in sports ranking
  coefs_rank = NULL
  
  for (i in seq(1:length(levels(datasetWithSites[[siteColumn]]))))
  {
    #rank absolute value excluding the intercept for each model
    coefs_rank = rbind(coefs_rank, rank(abs(coefs_extract[i, 1:lengthC]), ties.method = "min"))
  }
  
  # rank the mean ranks of each column across all models
  coefs_rank_mean = colMeans(coefs_rank)
  
  #Invert order of rank to identify top models
  coefs_order = rank(-coefs_rank_mean)
  
  #Get the column names (not the intercept)
  
  coef_names = dimnames(coefs[[1]])[[1]][1:lengthC]
  
  coefs_means = colMeans(coefs_extract)[1:lengthC]
  
  df = data.frame(coef_names,
                  coefs_order,
                  coefs_means,
                  colMeans(coefs_presenceint))
}

############################################################
#sample size
pmsampsize(
  type = "b",
  rsquared = 0.25,
  parameters = 95,
  shrinkage = 0.9,
  #lower value = smaller sample required
  prevalence = 0.5
)
#Need larger sample so did lasso aware of limitations

#Outlookall = read_spss("OutlookData/PSYGRID 1E.sav")
#Outlookall_labels = makeVlist(Outlookall)
#write_csv(all_to_factor(Outlookall), path = "Outlookall.csv") #Have to manually remove "NA"s from csv file
#write_csv(format.data.frame(Outlookall_labels), path = "Outlookall_labels.csv")

#################################################################
#Impute 12 month outcome for EDEN and Outlook datasets separately, using only PANSS at 6 months and 12 months.
#EDEN PANSS data
eden_6_12 = read_csv("EDEN_6_12_PANSS.csv")
eden_6_12 = csv_to_factor(eden_6_12)
eden_6_12$row = rownames(eden_6_12)
#remove empty rows (not including coloumn 64 which is the row number and never empty)
eden_6_12 = eden_6_12[rowSums(is.na(eden_6_12[,-64])) != ncol(eden_6_12[,-64]),]
#Impute the missing data (not using row number column 64)
tempData_eden_6_12 = mice(eden_6_12[,-64], m = 1, seed = 987)
#Take the imputed outcome
complete_eden_6_12 = complete(tempData_eden_6_12, 1)
#preserve original rownames
rownames(complete_eden_6_12) = eden_6_12$row

#Outlook PANSS data
outlook_6_12 = read_csv("Outlook_6_12_PANSS.csv")
outlook_6_12 = csv_to_factor(outlook_6_12)
outlook_6_12$row = rownames(outlook_6_12)
#remove empty rows (not including coloumn 64 which is the row number and never empty)
outlook_6_12 = outlook_6_12[rowSums(is.na(outlook_6_12[,-64])) != ncol(outlook_6_12[,-64]),]
#Impute the missing data (not using row number column 64)
tempData_outlook_6_12 = mice(outlook_6_12[,-64], m = 1, seed = 987)
#Take the imputed outcome
complete_outlook_6_12 = complete(tempData_outlook_6_12, 1)
#preserve original rownames
rownames(complete_outlook_6_12) = outlook_6_12$row

###############################################################
#Load study data
eden_outlook_preproc = read_csv("outlook_eden_preproc_rem.csv")
eden_outlook_preproc_factor = csv_to_factor(eden_outlook_preproc)
eden_all = eden_outlook_preproc_factor[which(eden_outlook_preproc_factor$Study == "EDEN"),!(colnames(eden_outlook_preproc_factor) %in% c("Study"))]
outlook_all = eden_outlook_preproc_factor[which(eden_outlook_preproc_factor$Study == "Outlook"),!(colnames(eden_outlook_preproc_factor) %in% c("Study"))]

#Choose Final predictors based on expert knowledge
eden_all_final = eden_all[rownames(complete_eden_6_12),-c(3:9, 17, 19:24, 113:119)]
outlook_all_final = outlook_all[rownames(complete_outlook_6_12),-c(3:9, 17, 19:24, 113:119)]

#Replace missing outcome column with imputed outcome column from above
eden_all_final$M12_PANSS_Period_Rem = complete_eden_6_12$M12_PANSS_Period_Rem
outlook_all_final$M12_PANSS_Period_Rem = complete_outlook_6_12$M12_PANSS_Period_Rem

#do Lasso as suggested by reviewer
tuneGridLasso = expand.grid(alpha = 1, lambda = 10 ^ seq(-0.3, -5, length = 100))

#Remission Data set up
eden_all_final_Rem = processData(
  dataset = eden_all_final,
  outcomeVariable = "M12_PANSS_Period_Rem",
  siteOutcomeColumnNames = c("Site", "M12_PANSS_Period_Rem")
)
#Don't have 25% cut off for test dataset predictors
outlook_all_final_Rem = processData(
  dataset = outlook_all_final,
  outcomeVariable = "M12_PANSS_Period_Rem",
  removeMissing = F,
  siteOutcomeColumnNames = c("Site", "M12_PANSS_Period_Rem")
)
#same predictor columns in both datasets
eden_all_final_Rem = matchData(otherData = outlook_all_final_Rem, yourData = eden_all_final_Rem)
outlook_all_final_Rem = matchData(otherData = eden_all_final_Rem, yourData = outlook_all_final_Rem)

#########################################################################
#internal validation
#nested cv
results_mods_final_RemLassoNested = nestedSiteCVUpdate(
  datasetWithSites = eden_all_final_Rem,
  outcomeVariable = "M12_PANSS_Period_Rem",
  tuneGrid = tuneGridLasso,
  nestedSite = F #nested cv
)

results_seq_final_RemLassoNested = combineResultsInSequence(results_mods_final_RemLassoNested,
                                                            predictYes = F,
                                                            nestedSite = F) #nested cv 10fold repeated 10 times default
results_seq_final_RemLassoNested_roc = roc(
  predictor = results_seq_final_RemLassoNested$pred,
  response = results_seq_final_RemLassoNested$obs,
  ci = T,
  levels = c("No", "Yes"),
  direction = ">"
)
results_seq_final_RemLassoNested_auc = results_seq_final_RemLassoNested_roc$auc
results_seq_final_RemLassoNested_roc #AUC 0.699
df_coefs_final_RemLassoNested = coefEvaluation(
  datasetWithSites = eden_all_final_Rem,
  results_mods = results_mods_final_RemLassoNested,
  isNestedCV = T
) #stability 0.618
View(df_coefs_final_RemLassoNested) 
results_seq_final_RemLassoNested_auc_p = permutationPValue(results_seq = results_seq_final_RemLassoNested,
                                                           auc_seq = results_seq_final_RemLassoNested_auc,
                                                           predictYes = F)
results_seq_final_RemLassoNested_auc_p # p <0.0001
#Final Model on whole dataset
View(coef(
  results_mods_final_RemLassoNested$finalMod$finalModel,
  results_mods_final_RemLassoNested$finalMod$bestTune$lambda
)[,1])

val.prob.ci.2(
  p = results_seq_final_RemLassoNested$pred,
  y = results_seq_final_RemLassoNested$obs == "No",
  g = 10,
  logistic.cal = T,
  lty.log = 9,
  col.log = "red",
  lwd.log = 1.5,
  col.ideal = "blue",
  lwd.ideal = 0.5,
  dostats = T
)

###############################################################
#losocv
#internal-external validation
results_mods_final_RemLasso = nestedSiteCVUpdate(
  datasetWithSites = eden_all_final_Rem,
  outcomeVariable = "M12_PANSS_Period_Rem",
  tuneGrid = tuneGridLasso,
  nestedSite = T #across 14 sites
)

results_seq_final_RemLasso = combineResultsInSequence(results_mods_final_RemLasso, predictYes = F)

results_seq_final_RemLasso_roc = roc(
  predictor = results_seq_final_RemLasso$pred,
  response = results_seq_final_RemLasso$obs,
  ci = T,
  levels = c("No", "Yes"),
  direction = ">"
)
results_seq_final_RemLasso_auc = results_seq_final_RemLasso_roc$auc
results_seq_final_RemLasso_roc #AUC 0.67
df_coefs_final_RemLasso = coefEvaluation(datasetWithSites = eden_all_final_Rem, results_mods = results_mods_final_RemLasso)
View(df_coefs_final_RemLasso) #stability 0.61
results_seq_final_RemLasso_auc_p = permutationPValue(results_seq = results_seq_final_RemLasso,
                                                     auc_seq = results_seq_final_RemLasso_auc,
                                                     predictYes = F)
results_seq_final_RemLasso_auc_p
#Final Model on whole dataset
View(coef(
  results_mods_final_RemLasso$finalMod$finalModel,
  results_mods_final_RemLasso$finalMod$bestTune$lambda
)[,1])

val.prob.ci.2(
  p = results_seq_final_RemLasso$pred,
  y = results_seq_final_RemLasso$obs == "No",
  g = 10,
  logistic.cal = T,
  lty.log = 9,
  col.log = "red",
  lwd.log = 1.5,
  col.ideal = "blue",
  lwd.ideal = 0.5,
  dostats = T
)

##############################################################
#external  model
outlook_all_final_Rem_glm_result_LASSO = predict(
  results_mods_final_RemLasso$finalMod,
  outlook_all_final_Rem,
  type = "prob",
  na.action = na.pass
)
val.prob.ci.2(
  p = outlook_all_final_Rem_glm_result_LASSO$No,
  y = outlook_all_final_Rem$M12_PANSS_Period_Rem == "No",
  g = 10,
  logistic.cal = T,
  lty.log = 9,
  col.log = "red",
  lwd.log = 1.5,
  col.ideal = "blue",
  lwd.ideal = 0.5
)
outlook_all_final_Rem_glm_result_LASSO_roc = roc(
  predictor = outlook_all_final_Rem_glm_result_LASSO$No,
  response = outlook_all_final_Rem$M12_PANSS_Period_Rem,
  ci = T,
  levels = c("No", "Yes"),
  direction = ">"
)
outlook_all_final_Rem_glm_result_LASSO_auc = outlook_all_final_Rem_glm_result_LASSO_roc$auc
outlook_all_final_Rem_glm_result_LASSO_roc #0.721
outlook_all_final_Rem_glm_result_LASSO_roc_p = permutationPValue(
  results_seq = outlook_all_final_Rem_glm_result_LASSO,
  seq = F,
  dataSet = outlook_all_final_Rem,
  auc_seq = outlook_all_final_Rem_glm_result_LASSO_auc,
  predictYes = F,
  outcomeVariable = "M12_PANSS_Period_Rem"
)
outlook_all_final_Rem_glm_result_LASSO_roc_p

#dca
#Decision Curve Analysis External
dca_glm_Rem_ext = NULL
dca_glm_Rem_ext$M12_PANSS_Period_Rem = as.integer(as.character(outlook_all_final_Rem$M12_PANSS_Period_Rem) ==
                                                    "No")
dca_glm_Rem_ext$Model = outlook_all_final_Rem_glm_result_LASSO$No
dca_glm_Rem_ext$ADJ_DUP = outlook_all_final_Rem$ADJ_DUP
#get data across all thresholds
dca(
  data = as.data.frame(dca_glm_Rem_ext),
  outcome = "M12_PANSS_Period_Rem",
  predictors = c("Model", "ADJ_DUP"),
  smooth = "TRUE",
  loess.span = 0.5,
  probability = c(TRUE, FALSE),
  graph = T,
  xstart = 0.35,
  xstop = 0.74
)


# 31 predictor variables in final LASSO model = 1,2,3,4,8,12,13,19,23,25,27,28,29,32,34,35,36,38,43,49,51,52,55,62,70,74,77,78,86,88,93
################################################################################################################################
#make model for app
#construct glm then replace the coefficients with the LASSO ones
eden_all_final_Rem_reduced = eden_all_final_Rem[,c(1,2,3,4,8,12,13,19,23,25,27,28,29,32,34,35,36,38,43,49,51,52,55,62,70,74,77,78,86,88,93,97)] 
write_csv(eden_all_final_Rem_reduced,"unstandardised_dep.csv", na = "")
#for online model, refit with prestandardised local concentration on based on whole range of values from IMD 2007 just to get correct preprocess info
#replace local concentration with population standardised
eden_all_final_Rem_reduced = read_csv("standardised_dep.csv")
eden_all_final_Rem_reduced$M12_PANSS_Period_Rem = as.factor(eden_all_final_Rem_reduced$M12_PANSS_Period_Rem)

set.seed(987)
#fit a model to get imputation model and replace coefficients with those from final LASSO above
mod_eden_all_final_Rem_outlook_glm = train(M12_PANSS_Period_Rem ~ ., data=eden_all_final_Rem_reduced, method="glm", metric="ROC", preProc = c("center", "scale","knnImpute"), 
                                           trControl = trainControl(method="none", classProbs=TRUE, summaryFunction=twoClassSummary), na.action = na.pass)
mod_eden_all_final_Rem_outlook_LASSO_reduced = mod_eden_all_final_Rem_outlook_glm
#get coefficients from lasso
View(coef(
  results_mods_final_RemLasso$finalMod$finalModel,
  results_mods_final_RemLasso$finalMod$bestTune$lambda
)[,1])
#replace glm coefficients with lasso coefficients
mod_eden_all_final_Rem_outlook_LASSO_reduced$finalModel$coefficients = c(-0.118644, #intercept
                                                                         -0.066973, #Sex.Male
                                                                         0.128178, #Qualification_Level_Ordinal
                                                                         0.065110, #BL_Paid_Employ.Yes
                                                                         -0.015504, #BL_Education.Yes
                                                                         -0.215642, #ADJ_DUP
                                                                         -0.005444, #Client_Scholastic_Performance_Childhood
                                                                         -0.050767, #Client_Adaption_School_Childhood
                                                                         -0.176558, #Client_Sociability_Withdrawal_Late_Adolescence
                                                                         -0.009559, #Client_Social_Sexual_Aspects_Late_Adolescence
                                                                         -0.037743, #Client_Employed_At_School
                                                                         -0.056620, #Client_Job_Change_Interrupted_School_Attendance
                                                                         -0.016662, #Client_Establishment_Independence
                                                                         -0.103722, #Client_Highest_Functioning_Achieved_Life
                                                                         -0.108008, #Client_Energy_Level
                                                                         0.155397, #BL_PANSS_P2_Conceptual_Disorganization
                                                                         -0.353536, #BL_PANSS_P3_Hallucinatory_Behaviour
                                                                         0.075174, #BL_PANSS_P4_Excitement
                                                                         -0.016680, #BL_PANSS_P6_Suspiciousness
                                                                         -0.108625, #BL_PANSS_N4_Passive_Social_Withdrawal
                                                                         -0.034999, #BL_PANSS_G3_Guilt
                                                                         0.002686, #BL_PANSS_G5_Mannerisms
                                                                         0.043360, #BL_PANSS_G6_Depression
                                                                         -0.044998, #BL_PANSS_G9_Unusual_Thought_Content
                                                                         -0.012430, #BL_PANSS_G16_Active_Social_Avoidance
                                                                         -0.042242, #YMRS_Appearance
                                                                         0.016379, #IS_Do_Not_Need_Med
                                                                         0.003058, #IS_Nervous_Mental_Illness
                                                                         0.054798, #IS_Not_Due_Illness
                                                                         -0.034980, #CDSS_Suicide
                                                                         0.065991, #BL_GAF_Symptoms
                                                                         -0.104312 #PCT_Local_Concentration_2007
)

#remove training data for sharing model
mod_eden_all_final_Rem_outlook_LASSO_reduced$trainingData = NULL
save(mod_eden_all_final_Rem_outlook_LASSO_reduced, file = "mod_eden_all_final_Rem_outlook_LASSO_reduced.rda")
