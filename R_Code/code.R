library(doParallel)
library(mice)
library(readr)
library(caret)
library(CalibrationCurves)
library(pmsampsize)
library(pROC)
library(gtools)
library(dcurves)

#enable multicore (windows) which roughly halfs time for analysis runs
cl <- makeCluster(detectCores(), type = 'PSOCK')
registerDoParallel(cl)

options(max.print = 1000000)

#don't use scientific notation (revert back with options(scipen=0)
options(scipen = 999)
options(digits = 4)

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

#sample size with 14 expert chosen predictors
pmsampsize(
  type = "b",
  rsquared = 0.25,
  parameters = 14,
  shrinkage = 0.9,
  prevalence = 0.5
)

#Load study data
#EDEN
eden = read_csv("eden_all.csv")
eden$Study = NULL
eden = csv_to_factor(eden)
#
tempData <- mice(eden,m=10,seed=987)

control <- trainControl(## 10-fold CV repeated 5 times
  method = "repeatedcv",
  number = 10,
  repeats = 5,
  classProbs=TRUE, 
  summaryFunction=twoClassSummary,
  savePredictions = T)

finalModels = list()
crossValModels = list()
for (i in seq(1:tempData$m))
{
  #Get imputed data
  eden_imp = complete(tempData,i)
  #just take the columns we are using except outcome as standardising first
  eden_imp_exp = eden_imp[,c(2,15,18,25,49,53,54,55,63,72,100,111,112,120)]
  #standardise the columns before building model
  preProcValues = preProcess(eden_imp_exp, method = c("center", "scale"))
  eden_imp_exp_stand = predict(preProcValues, eden_imp_exp)
  #Add factor outcome back in
  eden_imp_exp_stand$M12_PANSS_Period_Rem = eden$M12_PANSS_Period_Rem
  #Remove rows with missing outcomes
  eden_imp_exp_stand_MID = eden_imp_exp_stand[complete.cases(eden_imp_exp_stand), ]
  #need to return design matrix to reestimate intercept
  finalModels[[i]] = glm(M12_PANSS_Period_Rem ~ ., data = eden_imp_exp_stand_MID, family = "binomial", x = T)
  crossValModels[[i]] = train(M12_PANSS_Period_Rem ~ ., data = eden_imp_exp_stand_MID, method = "glm", metric = "ROC", trControl=control, na.action=na.pass)
}

internalPreds = NULL
internalOutcomes = NULL
for(i in seq(1:tempData$m))
{
  internalPreds = c(internalPreds,crossValModels[[i]]$pred$No)
  internalOutcomes = c(internalOutcomes, as.character(crossValModels[[i]]$pred$obs))
}

#C-statistic
results_internal_roc = roc(
  predictor = internalPreds,
  response = as.factor(internalOutcomes),
  ci = T,
  levels = c("No", "Yes"),
  direction = ">"
)
results_internal_roc
results_internal_auc = results_internal_roc$auc

#permutation p value
set.seed(987)
auc_null = NULL
for(i in seq (1:10001))
{
  perm = permute(as.factor(internalOutcomes))
  auc_null = c(auc_null, roc(predictor = internalPreds, response = perm, levels=c("No", "Yes"), direction=">")$auc)
}
pPermInternal = (1+sum(auc_null >= results_internal_auc))/10001
pPermInternal

#calibration plot
#increase memory allocated to R
memory.limit(size=56000)
intval_cal = val.prob.ci.2(p=internalPreds, y=internalOutcomes=="No", g=10, logistic.cal = T, lty.log=9,
              col.log="red", lwd.log=1.5, col.ideal="blue", lwd.ideal=0.5, dostats = T, smooth = "rcs")
#shrinkage factor = 0.8381
shrinkage = intval_cal["Slope"]

#Pool results to get predictor estimates based on Rubin's rule (coefficients reversed as using to predict non-remission)
View(-summary(pool(finalModels), conf.int = T, exponentiate = F, conf.level = 0.95)[,c(2,7,8)])
View(exp(-summary(pool(finalModels), conf.int = T, exponentiate = F, conf.level = 0.95)[,c(2,7,8)]))
#With shrinkage applied (N.B. intercept still to be re-estimated)
View(-summary(pool(finalModels), conf.int = T, exponentiate = F, conf.level = 0.95)[,c(2,7,8)]*shrinkage)
View(exp(-summary(pool(finalModels), conf.int = T, exponentiate = F, conf.level = 0.95)[,c(2,7,8)]*shrinkage))

#for each imputed dataset recalculate intercept with pooled shrunk final coefficients
#then average the intercepts to get final intercept
shrunk_coefs = summary(pool(finalModels))[[2]]*shrinkage
models = list()
intercepts = NULL
for(i in seq(1:tempData$m))
{
  shrunk_LP = finalModels[[i]]$x[,2:15] %*% shrunk_coefs[2:15]
  models[[i]] = glm(finalModels[[i]]$y ~ offset(shrunk_LP), family = "binomial")
}
#pooled recalculated intercept
recalc_intercept = summary(pool(models))[[2]]

#Get a GLM model
finalModel = finalModels[[1]]
#replace coefficients with pooled ones shrunk coefs
finalModel$coefficients = shrunk_coefs
#replace intercept with pooled recalculated intercept
finalModel$coefficients[1] = recalc_intercept

############################################################
#Load study data
#outlook
outlook = read_csv("outlook_all.csv")
outlook$Study = NULL
#remove correlated variables - identified as problem in MICE
outlook$PCT_Average_Rank_2007 = NULL
outlook$PCT_Employment_Scale_2007 = NULL
outlook$PCT_Extent_2007 = NULL
outlook$PCT_Local_Concentration_2007 = NULL
outlook$PCT_Income_Scale_2007 = NULL
outlook = csv_to_factor(outlook)

tempData2 <- mice(outlook,m=10,seed=987)

Results = list()
Obs = list()
DUP = list()
for (i in seq(1:tempData2$m))
{
  #Get imputed data
  outlook_imp = complete(tempData2,i)
  #just take the columns we are using except outcome as standardising first
  outlook_imp_exp = outlook_imp[,c(2,15,18,25,49,53,54,55,63,72,100,111,112,120)]
  #standardise the columns before building model
  preProcValues = preProcess(outlook_imp_exp, method = c("center", "scale"))
  outlook_imp_exp_stand = predict(preProcValues, outlook_imp_exp)
  #Add factor outcome back in
  outlook_imp_exp_stand$M12_PANSS_Period_Rem = outlook$M12_PANSS_Period_Rem
  #Remove rows with missing outcomes
  outlook_imp_exp_stand_MID = outlook_imp_exp_stand[complete.cases(outlook_imp_exp_stand), ]
  Results[[i]] = predict(finalModel, outlook_imp_exp_stand_MID, type = "response", na.action = na.pass)
  #Predicting No
  Results[[i]] = 1-Results[[i]]
  Obs[[i]] =  outlook_imp_exp_stand_MID$M12_PANSS_Period_Rem
  DUP[[i]] = outlook_imp_exp_stand_MID$ADJ_DUP
}

externalPreds = NULL
externalOutcomes = NULL
externalDUP = NULL
for(i in seq(1:tempData2$m))
{
  externalPreds = c(externalPreds,Results[[i]])
  externalOutcomes = c(externalOutcomes, as.character(Obs[[i]]))
  externalDUP = c(externalDUP, DUP[[i]])
}

#C-statistic
results_external_roc = roc(
  predictor = externalPreds,
  response = as.factor(externalOutcomes),
  ci = T,
  levels = c("No", "Yes"),
  direction = ">"
)
results_external_roc
results_external_auc = results_external_roc$auc

#permutation p value
set.seed(987)
auc_null = NULL
for(i in seq (1:10001))
{
  perm = permute(as.factor(externalOutcomes))
  auc_null = c(auc_null, roc(predictor = externalPreds, response = perm, levels=c("No", "Yes"), direction=">")$auc)
}
pPermExternal = (1+sum(auc_null >= results_external_auc))/10001
pPermExternal

#calibration plot
pdf("figure2.pdf", width = 7, height = 7)
val.prob.ci.2(p=externalPreds, y=externalOutcomes=="No", g=5, logistic.cal = T, lty.log=9,
              col.log="red", lwd.log=1.5, col.ideal="blue", lwd.ideal=0.5)
dev.off()

#dca
dca_ext = NULL
dca_ext$M12_PANSS_Period_Rem = as.integer(externalOutcomes=="No")
dca_ext$Model = externalPreds
dca_ext$DUP = externalDUP
#get data across all thresholds
dca_ext_calc = dca(M12_PANSS_Period_Rem ~ Model + DUP, data = as.data.frame(dca_ext), as_probability =  "DUP",
                   thresholds = seq(0.3, 0.77, by = 0.01))
pdf("figure3.pdf", width = 7, height = 7)
dca_ext_calc %>%   
  plot(smooth = TRUE)
dev.off()

#######################################################################
#Demographic comparisons
#######################################################################
#summary statistics
#Age
eden_final_demographics = eden[complete.cases(eden[ , "M12_PANSS_Period_Rem"]),]
outlook_final_demographics = outlook[complete.cases(outlook[ , "M12_PANSS_Period_Rem"]),]

summary(eden$Age_Entry)
mean(eden$Age_Entry, na.rm = T)
sd(eden$Age_Entry, na.rm = T)

summary(eden_final_demographics$Age_Entry)
mean(eden_final_demographics$Age_Entry, na.rm = T)
sd(eden_final_demographics$Age_Entry, na.rm = T)

summary(outlook$Age_Entry)
mean(outlook$Age_Entry, na.rm = T)
sd(outlook$Age_Entry, na.rm = T)

summary(outlook_final_demographics$Age_Entry)
mean(outlook_final_demographics$Age_Entry, na.rm = T)
sd(outlook_final_demographics$Age_Entry, na.rm = T)

summary(aov(y~group, data = data.frame(group=factor(rep(1:4, c(1027,673,399,191))),y=c(eden$Age_Entry,eden_final_demographics$Age_Entry,
                                                                                       outlook$Age_Entry,outlook_final_demographics$Age_Entry))))

#Sex
summary(eden$Sex)
summary(eden_final_demographics$Sex)
summary(outlook$Sex)
summary(outlook_final_demographics$Sex)

chisq.test(as.table(rbind(c(318,709),c(210,463),c(153,246),c(73,118))), correct = F)

#EET
summary(eden$BL_EET)
summary(eden_final_demographics$BL_EET)
summary(outlook$BL_EET)
summary(outlook_final_demographics$BL_EET)

chisq.test(as.table(rbind(c(589,284),c(383,190),c(225,174),c(106,85))), correct = F)

#Qualification
summary(as.factor(eden$Qualification_Level_Ordinal))
summary(as.factor(eden_final_demographics$Qualification_Level_Ordinal))
summary(as.factor(outlook$Qualification_Level_Ordinal))
summary(as.factor(outlook_final_demographics$Qualification_Level_Ordinal))

chisq.test(as.table(rbind(c(245,399,262,98),c(156,255,173,74),c(89,130,92,69),c(40,67,46,35))), correct = F)

#DUP
summary(eden$ADJ_DUP)
mean(eden$ADJ_DUP, na.rm = T)
sd(eden$ADJ_DUP, na.rm = T)

summary(eden_final_demographics$ADJ_DUP)
mean(eden_final_demographics$ADJ_DUP, na.rm = T)
sd(eden_final_demographics$ADJ_DUP, na.rm = T)

summary(outlook$ADJ_DUP)
mean(outlook$ADJ_DUP, na.rm = T)
sd(outlook$ADJ_DUP, na.rm = T)

summary(outlook_final_demographics$ADJ_DUP)
mean(outlook_final_demographics$ADJ_DUP, na.rm = T)
sd(outlook_final_demographics$ADJ_DUP, na.rm = T)

summary(aov(y~group, data = data.frame(group=factor(rep(1:4, c(1027,673,399,191))),y=c(eden$ADJ_DUP,
                                                                                       eden_final_demographics$ADJ_DUP,
                                                                                       outlook$ADJ_DUP,
                                                                                       outlook_final_demographics$ADJ_DUP))))

#Deprivation
summary(eden$PCT_Average_Scrore_2007)
mean(eden$PCT_Average_Scrore_2007, na.rm = T)
sd(eden$PCT_Average_Scrore_2007, na.rm = T)

summary(eden_final_demographics$PCT_Average_Scrore_2007)
mean(eden_final_demographics$PCT_Average_Scrore_2007, na.rm = T)
sd(eden_final_demographics$PCT_Average_Scrore_2007, na.rm = T)

summary(outlook$PCT_Average_Scrore_2007)
mean(outlook$PCT_Average_Scrore_2007, na.rm = T)
sd(outlook$PCT_Average_Scrore_2007, na.rm = T)

summary(outlook_final_demographics$PCT_Average_Scrore_2007)
mean(outlook_final_demographics$PCT_Average_Scrore_2007, na.rm = T)
sd(outlook_final_demographics$PCT_Average_Scrore_2007, na.rm = T)

summary(aov(y~group, data = data.frame(group=factor(rep(1:4, c(1027,673,399,191))),y=c(eden$PCT_Average_Scrore_2007,
                                                                                       eden_final_demographics$PCT_Average_Scrore_2007,
                                                                                       outlook$PCT_Average_Scrore_2007,
                                                                                       outlook_final_demographics$PCT_Average_Scrore_2007))))

#PANSS Totals
eden$BL_PANSS_Total = rowSums(eden[,53:82])
outlook$BL_PANSS_Total = rowSums(outlook[,53:82])
eden_final_demographics$BL_PANSS_Total = rowSums(eden_final_demographics[,53:82])
outlook_final_demographics$BL_PANSS_Total = rowSums(outlook_final_demographics[,53:82])

summary(eden$BL_PANSS_Total)
mean(eden$BL_PANSS_Total, na.rm = T)
sd(eden$BL_PANSS_Total, na.rm = T)

summary(eden_final_demographics$BL_PANSS_Total)
mean(eden_final_demographics$BL_PANSS_Total, na.rm = T)
sd(eden_final_demographics$BL_PANSS_Total, na.rm = T)

summary(outlook$BL_PANSS_Total)
mean(outlook$BL_PANSS_Total, na.rm = T)
sd(outlook$BL_PANSS_Total, na.rm = T)

summary(outlook_final_demographics$BL_PANSS_Total)
mean(outlook_final_demographics$BL_PANSS_Total, na.rm = T)
sd(outlook_final_demographics$BL_PANSS_Total, na.rm = T)

summary(aov(y~group, data = data.frame(group=factor(rep(1:4, c(1027,673,399,191))),y=c(eden$BL_PANSS_Total,
                                                                                       eden_final_demographics$BL_PANSS_Total,
                                                                                       outlook$BL_PANSS_Total,
                                                                                       outlook_final_demographics$BL_PANSS_Total))))
