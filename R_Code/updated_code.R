library(doParallel)
library(mice)
library(readr)
library(caret)
library(CalibrationCurves)
library(pmsampsize)
library(pROC)
library(gtools)
library(dcurves)
library(psfmi)

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

# 10-fold CV repeated 5 times
control <- trainControl(
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

#pool AUCs using Rubin's Rules
internalROCs = list()
internalROCsValues = list()
internalROCsSEs= list()
for(i in seq(1:tempData$m))
{
  internalROCs[[i]] = roc(
    predictor = crossValModels[[i]]$pred$No,
    response = crossValModels[[i]]$pred$obs,
    ci = T,
    levels = c("No", "Yes"),
    direction = ">"
  )
  internalROCsValues[[i]] = internalROCs[[i]]$auc
  internalROCsSEs[[i]] = (internalROCs[[i]]$auc - internalROCs[[i]]$ci[1])/1.96
}
#Correctly pooled C-statistic and 95% CI using Rubin's Rules with logit transformation
pool_auc(internalROCsValues, internalROCsSEs, nimp = 10, log_auc = T)

#Individual permutation tests for each multiple imputation dataset
psPermInternal = list()
set.seed(987)
for(i in seq(1:tempData$m))
{
  #permutation p value
  auc_null = NULL
  for(j in seq (1:10001))
  {
    perm = permute(crossValModels[[i]]$pred$obs)
    auc_null = c(auc_null, roc(predictor = crossValModels[[i]]$pred$No, response = perm, levels=c("No", "Yes"), direction=">")$auc)
  }
  psPermInternal[[i]] = (1+sum(auc_null >= internalROCsValues[[i]]))/10001
}
psPermInternal

internalCalIntValues = list()
internalCalIntSEs = list()
internalCalSlopeValues = list()
internalCalSlopeSEs = list()
for(i in seq(1:tempData$m))
{
  internalCal = val.prob.ci.3(p=crossValModels[[i]]$pred$No, y=as.character(crossValModels[[i]]$pred$obs)=="No", g=5, logistic.cal = T, lty.log=9,
                              col.log="red", lwd.log=1.5, col.ideal="blue", lwd.ideal=0.5)
  
  internalCalIntValues[[i]] = internalCal$stats[["Intercept"]]
  internalCalIntSEs[[i]] = (internalCal$stats[["Intercept"]] - internalCal$cl.interc[1])/1.96
  internalCalSlopeValues[[i]] = internalCal$stats[["Slope"]]
  internalCalSlopeSEs[[i]] = (internalCal$stats[["Slope"]] - internalCal$cl.slope[[1]])/1.96
}

#Correctly pooled internal Calibration intercept and SE - should be zero for internal validation
rubin.rules(unlist(internalCalIntValues), unlist(internalCalIntSEs))
#Correctly pooled internal Calibration slope and SE
rubin.rules(unlist(internalCalSlopeValues), unlist(internalCalSlopeSEs))

#combined data 
internalPreds = NULL
internalOutcomes = NULL
for(i in seq(1:tempData$m))
{
  internalPreds = c(internalPreds,crossValModels[[i]]$pred$No)
  internalOutcomes = c(internalOutcomes, as.character(crossValModels[[i]]$pred$obs))
}

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
#alternatively using rms package as per Steyerberg - see commented code
#library(rms)
intercepts = NULL
for(i in seq(1:tempData$m))
{
  shrunk_LP = finalModels[[i]]$x[,2:15] %*% shrunk_coefs[2:15]
  models[[i]] = glm(finalModels[[i]]$y ~ offset(shrunk_LP), family = "binomial")
  #intercepts = c(intercepts, lrm.fit(y = finalModels[[i]]$y, offset= finalModels[[i]]$x[,2:15] %*% shrunk_coefs[2:15])$coef[1])
}
#pooled recalculated intercept
recalc_intercept = summary(pool(models))[[2]]
#mean(intercepts)

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

#pool AUCs using Rubin's Rules
externalROCs = list()
externalROCsValues = list()
externalROCsSEs= list()
for(i in seq(1:tempData2$m))
{
  externalROCs[[i]] = roc(
    predictor = Results[[i]],
    response = Obs[[i]],
    ci = T,
    levels = c("No", "Yes"),
    direction = ">"
  )
  externalROCsValues[[i]] = externalROCs[[i]]$auc
  externalROCsSEs[[i]] = (externalROCs[[i]]$auc - externalROCs[[i]]$ci[1])/1.96
}
#Correctly pooled C-statistic and 95% CI using Rubin's Rules with logit transformation
pool_auc(externalROCsValues, externalROCsSEs, nimp = 10, log_auc = T)

#Individual permutation tests for each multiple imputation dataset
psPermExternal = list()
set.seed(987)
for(i in seq(1:tempData2$m))
{
  #permutation p value
  auc_null = NULL
  for(j in seq (1:10001))
  {
    perm = permute(Obs[[i]])
    auc_null = c(auc_null, roc(predictor = Results[[i]], response = perm, levels=c("No", "Yes"), direction=">")$auc)
  }
  psPermExternal[[i]] = (1+sum(auc_null >= externalROCsValues[[i]]))/10001
}
psPermExternal

externalCalIntValues = list()
externalCalIntSEs = list()
externalCalSlopeValues = list()
externalCalSlopeSEs = list()
for(i in seq(1:tempData2$m))
{
  externalCal = val.prob.ci.3(p=Results[[i]], y=Obs[[i]]=="No", g=5, logistic.cal = T, lty.log=9,
                              col.log="red", lwd.log=1.5, col.ideal="blue", lwd.ideal=0.5)
  
  externalCalIntValues[[i]] = externalCal$stats[["Intercept"]]
  externalCalIntSEs[[i]] = (externalCal$stats[["Intercept"]] - externalCal$cl.interc[1])/1.96
  externalCalSlopeValues[[i]] = externalCal$stats[["Slope"]]
  externalCalSlopeSEs[[i]] = (externalCal$stats[["Slope"]] - externalCal$cl.slope[[1]])/1.96
}

#Correctly pooled Calibration intercept and SE
rubin.rules(unlist(externalCalIntValues), unlist(externalCalIntSEs))
#Correctly pooled Calibration slope and SE
rubin.rules(unlist(externalCalSlopeValues), unlist(externalCalSlopeSEs))

#combined data 
externalPreds = NULL
externalOutcomes = NULL
externalDUP = NULL
for(i in seq(1:tempData2$m))
{
  externalPreds = c(externalPreds,Results[[i]])
  externalOutcomes = c(externalOutcomes, as.character(Obs[[i]]))
  externalDUP = c(externalDUP, DUP[[i]])
}

#calibration plot - ignore confidence intervals, use above via Rubin's rules instead
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

#My custom functions
#####################################################################################
#Calibration Curves function also returning confidence intervals around intercept, slope and c-statistic
#Only confirmed to work with default options
val.prob.ci.3 <- function(p, y, logit, group, weights = rep(1, length(y)), normwt = F, pl = T,
                          smooth = c("loess","rcs",F), CL.smooth="fill",CL.BT=F,lty.smooth=1,col.smooth="black",lwd.smooth=1,
                          nr.knots=5,logistic.cal = F,lty.log=1,col.log="black",lwd.log=1, xlab = "Predicted probability", ylab =
                            "Observed proportion", xlim = c(-0.02, 1),ylim = c(-0.15,1), m, g, cuts, emax.lim = c(0, 1),
                          legendloc =  c(0.50 , 0.27), statloc = c(0,.85),dostats=T,cl.level=0.95,method.ci="pepe",roundstats=2,
                          riskdist = "predicted", cex=0.75,cex.leg = 0.75, connect.group =
                            F, connect.smooth = T, g.group = 4, evaluate = 100, nmin = 0, d0lab="0", d1lab="1", cex.d01=0.7,
                          dist.label=0.04, line.bins=-.05, dist.label2=.03, cutoff, las=1, length.seg=1,
                          y.intersp=1,lty.ideal=1,col.ideal="red",lwd.ideal=1,...)
{
  if(smooth[1]==F){smooth <- "F"}
  smooth <- match.arg(smooth)
  if(!missing(p))
    if(any(!(p>=0 | p<=1))){stop("Probabilities can not be > 1 or < 0.")}
  if(missing(p))
    p <- 1/(1 + exp( - logit))
  else logit <- log(p/(1 - p))
  if(!all(y%in%0:1)){stop("The vector with the binary outcome can only contain the values 0 and 1.")}
  if(length(p) != length(y))
    stop("lengths of p or logit and y do not agree")
  names(p) <- names(y) <- names(logit) <- NULL
  if(!missing(group)) {
    if(length(group) == 1 && is.logical(group) && group)
      group <- rep("", length(y))
    if(!is.factor(group))
      group <- if(is.logical(group) || is.character(group))
        as.factor(group) else cut2(group, g =
                                     g.group)
    names(group) <- NULL
    nma <- !(is.na(p + y + weights) | is.na(group))
    ng <- length(levels(group))
  }
  else {
    nma <- !is.na(p + y + weights)
    ng <- 0
  }
  logit <- logit[nma]
  y <- y[nma]
  p <- p[nma]
  if(ng > 0) {
    group <- group[nma]
    weights <- weights[nma]
    return(val.probg(p, y, group, evaluate, weights, normwt, nmin)
    )
  }
  
  # Sort vector with probabilities
  y     <- y[order(p)]
  logit <- logit[order(p)]
  p     <- p[order(p)]
  
  
  if(length(p)>5000 & smooth=="loess"){warning("Number of observations > 5000, RCS is recommended.",immediate. = T)}
  if(length(p)>1000 & CL.BT==T){warning("Number of observations is > 1000, this could take a while...",immediate. = T)}
  
  
  if(length(unique(p)) == 1) {
    #22Sep94
    P <- mean(y)
    Intc <- log(P/(1 - P))
    n <- length(y)
    D <- -1/n
    L01 <- -2 * sum(y * logit - log(1 + exp(logit)), na.rm = T)
    L.cal <- -2 * sum(y * Intc - log(1 + exp(Intc)), na.rm = T)
    U.chisq <- L01 - L.cal
    U.p <- 1 - pchisq(U.chisq, 1)
    U <- (U.chisq - 1)/n
    Q <- D - U
    
    stats <- c(0, 0.5, 0, D, 0, 1, U, U.chisq, U.p, Q, mean((y - p[
      1])^2), Intc, 0, rep(abs(p[1] - P), 2))
    names(stats) <- c("Dxy", "C (ROC)", "R2", "D", "D:Chi-sq",
                      "D:p", "U", "U:Chi-sq", "U:p", "Q", "Brier",
                      "Intercept", "Slope", "Emax", "Eavg", "ECI")
    return(stats)
  }
  i <- !is.infinite(logit)
  nm <- sum(!i)
  if(nm > 0)
    warning(paste(nm, "observations deleted from logistic calibration due to probs. of 0 or 1"))
  i.2 <- i
  f.or <- lrm(y[i]~logit[i])
  f <- lrm.fit(logit[i], y[i])
  cl.slope <- confint(f,level=cl.level)[2,]
  f2 <-	lrm.fit(offset=logit[i], y=y[i])
  if(f2$fail){
    warning("The lrm function did not converge when computing the calibration intercept!",immediate.=T)
    f2 <- list()
    f2$coef <- NA
    cl.interc <- rep(NA,2)
  }else{
    cl.interc <- confint(f2,level=cl.level)
  }
  stats <- f$stats
  cl.auc <- CalibrationCurves:::ci.auc(y,p,cl.level,method.ci)
  
  n <- stats["Obs"]
  predprob <- seq(emax.lim[1], emax.lim[2], by = 0.0005)
  lt <- f$coef[1] + f$coef[2] * log(predprob/(1 - predprob))
  calp <- 1/(1 + exp( - lt))
  emax <- max(abs(predprob - calp))
  if (pl) {
    plot(0.5, 0.5, xlim = xlim, ylim = ylim, type = "n", xlab = xlab,
         ylab = ylab, las=las,...)
    clip(0,1,0,1)
    abline(0, 1, lty = lty.ideal,col=col.ideal,lwd=lwd.ideal)
    do.call("clip", as.list(par()$usr))
    
    
    lt <- lty.ideal
    lw.d <- lwd.ideal
    all.col <- col.ideal
    leg <- "Ideal"
    marks <- -1
    if (logistic.cal) {
      lt <- c(lt, lty.log)
      lw.d <- c(lw.d,lwd.log)
      all.col <- c(all.col,col.log)
      leg <- c(leg, "Logistic calibration")
      marks <- c(marks, -1)
    }
    if(smooth!="F"){all.col <- c(all.col,col.smooth)}
    if (smooth=="loess") {
      #Sm <- lowess(p,y,iter=0)
      Sm <- loess(y~p,degree=2)
      Sm <- data.frame(Sm$x,Sm$fitted); Sm.01 <- Sm
      
      if (connect.smooth==T & CL.smooth!="fill") {
        clip(0,1,0,1)
        lines(Sm, lty = lty.smooth,lwd=lwd.smooth,col=col.smooth)
        do.call("clip", as.list(par()$usr))
        lt <- c(lt, lty.smooth)
        lw.d <- c(lw.d,lwd.smooth)
        marks <- c(marks, -1)
      }else if(connect.smooth==F & CL.smooth!="fill"){
        clip(0,1,0,1)
        points(Sm,col=col.smooth)
        do.call("clip", as.list(par()$usr))
        lt <- c(lt, 0)
        lw.d <- c(lw.d,1)
        marks <- c(marks, 1)
      }
      if(CL.smooth==T | CL.smooth=="fill"){
        to.pred <- seq(min(p),max(p),length=200)
        if(CL.BT==T){
          cat("Bootstrap samples are being generated.\n\n\n")
          
          replicate(2000,CalibrationCurves:::BT.samples(y,p,to.pred)) -> res.BT
          apply(res.BT,1,quantile,c(0.025,0.975)) -> CL.BT
          colnames(CL.BT) <- to.pred
          
          if(CL.smooth=="fill"){
            clip(0,1,0,1)
            polygon(x = c(to.pred, rev(to.pred)), y = c(CL.BT[2,],
                                                        rev(CL.BT[1,])),
                    col = rgb(177, 177, 177, 177, maxColorValue = 255), border = NA)
            if (connect.smooth==T) {
              lines(Sm, lty = lty.smooth,lwd=lwd.smooth,col=col.smooth)
              lt <- c(lt, lty.smooth)
              lw.d <- c(lw.d,lwd.smooth)
              marks <- c(marks, -1)
            }else if(connect.smooth==F){
              points(Sm,col=col.smooth)
              lt <- c(lt, 0)
              lw.d <- c(lw.d,1)
              marks <- c(marks, 1)
            }
            do.call("clip", as.list(par()$usr))
            leg <- c(leg, "Flexible calibration (Loess)")
          }else{
            
            clip(0,1,0,1)
            lines(to.pred,CL.BT[1,],lty=2,lwd=1,col=col.smooth);clip(0,1,0,1);lines(to.pred,CL.BT[2,],lty=2,lwd=1,col=col.smooth)
            do.call("clip", as.list(par()$usr))
            leg <- c(leg,"Flexible calibration (Loess)","CL flexible")
            lt <- c(lt,2)
            lw.d <- c(lw.d,1)
            all.col <- c(all.col,col.smooth)
            marks <- c(marks,-1)
          }
          
        }else{
          Sm.0 <- loess(y~p,degree=2)
          predict(Sm.0,type="fitted",se=T) -> cl.loess
          clip(0,1,0,1)
          if(CL.smooth=="fill"){
            polygon(x = c(Sm.0$x, rev(Sm.0$x)), y = c(cl.loess$fit+cl.loess$se.fit*1.96,
                                                      rev(cl.loess$fit-cl.loess$se.fit*1.96)),
                    col = rgb(177, 177, 177, 177, maxColorValue = 255), border = NA)
            if (connect.smooth==T) {
              lines(Sm, lty = lty.smooth,lwd=lwd.smooth,col=col.smooth)
              lt <- c(lt, lty.smooth)
              lw.d <- c(lw.d,lwd.smooth)
              marks <- c(marks, -1)
            }else if(connect.smooth==F){
              points(Sm,col=col.smooth)
              lt <- c(lt, 0)
              lw.d <- c(lw.d,1)
              marks <- c(marks, 1)
            }
            do.call("clip", as.list(par()$usr))
            leg <- c(leg, "Flexible calibration (Loess)")
          }else{
            lines(Sm.0$x,cl.loess$fit+cl.loess$se.fit*1.96,lty=2,lwd=1,col=col.smooth)
            lines(Sm.0$x,cl.loess$fit-cl.loess$se.fit*1.96,lty=2,lwd=1,col=col.smooth)
            do.call("clip", as.list(par()$usr))
            leg <- c(leg,"Flexible calibration (Loess)","CL flexible")
            lt <- c(lt,2)
            lw.d <- c(lw.d,1)
            all.col <- c(all.col,col.smooth)
            marks <- c(marks,-1)
          }
          
        }
        
      }else{
        leg <- c(leg, "Flexible calibration (Loess)")}
      cal.smooth <- approx(Sm.01, xout = p)$y
      eavg <- mean(abs(p - cal.smooth))
      ECI <- mean((p-cal.smooth)^2)*100
    }
    if(smooth=="rcs"){
      par(lwd=lwd.smooth,bty="n",col=col.smooth)
      if(!is.numeric(nr.knots)){stop("Nr.knots must be numeric.")}
      if(nr.knots==5){
        tryCatch(CalibrationCurves:::rcspline.plot(p,y,model="logistic",nk=5,show="prob", statloc = "none"
                               ,add=T,showknots=F,xrange=c(min(na.omit(p)),max(na.omit(p))),lty=lty.smooth),error=function(e){
                                 warning("The number of knots led to estimation problems, nk will be set to 4.",immediate. = T)
                                 tryCatch(CalibrationCurves:::rcspline.plot(p,y,model="logistic",nk=4,show="prob", statloc = "none"
                                                        ,add=T,showknots=F,xrange=c(min(na.omit(p)),max(na.omit(p))),lty=lty.smooth)
                                          ,error=function(e){
                                            warning("Nk 4 also led to estimation problems, nk will be set to 3.",immediate.=T)
                                            CalibrationCurves:::rcspline.plot(p,y,model="logistic",nk=3,show="prob", statloc = "none"
                                                          ,add=T,showknots=F,xrange=c(min(na.omit(p)),max(na.omit(p)))
                                                          ,lty=lty.smooth)
                                          })
                               })
      }else if(nr.knots==4){
        tryCatch(CalibrationCurves:::rcspline.plot(p,y,model="logistic",nk=4,show="prob", statloc = "none"
                               ,add=T,showknots=F,xrange=c(min(na.omit(p)),max(na.omit(p))),lty=lty.smooth),error=function(e){
                                 warning("The number of knots led to estimation problems, nk will be set to 3.",immediate.=T)
                                 CalibrationCurves:::rcspline.plot(p,y,model="logistic",nk=3,show="prob", statloc = "none"
                                               ,add=T,showknots=F,xrange=c(min(na.omit(p)),max(na.omit(p))),lty=lty.smooth)
                               })
      }else if(nr.knots==3){
        tryCatch(CalibrationCurves:::rcspline.plot(p,y,model="logistic",nk=3,show="prob", statloc = "none"
                               ,add=T,showknots=F,xrange=c(min(na.omit(p)),max(na.omit(p))),lty=lty.smooth),
                 error=function(e){
                   stop("Nk = 3 led to estimation problems.")
                 })
      }else{stop(paste("Number of knots = ",nr.knots,sep="", ", only 5 >= nk >=3 is allowed."))}
      
      par(lwd=1,bty="o",col="black")
      leg <- c(leg,"Flexible calibration (RCS)","CL flexible")
      lt <- c(lt,lty.smooth,2)
      lw.d <- c(lw.d,rep(lwd.smooth,2))
      all.col <- c(all.col,col.smooth)
      marks <- c(marks,-1,-1)
    }
    if(!missing(m) | !missing(g) | !missing(cuts)) {
      if(!missing(m))
        q <- cut2(p, m = m, levels.mean = T, digits = 7)
      else if(!missing(g))
        q <- cut2(p, g = g, levels.mean = T, digits = 7)
      else if(!missing(cuts))
        q <- cut2(p, cuts = cuts, levels.mean = T, digits = 7)
      means <- as.single(levels(q))
      prop <- tapply(y, q, function(x)mean(x, na.rm = T))
      points(means, prop, pch = 2, cex=1)
      #18.11.02: CI triangles
      ng	<-tapply(y, q, length)
      og	<-tapply(y, q, sum)
      ob	<-og/ng
      se.ob	<-sqrt(ob*(1-ob)/ng)
      g		<- length(as.single(levels(q)))
      
      for (i in 1:g) lines(c(means[i], means[i]), c(prop[i],min(1,prop[i]+1.96*se.ob[i])), type="l")
      for (i in 1:g) lines(c(means[i], means[i]), c(prop[i],max(0,prop[i]-1.96*se.ob[i])), type="l")
      
      if(connect.group) {
        lines(means, prop)
        lt <- c(lt, 1)
        lw.d <- c(lw.d,1)
      }
      else {
        lt <- c(lt, 0)
        lw.d <- c(lw.d, 0)
      }
      leg <- c(leg, "Grouped observations")
      all.col <- c(all.col, col.smooth)
      marks <- c(marks, 2)
    }
  }
  lr <- stats["Model L.R."]
  p.lr <- stats["P"]
  D <- (lr - 1)/n
  L01 <- -2 * sum(y * logit - logb(1 + exp(logit)), na.rm = TRUE)
  U.chisq <- L01 - f$deviance[2]
  p.U <- 1 - pchisq(U.chisq, 2)
  U <- (U.chisq - 2)/n
  Q <- D - U
  Dxy <- stats["Dxy"]
  C <- stats["C"]
  R2 <- stats["R2"]
  B <- sum((p - y)^2)/n
  # ES 15dec08 add Brier scaled
  Bmax  <- mean(y) * (1-mean(y))^2 + (1-mean(y)) * mean(y)^2
  Bscaled <- 1 - B/Bmax
  stats <- c(Dxy, C, R2, D, lr, p.lr, U, U.chisq, p.U, Q, B,
             f2$coef[1], f$coef[2], emax, Bscaled)
  names(stats) <- c("Dxy", "C (ROC)", "R2", "D", "D:Chi-sq",
                    "D:p", "U", "U:Chi-sq", "U:p", "Q", "Brier", "Intercept",
                    "Slope", "Emax", "Brier scaled")
  if(smooth=="loess")
    stats <- c(stats, c(Eavg = eavg),c(ECI = ECI))
  
  # Cut off definition
  if(!missing(cutoff)) {
    arrows(x0=cutoff,y0=.1,x1=cutoff,y1=-0.025,length=.15)
  }
  if(pl) {
    if(min(p)>plogis(-7) | max(p)<plogis(7)){
      
      lrm(y[i.2]~qlogis(p[i.2]))-> lrm.fit.1
      if(logistic.cal)  lines(p[i.2],plogis(lrm.fit.1$linear.predictors),lwd=lwd.log,lty=lty.log,col=col.log)
      
    }else{logit <- seq(-7, 7, length = 200)
    prob <- 1/(1 + exp( - logit))
    pred.prob <- f$coef[1] + f$coef[2] * logit
    pred.prob <- 1/(1 + exp( - pred.prob))
    if(logistic.cal) lines(prob, pred.prob, lty=lty.log,lwd=lwd.log,col=col.log)
    }
    #	pc <- rep(" ", length(lt))
    #	pc[lt==0] <- "."
    lp <- legendloc
    if (!is.logical(lp)) {
      if (!is.list(lp))
        lp <- list(x = lp[1], y = lp[2])
      legend(lp, leg, lty = lt, pch = marks, cex = cex.leg, bty = "n",lwd=lw.d,
             col=all.col,y.intersp = y.intersp)
    }
    if(!is.logical(statloc)) {
      if(dostats[1]==T){
        stats.2 <- paste('Calibration\n',
                         '...intercept: '
                         , sprintf(paste("%.",roundstats,"f",sep=""), stats["Intercept"]), " (",
                         sprintf(paste("%.",roundstats,"f",sep=""),cl.interc[1])," to ",
                         sprintf(paste("%.",roundstats,"f",sep=""),cl.interc[2]),")",'\n',
                         '...slope: '
                         , sprintf(paste("%.",roundstats,"f",sep=""), stats["Slope"]), " (",
                         sprintf(paste("%.",roundstats,"f",sep=""),cl.slope[1])," to ",
                         sprintf(paste("%.",roundstats,"f",sep=""),cl.slope[2]),")",'\n',
                         'Discrimination\n',
                         '...c-statistic: '
                         , sprintf(paste("%.",roundstats,"f",sep=""), stats["C (ROC)"]), " (",
                         sprintf(paste("%.",roundstats,"f",sep=""),cl.auc[2])," to ",
                         sprintf(paste("%.",roundstats,"f",sep=""),cl.auc[3]),")"
                         , sep = '')
        text(statloc[1], statloc[2],stats.2,pos=4,cex=cex)
        
      }else{
        dostats <- dostats
        leg <- format(names(stats)[dostats])	#constant length
        leg <- paste(leg, ":", format(stats[dostats], digits=roundstats), sep =
                       "")
        if(!is.list(statloc))
          statloc <- list(x = statloc[1], y = statloc[2])
        text(statloc, paste(format(names(stats[dostats])),
                            collapse = "\n"), adj = 0, cex = cex)
        text(statloc$x + (xlim[2]-xlim[1])/3 , statloc$y, paste(
          format(round(stats[dostats], digits=roundstats)), collapse =
            "\n"), adj = 1, cex = cex)
      }
    }
    if(is.character(riskdist)) {
      if(riskdist == "calibrated") {
        x <- f$coef[1] + f$coef[2] * log(p/(1 - p))
        x <- 1/(1 + exp( - x))
        x[p == 0] <- 0
        x[p == 1] <- 1
      }
      else x <- p
      bins <- seq(0, min(1,max(xlim)), length = 101)
      x <- x[x >= 0 & x <= 1]
      #08.04.01,yvon: distribution of predicted prob according to outcome
      f0	<-table(cut(x[y==0],bins))
      f1	<-table(cut(x[y==1],bins))
      j0	<-f0 > 0
      j1	<-f1 > 0
      bins0 <-(bins[-101])[j0]
      bins1 <-(bins[-101])[j1]
      f0	<-f0[j0]
      f1	<-f1[j1]
      maxf <-max(f0,f1)
      f0	<-(0.1*f0)/maxf
      f1	<-(0.1*f1)/maxf
      
      segments(bins1,line.bins,bins1,length.seg*f1+line.bins)
      segments(bins0,line.bins,bins0,length.seg*-f0+line.bins)
      lines(c(min(bins0,bins1)-0.01,max(bins0,bins1)+0.01),c(line.bins,line.bins))
      text(max(bins0,bins1)+dist.label,line.bins+dist.label2,d1lab,cex=cex.d01)
      text(max(bins0,bins1)+dist.label,line.bins-dist.label2,d0lab,cex=cex.d01)
      
    }
  }
  if(dostats==T){
    cat(paste("\n\n A ",cl.level*100,
              "% confidence interval is given for the calibration intercept, calibration slope and c-statistic. \n\n",
              sep=""))}
  
  stats_ci <- list("stats" = stats, "cl.interc" = cl.interc, "cl.slope" = cl.slope, "cl.auc" = cl.auc)
  return(stats_ci)
}

#function to calculate rubin's rules
#https://bookdown.org/mwheymans/bookmi/rubins-rules.html
rubin.rules <- function(means, SEs)
{
  n = length(SEs)
  rubin_mean = mean(means)
  variance_within = (sum(SEs^2))/n
  variance_between = (sum((means-rubin_mean)^2))/(n-1)
  variance_total = variance_within + variance_between + variance_between/n
  rubin_se = sqrt(variance_total)
  
  rubins = list("rubin_mean" = rubin_mean, "rubin_se" = rubin_se)
  return(rubins)
}
