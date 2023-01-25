# Development and validation of a non-remission risk prediction model in First Episode Psychosis: An analysis of two longitudinal studies.

This is the R code for our manuscript.

Comments are inline in the code.R and should be easy to follow. Any issues contact the authors.

We also upload a corrected version of Table 1 from the manuscript detailing the final logistic regression nonremission prediction model specification. We now provide mean and standard deviation values to allow the transformation of the predictor variables to Z-scores for their use in the model. This was omitted from the original published paper but is required to apply the model to new patients.

We would like to acknowledge the authors of the following packages which were required for our analysis:
* https://CRAN.R-project.org/package=doParallel
* https://CRAN.R-project.org/package=readr
* https://CRAN.R-project.org/package=caret
* https://CRAN.R-project.org/package=pROC
* https://CRAN.R-project.org/package=gtools
* https://CRAN.R-project.org/package=mice
* https://CRAN.R-project.org/package=pmsampsize
* https://github.com/BavoDC/CalibrationCurves
* https://github.com/ddsjoberg/dcurves
* https://CRAN.R-project.org/package=psfmi
