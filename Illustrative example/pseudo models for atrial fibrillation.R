#############################################################################################
# Code for illustrative example analyses in the Framingham Heart Study

# Derives the pseudo-observations for AF and death without AF and
# fit logistic models for pseudo-observation.

# Requires 'Rcpp' package and 'pseudo_ltr_cpp' file of Rcpp functions for the lifetime risk, 
# 'pseudo.lr' function for pseudo-observations, 
# and 'geepack' package for 'geese' function.

# 'tidyverse' package only used for presenting results/subject predictions
#############################################################################################


library(geepack)
library(Rcpp)
library(tidyverse)

options(scipen=999)
source('\\\\bumc-y.bu.edu\\BUMC\\SPH\\Projects\\FHS_AF\\19SCR-lifetimerisk_varselection\\R programs\\pseudo left truncation\\pseudo_ltr_function.cpp')
sourceCpp('\\\\bumc-y.bu.edu\\BUMC\\SPH\\Projects\\FHS_AF\\19SCR-lifetimerisk_varselection\\R programs\\pseudo left truncation\\pseudo_ltr_cpp.cpp')
respath <- 'Y:\\19SCR-lifetimerisk_varselection\\FHS example\\'
af <- read.csv('Y:\\19SCR-lifetimerisk_varselection\\FHS example\\af.csv')



####### Derive pseudo-observations for AF and death without AF #######
af$pseudo.af <- pseudo.lr(af, 'startage_55', 'endage_55', 'event_55', 1, 95) 
af$pseudo.death <- pseudo.lr(af, 'startage_55', 'endage_55', 'event_55', 2, 95) 



####### Pseudo Models for AF ########

# Logit link 
logit.fit <- geese(pseudo.af ~ startage_55 + male + sbp_55 + dbp_55 + hrx_55 + currsmk_55 + 
                     alco_elevated + bmi_55 + hx_diab_55 + pchf_55 + pmi_55, 
                   data=af, id=id, jack = TRUE, scale.fix=TRUE, family=gaussian,
                   mean.link = "logit", corstr="independence")

logit.res <- cbind(mean = logit.fit$beta, 
                   san.se = sqrt(diag(logit.fit$vbeta)),
                   san.pval = 2-2*pnorm(abs(logit.fit$beta/sqrt(diag(logit.fit$vbeta)))),
                   or = exp(logit.fit$beta),
                   cil = exp(logit.fit$beta - 1.96*sqrt(diag(logit.fit$vbeta))),
                   ciu = exp(logit.fit$beta + 1.96*sqrt(diag(logit.fit$vbeta))))

# Arrange results, obtain ORs for 1 SD for continuous vars
logit.res <- logit.res[2:12,]
logit.res[1,4:6] <- logit.res[1,4:6]^sd(af$startage_55)
logit.res[3,4:6] <- logit.res[3,4:6]^sd(af$sbp_55)
logit.res[4,4:6] <- logit.res[4,4:6]^sd(af$dbp_55)
logit.res[8,4:6] <- logit.res[8,4:6]^sd(af$bmi_55)

logit.tab <- cbind(paste0(round(logit.res[,4], 2), " (", 
                          round(logit.res[,5], 2), ", ", 
                          round(logit.res[,6], 2), ")"),
                   round(logit.res[,3], 4))

## Calibration of AF model
# Subject-level predictions
beta <- unname(logit.fit$beta)
subjects <- as.matrix(unname(cbind(rep(1, nrow(af)), 
                  af[, c("startage_55", "male", "sbp_55", "dbp_55", "hrx_55", "currsmk_55", "alco_elevated", 
                         "bmi_55", "hx_diab_55", "pchf_55", "pmi_55")])))
lp <- subjects %*% beta
af$pred.ltr <- 1/(1+exp(-lp))
summary(af$pred.ltr)

# Divide into deciles
af <- mutate(af, quantile = ntile(pred.ltr, 10))





####### Pseudo models for death ####### 
logit.fit2 <- geese(pseudo.death ~ startage_55 + male + sbp_55 + dbp_55 + hrx_55 + currsmk_55 + 
                      alco_elevated + bmi_55 + hx_diab_55 + pchf_55 + pmi_55, 
                   data=af, id=id, jack = TRUE, scale.fix=TRUE, family=gaussian,
                   mean.link = "logit", corstr="independence")

logit.res2 <- cbind(mean = logit.fit2$beta, 
                   san.se = sqrt(diag(logit.fit2$vbeta)),
                   san.pval = 2-2*pnorm(abs(logit.fit2$beta/sqrt(diag(logit.fit2$vbeta)))),
                   or = exp(logit.fit2$beta),
                   cil = exp(logit.fit2$beta - 1.96*sqrt(diag(logit.fit2$vbeta))),
                   ciu = exp(logit.fit2$beta + 1.96*sqrt(diag(logit.fit2$vbeta))))

# Arrange results, obtain ORs for 1 SD for continuous vars
logit.res2 <- logit.res2[2:12,]
logit.res2[1,4:6] <- logit.res2[1,4:6]^sd(af$startage_55)
logit.res2[3,4:6] <- logit.res2[3,4:6]^sd(af$sbp_55)
logit.res2[4,4:6] <- logit.res2[4,4:6]^sd(af$dbp_55)
logit.res2[8,4:6] <- logit.res2[8,4:6]^sd(af$bmi_55)

logit.tab2 <- cbind(paste0(round(logit.res2[,4], 2), " (", 
                          round(logit.res2[,5], 2), ", ", 
                          round(logit.res2[,6], 2), ")"),
                   round(logit.res2[,3], 4))


