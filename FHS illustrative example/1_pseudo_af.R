# Pseudo observation models for AF with logit link & LSmeans-type predictions
# Sarah Conner - Nov 23 2020

library(sas7bdat)
library(geepack)
library(survival)
library(etm)
library(Rcpp)
library(tidyverse)

options(scipen=999)
af <- read.csv('af_pseudo.csv')
#af$entry_alc_elev <- ifelse(af$entry_alc>14 & af$male==1 | af$entry_alc>7 & af$male==0, 1, 0)

# ----- Fit Pseudo Model for AF with logit link -----

logit.fit <- geese(pseudo.af ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi,
                   data=af, id=fid, jack = TRUE, scale.fix=TRUE, family=gaussian,
                   mean.link = "logit", corstr="independence")
summary(logit.fit)

logit.res <- cbind(mean = logit.fit$beta, 
                   san.se = sqrt(diag(logit.fit$vbeta)),
                   san.pval = 2-2*pnorm(abs(logit.fit$beta/sqrt(diag(logit.fit$vbeta)))),
                   or = exp(logit.fit$beta),
                   cil = exp(logit.fit$beta - 1.96*sqrt(diag(logit.fit$vbeta))),
                   ciu = exp(logit.fit$beta + 1.96*sqrt(diag(logit.fit$vbeta))))

logit.res <- logit.res[2:12,]
logit.res[2,4:6] <- logit.res[2,4:6]^sd(af$entry_hgt)
logit.res[3,4:6] <- logit.res[3,4:6]^sd(af$entry_wgt)
logit.res[4,4:6] <- logit.res[4,4:6]^sd(af$entry_sbp)
logit.res[5,4:6] <- logit.res[5,4:6]^sd(af$entry_dbp)

#logit.fmt <- cbind(vars, ltrdiff=paste0(sprintf('%.2f', logit.res[,4]), " (", sprintf('%.2f', logit.res[,5]), ", ", sprintf('%.2f', logit.res[,6]), ")"))
#write.csv(logit.fmt, 'results//pseudologit_af_mod.csv', row.names=FALSE)



# ----- Predict lifetime risk and difference for each risk factor-----

logit.beta <- unname(logit.fit$beta)

vars <- c('male', 'entry_hgt', 'entry_wgt', 'entry_sbp', 'entry_dbp', 'entry_hrx', 'entry_smk', 'entry_alc_elev', 'entry_diab', 'entry_pchf', 'entry_pmi')
binary <- c(1, rep(0,4), rep(1,6))

average.pred <- c(1, unname(apply(af[, vars], 2, mean)))
sd.pred <- c(1, unname(apply(af[, vars], 2, sd)))
average1sd.pred <- average.pred + sd.pred

lp0 <- NULL
ltr0 <- NULL
grad0 <- NULL
ltr0.se <- NULL

lp1 <- NULL
ltr1 <- NULL
grad1 <- NULL
ltr1.se <- NULL

graddiff <- NULL
ltrdiff <- NULL
ltrdiff.se <- NULL


for(i in 2:12){
  
  # Covariate of interest = 0 if binary, average if continuous. All other covariates set to average
  average.pred0 <- average.pred
  if (binary[i-1]==1) average.pred0[i] <- 0
  lp0 <- logit.beta %*% average.pred0
  ltr0[i] <- 1/(1 + exp(-lp0))
  
  # Covariate of interest = 1 if binary, average+1 if continuous. All other covariates set to average
  average.pred1 <- average.pred
  if (binary[i-1]==1) {
    average.pred1[i] <- 1
  } else {
    average.pred1[i] <- average1sd.pred[i]
  }
  lp1 <- logit.beta %*% average.pred1
  ltr1[i] <- 1/(1 + exp(-lp1))
  
  # Difference per each covariate
  ltrdiff[i] <- ltr1[i] - ltr0[i]
  
  # Delta method for SEs
  
  # Gradient
  for(j in 1:12){
    grad0[j] <- (average.pred0[j]*exp(-lp0))/(1 + exp(-lp0))^2
    grad1[j] <- (average.pred1[j]*exp(-lp1))/(1 + exp(-lp1))^2
    graddiff[j] <- grad1[j] - grad0[j]
  }
  
  ltr0.se[i] <- sqrt(t(grad0) %*% logit.fit$vbeta %*% grad0)
  ltr1.se[i] <- sqrt(t(grad1) %*% logit.fit$vbeta %*% grad1)
  ltrdiff.se[i] <- sqrt(t(graddiff) %*% logit.fit$vbeta %*% graddiff)
}

ltrpreds <- cbind(ltrdiff, ltrdiff-1.96*ltrdiff.se, ltrdiff+1.96*ltrdiff.se)
ltrpreds.fmt <- cbind(vars, paste0(sprintf('%.2f', ltrpreds[2:12,1]), " (", sprintf('%.2f', ltrpreds[2:12,2]), ", ", sprintf('%.2f', ltrpreds[2:12,3]), ")"))

#write.csv(ltrpreds, 'results//pseudologit_af_pred.csv', row.names=FALSE)


#  ----- Predictions for each individual (calibration) ----- 

beta <- unname(logit.fit$beta)
subjects <- as.matrix(unname(cbind(rep(1, nrow(af)), af[, vars]))) # add vector of 1s for intercept 
lp <- subjects %*% beta
af$pred.ltr <- 1/(1+exp(-lp))


write.csv(af$pred.ltr, "results//pseudo_af_individual_ltr.csv", row.names=FALSE)
