# Pseudo observation models for Death without AF with logit link & LSmeans-type predictions
# Sarah Conner - Nov 23 2020

library(sas7bdat)
library(pseudo)
library(geepack)
library(survival)
library(etm)
library(Rcpp)
library(tidyverse)

options(scipen=999)
af <- read.csv('af_pseudo.csv')

# ----- Fit Pseudo Model for Death without AF, logit link -----

# Logit link 
logit.fit <- geese(pseudo.death ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi,
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

logit.fmt <- cbind(vars, ltrdiff=paste0(sprintf('%.2f', logit.res[,4]), " (", sprintf('%.2f', logit.res[,5]), ", ", sprintf('%.2f', logit.res[,6]), ")"))
write.csv(logit.fmt, 'results//pseudologit_death_mod.csv', row.names=FALSE)


