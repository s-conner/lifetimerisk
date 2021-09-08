# Pseudo observation models for Death without AF with logit link & LSmeans-type predictions
# Sarah Conner - updated 8/7/2021

library(sas7bdat)
library(pseudo)
library(geepack)
library(survival)
library(etm)
library(Rcpp)
library(tidyverse)

options(scipen=999)
af <- read.csv('af_pseudo.csv')


# ----- Check Pseudo Models for Death without AF with sex interactions -----

logit.fit <- geese(pseudo.death ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi,
                   data=af, id=fid, jack = TRUE, scale.fix=TRUE, family=gaussian,
                   mean.link = "logit", corstr="independence")
summary(logit.fit)

qic(logit.fit, af, af$pseudo.death) # QIC 5947, QICu 5951

int1 <- c('entry_hgt', 'entry_wgt', 'entry_sbp', 'entry_dbp', 'entry_hrx',
          'entry_smk', 'entry_alc_elev', 'entry_diab', 'entry_pchf', 'entry_pmi')
n.vars.int1 <- length(int1)
int.qic1 <- matrix(NA, nrow=n.vars.int1, ncol=3)
j <- 0

# Create sex interactions
for(i in int1){
  j <- j+1
  af[[18+j]] <- af[[i]] * af$male
}

intnames <- c('male_hgt', 'male_wgt', 'male_sbp', 'male_dbp', 'male_hrx',
              'male_smk', 'male_alc_elev', 'male_diab', 'male_pchf', 'male_pmi')
colnames(af)[19:28] <- intnames

logit.fit1 <- geese(pseudo.death ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male_hgt,
                     data=af, id=fid, jack = TRUE, scale.fix=TRUE, family=gaussian, mean.link = "logit", corstr="independence")
logit.fit2 <- geese(pseudo.death ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male_wgt,
                    data=af, id=fid, jack = TRUE, scale.fix=TRUE, family=gaussian, mean.link = "logit", corstr="independence")
logit.fit3 <- geese(pseudo.death ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male_sbp,
                    data=af, id=fid, jack = TRUE, scale.fix=TRUE, family=gaussian, mean.link = "logit", corstr="independence")
logit.fit4 <- geese(pseudo.death ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male_dbp,
                    data=af, id=fid, jack = TRUE, scale.fix=TRUE, family=gaussian, mean.link = "logit", corstr="independence")
logit.fit5 <- geese(pseudo.death ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male_hrx,
                    data=af, id=fid, jack = TRUE, scale.fix=TRUE, family=gaussian, mean.link = "logit", corstr="independence")
logit.fit6 <- geese(pseudo.death ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male_smk,
                    data=af, id=fid, jack = TRUE, scale.fix=TRUE, family=gaussian, mean.link = "logit", corstr="independence")
logit.fit7 <- geese(pseudo.death ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male_alc_elev,
                    data=af, id=fid, jack = TRUE, scale.fix=TRUE, family=gaussian, mean.link = "logit", corstr="independence")
logit.fit8 <- geese(pseudo.death ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male_diab,
                    data=af, id=fid, jack = TRUE, scale.fix=TRUE, family=gaussian, mean.link = "logit", corstr="independence")
logit.fit9 <- geese(pseudo.death ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male_pchf,
                    data=af, id=fid, jack = TRUE, scale.fix=TRUE, family=gaussian, mean.link = "logit", corstr="independence")
logit.fit10 <- geese(pseudo.death ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male_pmi,
                    data=af, id=fid, jack = TRUE, scale.fix=TRUE, family=gaussian, mean.link = "logit", corstr="independence")

# QICu (AIC-like) and QIC (BIC-like) increase
do.call("rbind", list(qic(logit.fit, af, af$pseudo.death),
                      qic(logit.fit1, af, af$pseudo.death),
                      qic(logit.fit2, af, af$pseudo.death),
                      qic(logit.fit3, af, af$pseudo.death),
                      qic(logit.fit4, af, af$pseudo.death),
                      qic(logit.fit5, af, af$pseudo.death),
                      qic(logit.fit6, af, af$pseudo.death),
                      qic(logit.fit8, af, af$pseudo.death),
                      qic(logit.fit9, af, af$pseudo.death)))

# MI decreased QIC by ~2, not substantial

# P-values all very non-significant
c(summary(logit.fit1)$mean[13,5],
  summary(logit.fit2)$mean[13,5],
  summary(logit.fit3)$mean[13,5],
  summary(logit.fit4)$mean[13,5],
  summary(logit.fit5)$mean[13,5],
  summary(logit.fit6)$mean[13,5],
  summary(logit.fit8)$mean[13,5],
  summary(logit.fit9)$mean[13,5])


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


