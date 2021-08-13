### Fine-Gray models for FHS - Lifetime risk of AF
### Mstate package to prepare data, no errors in R v. 4.0.2
### Include interactions with log(time), determine using AIC
### Updated analyses 1/12/2021

library(survival)
library(mstate)
options(scipen=999)

# ----- Data prep -----

af <- read.csv('af_pseudo.csv')
vars <- c('male', 'entry_hgt', 'entry_wgt', 'entry_sbp', 'entry_dbp', 'entry_hrx', 
          'entry_smk', 'entry_alc_elev', 'entry_diab', 'entry_pchf', 'entry_pmi')
means <- apply(af[, vars], 2, mean)
sds <- apply(af[, vars], 2, sd)

# Weighted dataset
af.fg <- crprep(Tstop="exitage", status="event", data=af, trans=1, cens=0, Tstart="entryage_55", id="fid", keep=vars)
af.fg$weight_ltrc <- af.fg$weight.cens*af.fg$weight.trunc
af.fg$ev1 <- ifelse(af.fg$status==1,1,0)


# ----- PH model ----- 

fit.fg <- coxph(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                  entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, weight=weight_ltrc, data=af.fg)
summary(fit.fg)

# AIC
2*length(coef(fit.fg)) - 2*fit.fg$loglik[2]


# ----- Round 1: add 1 log(time) interaction at a time ----- 

n.vars1 <- length(vars)
nonph.pval1 <- rep(NA, n.vars1)
nonph.aic1 <- rep(NA, n.vars1)
j <- 0

for(i in vars){
  j <- j+1
  fit.fg.int <- coxph(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                        entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + tt(get(i)), 
                      tt = function(x, t, ...) x * log(t-55),
                      weight=weight_ltrc, data=af.fg)
  nonph.pval1[j] <- summary(fit.fg.int)$coef[n.vars1+1,5]
  nonph.aic1[j] <- 2*length(coef(fit.fg.int)) - 2*fit.fg.int$loglik[2]
}

nonph1 <- data.frame(vars, pval=round(nonph.pval1,6), aic=nonph.aic1)
nonph1[order(nonph1$aic),]


# ----- Round 2: tt(entry_wgt) + add 1 log(time) interaction at a time 

vars2 <- vars[which(vars!='entry_wgt')]
n.vars2 <- length(vars2)
nonph.pval2 <- rep(NA, n.vars2)
nonph.aic2 <- rep(NA, n.vars2)
j <- 0

for(i in vars2){
  j <- j+1
  fit.fg.int <- coxph(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                        entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + tt(entry_wgt) + tt(get(i)), 
                      tt = function(x, t, ...) x * log(t-55),
                      weight=weight_ltrc, data=af.fg)
  nonph.pval2[j] <- summary(fit.fg.int)$coef[n.vars2+1,5]
  nonph.aic2[j] <- 2*length(coef(fit.fg.int)) - 2*fit.fg.int$loglik[2]
}

nonph2 <- data.frame(vars2, pval=round(nonph.pval2,6), aic=nonph.aic2)
nonph2[order(nonph2$aic),]


# ----- Round 3: tt(entry_wgt) + tt(entry_pmi) + add 1 log(time) interaction at a time 

vars3 <- vars2[which(vars2!='entry_pmi')]
n.vars3 <- length(vars3)
nonph.pval3 <- rep(NA, n.vars3)
nonph.aic3 <- rep(NA, n.vars3)
j <- 0

for(i in vars3){
  j <- j+1
  fit.fg.int <- coxph(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                        entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + tt(entry_wgt) + tt(entry_pmi) + tt(get(i)), 
                      tt = function(x, t, ...) x * log(t-55),
                      weight=weight_ltrc, data=af.fg)
  nonph.pval3[j] <- summary(fit.fg.int)$coef[n.vars3+1,5]
  nonph.aic3[j] <- 2*length(coef(fit.fg.int)) - 2*fit.fg.int$loglik[2]
}

nonph3 <- data.frame(vars3, pval=round(nonph.pval3,6), aic=nonph.aic3)
nonph3[order(nonph3$aic),]



# ----- Round 4: tt(entry_wgt) + tt(entry_pmi) + tt(entry_smk) + add 1 log(time) interaction at a time 

vars4 <- vars3[which(vars3!='entry_smk')]
n.vars4 <- length(vars4)
nonph.pval4 <- rep(NA, n.vars4)
nonph.aic4 <- rep(NA, n.vars4)
j <- 0

for(i in vars4){
  j <- j+1
  fit.fg.int <- coxph(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                        entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + 
                        tt(entry_wgt) + tt(entry_pmi) + tt(entry_smk) + tt(get(i)), 
                      tt = function(x, t, ...) x * log(t-55),
                      weight=weight_ltrc, data=af.fg)
  nonph.pval4[j] <- summary(fit.fg.int)$coef[n.vars4+1,5]
  nonph.aic4[j] <- 2*length(coef(fit.fg.int)) - 2*fit.fg.int$loglik[2]
}

nonph4 <- data.frame(vars4, pval=round(nonph.pval4,6), aic=nonph.aic4)
nonph4[order(nonph4$aic),]



# ----- Round 5:  tt(entry_wgt) + tt(entry_pmi) + tt(entry_smk) + tt(entry_hgt) + add 1 log(time) interaction at a time 

vars5 <- vars4[which(vars4!='entry_hgt')]
n.vars5 <- length(vars5)
nonph.pval5 <- rep(NA, n.vars5)
nonph.aic5 <- rep(NA, n.vars5)
j <- 0

for(i in vars5){
  j <- j+1
  fit.fg.int <- coxph(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                        entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + 
                        tt(entry_wgt) + tt(entry_pmi) + tt(entry_smk) + tt(entry_hgt) + tt(get(i)), 
                      tt = function(x, t, ...) x * log(t-55),
                      weight=weight_ltrc, data=af.fg)
  nonph.pval5[j] <- summary(fit.fg.int)$coef[n.vars5+1,5]
  nonph.aic5[j] <- 2*length(coef(fit.fg.int)) - 2*fit.fg.int$loglik[2]
}

nonph5 <- data.frame(vars5, pval=round(nonph.pval5,6), aic=nonph.aic5)
nonph5[order(nonph5$aic),]


# ----- Round 6: tt(entry_wgt) + tt(entry_pmi) + tt(entry_smk) + tt(entry_hgt) + tt(entry_diab) + add 1 log(time) interaction at a time 

vars6 <- vars5[which(vars5!='entry_diab')]
n.vars6 <- length(vars6)
nonph.pval6 <- rep(NA, n.vars6)
nonph.aic6 <- rep(NA, n.vars6)
j <- 0

for(i in vars6){
  j <- j+1
  fit.fg.int <- coxph(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                        entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + 
                        tt(entry_wgt) + tt(entry_pmi) + tt(entry_smk) + tt(entry_hgt) + tt(entry_diab) + tt(get(i)), 
                      tt = function(x, t, ...) x * log(t-55),
                      weight=weight_ltrc, data=af.fg)
  nonph.pval6[j] <- summary(fit.fg.int)$coef[n.vars6+1,5]
  nonph.aic6[j] <- 2*length(coef(fit.fg.int)) - 2*fit.fg.int$loglik[2]
}

nonph6 <- data.frame(vars6, pval=round(nonph.pval6,6), aic=nonph.aic6)
nonph6[order(nonph6$aic),]

# AIC increased, proceed with Round 5 model!


