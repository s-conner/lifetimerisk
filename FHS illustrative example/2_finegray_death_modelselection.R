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
af.fg <- crprep(Tstop="exitage", status="event", data=af, trans=2, cens=0, Tstart="entryage_55", id="fid", keep=vars, shorten=FALSE)
af.fg$weight_ltrc <- af.fg$weight.cens*af.fg$weight.trunc
af.fg$ev2 <- ifelse(af.fg$status==2,1,0)


# ----- PH model ----- 

fit.fg <- coxph(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                  entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, weight=weight_ltrc, data=af.fg)
summary(fit.fg)

# AIC 19324.94
full.aic0 <- AIC(fit.fg)


# ----- Round 1: add 1 log(time) interaction at a time ----- 

n.vars1 <- length(vars)
nonph.pval1 <- rep(NA, n.vars1)
nonph.aic1 <- rep(NA, n.vars1)
j <- 0

for(i in vars){
  j <- j+1
  fit.fg.int <- coxph(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                        entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + tt(get(i)), 
                      tt = function(x, t, ...) x * log(t-55),
                      weight=weight_ltrc, data=af.fg)
  nonph.aic1[j] <- 2*length(coef(fit.fg.int)) - 2*fit.fg.int$loglik[2]
}

nonph1 <- data.frame(vars, aic=nonph.aic1)
nonph1$delta <- full.aic0 - nonph1$aic
nonph1[order(nonph1$aic),]
# lowest AIC entry_pchf 19308.89, delta=16
full.aic1 <- min(nonph1$aic)


# ----- Round 2: tt(entry_pchf) + add 1 log(time) interaction at a time 

vars2 <- vars[which(vars!='entry_pchf')]
n.vars2 <- length(vars2)
nonph.pval2 <- rep(NA, n.vars2)
nonph.aic2 <- rep(NA, n.vars2)
j <- 0

for(i in vars2){
  j <- j+1
  fit.fg.int <- coxph(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                        entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + tt(entry_pchf) + tt(get(i)), 
                      tt = function(x, t, ...) x * log(t-55),
                      weight=weight_ltrc, data=af.fg)
  nonph.aic2[j] <- 2*length(coef(fit.fg.int)) - 2*fit.fg.int$loglik[2]
}

nonph2 <- data.frame(vars2, aic=nonph.aic2)
nonph2$delta <- full.aic1 - nonph2$aic
nonph2[order(nonph2$aic),]

# lowest AIC male 19296.07, delta=12.8
full.aic2 <- min(nonph2$aic)


# ----- Round 3: tt(entry_pchf) + tt(male) + add 1 log(time) interaction at a time 

vars3 <- vars2[which(vars2!='male')]
n.vars3 <- length(vars3)
nonph.pval3 <- rep(NA, n.vars3)
nonph.aic3 <- rep(NA, n.vars3)
j <- 0

for(i in vars3){
  j <- j+1
  fit.fg.int <- coxph(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                        entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + tt(entry_pchf) + tt(male) + tt(get(i)), 
                      tt = function(x, t, ...) x * log(t-55),
                      weight=weight_ltrc, data=af.fg)
  nonph.aic3[j] <- 2*length(coef(fit.fg.int)) - 2*fit.fg.int$loglik[2]
}

nonph3 <- data.frame(vars3, aic=nonph.aic3)
nonph3$delta <- full.aic2 - nonph3$aic
nonph3[order(nonph3$aic),]

# lowest AIC 19289.58 entry_smk, delta=6.5 
full.aic3 <- min(nonph3$aic)


# ----- Round 4: tt(entry_pchf) + tt(male) + tt(entry_smk) + add 1 log(time) interaction at a time 

vars4 <- vars3[which(vars3!='entry_smk')]
n.vars4 <- length(vars4)
nonph.pval4 <- rep(NA, n.vars4)
nonph.aic4 <- rep(NA, n.vars4)
j <- 0

for(i in vars4){
  j <- j+1
  fit.fg.int <- coxph(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                        entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + 
                        tt(entry_pchf) + tt(male) + tt(entry_smk) + tt(get(i)), 
                      tt = function(x, t, ...) x * log(t-55),
                      weight=weight_ltrc, data=af.fg)
  nonph.aic4[j] <- 2*length(coef(fit.fg.int)) - 2*fit.fg.int$loglik[2]
}

nonph4 <- data.frame(vars4, aic=nonph.aic4)
nonph4$delta <- full.aic3 - nonph4$aic
nonph4[order(nonph4$aic),]

# lowest AIC  entry_hrx 19283.23, delta=6
full.aic4 <- min(nonph4$aic)


# ----- Round 5:  tt(entry_pchf) + tt(male) + tt(entry_smk) + tt(entry_hrx) add 1 log(time) interaction at a time 

vars5 <- vars4[which(vars4!='entry_hrx')]
n.vars5 <- length(vars5)
nonph.pval5 <- rep(NA, n.vars5)
nonph.aic5 <- rep(NA, n.vars5)
j <- 0

for(i in vars5){
  j <- j+1
  fit.fg.int <- coxph(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                        entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + 
                        tt(entry_pchf) + tt(male) + tt(entry_smk) + tt(entry_hrx) + tt(get(i)), 
                      tt = function(x, t, ...) x * log(t-55),
                      weight=weight_ltrc, data=af.fg)
  nonph.aic5[j] <- 2*length(coef(fit.fg.int)) - 2*fit.fg.int$loglik[2]
}

nonph5 <- data.frame(vars5, aic=nonph.aic5)
nonph5$delta <- full.aic4 - nonph5$aic
nonph5[order(nonph5$aic),]


# AIC increased, proceed with Round 4 model!





# ----- Check interactions with sex ----- 

int1 <- c('entry_hgt', 'entry_wgt', 'entry_sbp', 'entry_dbp', 'entry_hrx', 
          'entry_smk', 'entry_alc_elev', 'entry_diab', 'entry_pchf', 'entry_pmi')
n.vars.int1 <- length(int1)
int.pval1 <- rep(NA, n.vars.int1)
int.aic1 <- rep(NA, n.vars.int1)
j <- 0

for(i in int1){
  j <- j+1
  fit.fg.int <- coxph(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                        entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:get(i) +
                        tt(entry_pchf) + tt(male) + tt(entry_smk) + tt(entry_hrx), 
                      tt = function(x, t, ...) x * log(t-55),
                      weight=weight_ltrc, data=af.fg)
  int.aic1[j] <- AIC(fit.fg.int)
}

intres1 <- data.frame(int1, aic=int.aic1)
intres1$delta <- full.aic4 - intres1$aic
intres1[order(intres1$aic),]

# Diabetes has delta=3; do not add any sex interactions


