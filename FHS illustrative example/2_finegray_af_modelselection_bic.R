### Fine-Gray models for FHS - Lifetime risk of AF
### Mstate package to prepare data, no errors in R v. 4.0.2
### Include interactions with log(time), determine using bic
### Updated analyses for revision 7/15/2021


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

# bic
BIC(fit.fg)
full.bic0 <- BIC(fit.fg)


# ----- Round 1: add 1 log(time) interaction at a time ----- 

n.vars1 <- length(vars)
nonph.pval1 <- rep(NA, n.vars1)
nonph.bic1 <- rep(NA, n.vars1)
j <- 0

for(i in vars){
  j <- j+1
  fit.fg.int <- coxph(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                        entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + tt(get(i)), 
                      tt = function(x, t, ...) x * log(t-55),
                      weight=weight_ltrc, data=af.fg)
  nonph.pval1[j] <- summary(fit.fg.int)$coef[n.vars1+1,5]
  nonph.bic1[j] <- BIC(fit.fg.int)
}

nonph1 <- data.frame(vars, pval=round(nonph.pval1,6), bic=nonph.bic1)
nonph1$delta <- full.bic0 - nonph1$bic
nonph1[order(nonph1$bic),]
full.bic1 <- min(nonph1$bic) # Weight delta 25.6


# ----- Round 2: tt(entry_wgt) + add 1 log(time) interaction at a time 

vars2 <- vars[which(vars!='entry_wgt')]
n.vars2 <- length(vars2)
nonph.pval2 <- rep(NA, n.vars2)
nonph.bic2 <- rep(NA, n.vars2)
j <- 0

for(i in vars2){
  j <- j+1
  fit.fg.int <- coxph(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                        entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + tt(entry_wgt) + tt(get(i)), 
                      tt = function(x, t, ...) x * log(t-55),
                      weight=weight_ltrc, data=af.fg)
  nonph.pval2[j] <- summary(fit.fg.int)$coef[n.vars2+1,5]
  nonph.bic2[j] <- BIC(fit.fg.int)
}

nonph2 <- data.frame(vars2, pval=round(nonph.pval2,6), bic=nonph.bic2)
nonph2$delta <- full.bic1 - nonph2$bic
nonph2[order(nonph2$bic),]
full.bic2 <- min(nonph2$bic) # Delta 15.2 for MI



# ----- Round 3: tt(entry_wgt) + tt(entry_pmi) + add 1 log(time) interaction at a time 

vars3 <- vars2[which(vars2!='entry_pmi')]
n.vars3 <- length(vars3)
nonph.pval3 <- rep(NA, n.vars3)
nonph.bic3 <- rep(NA, n.vars3)
j <- 0

for(i in vars3){
  j <- j+1
  fit.fg.int <- coxph(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                        entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + tt(entry_wgt) + tt(entry_pmi) + tt(get(i)), 
                      tt = function(x, t, ...) x * log(t-55),
                      weight=weight_ltrc, data=af.fg)
  nonph.pval3[j] <- summary(fit.fg.int)$coef[n.vars3+1,5]
  nonph.bic3[j] <- BIC(fit.fg.int)
}

nonph3 <- data.frame(vars3, pval=round(nonph.pval3,6), bic=nonph.bic3)
nonph3$delta <- full.bic2 - nonph3$bic
nonph3[order(nonph3$bic),]
full.bic3 <- min(nonph3$bic) # Delta 8.9 for smoking



# ----- Round 4: tt(entry_wgt) + tt(entry_pmi) + tt(entry_smk) + add 1 log(time) interaction at a time 

vars4 <- vars3[which(vars3!='entry_smk')]
n.vars4 <- length(vars4)
nonph.pval4 <- rep(NA, n.vars4)
nonph.bic4 <- rep(NA, n.vars4)
j <- 0

for(i in vars4){
  j <- j+1
  fit.fg.int <- coxph(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                        entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + 
                        tt(entry_wgt) + tt(entry_pmi) + tt(entry_smk) + tt(get(i)), 
                      tt = function(x, t, ...) x * log(t-55),
                      weight=weight_ltrc, data=af.fg)
  nonph.pval4[j] <- summary(fit.fg.int)$coef[n.vars4+1,5]
  nonph.bic4[j] <- BIC(fit.fg.int)
}

nonph4 <- data.frame(vars4, pval=round(nonph.pval4,6), bic=nonph.bic4)
nonph4$delta <- full.bic3 - nonph4$bic
nonph4[order(nonph4$bic),]


# Note: stop here  bic, increased




# ----- Check interactions with sex ----- 




# ----- Round 1: tt(entry_wgt) + tt(entry_pmi) + tt(entry_smk) + add 1 covariate-sex interaction at a time 

int1 <- c('entry_hgt', 'entry_wgt', 'entry_sbp', 'entry_dbp', 'entry_hrx', 
          'entry_smk', 'entry_alc_elev', 'entry_diab', 'entry_pchf', 'entry_pmi')
n.vars.int1 <- length(int1)
int.pval1 <- rep(NA, n.vars.int1)
int.bic1 <- rep(NA, n.vars.int1)
j <- 0

for(i in int1){
  j <- j+1
  fit.fg.int <- coxph(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                        entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + 
                        tt(entry_wgt) + tt(entry_pmi) + tt(entry_smk) + male:get(i), 
                      tt = function(x, t, ...) x * log(t-55),
                      weight=weight_ltrc, data=af.fg)
  int.bic1[j] <- BIC(fit.fg.int)
}

intres1 <- data.frame(int1, bic=int.bic1)
intres1$delta <- full.bic3 - intres1$bic
intres1[order(intres1$bic),] # BIC increased, no interactions 


