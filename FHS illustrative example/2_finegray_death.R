### Fine-Gray models for FHS - Lifetime risk of AF
### Mstate package to prepare data, no errors in R v. 4.0.2
### Include interactions with log(time), determine using AIC
### Updated analyses 1/12/2021

library(survival)
library(mstate)
options(scipen=999)

# ----- Data prep -----

af <- read.csv('af_pseudo.csv')

af$entry_hgt_z <- scale(af$entry_hgt)
af$entry_wgt_z <- scale(af$entry_wgt)
af$entry_sbp_z <- scale(af$entry_sbp)
af$entry_dbp_z <- scale(af$entry_dbp)

vars <- c('male', 'entry_hgt_z', 'entry_wgt_z', 'entry_sbp_z', 'entry_dbp_z', 'entry_hrx', 
          'entry_smk', 'entry_alc_elev', 'entry_diab', 'entry_pchf', 'entry_pmi')
means <- apply(af[, vars], 2, mean)
sds <- apply(af[, vars], 2, sd)

# Weighted dataset
af.ltrc <- crprep(Tstop="exitage", status="event", data=af, trans=2, cens=0, Tstart="entryage_55", id="fid", keep=vars, shorten=FALSE)
af.ltrc$weight_ltrc <- af.ltrc$weight.cens*af.ltrc$weight.trunc
af.ltrc$ev2 <- ifelse(af.ltrc$status==2,1,0)


# ----- Long dataset ----- 

evtimes <- sort(unique(af.ltrc$Tstop[af.ltrc$ev2==1]))
af.long <- survSplit(Surv(Tstart, Tstop, ev2) ~., af.ltrc, cut=evtimes, episode="timegroup")

af.long$logt <- log(af.long$Tstop-55)
af.long$entry_pchf.logt <- af.long$entry_pchf * af.long$logt
af.long$male.logt <- af.long$male * af.long$logt
af.long$entry_smk.logt <- af.long$entry_smk * af.long$logt
af.long$entry_hrx.logt <- af.long$entry_hrx * af.long$logt


# ----- Fit model ----- 

final.long <- coxph(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt_z + entry_wgt_z + entry_sbp_z + entry_dbp_z + entry_hrx + 
                      entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + 
                      entry_pchf.logt + male.logt + entry_smk.logt + entry_hrx.logt,
                    weight=weight_ltrc, data=af.long)
s.final.long <- summary(final.long)
s.hrs <- s.final.long$conf.int[,c(1,3,4)]

fg.fmt <- cbind(rownames(s.hrs), paste0(sprintf('%.2f', s.hrs[,1]), " (", sprintf('%.2f', s.hrs[,2]), ", ", sprintf('%.2f', s.hrs[,3]), ")"))
colnames(fg.fmt) <- c('var', 'est')
write.csv(fg.fmt, 'results//finegray_death_mod.csv', row.names=FALSE)

