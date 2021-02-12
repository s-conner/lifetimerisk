### Flexible parametric models for FHS - Lifetime risk of AF
### Choose # of knots and include interactions with log(time), determine using AIC
### Updated analyses 12/14/2020

library(mstate)
library(rstpm2) 

# ----- Data prep -----

af <- read.csv('af_pseudo.csv')

# Convert continuous vars to SDUs to simplify later
af$entry_hgt_z <- scale(af$entry_hgt)
af$entry_wgt_z <- scale(af$entry_wgt)
af$entry_sbp_z <- scale(af$entry_sbp)
af$entry_dbp_z <- scale(af$entry_dbp)

vars <- c('male', 'entry_hgt_z', 'entry_wgt_z', 'entry_sbp_z', 'entry_dbp_z', 'entry_hrx', 
          'entry_smk', 'entry_alc_elev', 'entry_diab', 'entry_pchf', 'entry_pmi')
means <- apply(af[, vars], 2, mean)

# Create LTRC weighted dataset
af.fg <- crprep(Tstop="exitage", status="event", data=af, trans=2, cens=0, Tstart="entryage_55", id="fid", keep=vars, shorten=FALSE)
af.fg$weight_ltrc <- af.fg$weight.cens*af.fg$weight.trunc
af.fg$ev2 <- ifelse(af.fg$status==2,1,0)


# ----- Final model ----- 

fp.final <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt_z + entry_wgt_z + entry_sbp_z + entry_dbp_z + entry_hrx + 
                    entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                  weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pchf=1, entry_wgt_z=1))
fp.s <- summary(fp.final)

fp.final2 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt_z + entry_wgt_z + entry_sbp_z + entry_dbp_z + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   link='PH', weight=weight_ltrc, data=af.fg, 
                   smooth.formula=~nsx(log(Tstop),df=3) + male:nsx(log(Tstop-53),df=1))
summary(fp.final2)


fp.final2 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt_z + entry_wgt_z + entry_sbp_z + entry_dbp_z + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   link='PH', weight=weight_ltrc, data=af.fg, 
                   smooth.formula=~nsx(log(Tstop-55),df=3) + male:nsx(log(Tstop-55),df=1) + entry_hrx:nsx(log(Tstop-55),df=1) + entry_diab:nsx(log(Tstop-55),df=1) +
                     entry_smk:nsx(log(Tstop-55),df=1) + entry_pchf:nsx(log(Tstop-55),df=1) + entry_wgt_z:nsx(log(Tstop-55),df=1))



fp.final2 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt_z + entry_wgt_z + entry_sbp_z + entry_dbp_z + entry_hrx + 
                    entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                  weight=weight_ltrc, data=af.fg, df=3, 
                  smooth.formula=~male:log(Tstop-55) + entry_hrx:log(Tstop-55) + entry_diab:log(Tstop-55) + 
                    entry_pchf:log(Tstop-55) + entry_wgt_z:log(Tstop-55))
fp.s2 <- summary(fp.final2)




fp.res <- cbind(coef(fp.final), confint(fp.final)) # slow
fp.res2 <- exp(fp.res)  

fp.fmt <- cbind(rownames(fp.res2), paste0(sprintf('%.2f', fp.res2[,1]), " (", sprintf('%.2f', fp.res2[,2]), ", ", sprintf('%.2f', fp.res2[,3]), ")"))
write.csv(fp.fmt[c(2:12, 16:21),], 'results//flex_death_mod.csv', row.names=FALSE)

