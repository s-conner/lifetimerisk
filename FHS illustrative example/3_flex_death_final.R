### Flexible parametric models for FHS - Lifetime risk of AF
### Choose # of knots and include interactions with log(time), determine using AIC
### Updated analyses 12/14/2020

library(mstate)
library(rstpm2) 

# ----- Data prep -----

af <- read.csv('af_pseudo.csv')

vars <- c('male', 'entry_hgt', 'entry_wgt', 'entry_sbp', 'entry_dbp', 'entry_hrx', 
          'entry_smk', 'entry_alc_elev', 'entry_diab', 'entry_pchf', 'entry_pmi')
means <- apply(af[, vars], 2, mean)

# Create LTRC weighted dataset
af.fg <- crprep(Tstop="exitage", status="event", data=af, trans=2, cens=0, Tstart="entryage_55", id="fid", keep=vars, shorten=FALSE)
af.fg$weight_ltrc <- af.fg$weight.cens*af.fg$weight.trunc
af.fg$ev2 <- ifelse(af.fg$status==2,1,0)


# ----- Final model ----- 

fp.final <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi , 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pchf=1))
fp.s <- summary(fp.final)

# Prep for export
fp.res <- cbind(coef(fp.s)[,1], coef(fp.s)[,1]+1.96*coef(fp.s)[,2], coef(fp.s)[,1]-1.96*coef(fp.s)[,2]) 

# convert to 1 SDU increase
fp.res <- exp(fp.res)
fp.res[3,] <- fp.res[3,] ^ sds[2]
fp.res[4,] <- fp.res[4,] ^ sds[3]
fp.res[5,] <- fp.res[5,] ^ sds[4]
fp.res[6,] <- fp.res[6,] ^ sds[5]

fp.fmt <- cbind(rownames(fp.res), paste0(sprintf('%.2f', fp.res[,1]), " (", sprintf('%.2f', fp.res[,2]), ", ", sprintf('%.2f', fp.res[,3]), ")"))
write.csv(fp.fmt[c(2:12,16:20),], 'results//flex_death_mod_revision.csv', row.names=FALSE)

