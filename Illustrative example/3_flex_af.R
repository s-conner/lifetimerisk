### Flexible parametric models for FHS - Lifetime risk of AF
### Choose # of knots and include interactions with log(time), determine using AIC
### Updated analyses 12/14/2020

library(mstate)
library(rstpm2) 

# ----- Data prep -----

af <- read.csv('af_pseudo.csv')

# Convert continuous vars to SDUs to simplify later
af$entry_hgt_z <- as.numeric(scale(af$entry_hgt))
af$entry_wgt_z <- as.numeric(scale(af$entry_wgt))
af$entry_sbp_z <- as.numeric(scale(af$entry_sbp))
af$entry_dbp_z <- as.numeric(scale(af$entry_dbp))

vars <- c('male', 'entry_hgt_z', 'entry_wgt_z', 'entry_sbp_z', 'entry_dbp_z', 'entry_hrx', 
          'entry_smk', 'entry_alc_elev', 'entry_diab', 'entry_pchf', 'entry_pmi')
means <- apply(af[, vars], 2, mean)

# Weighted dataset
af.fg <- crprep(Tstop="exitage", status="event", data=af, trans=1, cens=0, Tstart="entryage_55", id="fid", keep=vars, shorten=FALSE)
af.fg$weight_ltrc <- af.fg$weight.cens*af.fg$weight.trunc
af.fg$ev1 <- ifelse(af.fg$status==1,1,0)


# ----- Final model ----- 

fp.final <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt_z + entry_wgt_z + entry_sbp_z + entry_dbp_z + entry_hrx + 
                    entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                  weight=weight_ltrc, data=af.fg, df=5, tvc=list(entry_wgt_z=1))
fp.s <- summary(fp.final)

# export
#fp.res <- cbind(coef(fp.final), confint(fp.final)) # slow
#fp.res2 <- exp(fp.res)  
#fp.fmt <- cbind(rownames(fp.res2), paste0(sprintf('%.2f', fp.res2[,1]), " (", sprintf('%.2f', fp.res2[,2]), ", ", sprintf('%.2f', fp.res2[,3]), ")"))
#write.csv(fp.fmt[c(2:12,18),], 'results//flex_af_mod.csv', row.names=FALSE)


#  ----- Predictions for each variable ----- 

# Get time intervals
min.death.id <- af$fid[af$exitage == min(af$exitage[af$event==2])]
min.death.id <- min.death.id[1]
intervals <- af.fg[af.fg$fid==min.death.id, c('Tstart', 'Tstop', 'ev1')]
intervals$ev1 <- 1

# Add LSMeans type variables
means.long <- do.call("rbind", replicate(nrow(intervals), means, simplify = FALSE))
preddat <- cbind(means.long, intervals)

# Compare variables - 1 - 0
diffs <- matrix(NA, 11, 3)

# Predict difference in lifetime risks
j <- 0

for(i in vars){
  j <- j+1
  preddat0 <- preddat
  preddat0[,j] <- 0
  diffs[j,] <- as.matrix((-1)*tail(predict(fp.final, newdata=preddat0, var=i, type='sdiff', se.fit=TRUE), n=1))
}

diffs <- diffs*100

cbind(vars, diffs)

diffs.fmt <- cbind(vars, paste0(sprintf('%.2f', diffs[,1]), " (", sprintf('%.2f', diffs[,3]), ", ", sprintf('%.2f', diffs[,2]), ")"))
colnames(diffs.fmt) <- c('var', 'ltrdiff')
#write.csv(diffs.fmt, "results//flex_af_diffltr_se.csv", row.names=FALSE)


#  ----- Predictions for each individual (calibration) ----- 


# Get time intervals
min.death.id <- af$fid[af$exitage == min(af$exitage[af$event==2])]
min.death.id <- min.death.id[1]
intervals <- af.fg[af.fg$fid==min.death.id, c('Tstart', 'Tstop', 'ev1')]
intervals$ev1 <- 1 

n <- nrow(af)
pred.ltr <- rep(NA, nrow(af))


for(i in 1:n){
  
  # Covariate profile for individual i
  cov.i <- do.call("rbind", replicate(nrow(intervals), af[i,vars], simplify=FALSE))
  preddat.i <- cbind(cov.i, intervals)
  
  # Predict lifetime risks
  pred.ltr[i] <- 1-tail(predict(fp.final, newdata=preddat.i, type='surv', se.fit=FALSE), n=1)
}

write.csv(pred.ltr, "results//flex_af_individual_ltr.csv", row.names=FALSE)



