### Flexible parametric models for FHS - Lifetime risk of AF
### Choose # of knots and include interactions with log(time), determine using AIC
### Updated analyses for revision 7/16/2021

library(mstate)
library(rstpm2) 
options(scipen=999)

# ----- Data prep -----

af <- read.csv('af_pseudo.csv')
vars <- c('male', 'entry_hgt', 'entry_wgt', 'entry_sbp', 'entry_dbp', 'entry_hrx', 
          'entry_smk', 'entry_alc_elev', 'entry_diab', 'entry_pchf', 'entry_pmi')
means <- apply(af[, vars], 2, mean)
sds <- apply(af[, vars], 2, sd)

start_time <- Sys.time()

# Weighted dataset
af.fg <- crprep(Tstop="exitage", status="event", data=af, trans=1, cens=0, Tstart="entryage_55", id="fid", keep=vars, shorten=FALSE)
af.fg$weight_ltrc <- af.fg$weight.cens*af.fg$weight.trunc
af.fg$ev1 <- ifelse(af.fg$status==1,1,0)

end_time <- Sys.time()
end_time - start_time

# ----- Final model ----- 

start_time2 <- Sys.time()
fp.final <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                    entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_sbp + male:entry_pchf, 
                  weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.s <- summary(fp.final)
fp.s

# Prep for export
fp.res <- cbind(coef(fp.s)[,1], coef(fp.s)[,1]-1.96*coef(fp.s)[,2], coef(fp.s)[,1]+1.96*coef(fp.s)[,2]) 

# convert to 1 SDU increase
fp.res <- exp(fp.res)
fp.res[c(3,21),] <- fp.res[c(3,21),] ^ sds[2]
fp.res[c(4,19),] <- fp.res[c(4,19),] ^ sds[3]
fp.res[c(5,17),] <- fp.res[c(5,17),] ^ sds[4]
fp.res[6,] <- fp.res[6,] ^ sds[5]

fp.fmt <- cbind(rownames(fp.res), paste0(sprintf('%.2f', fp.res[,1]), " (", sprintf('%.2f', fp.res[,2]), ", ", sprintf('%.2f', fp.res[,3]), ")"))
write.csv(fp.fmt[c(2:12,17:21),], 'results//flex_af_mod_revision.csv', row.names=FALSE)


#  ----- Predictions for each variable ----- 

# Get time intervals
min.death.id <- af$fid[af$exitage == min(af$exitage[af$event==2])]
min.death.id <- min.death.id[1]
intervals <- af.fg[af.fg$fid==min.death.id, c('Tstart', 'Tstop', 'ev1')]
intervals$ev1 <- 1

# Add LSMeans type variables
means.long <- do.call("rbind", replicate(nrow(intervals), means, simplify = FALSE))
preddat <- cbind(means.long, intervals)

# Binary variables - 1 - 0
b.vars <- c('male',  'entry_hrx', 'entry_smk', 'entry_alc_elev', 'entry_diab', 'entry_pchf', 'entry_pmi')
b.diffs <- matrix(NA, 7, 3)
j <- 0

for(i in b.vars){
  j <- j+1
  preddat0 <- preddat
  preddat0[,j] <- 0
  
  # Default prediction is increment of 1, which works for binary variables
  b.diffs[j,] <- as.matrix((-1)*tail(predict(fp.final, newdata=preddat0, var=i, type='sdiff', se.fit=TRUE), n=1))
}


# Continuous variables - difference per 1 SDU increase from mean
c.vars <- c('entry_hgt', 'entry_wgt', 'entry_sbp', 'entry_dbp')
c.diffs <- matrix(NA, 4, 3)
j <- 0

for(i in c.vars){
  j <- j+1
  
  # Specify we want an SD increase with exposed argument
  c.diffs[j,] <-  as.matrix((-1)*tail(predict(fp.final, newdata=preddat, var=i, exposed=incrVar(i, increment=sds[i]), type='sdiff', se.fit=TRUE), n=1))
}

end_time2 <- Sys.time()
end_time2 - start_time2

diffs <- 100*rbind(b.diffs, c.diffs)
diffs.fmt <- cbind(c(b.vars, c.vars), paste0(sprintf('%.2f', diffs[,1]), " (", sprintf('%.2f', diffs[,3]), ", ", sprintf('%.2f', diffs[,2]), ")"))
colnames(diffs.fmt) <- c('var', 'ltrdiff')
diffs.fmt <- diffs.fmt[c(1,8:11,2:7),]
write.csv(diffs.fmt, "results//flex_af_diffltr_revision.csv", row.names=FALSE)

