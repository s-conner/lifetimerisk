### Fine-Gray models for FHS - Lifetime risk of AF
### Mstate package to prepare data, no errors in R v. 4.0.2
### Include interactions with log(time), determine using AIC
### Updated analyses for revision 7/16/2021

library(survival)
library(mstate)
options(scipen=999)

# ----- Data prep -----

af <- read.csv('af_pseudo.csv')
vars <- c('male', 'entry_hgt', 'entry_wgt', 'entry_sbp', 'entry_dbp', 'entry_hrx', 
          'entry_smk', 'entry_alc_elev', 'entry_diab', 'entry_pchf', 'entry_pmi')
means <- apply(af[, vars], 2, mean)
sds <- apply(af[, vars], 2, sd)

start_time <- Sys.time()

# Create LTRC weighted dataset
af.fg <- crprep(Tstop="exitage", status="event", data=af, trans=1, cens=0, Tstart="entryage_55", id="fid", keep=vars, shorten=FALSE)
af.fg$weight_ltrc <- af.fg$weight.cens*af.fg$weight.trunc
af.fg$ev1 <- ifelse(af.fg$status==1,1,0)


# ----- Final model ----- 

final <- coxph(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx +
                  entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi +
                  tt(entry_wgt) + tt(entry_pmi) + tt(entry_smk),
                tt = function(x, t, ...) x * log(t-55),
                weight=weight_ltrc, data=af.fg)
s.final <- summary(final)

s.hrs <- s.final$conf.int[,c(1,3,4)]
s.hrs[c(2:5),] <- s.hrs[c(2:5),]^sds[2:5]
s.hrs[12,] <- s.hrs[12,]^sds[3]

fg.fmt <- cbind(rownames(s.hrs), paste0(sprintf('%.2f', s.hrs[,1]), " (", sprintf('%.2f', s.hrs[,2]), ", ", sprintf('%.2f', s.hrs[,3]), ")"))
colnames(fg.fmt) <- c('var', 'est')
write.csv(fg.fmt, 'results//finegray_af_mod_revision.csv', row.names=FALSE)


# ----- Selected model in long dataset form ------

# While anyone who died already has long dataset form, need to additionally split those with AF or censored
evtimes <- sort(unique(af.fg$Tstop[af.fg$ev1==1]))
af.long <- survSplit(Surv(Tstart, Tstop, ev1) ~., af.fg, cut=evtimes, episode="timegroup")

af.long$logt <- log(af.long$Tstop-55)
af.long$entry_wgt.logt <- af.long$entry_wgt * af.long$logt
af.long$entry_pmi.logt <- af.long$entry_pmi * af.long$logt
af.long$entry_smk.logt <- af.long$entry_smk * af.long$logt

end_time <- Sys.time()
end_time - start_time

start_time2 <- Sys.time()
final.long <- coxph(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                      entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + 
                      entry_wgt.logt + entry_pmi.logt + entry_smk.logt,
                    weight=weight_ltrc, data=af.long)


#  ----- Predictions for each variable -----

# We need all intervals through tau to make a prediction at tau
# This could be anyone who died or survived event-free until tau, but who did not enter late
# If all individuals enter delayed, this could pose a complication - then grab all the intervals
# Maybe there is an easier way

min.death.id <- af$fid[af$exitage == min(af$exitage[af$event==2])]
intervals <- af.long[af.long$fid==min.death.id, c('Tstart', 'Tstop', 'ev1', 'logt')]
intervals$ev1 <- 1

# Derive LSMeans type variables
means.long <- do.call("rbind", replicate(nrow(intervals), means, simplify = FALSE))
preddat <- cbind(means.long, intervals)

preddat$entry_wgt.logt <- preddat$entry_wgt * preddat$logt
preddat$entry_pmi.logt <- preddat$entry_pmi * preddat$logt
preddat$entry_smk.logt <- preddat$entry_smk * preddat$logt


# Binary variables - 1 - 0
b.vars <- c('male',  'entry_hrx', 'entry_smk', 'entry_alc_elev', 'entry_diab', 'entry_pchf', 'entry_pmi')
b.diffs <- matrix(NA, 7, 4)
j <- 0

for(i in b.vars){
  j <- j+1

  # Exclude var of interest
  tmp <- preddat[, !(colnames(preddat) %in% c(i, paste0(i, '.logt')))]

  # Add contrast variable
  preddat1 <- cbind(tmp, rep(1, nrow(tmp)), tmp$logt)
  colnames(preddat1)[c(ncol(preddat1)-1,ncol(preddat1))] <- c(i, paste0(i, '.logt'))
  preddat0 <- cbind(tmp, rep(0, nrow(tmp)), rep(0, nrow(tmp)))
  colnames(preddat0)[c(ncol(preddat0)-1,ncol(preddat0))] <- c(i, paste0(i, '.logt'))

  # Predict lifetime risks
  pred1 <- survfit(final.long, newdata=preddat1, se.fit=FALSE, individual=TRUE)
  pred0 <- survfit(final.long, newdata=preddat0, se.fit=FALSE, individual=TRUE)

  ltr1 <- 1-tail(pred1$surv, n=1)
  ltr0 <- 1-tail(pred0$surv, n=1)
  b.diffs[j,] <- c(i, ltr1, ltr0, ltr1-ltr0)
}


# Continuous variables - difference per 1 SDU increase from mean
c.vars <- c('entry_hgt', 'entry_wgt', 'entry_sbp', 'entry_dbp' )
c.diffs <- matrix(NA, 4, 4)
j <- 0

for(i in c.vars){
  j <- j+1

  # Exclude var of interest
  tmp <- preddat[, !(colnames(preddat) %in% c(i, paste0(i, '.logt')))]

  # Derive contrast (mean vs 1 SD above mean)
  c.avg <- mean(af[,i])
  c.avg1sd <- mean(af[,i]) + sd(af[,i])

  # Add contrast variable
  preddat1 <- cbind(tmp, rep(c.avg1sd, nrow(tmp)), c.avg1sd*tmp$logt)
  colnames(preddat1)[c(ncol(preddat1)-1,ncol(preddat1))] <- c(i, paste0(i, '.logt'))

  preddat0 <- cbind(tmp, rep(c.avg, nrow(tmp)), c.avg*tmp$logt)
  colnames(preddat0)[c(ncol(preddat0)-1,ncol(preddat0))] <- c(i, paste0(i, '.logt'))

  # Predict lifetime risks
  pred1 <- survfit(final.long, newdata=preddat1, se.fit=FALSE, individual=TRUE)
  pred0 <- survfit(final.long, newdata=preddat0, se.fit=FALSE, individual=TRUE)

  ltr1 <- 1-tail(pred1$surv, n=1)
  ltr0 <- 1-tail(pred0$surv, n=1)
  c.diffs[j,] <- c(i, ltr1, ltr0, ltr1-ltr0)
}

end_time2 <- Sys.time()
end_time2 - start_time2

diffs <- rbind(b.diffs, c.diffs)
colnames(diffs) <- c('var', 'ltr1', 'ltr0', 'ltrdiff')
write.csv(diffs, "results//finegray_af_diffltr_revision2.csv", row.names=FALSE)



