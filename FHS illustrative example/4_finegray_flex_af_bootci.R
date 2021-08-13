### Fine-Gray models for FHS - Lifetime risk of AF
### Bootstrap CI (use SCC to reduce computation)
### Updated analyses 1/12/2021

library(survival)
library(mstate)
library(rstpm2)
options(scipen=999)

a <- as.numeric(Sys.getenv("SGE_TASK_ID")) 
nboot <- 1000
nboot.a <- nboot/20
start <- (a-1)*nboot.a + 1
end <- a*nboot.a


# ----- Data prep -----

af <- read.csv('af_pseudo.csv')
n <- nrow(af)
af$id <- 1:n
vars <- c('male', 'entry_hgt', 'entry_wgt', 'entry_sbp', 'entry_dbp', 'entry_hrx', 
          'entry_smk', 'entry_alc_elev', 'entry_diab', 'entry_pchf', 'entry_pmi')
means <- apply(af[, vars], 2, mean)
sds <- apply(af[, vars], 2, sd)


# ----- Bootstrap differences in LTR for each variable ------

boot.fp.ltr1 <- matrix(NA, nboot, length(vars))
boot.fp.ltr0 <- matrix(NA, nboot, length(vars))
boot.fp.ltr.diff <- matrix(NA, nboot, length(vars))
boot.fg.ltr1 <- matrix(NA, nboot, length(vars))
boot.fg.ltr0 <- matrix(NA, nboot, length(vars))
boot.fg.ltr.diff <- matrix(NA, nboot, length(vars))
b.vars <- c('male',  'entry_hrx', 'entry_smk', 'entry_alc_elev', 'entry_diab', 'entry_pchf', 'entry_pmi')
c.vars <- c('entry_hgt', 'entry_wgt', 'entry_sbp', 'entry_dbp')

for(i in start:end){
  print(i)
  set.seed(i)
  
  # Bootstrap sample
  ids <- sample(af$id, n, replace=TRUE)
  boot.dat <- af[ids,]
  boot.dat$bootid <- 1:n
  
  # Create LTRC-weighted data on bootstrap sample
  boot.dat.fg <- crprep(Tstop="exitage", status="event", data=boot.dat, trans=1, cens=0, Tstart="entryage_55", id="bootid", keep=vars, shorten=FALSE)
  boot.dat.fg$weight_ltrc <- boot.dat.fg$weight.cens * boot.dat.fg$weight.trunc
  boot.dat.fg$ev1 <- ifelse(boot.dat.fg$status==1,1,0)
  
  
  # ----- Flexible parametric ------
  
  if(exists('boot.fp')){rm(boot.fp)} # clear model from previous iteration
  if(exists('s')){rm(s)} # clear model from previous iteration
  
  # Fit flexible parametric model
  tryCatch({boot.fp <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=boot.dat.fg, df=5, tvc=list(entry_wgt=1))
  }, error=function(e) {}, warning=function(w) {})
  
  # Predictions are bad when model doesn't fully converge, don't use
  tryCatch({s <-summary(boot.fp)}, error=function(e) {}, warning=function(w) {})
  
  if(exists('s')){
    # Prepare prediction data
    min.death.id <- boot.dat$bootid[boot.dat$exitage == min(boot.dat$exitage[boot.dat$event==2])]
    min.death.id <- min.death.id[1] # due to bootstrap, this person could appear in data more than once
    intervals <- boot.dat.fg[boot.dat.fg$bootid==min.death.id, c('Tstart', 'Tstop', 'ev1')]
    intervals$ev1 <- 1 
    
    # Add LSMeans type variables
    means.long <- do.call("rbind", replicate(nrow(intervals), means, simplify = FALSE))
    preddat <- cbind(means.long, intervals)
    
    # Predict for binary variables - 1 - 0
    j <- 0
    
    for(k in b.vars){
      j <- j+1
      
      # Exclude var of interest
      tmp <- preddat[, !(colnames(preddat)==k)]
      
      # Add contrast variable
      preddat1 <- cbind(tmp, rep(1, nrow(tmp)))
      colnames(preddat1)[ncol(preddat1)] <- k
      preddat0 <- cbind(tmp, rep(0, nrow(tmp)))
      colnames(preddat0)[ncol(preddat0)] <- k
      
      # Predict lifetime risks
      pred1 <- predict(boot.fp, newdata=preddat1, type='fail', se.fit=FALSE)
      pred0 <- predict(boot.fp, newdata=preddat0, type='fail', se.fit=FALSE)
      
      boot.fp.ltr1[i,j] <- tail(pred1, n=1)
      boot.fp.ltr0[i,j] <- tail(pred0, n=1)
      boot.fp.ltr.diff[i,j] <- boot.fp.ltr1[i,j] - boot.fp.ltr0[i,j]
    }
    
    # Predict for continuous variables - difference per 1 SDU increase from mean
    
    for(k in c.vars){
      j <- j+1
      
      # Exclude var of interest
      tmp <- preddat[, !(colnames(preddat)==k)]
      
      # Derive contrast (mean vs 1 SD above mean)
      c.avg <- mean(af[,k])
      c.avg1sd <- mean(af[,k]) + sd(af[,k])
      
      # Add contrast variable
      preddat1 <- cbind(tmp, rep(c.avg1sd, nrow(tmp)))
      colnames(preddat1)[ncol(preddat1)] <- k
      preddat0 <- cbind(tmp, rep(c.avg, nrow(tmp)))
      colnames(preddat0)[ncol(preddat0)] <- k
      
      # Predict lifetime risks
      pred1 <- predict(boot.fp, newdata=preddat1, type='fail', se.fit=FALSE)
      pred0 <- predict(boot.fp, newdata=preddat0, type='fail', se.fit=FALSE)
      
      boot.fp.ltr1[i,j] <- tail(pred1, n=1)
      boot.fp.ltr0[i,j] <- tail(pred0, n=1)
      boot.fp.ltr.diff[i,j] <- boot.fp.ltr1[i,j] - boot.fp.ltr0[i,j]
    }
  }

    
  # ----- Fine-Gray ------
  
  if(exists('boot.fg')){rm(boot.fg)} # clear model from previous iteration
  
  # Further split LTRC data at event times
  boot.evtimes <- sort(unique(boot.dat$exitage[boot.dat$event==1]))
  boot.dat.long <- survSplit(Surv(Tstart, Tstop, ev1) ~., boot.dat.fg, cut=boot.evtimes, episode ="timegroup")
  
  # Create covariate-log(time) interaction variables
  boot.dat.long$logt <- log(boot.dat.long$Tstop-55)
  boot.dat.long$entry_wgt.logt <- boot.dat.long$entry_wgt * boot.dat.long$logt
  boot.dat.long$entry_pmi.logt <- boot.dat.long$entry_pmi * boot.dat.long$logt
  boot.dat.long$entry_smk.logt <- boot.dat.long$entry_smk * boot.dat.long$logt
  boot.dat.long$entry_hgt.logt <- boot.dat.long$entry_hgt * boot.dat.long$logt
  boot.dat.long$entry_diab.logt <- boot.dat.long$entry_diab * boot.dat.long$logt
  
  # Fit fine-gray model with long data
  tryCatch({boot.fg <- coxph(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                        entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + 
                        entry_wgt.logt + entry_pmi.logt + entry_smk.logt + entry_hgt.logt + entry_diab.logt,
                      weight=weight_ltrc, data=boot.dat.long, control=coxph.control(timefix = FALSE))
  }, error = function(e) {})
  
  summary(boot.fg)
  
  if(exists('boot.fg')){
    # Prepare prediction data
    min.death.id <- boot.dat$bootid[boot.dat$exitage == min(boot.dat$exitage[boot.dat$event==2])]
    intervals <- boot.dat.long[boot.dat.long$bootid==min.death.id, c('Tstart', 'Tstop', 'ev1', 'logt')]
    intervals$ev1 <- 1 
    
    # Derive LSMeans type variables
    means.long <- do.call("rbind", replicate(nrow(intervals), means, simplify = FALSE))
    preddat <- cbind(means.long, intervals)
    
    preddat$entry_wgt.logt <- preddat$entry_wgt * preddat$logt
    preddat$entry_pmi.logt <- preddat$entry_pmi * preddat$logt
    preddat$entry_smk.logt <- preddat$entry_smk * preddat$logt
    preddat$entry_hgt.logt <- preddat$entry_hgt * preddat$logt
    preddat$entry_diab.logt <- preddat$entry_diab * preddat$logt
    
    # Predict for binary variables - 1 - 0
    l <- 0
    
    for(m in b.vars){
      l <- l+1
      
      # Exclude var of interest
      tmp <- preddat[, !(colnames(preddat) %in% c(m, paste0(m, '.logt')))]
      
      # Add contrast variable
      preddat1 <- cbind(tmp, rep(1, nrow(tmp)), tmp$logt)
      colnames(preddat1)[c(ncol(preddat1)-1,ncol(preddat1))] <- c(m, paste0(m, '.logt'))
      preddat0 <- cbind(tmp, rep(0, nrow(tmp)), rep(0, nrow(tmp)))
      colnames(preddat0)[c(ncol(preddat0)-1,ncol(preddat0))] <- c(m, paste0(m, '.logt'))
      
      # Predict lifetime risks
      pred1 <- survfit(boot.fg, newdata=preddat1, se.fit=FALSE, individual=TRUE)
      pred0 <- survfit(boot.fg, newdata=preddat0, se.fit=FALSE, individual=TRUE)
      
      boot.fg.ltr1[i,l] <- 1-tail(pred1$surv, n=1)
      boot.fg.ltr0[i,l] <- 1-tail(pred0$surv, n=1)
      boot.fg.ltr.diff[i,l] <- boot.fg.ltr1[i,l] - boot.fg.ltr0[i,l]
    }
    
    # Predict for continuous variables - difference per 1 SDU increase from mean
    
    for(m in c.vars){
      l <- l+1
      
      # Exclude var of interest
      tmp <- preddat[, !(colnames(preddat) %in% c(m, paste0(m, '.logt')))]
      
      # Derive contrast (mean vs 1 SD above mean)
      c.avg <- mean(af[,m])
      c.avg1sd <- mean(af[,m]) + sd(af[,m])
      
      # Add contrast variable
      preddat1 <- cbind(tmp, rep(c.avg1sd, nrow(tmp)), c.avg1sd*tmp$logt)
      colnames(preddat1)[c(ncol(preddat1)-1,ncol(preddat1))] <- c(m, paste0(m, '.logt'))
      preddat0 <- cbind(tmp, rep(c.avg, nrow(tmp)), c.avg*tmp$logt)
      colnames(preddat0)[c(ncol(preddat0)-1,ncol(preddat0))] <- c(m, paste0(m, '.logt'))
      
      # Predict lifetime risks
      pred1 <- survfit(boot.fg, newdata=preddat1, se.fit=FALSE, individual=TRUE)
      pred0 <- survfit(boot.fg, newdata=preddat0, se.fit=FALSE, individual=TRUE)
      
      boot.fg.ltr1[i,l] <- 1-tail(pred1$surv, n=1)
      boot.fg.ltr0[i,l] <- 1-tail(pred0$surv, n=1)
      boot.fg.ltr.diff[i,l] <- boot.fg.ltr1[i,l] - boot.fg.ltr0[i,l]
    }
  }
}


# ----- Export ------

write.csv(boot.fp.ltr1, paste0("bootci_af//fp_ltr1_ci_", a, ".csv"), row.names=FALSE)
write.csv(boot.fp.ltr0, paste0("bootci_af//fp_ltr0_ci_", a, ".csv"), row.names=FALSE)
write.csv(boot.fp.ltr.diff, paste0("bootci_af//fp_ltrdiff_ci_", a, ".csv"), row.names=FALSE)

write.csv(boot.fg.ltr1, paste0("bootci_af//fg_ltr1_ci_", a, ".csv"), row.names=FALSE)
write.csv(boot.fg.ltr0, paste0("bootci_af//fg_ltr0_ci_", a, ".csv"), row.names=FALSE)
write.csv(boot.fg.ltr.diff, paste0("bootci_af//fg_ltrdiff_ci_", a, ".csv"), row.names=FALSE)


