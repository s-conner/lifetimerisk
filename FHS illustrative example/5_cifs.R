# Plot predicted CIFs

library(mstate)
library(geepack)
library(rstpm2) 
options(scipen=999)

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
af.ltrc <- crprep(Tstop="exitage", status="event", data=af, trans=1, cens=0, Tstart="entryage_55", id="fid", keep=vars, shorten=FALSE)
af.ltrc$weight_ltrc <- af.ltrc$weight.cens*af.ltrc$weight.trunc
af.ltrc$ev1 <- ifelse(af.ltrc$status==1,1,0)

# Long dataset for Fine Gray
evtimes <- sort(unique(af.ltrc$Tstop[af.ltrc$ev1==1]))
af.long <- survSplit(Surv(Tstart, Tstop, ev1) ~., af.ltrc, cut=evtimes, episode="timegroup")

af.long$logt <- log(af.long$Tstop-55)
af.long$entry_wgt_z.logt <- af.long$entry_wgt_z * af.long$logt
af.long$entry_pmi.logt <- af.long$entry_pmi * af.long$logt
af.long$entry_smk.logt <- af.long$entry_smk * af.long$logt
af.long$entry_hgt_z.logt <- af.long$entry_hgt_z * af.long$logt
af.long$entry_diab.logt <- af.long$entry_diab * af.long$logt



# ----- Fit models -----


# Pseudo
pseudo <- geese(pseudo.af ~ male + entry_hgt_z + entry_wgt_z + entry_sbp_z + entry_dbp_z + entry_hrx + 
                  entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                data=af, id=fid, jack = TRUE, scale.fix=TRUE, family=gaussian,
                mean.link = "logit", corstr="independence")
pseudo.beta <- unname(pseudo$beta)

# Fine Gray
fg <- coxph(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt_z + entry_wgt_z + entry_sbp_z + entry_dbp_z + entry_hrx + 
              entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + 
              entry_wgt_z.logt + entry_pmi.logt + entry_smk.logt + entry_hgt_z.logt + entry_diab.logt, 
            weight=weight_ltrc, data=af.long)

# Flexible parametric
fp <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt_z + entry_wgt_z + entry_sbp_z + entry_dbp_z + entry_hrx + 
              entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
            weight=weight_ltrc, data=af.ltrc, df=5, tvc=list(entry_wgt_z=1))



# ----- Plot CIFs -----

# Get time intervals
min.death.id <- af$fid[af$exitage == min(af$exitage[af$event==2])]
min.death.id <- min.death.id[1]
intervals <- af.long[af.long$fid==min.death.id, c('Tstart', 'Tstop', 'ev1', 'logt')]
intervals$ev1 <- 1

# Add LSMeans type variables
means[2:5] <- 0 
means.long <- do.call("rbind", replicate(nrow(intervals), means, simplify = FALSE))
preddat <- cbind(means.long, intervals)

# Add time-covariate interactions (note - continuous vars are set to value 0)
preddat$entry_wgt_z.logt <- preddat$entry_wgt_z * preddat$logt
preddat$entry_pmi.logt <- preddat$entry_pmi * preddat$logt
preddat$entry_smk.logt <- preddat$entry_smk * preddat$logt
preddat$entry_hgt_z.logt <- preddat$entry_hgt_z * preddat$logt
preddat$entry_diab.logt <- preddat$entry_diab * preddat$logt


# Predict difference in lifetime risks

pdf("results//af_cifs.pdf", width=10, height=10)
layout(matrix(c(1:12), nrow=4, ncol=3, byrow=TRUE),
       heights=rep(3,4), widths=c(rep(3,3)))
par(mai=c(0,0,0,0), oma=c(6,6,2,2))

varnames <- c('Sex, male vs. female', 'Height, 70.0 vs. 66.2 in.', 'Weight, 213.3 vs. 173.8 lb.', 
              'Systolic blood pressure, 142.8 vs. 125.7 mmHg', 'Diastolic blood pressure, 88.2 vs. 78.4 mmHg', 
              'Antihypertensive use, yes vs. no', 'Current smoker, yes vs. no', 'Elevated alcohol use, yes vs. no', 
              'Prior diabetes, yes vs. no', 'Prior heart failure, yes vs. no', 'Prior myocardial infarction, yes vs. no')
j <- 0

for(i in vars){
  
  j <- j+1
  
  # Recode var to 1 vs. 0
  preddat1 <- preddat0 <- preddat
  preddat0[, colnames(preddat0) == i] <- 0
  preddat1[, colnames(preddat1) == i] <- 1
  preddat1[, colnames(preddat1) == paste0(i, '.logt')] <- preddat1$logt
  
  
  # Predict CIFs from Fine-Gray and Flexible Parametric model 
  fg.cif0 <- survfit(fg, newdata=preddat0, se.fit=FALSE, individual=TRUE)
  fg.cif1 <- survfit(fg, newdata=preddat1, se.fit=FALSE, individual=TRUE)
 
  fp.cif0 <- 1-predict(fp, newdata=preddat0, type='surv', se.fit=FALSE)
  fp.cif1 <- 1-predict(fp, newdata=preddat1, type='surv', se.fit=FALSE)
  
  
  # Predict lifetime risk (pseudo)
  lp0 <- pseudo.beta %*% c(1, as.matrix(preddat0[1,1:11])) # add 1 for intercept
  lp1 <- pseudo.beta %*% c(1, as.matrix(preddat1[1,1:11]))
  pseudo.ltr0 <- 1/(1 + exp(-lp0))
  pseudo.ltr1 <- 1/(1 + exp(-lp1))
  
  
  # Plot CIF
  plot(evtimes, unique(1-fg.cif0$surv), lwd=1, type='s', lty=1, col=1, ylim=c(0,.6), 
       xlab='', ylab='Cumulative incidence of AF', axes = FALSE)
  lines(evtimes, unique(1-fg.cif1$surv), lwd=1,type='s', lty=1, col=2)
  
  lines(evtimes, fp.cif0[1:829], lwd=1, type='s', lty=2, col=1)
  lines(evtimes, fp.cif1[1:829], lwd=1, type='s', lty=2, col=2)
  
  points(max(evtimes), pseudo.ltr0, pch=23, col=1)
  points(max(evtimes), pseudo.ltr1, pch=23, col=2)
  
  # Add labels
  mtext(varnames[j], side=3, adj=0.05, line=-1.5, cex=.75, font=2)
  if(j %in% c(9,10,11)){axis(1L, at=c(55, 65, 75, 85, 95))}
  if(j %in% c(1,4,7,10)){axis(2L, at=c(0, 0.2, 0.4, 0.6))}
  if(j %in% c(3,6,9)){axis(4L, at=c(0, 0.2, 0.4, 0.6))}
  if(j==11){axis(4L, at=c(0, 0.2, 0.4))}
  
  box()
}

mtext("Cumulative incidence of atrial fibrillation", side=2, line=2.5, outer=TRUE, cex=.8, at=.54)
mtext("Age (years)", side=1, line=3, outer=TRUE, cex=.8)

plot.new()
legend(x="center", ncol=1, lwd=1, col=c(2,2,2,1,1,1), lty=c(NA,1,2,NA,1,2), pch=c(23,NA,NA,23,NA,NA),
       legend=c('Pseudo-observation Z=1', 
                'Fine-Gray Z=1', 
                'Flexible parametric Z=1', 
                'Pseudo-observation Z=0', 
                'Fine-Gray Z=0',
                'Flexible parametric Z=0'), cex=.9)


dev.off()

