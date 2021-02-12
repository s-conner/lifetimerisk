# Simulation study to assess Fine-Gray PH Model with log(time) interactions using Geskus (2011) method
# CSH Data
# Sarah Conner - Updated 1/20/2021

library(survival)
library(mstate)
options(scipen = 999)

a <- as.numeric(Sys.getenv("SGE_TASK_ID")) #Identify the job task ID from the t variables in the shell, see below
sim <- read.csv("//restricted//projectnb//conner-thesis//lifetimerisk//simulation//fg_splitjobs.csv", header=TRUE, sep=",")
sim <- sim[a, ]

l <- sim$simid
nstart <- sim$start
nend <- sim$end

ltr1 <- rep(NA, 500)
ltr0 <- rep(NA, 500)
ltr.diff <- rep(NA, 500)

# For each simulated dataset, fit the model, predict LTRs
j <- 0

for(k in nstart:nend){
  
  j <- j + 1

  
  # ----- Read in data ----- 
  dat <- read.csv(paste0("//restricted//projectnb//conner-thesis//lifetimerisk//simulation//data_csh//data", l, "_", k, ".csv"), header=TRUE, sep=",")

  dat$cause[dat$times>40] <- 0
  dat$times[dat$times>40] <- 40
  
  dat$entry <- 55+round(dat$entry, 4)
  dat$times <- 55+round(dat$times, 4)
  
  dat <- dat[dat$entry!=dat$times,] # In case the rounding means entry and exit times are close together
  
  
  # ----- LTRC weighted dataset in long form to accommodate log(time) interaction  ----- 
  
  # LTRC weights (Geskus 2011)
  dat.ltrc <- crprep(Tstop='times', status='cause', data=dat, trans=1, cens=0, Tstart='entry', id='id', keep='arm', shorten=FALSE)
  dat.ltrc$weight.ltrc <- dat.ltrc$weight.cens*dat.ltrc$weight.trunc
  dat.ltrc$ev1 <- ifelse(dat.ltrc$status==1,1,0)
  
  # long dataset
  evtimes <- sort(unique(dat$times[dat$cause==1]))
  dat.long <- survSplit(Surv(Tstart, Tstop, ev1) ~., dat.ltrc, cut=evtimes, episode ="timegroup")
  dat.long$logt <- log(dat.long$Tstop)
  dat.long$arm.logt <- dat.long$arm * dat.long$logt
  
  # fit model
  nonph <- coxph(Surv(Tstart, Tstop, ev1) ~ arm + arm.logt, weight=weight.ltrc, data=dat.long, control=coxph.control(timefix = FALSE))

  
  # ----- Predict CIF -----
  
  min.death.id <- dat$id[dat$times == min(dat$times[dat$cause==2])]
  min.death.id <- min.death.id[1]
  intervals <- dat.long[dat.long$id==min.death.id, c('Tstart', 'Tstop', 'ev1', 'logt')]
  intervals$ev1 <- 1
  
  arm1 <- data.frame(arm=1, intervals)
  arm1$arm.logt <- arm1$logt
  arm0 <- data.frame(arm=0, intervals)
  arm0$arm.logt <- 0
  
  pred.arm1 <- survfit(nonph, newdata=arm1, se.fit=FALSE, individual=TRUE)
  pred.arm0 <- survfit(nonph, newdata=arm0, se.fit=FALSE, individual=TRUE)
  
  # Lifetime risks
  ltr1[j] <- 1-tail(pred.arm1$surv, n=1)
  ltr0[j] <- 1-tail(pred.arm0$surv, n=1)
  ltr.diff[j] <- ltr1[j] - ltr0[j]

}  


# ----- Output results ----- 
res <- cbind(ltr1, ltr0, ltr.diff)
colnames(res) <- c('ltr1', 'ltr0', 'ltr.diff')
write.csv(res, paste0("//restricted//projectnb//conner-thesis//lifetimerisk//simulation//finegray_logtime//cshv2//res_", l, "_", nend, ".csv"), row.names=FALSE)





