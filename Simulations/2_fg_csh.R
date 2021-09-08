# Simulation study to assess Fine-Gray PH Model using Geskus (2011) method
# CSH Data
# Sarah Conner - Updated 1/20/2021

library(survival)
library(mstate)
options(scipen = 999)

a <- as.numeric(Sys.getenv("SGE_TASK_ID")) #Identify the job task ID from the t variables in the shell, see below

ltr1 <- rep(NA, 2000)
ltr0 <- rep(NA, 2000)
ltr.diff <- rep(NA, 2000)

# For each simulated dataset, fit the model, predict LTRs

for(j in 1:2000){
  

  # ----- Read in data ----- 
  dat <- read.csv(paste0("//restricted//projectnb//conner-thesis//lifetimerisk//simulation//data_csh//data", a, "_", j, ".csv"), header=TRUE, sep=",")

  dat$cause[dat$times>40] <- 0
  dat$times[dat$times>40] <- 40
  
  # Shift data to be similar to FHS; lifetime risk at age 95
  dat$entry <- 55+round(dat$entry, 4)
  dat$times <- 55+round(dat$times, 4)
  
  dat <- dat[dat$entry!=dat$times,] # In case the rounding means entry and exit times are close together
  
  
  # ----- PH FG with left truncation ----- 
  dat.ltrc <- crprep(Tstop='times', status='cause', data=dat, trans=1, cens=0, Tstart='entry', id='id', keep='arm', shorten=FALSE)
  dat.ltrc$weight.ltrc <- dat.ltrc$weight.cens*dat.ltrc$weight.trunc
  dat.ltrc$ev1 <- ifelse(dat.ltrc$status==1,1,0)
  
  ph <- coxph(Surv(Tstart, Tstop, ev1) ~ arm, weight=weight.ltrc, data=dat.ltrc)

  
  # ----- Predict CIF -----
  arm1 <- survfit(ph, data.frame(arm=1), se.fit=FALSE)
  ltr1[j] <- 1-tail(arm1$surv, n=1)
  
  arm0 <- survfit(ph, data.frame(arm=0), se.fit=FALSE)
  ltr0[j] <- 1-tail(arm0$surv, n=1)
  
  ltr.diff[j] <- ltr1[j]-ltr0[j]
  
}  


# ----- Output results ----- 
res <- cbind(ltr1, ltr0, ltr.diff)
colnames(res) <- c('ltr1', 'ltr0', 'ltr.diff')
write.csv(res, paste0("//restricted//projectnb//conner-thesis//lifetimerisk//simulation//finegray//cshv2//res_", a, ".csv"), row.names=FALSE)





