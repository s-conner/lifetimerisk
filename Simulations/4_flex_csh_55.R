# Simulation study to assess flexible parametric model using Geskus (2011) weighted data
# CSH Data
# Sarah Conner - Updated 1/20/2021

library(survival)
library(mstate)
library(rstpm2) 
options(scipen = 999)

a <- as.numeric(Sys.getenv("SGE_TASK_ID")) #Identify the job task ID from the t variables in the shell, see below
sim <- read.csv("//restricted//projectnb//conner-thesis//lifetimerisk//simulation//fg_splitjobs.csv", header=TRUE, sep=",")
sim <- sim[a, ]

l <- sim$simid
nstart <- sim$start
nend <- sim$end

ltr1 <- matrix(data=NA, nrow=500, ncol=3)
ltr0 <- matrix(data=NA, nrow=500, ncol=3)
ltr.diff <- matrix(data=NA, nrow=500, ncol=3)

# For each simulated dataset, fit the model, predict LTRs
j <- 0

for(k in nstart:nend){
  print(k)
  
  j <- j + 1
  skip.e <- skip.w <- FALSE
  if(exists('fp')){rm(fp)} # clear model from previous iteration
  if(exists('s')){rm(s)} # clear model from previous iteration
  
  
  # ----- Read in data ----- 
  dat <- read.csv(paste0("//restricted//projectnb//conner-thesis//lifetimerisk//simulation//data_csh//data", l, "_", k, ".csv"), header=TRUE, sep=",")

  dat$cause[dat$times>40] <- 0
  dat$times[dat$times>40] <- 40
  
  # Shift data to be similar to FHS; lifetime risk at age 95
  dat$entry <- 55+round(dat$entry, 4)
  dat$times <- 55+round(dat$times, 4)
  
  dat <- dat[dat$entry!=dat$times,] # In case the rounding means entry and exit times are close together
  
  
  # ----- LTRC weighted dataset in long form to accommodate log(time) interaction  ----- 
  
  # LTRC weights (Geskus 2011)
  dat.ltrc <- crprep(Tstop='times', status='cause', data=dat, trans=1, cens=0, Tstart='entry', id='id', keep='arm', shorten=FALSE)
  dat.ltrc$weight.ltrc <- dat.ltrc$weight.cens*dat.ltrc$weight.trunc
  dat.ltrc$ev1 <- ifelse(dat.ltrc$status==1,1,0)
  
  # Fit model with CLL link
  tryCatch({fp <- stpm2(Surv(Tstart-55, Tstop-55, ev1) ~ arm, weight=weight.ltrc, data=dat.ltrc, df=3, tvc=list(arm=1), link='PH')
  }, error=function(e) {skip.e <<- TRUE}, warning=function(w) {skip.w <<- TRUE})
  if(skip.e | skip.w) { next }     
  if(exists('fp')){if(fp@details$convergence!=0) { next }}     
  
  # Predictions are bad when model doesn't fully converge, don't use
  tryCatch({s <- summary(fp)
  }, error=function(e) {skip.e <<- TRUE}, warning=function(w) {skip.w <<- TRUE})
  if(skip.e | skip.w) { next }
  
  
  # ----- Predict CIF and Lifetime Risks -----
  
  min.death.id <- dat$id[dat$times == min(dat$times[dat$cause==2])]
  min.death.id <- min.death.id[1]
  intervals <- dat.ltrc[dat.ltrc$id==min.death.id, c('Tstart', 'Tstop')]
  intervals$ev1 <- 1
  arm1 <- data.frame(arm=1, intervals)
  arm0 <- data.frame(arm=0, intervals)
  
  # Could also use type='fail', but there is no type='faildiff' to obtain SE and CI of the contrast/difference
  # Instead, we take the complement of survival probability (works here b/c we used FG-weighted Cox model)
  tryCatch({pred1 <- tail(predict(fp, newdata=arm1, type='surv', se.fit=TRUE), n=1)
  }, error=function(e) {skip.e <<- TRUE}, warning=function(w) {skip.w <<- TRUE})
  if(skip.e | skip.w) { next }     
  
  tryCatch({pred0 <- tail(predict(fp, newdata=arm0, type='surv', se.fit=TRUE), n=1)
  }, error=function(e) {skip.e <<- TRUE}, warning=function(w) {skip.w <<- TRUE})
  if(skip.e | skip.w) { next }     
  
  tryCatch({pred.diff <- tail(predict(fp, newdata=arm0, var='arm', type='sdiff', se.fit=TRUE), n=1)
  }, error=function(e) {skip.e <<- TRUE}, warning=function(w) {skip.w <<- TRUE})
  if(skip.e | skip.w) { next }    
  
  # Alternative syntax using plot(), allows specifying other contrasts than 1 vs. 0
  # fig <- plot(fp, newdata=arm0, type='sdiff', exposed=function(data) transform(data, arm=1))
  
  # Take complement F(t)=1-S(t), flip sign of the difference, flip CIL and CIU
  ltr1[j,] <- as.matrix(unname(1-pred1[c(1,3,2)]))
  ltr0[j,] <- as.matrix(unname(1-pred0[c(1,3,2)]))
  ltr.diff[j,] <- as.matrix(unname((-1) * pred.diff[c(1,3,2)])) 
}  


# ----- Output results ----- 
res <- cbind(ltr1, ltr0, ltr.diff)
colnames(res) <- c('ltr1', 'ltr1.cil', 'ltr1.ciu',
                   'ltr0', 'ltr0.cil', 'ltr0.ciu', 
                   'ltr.diff', 'ltr.diff.cil', 'ltr.diff.ciu')
write.csv(res, paste0("//restricted//projectnb//conner-thesis//lifetimerisk//simulation//flex//cshv2_sub55//res_", l, "_", nend, ".csv"), row.names=FALSE)





