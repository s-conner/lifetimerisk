# Generate Non-PSH according to Li 2015

library(survival)
#library(etm)

l <- as.numeric(Sys.getenv("SGE_TASK_ID")) 
if (is.na(l)) l <- 1 
scenarios <- read.csv("//restricted//projectnb//conner-thesis//lifetimerisk//sim2nov11//table_scenarios_sdh.csv", header=TRUE, sep=",")
#scenarios <- read.csv("Y:\\19SCR-lifetimerisk_varselection\\simulation\\objective 2\\table_scenarios_sdh.csv", header=TRUE, sep=",")

nsim <- 2000

# sampsize=10;
# n=sampsize/2; z=1; gamma=1.2; rho=-2; psi1=.5; theta=-2.2; psi2=.5
# censmin=0; censmax=1.2; ptrunc=.5

gen <- function(n, z, gamma, rho, psi1, theta, psi2){
  
  p2 <- exp((gamma*exp(psi1*z))/(rho + theta*z))
  cause <- 1 + rbinom(n, 1, p2)
  
  u <- runif(n, 0, 1)
  t <- rep(NA, n)
  
  for(i in 1:n){
    if(cause[i]==1){ # Conditional CIF for event 1
      a <- 1 - exp(gamma*exp(psi1*z)/(rho + theta*z))
      b <- (log(1 - u[i]*a) * (rho + theta*z))/(gamma * exp(psi1*z))
      t[i] <- (log(1 - b))/(rho + theta*z)
    } else { # Conditional CIF for event 2
      t[i] <- -(log(1-u[i]))/exp(psi2*z)
    }
  }
  return(cbind(cause, t))
}

simulate <- function(censmin, censmax, ptrunc, sampsize, gamma, rho, psi1, theta, psi2){
  
  # Generate causes and corresponding event times
  arm <- c(rep(0, (sampsize/2)), rep(1,(sampsize/2)))
  out0 <- gen(n=(sampsize/2), z=0, gamma=gamma, rho=rho, psi1=psi1, theta=theta, psi2=psi2)
  out1 <- gen(n=(sampsize/2), z=1, gamma=gamma, rho=rho, psi1=psi1, theta=theta, psi2=psi2)
  times <- rbind(out0, out1)
  times0 <- times[,2]
  cause0 <- times[,1]
  
  # Introduce censoring and truncation
  cens <- runif(sampsize, censmin, censmax)
  delayed <- rbinom(sampsize, 1, ptrunc)
  left <- delayed*runif(sampsize, 0, .125)
  
  times1 <- pmin(cens, times0) # Determine if censoring occurs before the event time
  cause1 <-  c(times1 < cens) * cause0 # Update cause to what we observe w/ censoring
  # Determine if entry occurs after event time
  # Rounding to prevent errors, i.e. in survival package
  include <- ifelse(round(times1 - left, 5) <= 0, 0, 1) 
  
  #obs.time <- pmin(rep(40, sampsize), times1) # Administratively censor at tau=40
  #obs.cause <-  c(obs.time==times1) * cause1 # Update cause for administrative censoring
  
  # Export non-truncated individuals as data frame
  #all.data <- data.frame(include=include, entry=left, times=obs.time, cause=obs.cause, arm, censtime=cens, origtime=times0)
  all.data <- data.frame(include=include, entry=left, times=times1, 
                         cause=cause1, arm, censtime=cens, origtime=times0)
  nottrunc.data <- all.data[all.data$include==1, ]
  nottrunc.data$id <- 1:nrow(nottrunc.data)
  return(nottrunc.data)
}


j <- 1
k <- 1
while(j<=nsim){
  
  dat <- simulate(scenarios$censmin[l], scenarios$censmax[l], scenarios$ptrunc[l], scenarios$sampsize[l],
                  scenarios$gamma[l], scenarios$rho[l], scenarios$psi1[l], scenarios$theta[l], scenarios$psi2[l])
  
  data <- dat[dat$include==1, ]
  
  #plot(etmCIF(Surv(entry, times, cause != 0) ~ arm, data=data, etype=cause, failcode=1), lty=c(2,1), xlim=c(0,1.5))

  if((sum(data$times[data$arm==0]>=1)>2) & (sum(data$times[data$arm==1]>=1)>2)){
    data$entry <- data$entry*40
    data$times <- data$times*40
    data$censtime <- data$censtime*40
    data$times <- data$origtimes*40
    
    write.csv(data, paste0("//restricted//projectnb//conner-thesis//lifetimerisk//sim2nov11//datasets//data", l, "_", j, ".csv"), row.names=FALSE)
    j <- j+1
  }
  k <- k+1
}

write.csv(k, paste0("//restricted//projectnb//conner-thesis//lifetimerisk//sim2nov11//datasets//numiter", l, ".csv"), row.names=FALSE)

