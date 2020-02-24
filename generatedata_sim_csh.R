# Simulate 2 groups with 'direct' CSH method
# Save datasets for further analyses

library(survival)

l <- as.numeric(Sys.getenv("SGE_TASK_ID")) 
if (is.na(l)) l <- 1 
scenarios <- read.csv("//restricted//projectnb//conner-thesis//lifetimerisk//sim1nov4//table_scenarios_csh.csv", header=TRUE, sep=",")
#scenarios <- read.csv("Y:\\19SCR-lifetimerisk_varselection\\simulation\\objective 1 nov4\\table_scenarios_csh.csv", header=TRUE, sep=",")
nsim <- 2000

# First, define function f(x)=0 for which we seek the root. 
# F(x)=cumulative all-cause hazard + y, where y=(ln(1-u)).
fxn.tosolve <- function(x, y, a1, b1, a2, b2) { return( ((x/b1)^a1) + ((x/b2)^a2) + y ) }


# Function to generate event times by solving inverse of cumulative all-cause hazard + y
generate.times <- function(n, max.int, a1, b1, a2, b2) {
  
  stime <- NULL
  i <- 1
  
  while(length(stime) < n) {
    u <- runif(1)
    # Bisection method of uniroot() requires endpoints to be opposite signs
    if (fxn.tosolve(0, log(1 - u), a1, b1, a2, b2) *
        fxn.tosolve(max.int, log(1 - u), a1, b1, a2, b2) < 0) {
      soln <- uniroot(fxn.tosolve, c(0, max.int), tol = 0.0001, y = log(1 - u), a1=a1, b1=b1, a2=a2, b2=b2)
      stime[i] <- soln$root
      i <- i + 1
    }
    else cat("Values at endpoints not of opposite signs.")
  }
  return(stime)
}


# Function to determine event type
cause <- function(x, a1, b1, a2, b2) {
  out <- rbinom(length(x), 1, 
                prob = (a1/(b1^a1))*(x^(a1-1)) / ((a1/(b1^a1))*(x^(a1-1)) + (a2/(b2^a2))*(x^(a2-1))) )
  ifelse(out == 0, 2, 1)
}


# Function to generate dataset
simulate <- function(censmin, censmax, ptrunc, sampsize, a01, b01, a11, b11, a02, b02, a12, b12){
  
  # Generate event times for each group
  # l=9
  # censmin=scenarios$censmin[l]; censmax=scenarios$censmax[l]; ptrunc=scenarios$ptrunc[l]; sampsize=scenarios$sampsize[l]
  # a01=scenarios$a01[l]; b01=scenarios$b01[l]; a02=scenarios$a02[l]; b02=scenarios$b02[l]
  # a11=scenarios$a11[l]; b11=scenarios$b11[l]; a12=scenarios$a12[l]; b12=scenarios$b12[l]

  times.z0 <- generate.times(n=(sampsize/2), max.int=99, a1=a01, b1=b01, a2=a02, b2=b02)
  times.z1 <- generate.times(n=(sampsize/2), max.int=99, a1=a11, b1=b11, a2=a12, b2=b12)
  
  # Assign cause for each group
  cause.z0 <- cause(times.z0, a1=a01, b1=b01, a2=a02, b2=b02)
  cause.z1 <- cause(times.z1, a1=a11, b1=b11, a2=a12, b2=b12)
  
  # Merge together the separate arms
  arm <- c(rep(0, (sampsize/2)), rep(1,(sampsize/2)))
  times0 <- c(times.z0, times.z1)
  cause0 <- c(cause.z0, cause.z1)
  
  # Introduce censoring and truncation
  cens <- runif(sampsize, censmin, censmax)
  delayed <- rbinom(sampsize, 1, ptrunc)
  left <- delayed*runif(sampsize, 0, 5)
  
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


# Call function to simulate data
j <- 1
k <- 1

while(j<=nsim){

  data <- simulate(censmin=scenarios$censmin[l], censmax=scenarios$censmax[l], ptrunc=scenarios$ptrunc[l], sampsize=scenarios$sampsize[l],
                  a01=scenarios$a01[l], b01=scenarios$b01[l], a11=scenarios$a11[l], b11=scenarios$b11[l], 
                  a02=scenarios$a02[l], b02=scenarios$b02[l], a12=scenarios$a12[l], b12=scenarios$b12[l])
  
  if((sum(data$times[data$arm==0]>=40)>2) & (sum(data$times[data$arm==1]>=40)>2)){
    write.csv(data, paste0("//restricted//projectnb//conner-thesis//lifetimerisk//sim1nov4//data//data", l, "_", j, ".csv"), row.names=FALSE)
    j <- j+1
  }
  k <- k+1
}

write.csv(k, paste0("//restricted//projectnb//conner-thesis//lifetimerisk//sim1nov4//data//numiter", l, "_", j, ".csv"), row.names=FALSE)


