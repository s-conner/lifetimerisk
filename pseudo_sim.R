# Simulate 2 groups according to Beyersmann Stat Med 2009

library(survival)
library(geepack)

# Load simulation parameters from shared computing cluster
l <- as.numeric(Sys.getenv("SGE_TASK_ID")) #Identify the job task ID from the t variables in the shell, see below
if (is.na(l)) l <- 1 

nsim <- 2000
res_logit <- matrix(data=NA, nrow=nsim, ncol=18)
res_cll <- matrix(data=NA, nrow=nsim, ncol=18)


for(j in 1:nsim){
  
  #### Step 1 - Load data ####
  data <- read.csv(paste0("//restricted//projectnb//conner-thesis//lifetimerisk//sim1nov4//datacsh//data", l, "_", j, ".csv"), header=TRUE, sep=",")
  
  #### Step 2 - Derive pseudo-observations ####  
  data$pseudo <- pseudo.lr(data, 'entry', 'times', 'cause', 1, 40) 
  
  #### Step 3 - Fit models ####  
  cll.fit <- geese(pseudo ~ arm, 
                   data=data, id=id, jack = TRUE, scale.fix=TRUE, family=gaussian,
                   mean.link = "cloglog", corstr="independence")
  
  cll.res <- cbind(mean = cll.fit$beta, 
                   san.se = sqrt(diag(cll.fit$vbeta)),
                   san.pval = 2-2*pnorm(abs(cll.fit$beta/sqrt(diag(cll.fit$vbeta)))))
  
  logit.fit <- geese(pseudo ~ arm, 
                   data=data, id=id, jack = TRUE, scale.fix=TRUE, family=gaussian,
                   mean.link = "logit", corstr="independence")
  
  logit.res <- cbind(mean = logit.fit$beta, 
                   san.se = sqrt(diag(logit.fit$vbeta)),
                   san.pval = 2-2*pnorm(abs(logit.fit$beta/sqrt(diag(logit.fit$vbeta)))))
                  

  #### Step 4 - Derive predicted LTR and SE for each group ####
  ### Logit ###
  logit.beta0 <- unname(logit.fit$beta[1])
  logit.beta1 <- unname(logit.fit$beta[2])
  logit.ltr0 <- 1/(1 + exp(-logit.beta0))
  logit.ltr1 <- 1/(1 + exp(-logit.beta0 - logit.beta1))
  logit.diff <- logit.ltr1 - logit.ltr0
  
  # ltr1 se
  logit.ltr1.dbeta1 <- exp(-(logit.beta0 + logit.beta1))/(1 + exp(-(logit.beta0 + logit.beta1)))^2
  logit.ltr1.dbeta0 <- logit.ltr1.dbeta1
  logit.ltr1.grad <- c(logit.ltr1.dbeta0, logit.ltr1.dbeta1)
  logit.ltr1.se <- sqrt(t(logit.ltr1.grad) %*% logit.fit$vbeta %*% logit.ltr1.grad)
  #deltamethod( ~ (1/(1 + exp(-(x1 + x2)))), c(logit.beta0, logit.beta1), logit.fit$vbeta)
  
  # ltr0 se
  logit.ltr0.dbeta1 <- 0
  logit.ltr0.dbeta0 <- exp(-logit.beta0)/(1 + exp(-logit.beta0))^2
  logit.ltr0.grad <- c(logit.ltr0.dbeta0, logit.ltr0.dbeta1)
  logit.ltr0.se <- sqrt(t(logit.ltr0.grad) %*% logit.fit$vbeta %*% logit.ltr0.grad)
  #deltamethod( ~ (1/(1 + exp(-x1))), c(logit.beta0, logit.beta1), logit.fit$vbeta)
  
  # diff se
  logit.diff.dbeta1 <- logit.ltr1.dbeta1 - logit.ltr0.dbeta1
  logit.diff.dbeta0 <- logit.ltr1.dbeta0 - logit.ltr0.dbeta0
  logit.diff.grad <- c(logit.diff.dbeta0, logit.diff.dbeta1)
  logit.diff.se <- sqrt(t(logit.diff.grad) %*% logit.fit$vbeta %*% logit.diff.grad)
  #deltamethod( ~ (1/(1 + exp(-(x1 + x2)))) - (1/(1 + exp(-x1))), c(logit.beta0, logit.beta1), logit.fit$vbeta)
  

  ### CLL ###
  cll.beta0 <- unname(cll.fit$beta[1])
  cll.beta1 <- unname(cll.fit$beta[2])
  cll.ltr0 <- 1-exp(-exp(cll.beta0))
  cll.ltr1 <- 1-exp(-exp(cll.beta0 + cll.beta1))
  cll.diff <- cll.ltr1 - cll.ltr0

  # ltr1 se
  cll.ltr1.dbeta1 <- exp(-exp(cll.beta0 + cll.beta1) + cll.beta0 + cll.beta1)
  cll.ltr1.dbeta0 <- cll.ltr1.dbeta1
  cll.ltr1.grad <- c(cll.ltr1.dbeta0, cll.ltr1.dbeta1)
  cll.ltr1.se <- sqrt(t(cll.ltr1.grad) %*% cll.fit$vbeta %*% cll.ltr1.grad)
  #deltamethod( ~ (1 - exp(-exp(x1 + x2))), c(cll.beta0, cll.beta1), cll.fit$vbeta)
  
  # ltr0 se
  cll.ltr0.dbeta1 <- 0
  cll.ltr0.dbeta0 <- exp(-exp(cll.beta0) + cll.beta0)
  cll.ltr0.grad <- c(cll.ltr0.dbeta0, cll.ltr0.dbeta1)
  cll.ltr0.se <- sqrt(t(cll.ltr0.grad) %*% cll.fit$vbeta %*% cll.ltr0.grad)
  #deltamethod( ~ (1 - exp(-exp(x1))), c(cll.beta0, cll.beta1), cll.fit$vbeta)
  
  # diff se
  cll.diff.dbeta1 <- cll.ltr1.dbeta1 - cll.ltr0.dbeta1
  cll.diff.dbeta0 <- cll.ltr1.dbeta0 - cll.ltr0.dbeta0
  cll.diff.grad <- c(cll.diff.dbeta0, cll.diff.dbeta1)
  cll.diff.se <- sqrt(t(cll.diff.grad) %*% cll.fit$vbeta %*% cll.diff.grad)
  #deltamethod( ~ (1 - exp(-exp(x1 + x2))) - (1 - exp(-exp(x1))), c(cll.beta0, cll.beta1), cll.fit$vbeta)
  

  #### Step 5 - Export!
  res_logit[j,] <- c(logit.res[1,1], logit.res[1,2], logit.res[1,3], logit.res[2,1], logit.res[2,2], logit.res[2,3], 
                    logit.ltr0, logit.ltr0.se, logit.ltr0 - 1.96*logit.ltr0.se, logit.ltr0 + 1.96*logit.ltr0.se,
                    logit.ltr1, logit.ltr1.se, logit.ltr1 - 1.96*logit.ltr1.se, logit.ltr1 + 1.96*logit.ltr1.se,
                    logit.diff, logit.diff.se, logit.diff - 1.96*logit.diff.se, logit.diff + 1.96*logit.diff.se)
  
  res_cll[j,] <- c(cll.res[1,1], cll.res[1,2], cll.res[1,3], cll.res[2,1], cll.res[2,2], cll.res[2,3], 
                  cll.ltr0, cll.ltr0.se, cll.ltr0 - 1.96*cll.ltr0.se, cll.ltr0 + 1.96*cll.ltr0.se,
                  cll.ltr1, cll.ltr1.se, cll.ltr1 - 1.96*cll.ltr1.se, cll.ltr1 + 1.96*cll.ltr1.se,
                  cll.diff, cll.diff.se, cll.diff - 1.96*cll.diff.se, cll.diff + 1.96*cll.diff.se)

}  

# Save results
colnames(res_logit) <- c('beta0', 'beta0_SE', 'beta0_pval', 'beta1', 'beta1_SE', 'beta1_pval', 
                         'ltr0', 'ltr0_SE', 'ltr0_CIL', 'ltr0_CIU', 
                         'ltr1', 'ltr1_SE', 'ltr1_CIL', 'ltr1_CIU',
                         'ltrdiff', 'ltrdiff_SE', 'ltrdiff_CIL', 'ltrdiff_CIU')
colnames(res_cll) <- c('beta0', 'beta0_SE', 'beta0_pval', 'beta1', 'beta1_SE', 'beta1_pval', 
                         'ltr0', 'ltr0_SE', 'ltr0_CIL', 'ltr0_CIU', 
                         'ltr1', 'ltr1_SE', 'ltr1_CIL', 'ltr1_CIU',
                         'ltrdiff', 'ltrdiff_SE', 'ltrdiff_CIL', 'ltrdiff_CIU')

write.csv(res_logit, paste0("//restricted//projectnb//conner-thesis//lifetimerisk//sim1nov4//logitres//logitres", l, ".csv"), row.names=FALSE)
write.csv(res_cll, paste0("//restricted//projectnb//conner-thesis//lifetimerisk//sim1nov4//cllres//cllres", l, ".csv"), row.names=FALSE)

