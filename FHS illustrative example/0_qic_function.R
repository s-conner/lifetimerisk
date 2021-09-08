
# QIC function designed for geese() model with normal distribution and cloglog, logit, or identity link
# model = geese model
# X = covariate matrix including vector of 1's for intercept

#model <- test; data <- af;

qic <- function(model, data, outc){
  require(MASS)
  model.indep <- update(model, corstr = "independence")
  
  # Predicted values
  names(model$beta) <- NULL
  if(length(model$xnames)>1){ 
    X <- as.matrix(cbind(rep(1, nrow(data)), data[model$xnames[2:length(model$xnames)]]))
    mu.R <- switch(model$model$mean.link, 
                   cloglog =  1 - exp(-exp(X %*% model$beta)),
                   logit = 1/(1 + exp(-(X %*% model$beta))),
                   identity = X %*% model$beta)
    
  } else{ 
    X <- rep(1, nrow(data))
    mu.R <- switch(model$model$mean.link, 
                   cloglog =  1 - exp(-exp(X*model$beta)),
                   logit = 1/(1 + exp(-(X*model$beta))),
                   identity = X*model$beta)
  }
  
  # Quasi-likelihood
  y <- outc
  names(model$gamma) <- NULL
  quasi.R <- (sum(((y - mu.R)^2)/-2))/model$gamma
  
  # Trace 
  omegaI <- ginv(model.indep$vbeta.naiv) 
  Vr <- model$vbeta
  trace.R <- sum(diag(omegaI %*% Vr))
  px <- length(model.indep$xnames) 
  
  # QIC
  QIC <- -2*quasi.R  + 2*trace.R  # Not really the same as BIC. trace is not the same as log(n)*px
  QICu <- -2*quasi.R + 2*px   # Matches AIC
  output <- cbind(QIC, QICu, quasi.R)
  colnames(output) <- c('QIC', 'QICu', 'LogQuasiLik')
  return(output)
}
