# Assess simulation results and compile into a dataset for nested loop plot
# Updated for revision on 7/9/2021

options(scipen=999)

scenarios.csh <- read.csv('table_scenarios_csh.csv')
scenarios.sdh <- read.csv('table_scenarios_sdh.csv')

pseudocsh <- 'pseudo//cshv2//'
pseudosdh <- 'pseudo//sdhv2//'

fg.csh <- 'finegray//cshv2//'
fg.sdh <- 'finegray//sdhv2//'

fg.logt.csh <- 'finegray_logtime//cshv2//'
fg.logt.sdh <- 'finegray_logtime//sdhv2//'

fp.csh <- 'flex//cshv2_sub55//'
fp.sdh <- 'flex//sdhv2_sub55//'

fp.2df.csh <- 'flex//cshv2_sub55_2df//'
fp.2df.sdh <- 'flex//sdhv2_sub55_2df//'

fp.3df.csh <- 'flex//cshv2_sub55_3df//'
fp.3df.sdh <- 'flex//sdhv2_sub55_3df//'

scenarios.sdh$setting <- scenarios.sdh$setting+4

for(i in 1:32){
  # CSH Truth
  a01=scenarios.csh$a01[i]; b01=scenarios.csh$b01[i]; a02=scenarios.csh$a02[i]; b02=scenarios.csh$b02[i];
  a11=scenarios.csh$a11[i]; b11=scenarios.csh$b11[i]; a12=scenarios.csh$a12[i]; b12=scenarios.csh$b12[i];
  
  integrandz0 <- function(x) {((a01/(b01^a01))*(x^(a01-1))) * exp(-(((x/b01)^a01) + ((x/b02)^a02)))}
  integrandz1 <- function(x) {((a11/(b11^a11))*(x^(a11-1))) * exp(-(((x/b11)^a11) + ((x/b12)^a12)))}
  
  scenarios.csh$ltr0[i] <- integrate(integrandz0, lower = 0 , upper = 40)$value
  scenarios.csh$ltr1[i] <- integrate(integrandz1, lower = 0 , upper = 40)$value
  scenarios.csh$ltrdiff[i] <- scenarios.csh$ltr1[i] - scenarios.csh$ltr0[i]
  
  # SDH Truth
  gamma=scenarios.sdh$gamma[i]; rho=scenarios.sdh$rho[i]; psi1=scenarios.sdh$psi1[i]; theta=scenarios.sdh$theta[i];
  scenarios.sdh$ltr1[i] <- 1 - exp(gamma*exp(psi1)*(1-exp(rho + theta))/(rho + theta))
  scenarios.sdh$ltr0[i] <- 1 - exp(gamma*(1-exp(rho))/rho)
  scenarios.sdh$ltrdiff[i] <- scenarios.sdh$ltr1[i] - scenarios.sdh$ltr0[i]
}

scenarios <- rbind(scenarios.csh[, c('setting', 'sampsize', 'probcens', 'ptrunc', 'ltr0', 'ltr1', 'ltrdiff')],
                   scenarios.sdh[, c('setting', 'sampsize', 'probcens', 'ptrunc', 'ltr0', 'ltr1', 'ltrdiff')])

scenarios <- scenarios[order(scenarios$setting, scenarios$sampsize, scenarios$probcens, scenarios$ptrunc), ]

bias.ltr0 <- matrix(NA,64,ncol=6)
bias.ltr1 <- matrix(NA,64,ncol=6)
bias.ltrdiff <- matrix(NA,64,ncol=6)
relbias.ltr0 <- matrix(NA,64,ncol=6)
relbias.ltr1 <- matrix(NA,64,ncol=6)
relbias.ltrdiff <- matrix(NA,64,ncol=6)
rmse.ltr0 <- matrix(NA,64,ncol=6)
rmse.ltr1 <- matrix(NA,64,ncol=6)
rmse.ltrdiff <- matrix(NA,64,ncol=6)
cov.ltr0 <- matrix(NA,64,ncol=4)
cov.ltr1 <- matrix(NA,64,ncol=4)
cov.ltrdiff <- matrix(NA,64,ncol=4)
power <- rep(NA,64)
effsampsize <- matrix(NA,64,ncol=6)


for(i in 1:64){
  
  # Truth
  ltr0 <- scenarios$ltr0[i]
  ltr1 <- scenarios$ltr1[i]
  ltrdiff <- ltr1-ltr0
  
  # Logit results
  if(i<33){
    res.p <- read.csv(paste0(pseudocsh, 'logitres', i, '.csv'))
  } else {
    res.p <- read.csv(paste0(pseudosdh, 'logitres', i-32, '.csv'))
  }
  
  res.p <- res.p[!is.na(res.p$ltrdiff),]
  effsampsize[i,1] <- nrow(res.p)
  
  avgest.ltr0 <- mean(res.p$ltr0)
  avgest.ltr1 <- mean(res.p$ltr1)
  avgest.ltrdiff <- mean(res.p$ltrdiff)
  
  bias.ltr0[i,1] <- avgest.ltr0-ltr0
  bias.ltr1[i,1] <- avgest.ltr1-ltr1
  bias.ltrdiff[i,1] <- avgest.ltrdiff-ltrdiff
  
  relbias.ltr0[i,1] <- bias.ltr0[i,1]/ltr0
  relbias.ltr1[i,1] <- bias.ltr1[i,1]/ltr1
  relbias.ltrdiff[i,1] <- bias.ltrdiff[i,1]/ltrdiff
  
  rmse.ltr0[i,1] <- sqrt(mean((rep(ltr0, effsampsize[i,1])-res.p$ltr0)^2))
  rmse.ltr1[i,1] <- sqrt(mean((rep(ltr1, effsampsize[i,1])-res.p$ltr1)^2))
  rmse.ltrdiff[i,1] <- sqrt(mean((rep(ltrdiff, effsampsize[i,1])-res.p$ltrdiff)^2))
  
  cov.ltr0[i,1] <- sum((res.p$ltr0_CIL <= ltr0) & (res.p$ltr0_CIU >= ltr0))/effsampsize[i,1]
  cov.ltr1[i,1] <- sum((res.p$ltr1_CIL <= ltr1) & (res.p$ltr1_CIU >= ltr1))/effsampsize[i,1]
  cov.ltrdiff[i,1] <- sum((res.p$ltrdiff_CIL <= ltrdiff) & (res.p$ltrdiff_CIU >= ltrdiff))/effsampsize[i,1]
  
  power[i] <- sum(res.p$beta1_pval <= .05)/effsampsize[i,1]
  
  
  # Fine-Gray results
  if(i<33){
    res.fg <- read.csv(paste0(fg.csh, 'res_', i, '.csv'))
  } else {
    res.fg <- read.csv(paste0(fg.sdh, 'res_', i-32, '.csv'))
  }
  
  res.fg <- res.fg[!is.na(res.fg$ltr.diff),]
  effsampsize[i,2] <- nrow(res.fg)
  
  avgest.ltr0 <- mean(res.fg$ltr0)
  avgest.ltr1 <- mean(res.fg$ltr1)
  avgest.ltrdiff <- mean(res.fg$ltr.diff)
  
  bias.ltr0[i,2] <- avgest.ltr0-ltr0
  bias.ltr1[i,2] <- avgest.ltr1-ltr1
  bias.ltrdiff[i,2] <- avgest.ltrdiff-ltrdiff
  
  relbias.ltr0[i,2] <- bias.ltr0[i,2]/ltr0
  relbias.ltr1[i,2] <- bias.ltr1[i,2]/ltr1
  relbias.ltrdiff[i,2] <- bias.ltrdiff[i,2]/ltrdiff
  
  rmse.ltr0[i,2] <- sqrt(mean((rep(ltr0, effsampsize[i,2])-res.fg$ltr0)^2))
  rmse.ltr1[i,2] <- sqrt(mean((rep(ltr1, effsampsize[i,2])-res.fg$ltr1)^2))
  rmse.ltrdiff[i,2] <- sqrt(mean((rep(ltrdiff, effsampsize[i,2])-res.fg$ltr.diff)^2))

  
  
  # Fine-Gray with logtime results
  if(i<33){
    res.fg.logt.path1 <- read.csv(paste0(fg.logt.csh, 'res_', i, '_500.csv'))
    res.fg.logt.path2 <- read.csv(paste0(fg.logt.csh, 'res_', i, '_1000.csv'))
    res.fg.logt.path3 <- read.csv(paste0(fg.logt.csh, 'res_', i, '_1500.csv'))
    res.fg.logt.path4 <- read.csv(paste0(fg.logt.csh, 'res_', i, '_2000.csv'))

  } else {
    res.fg.logt.path1 <- read.csv(paste0(fg.logt.sdh, 'res_', i-32, '_500.csv'))
    res.fg.logt.path2 <- read.csv(paste0(fg.logt.sdh, 'res_', i-32, '_1000.csv'))
    res.fg.logt.path3 <- read.csv(paste0(fg.logt.sdh, 'res_', i-32, '_1500.csv'))
    res.fg.logt.path4 <- read.csv(paste0(fg.logt.sdh, 'res_', i-32, '_2000.csv'))
  }
  
  res.fg.logt <- rbind(res.fg.logt.path1, res.fg.logt.path2, res.fg.logt.path3, res.fg.logt.path4)
  res.fg.logt <- res.fg.logt[!is.na(res.fg.logt$ltr.diff),]
  effsampsize[i,3] <- nrow(res.fg.logt)
  
  avgest.ltr0 <- mean(res.fg.logt$ltr0)
  avgest.ltr1 <- mean(res.fg.logt$ltr1)
  avgest.ltrdiff <- mean(res.fg.logt$ltr.diff)
  
  bias.ltr0[i,3] <- avgest.ltr0-ltr0
  bias.ltr1[i,3] <- avgest.ltr1-ltr1
  bias.ltrdiff[i,3] <- avgest.ltrdiff-ltrdiff
  
  relbias.ltr0[i,3] <- bias.ltr0[i,3]/ltr0
  relbias.ltr1[i,3] <- bias.ltr1[i,3]/ltr1
  relbias.ltrdiff[i,3] <- bias.ltrdiff[i,3]/ltrdiff
  
  rmse.ltr0[i,3] <- sqrt(mean((rep(ltr0, effsampsize[i,3])-res.fg.logt$ltr0)^2))
  rmse.ltr1[i,3] <- sqrt(mean((rep(ltr1, effsampsize[i,3])-res.fg.logt$ltr1)^2))
  rmse.ltrdiff[i,3] <- sqrt(mean((rep(ltrdiff, effsampsize[i,3])-res.fg.logt$ltr.diff)^2))
  
  
  
  # Flexible parametric results
  if(i<33){
    res.fp.path1 <- read.csv(paste0(fp.csh, 'res_', i, '_500.csv'))
    res.fp.path2 <- read.csv(paste0(fp.csh, 'res_', i, '_1000.csv'))
    res.fp.path3 <- read.csv(paste0(fp.csh, 'res_', i, '_1500.csv'))
    res.fp.path4 <- read.csv(paste0(fp.csh, 'res_', i, '_2000.csv'))
    
  } else {
    res.fp.path1 <- read.csv(paste0(fp.sdh, 'res_', i-32, '_500.csv'))
    res.fp.path2 <- read.csv(paste0(fp.sdh, 'res_', i-32, '_1000.csv'))
    res.fp.path3 <- read.csv(paste0(fp.sdh, 'res_', i-32, '_1500.csv'))
    res.fp.path4 <- read.csv(paste0(fp.sdh, 'res_', i-32, '_2000.csv'))
  }
  
  res.fp <- rbind(res.fp.path1, res.fp.path2, res.fp.path3, res.fp.path4)
  res.fp <- res.fp[!is.na(res.fp$ltr.diff),]
  effsampsize[i,4] <- nrow(res.fp)
  
  avgest.ltr0 <- mean(res.fp$ltr0)
  avgest.ltr1 <- mean(res.fp$ltr1)
  avgest.ltrdiff <- mean(res.fp$ltr.diff)
    
  bias.ltr0[i,4] <- avgest.ltr0-ltr0
  bias.ltr1[i,4] <- avgest.ltr1-ltr1
  bias.ltrdiff[i,4] <- avgest.ltrdiff-ltrdiff
  
  relbias.ltr0[i,4] <- bias.ltr0[i,4]/ltr0
  relbias.ltr1[i,4] <- bias.ltr1[i,4]/ltr1
  relbias.ltrdiff[i,4] <- bias.ltrdiff[i,4]/ltrdiff
  
  rmse.ltr0[i,4] <- sqrt(mean((rep(ltr0, effsampsize[i,4])-res.fp$ltr0)^2))
  rmse.ltr1[i,4] <- sqrt(mean((rep(ltr1, effsampsize[i,4])-res.fp$ltr1)^2))
  rmse.ltrdiff[i,4] <- sqrt(mean((rep(ltrdiff, effsampsize[i,4])-res.fp$ltr.diff)^2))
  
  cov.ltr0[i,2] <- sum((res.fp$ltr0.cil <= ltr0) & (res.fp$ltr0.ciu >= ltr0))/effsampsize[i,4]
  cov.ltr1[i,2] <- sum((res.fp$ltr1.cil <= ltr1) & (res.fp$ltr1.ciu >= ltr1))/effsampsize[i,4]
  cov.ltrdiff[i,2] <- sum((res.fp$ltr.diff.cil <= ltrdiff) & (res.fp$ltr.diff.ciu >= ltrdiff))/effsampsize[i,4]
  
  
  # Flexible parametric w/ 2df results
  if(i<33){
    res.fp.2df.path1 <- read.csv(paste0(fp.2df.csh, 'res_', i, '_500.csv'))
    res.fp.2df.path2 <- read.csv(paste0(fp.2df.csh, 'res_', i, '_1000.csv'))
    res.fp.2df.path3 <- read.csv(paste0(fp.2df.csh, 'res_', i, '_1500.csv'))
    res.fp.2df.path4 <- read.csv(paste0(fp.2df.csh, 'res_', i, '_2000.csv'))
    
  } else {
    res.fp.2df.path1 <- read.csv(paste0(fp.2df.sdh, 'res_', i-32, '_500.csv'))
    res.fp.2df.path2 <- read.csv(paste0(fp.2df.sdh, 'res_', i-32, '_1000.csv'))
    res.fp.2df.path3 <- read.csv(paste0(fp.2df.sdh, 'res_', i-32, '_1500.csv'))
    res.fp.2df.path4 <- read.csv(paste0(fp.2df.sdh, 'res_', i-32, '_2000.csv'))
  }
  
  res.fp.2df <- rbind(res.fp.2df.path1, res.fp.2df.path2, res.fp.2df.path3, res.fp.2df.path4)
  res.fp.2df <- res.fp.2df[!is.na(res.fp.2df$ltr.diff),]
  effsampsize[i,5] <- nrow(res.fp.2df)
  
  avgest.ltr0 <- mean(res.fp.2df$ltr0)
  avgest.ltr1 <- mean(res.fp.2df$ltr1)
  avgest.ltrdiff <- mean(res.fp.2df$ltr.diff)
  
  bias.ltr0[i,5] <- avgest.ltr0-ltr0
  bias.ltr1[i,5] <- avgest.ltr1-ltr1
  bias.ltrdiff[i,5] <- avgest.ltrdiff-ltrdiff
  
  relbias.ltr0[i,5] <- bias.ltr0[i,5]/ltr0
  relbias.ltr1[i,5] <- bias.ltr1[i,5]/ltr1
  relbias.ltrdiff[i,5] <- bias.ltrdiff[i,5]/ltrdiff
  
  rmse.ltr0[i,5] <- sqrt(mean((rep(ltr0, effsampsize[i,5])-res.fp.2df$ltr0)^2))
  rmse.ltr1[i,5] <- sqrt(mean((rep(ltr1, effsampsize[i,5])-res.fp.2df$ltr1)^2))
  rmse.ltrdiff[i,5] <- sqrt(mean((rep(ltrdiff, effsampsize[i,5])-res.fp.2df$ltr.diff)^2))
  
  cov.ltr0[i,3] <- sum((res.fp.2df$ltr0.cil <= ltr0) & (res.fp.2df$ltr0.ciu >= ltr0))/effsampsize[i,5]
  cov.ltr1[i,3] <- sum((res.fp.2df$ltr1.cil <= ltr1) & (res.fp.2df$ltr1.ciu >= ltr1))/effsampsize[i,5]
  cov.ltrdiff[i,3] <- sum((res.fp.2df$ltr.diff.cil <= ltrdiff) & (res.fp.2df$ltr.diff.ciu >= ltrdiff))/effsampsize[i,5]
  
  

 # Flexible parametric w/ 3df results
  if(i<33){
    res.fp.3df.path1 <- read.csv(paste0(fp.3df.csh, 'res_', i, '_500.csv'))
    res.fp.3df.path2 <- read.csv(paste0(fp.3df.csh, 'res_', i, '_1000.csv'))
    res.fp.3df.path3 <- read.csv(paste0(fp.3df.csh, 'res_', i, '_1500.csv'))
    res.fp.3df.path4 <- read.csv(paste0(fp.3df.csh, 'res_', i, '_2000.csv'))

  } else {
    res.fp.3df.path1 <- read.csv(paste0(fp.3df.sdh, 'res_', i-32, '_500.csv'))
    res.fp.3df.path2 <- read.csv(paste0(fp.3df.sdh, 'res_', i-32, '_1000.csv'))
    res.fp.3df.path3 <- read.csv(paste0(fp.3df.sdh, 'res_', i-32, '_1500.csv'))
    res.fp.3df.path4 <- read.csv(paste0(fp.3df.sdh, 'res_', i-32, '_2000.csv'))
  }

  res.fp.3df <- rbind(res.fp.3df.path1, res.fp.3df.path2, res.fp.3df.path3, res.fp.3df.path4)
  res.fp.3df <- res.fp.3df[!is.na(res.fp.3df$ltr.diff),]
  effsampsize[i,6] <- nrow(res.fp.3df)

  avgest.ltr0 <- mean(res.fp.3df$ltr0)
  avgest.ltr1 <- mean(res.fp.3df$ltr1)
  avgest.ltrdiff <- mean(res.fp.3df$ltr.diff)

  bias.ltr0[i,6] <- avgest.ltr0-ltr0
  bias.ltr1[i,6] <- avgest.ltr1-ltr1
  bias.ltrdiff[i,6] <- avgest.ltrdiff-ltrdiff

  relbias.ltr0[i,6] <- bias.ltr0[i,6]/ltr0
  relbias.ltr1[i,6] <- bias.ltr1[i,6]/ltr1
  relbias.ltrdiff[i,6] <- bias.ltrdiff[i,6]/ltrdiff

  rmse.ltr0[i,6] <- sqrt(mean((rep(ltr0, effsampsize[i,6])-res.fp.3df$ltr0)^2))
  rmse.ltr1[i,6] <- sqrt(mean((rep(ltr1, effsampsize[i,6])-res.fp.3df$ltr1)^2))
  rmse.ltrdiff[i,6] <- sqrt(mean((rep(ltrdiff, effsampsize[i,6])-res.fp.3df$ltr.diff)^2))

  cov.ltr0[i,4] <- sum((res.fp.3df$ltr0.cil <= ltr0) & (res.fp.3df$ltr0.ciu >= ltr0))/effsampsize[i,6]
  cov.ltr1[i,4] <- sum((res.fp.3df$ltr1.cil <= ltr1) & (res.fp.3df$ltr1.ciu >= ltr1))/effsampsize[i,6]
  cov.ltrdiff[i,4] <- sum((res.fp.3df$ltr.diff.cil <= ltrdiff) & (res.fp.3df$ltr.diff.ciu >= ltrdiff))/effsampsize[i,6]

}


results <- cbind(scenarios,
                 bias.ltr0, bias.ltr1, bias.ltrdiff,
                 relbias.ltr0, relbias.ltr1, relbias.ltrdiff,
                 rmse.ltr0, rmse.ltr1, rmse.ltrdiff,
                 cov.ltr0, cov.ltr1, cov.ltrdiff,
                 effsampsize, power
                 )

results.cols <- c("bias.ltr0.pseudo", "bias.ltr0.fg", "bias.ltr0.fglogt", "bias.ltr0.flex", "bias.ltr0.flex2df", "bias.ltr0.flex3df",
                  "bias.ltr1.pseudo", "bias.ltr1.fg", "bias.ltr1.fglogt", "bias.ltr1.flex", "bias.ltr1.flex2df", "bias.ltr1.flex3df",
                  "bias.ltrdiff.pseudo", "bias.ltrdiff.fg", "bias.ltrdiff.fglogt", "bias.ltrdiff.flex", "bias.ltrdiff.flex2df", "bias.ltrdiff.flex3df",
                  "relbias.ltr0.pseudo", "relbias.ltr0.fg", "relbias.ltr0.fglogt", "relbias.ltr0.flex", "relbias.ltr0.flex2df", "relbias.ltr0.flex3df",
                  "relbias.ltr1.pseudo", "relbias.ltr1.fg", "relbias.ltr1.fglogt", "relbias.ltr1.flex", "relbias.ltr1.flex2df", "relbias.ltr1.flex3df",
                  "relbias.ltrdiff.pseudo", "relbias.ltrdiff.fg", "relbias.ltrdiff.fglogt", "relbias.ltrdiff.flex", "relbias.ltrdiff.flex2df", "relbias.ltrdiff.flex3df",
                  "rmse.ltr0.pseudo", "rmse.ltr0.fg", "rmse.ltr0.fglogt", "rmse.ltr0.flex",  "rmse.ltr0.flex2df",  "rmse.ltr0.flex3df",
                  "rmse.ltr1.pseudo", "rmse.ltr1.fg", "rmse.ltr1.fglogt", "rmse.ltr1.flex", "rmse.ltr1.flex2df", "rmse.ltr1.flex3df",
                  "rmse.ltrdiff.pseudo", "rmse.ltrdiff.fg", "rmse.ltrdiff.fglogt", "rmse.ltrdiff.flex", "rmse.ltrdiff.flex2df", "rmse.ltrdiff.flex3df",
                  "cov.ltr0.pseudo", "cov.ltr0.flex", "cov.ltr0.flex2df", "cov.ltr0.flex3df",
                  "cov.ltr1.pseudo", "cov.ltr1.flex", "cov.ltr1.flex2df", "cov.ltr1.flex3df",
                  "cov.ltrdiff.pseudo", "cov.ltrdiff.flex", "cov.ltrdiff.flex2df", "cov.ltrdiff.flex3df",
                  "effsampsize.pseudo", "effsampsize.fg", "effsampsize.fglogt", "effsampsize.flex", "effsampsize.flex2df", "effsampsize.flex3df",
                  "power.pseudo")

colnames(results) <- c(colnames(scenarios), results.cols)

write.csv(results, "results_revision.csv", row.names=FALSE)

# Try reordering
# results2 <- results[order(results$setting, results$probcens, results$sampsize, results$ptrunc), ]
# results2 <- results2[, c("setting", "probcens", "sampsize", "ptrunc", results.cols)]
# write.csv(results2, "Y:\\19SCR-lifetimerisk_varselection\\simulation\\sim results\\results order 2.csv", row.names=FALSE)
# 



