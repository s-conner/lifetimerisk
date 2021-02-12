# Compile Fine-Gray bootstrap CIs for AF
# 1/20/2021


# Bootstrap CIs

ltrdiff <- matrix(NA, nrow=1000, ncol=11)

for(i in 1:20){
  start <- (50*i)-49
  end <- 50*i
  ltrdiff[start:end, ]  <- as.matrix(read.csv(paste0('bootci_af//fg_ltrdiff_ci_', i, '.csv')))[start:end, ]
}

ltr.diff.pil <- apply(ltrdiff, 2, function(x) quantile(x, .025))
ltr.diff.piu <- apply(ltrdiff, 2, function(x) quantile(x, .975))

cis <- cbind(ltr.diff.pil, ltr.diff.piu)


# Merge with predictions

ests <- read.csv('results//finegray_af_diffltr.csv')
res <- data.frame(cbind(ests$var, ests$ltrdiff, cis))
colnames(res) <- c('var', 'ltrdiff', 'cil', 'ciu')
res$var <- factor(res$var, levels = c('male', 'entry_hgt', 'entry_wgt', 'entry_sbp', 'entry_dbp', 'entry_smk', 
                                      'entry_alc_elev', 'entry_hrx', 'entry_diab', 'entry_pchf', 'entry_pmi'),
                  labels= c('Male', 'Height, SDU', 'Weight, SDU', 'Systolic blood pressure, SDU', 'Diastolic blood pressure, SDU',
                            'Current smoker', 'Elevated alcohol use', 'Anti-hypertensive use', 'Prior diabetes', 
                            'Prior heart failure', 'Prior myocardial infarction'))
res[order(res$var), ]


# Export

write.csv(cis, "results//finegray_af_diffltr_ci.csv", row.names=FALSE)



