# Compile Fine-Gray bootstrap CIs for AF
# 8/3/2021


# Bootstrap CIs

ltrdiff <- matrix(NA, nrow=1000, ncol=11)

for(i in 1:20){
  start <- (50*i)-49
  end <- 50*i
  ltrdiff[start:end, ]  <- as.matrix(read.csv(paste0('bootci_af_revision//fg_ltrdiff_ci_', i, '.csv')))[start:end, ]
}

ltr.diff.pil <- apply(ltrdiff, 2, function(x) quantile(x, .025))
ltr.diff.piu <- apply(ltrdiff, 2, function(x) quantile(x, .975))

cis <- cbind(ltr.diff.pil, ltr.diff.piu)


# Merge with predictions

ests <- read.csv('results//finegray_af_diffltr_revision.csv')
res <- data.frame(cbind(ests$var, ests$ltrdiff, cis))
colnames(res) <- c('var', 'ltrdiff', 'cil', 'ciu')
res$var <- factor(res$var, levels = c('male', 'entry_hgt', 'entry_wgt', 'entry_sbp', 'entry_dbp', 'entry_hrx', 
                                      'entry_smk', 'entry_alc_elev',  'entry_diab', 'entry_pchf', 'entry_pmi'),
                  labels= c('Male', 'Height, SDU', 'Weight, SDU', 'Systolic blood pressure, SDU', 'Diastolic blood pressure, SDU',
                            'Anti-hypertensive use', 'Current smoker', 'Elevated alcohol use', 'Prior diabetes', 
                            'Prior heart failure', 'Prior myocardial infarction'))
res2 <- res[order(res$var), ]


# Export

write.csv(res2, "results//finegray_af_diffltr_ci_revision.csv", row.names=FALSE)



