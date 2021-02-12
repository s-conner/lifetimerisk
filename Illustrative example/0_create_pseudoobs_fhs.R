# Derive pseudo-observations
# Sarah Conner - Nov 18 2020

library(sas7bdat)
library(Rcpp)
sourceCpp('0_cpp.cpp')

# Subset data to non-missing covariates
af0 <- read.sas7bdat('allcohorts_55.sas7bdat')
af0$male <- 2-af0$SEX
covs <- c('male', 'entryage_55', 'entry_sbp', 'entry_dbp', 'entry_smk', 'entry_alc', 'entry_hgt', 'entry_wgt', 'entry_hrx', 'entry_diab', 'entry_pchf', 'entry_pmi')
af <- na.omit(af0[, c('fid', 'exitage', 'event', covs)])

# Pseudo-observation function
pseudo.ltr <- function(dataset, entry, times, cause, eoi){
  data <- data.frame(entry=dataset[[entry]], times=dataset[[times]], cause=dataset[[cause]])
  jk.lr <- rep(NA, nrow(data))
  pseudo <- rep(NA, nrow(data))
  n <- nrow(data)
  data$jid <- 1:n
  
  event.ages <- data$times[data$cause!=0]
  event.ages <- unique(event.ages[order(event.ages)])
  overall.lr <- cif(start=data$entry, stop=data$times, event=data$cause, eventAges=event.ages, failcode=eoi)
  
  # Jackknife, my function
  for(i in 1:n){
    samp <- data[data$jid != i, ]
    event.ages <- samp$times[samp$cause!=0]
    event.ages <- unique(event.ages[order(event.ages)])
    jk.lr[i] <- cif(start=samp$entry, stop=samp$times, event=samp$cause, eventAges=event.ages, failcode=eoi)
    pseudo[i] <- n*overall.lr - (n-1)*jk.lr[i]
  }
  return(pseudo)
}

# Create pseudo-observations for AF and Death
af$pseudo.af <- pseudo.ltr(af, 'entryage_55', 'exitage', 'event', 1)
af$pseudo.death <- pseudo.ltr(af, 'entryage_55', 'exitage', 'event', 2)

# Export for future use
write.csv(af, "af_pseudo.csv", row.names=FALSE)


