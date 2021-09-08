# Derive pseudo-observations
# Sarah Conner - Nov 18 2020

library(sas7bdat)
library(Rcpp)
sourceCpp('0_cpp.cpp')


# Pseudo-observation function
pseudo.ltr <- function(dataset, entry, times, cause, eoi=NULL, tau=NULL){
  
  data <- data.frame(entry=dataset[[entry]], times=dataset[[times]], cause=dataset[[cause]])

  # Set tau
  if(length(data$times[data$times>=tau]) <=1){
    print('Please choose another tau. Event times do not reach tau in all jackknife samples.')
  }else{
    if(is.null(tau)){
      tau <- max(data$times[data$cause==eoi])
    }
  }

  # Censor any time after tau & drop any entry times after event times
  data$cause[data$times>tau] <- 0
  data$times[data$times>tau] <- tau
  data <- data[data$entry<data$times, ]
  n <- nrow(data)
  data$jid <- 1:n
  
  # Prepare pseudo-obs for 1 event, or all events
  if(is.null(eoi)){
    jk.lr <- matrix(NA, nrow=n, ncol=max(data$cause))
    pseudo <- matrix(NA, nrow=n, ncol=max(data$cause))
    eoi <- c(1:max(data$cause))
  }else{
    jk.lr <- matrix(NA, nrow=n, ncol=1)
    pseudo <- matrix(NA, nrow=n, ncol=1)
    }
    
  tj <- data$times[data$cause!=0]
  tj <- unique(tj[order(tj)])
  overall.lr <- rep(NA, length(eoi))
    
    # Overall LTRs
  for(k in eoi){
    overall.lr[k] <- cif(start=data$entry, stop=data$times, event=data$cause, eventAges=tj, failcode=k)
  }
    
    # Jackknife & derive pseudo-observations
  for(i in 1:n){
      
    samp <- data[data$jid != i, ]
    tj <- samp$times[samp$cause!=0]
    tj <- unique(tj[order(tj)])
      
    for(k in eoi){
      jk.lr[i,k] <- cif(start=samp$entry, stop=samp$times, event=samp$cause, eventAges=tj, failcode=k)
      pseudo[i,k] <- n*overall.lr[k] - (n-1)*jk.lr[i,k]
      }
    }

  return(pseudo)
}

# Create pseudo-observations for AF and Death

af0 <- read.sas7bdat('allcohorts_55.sas7bdat')
af0$male <- 2-af0$SEX
covs <- c('male', 'entryage_55', 'entry_sbp', 'entry_dbp', 'entry_smk', 'entry_alc', 'entry_hgt', 'entry_wgt', 'entry_hrx', 'entry_diab', 'entry_pchf', 'entry_pmi')
af <- na.omit(af0[, c('fid', 'exitage', 'event', covs)])

# pseudo1.95 <- pseudo.ltr(af, 'entryage_55', 'exitage', 'event', eoi=1, tau=95)
# pseudoboth.95 <- pseudo.ltr(af, 'entryage_55', 'exitage', 'event', tau=95)
# pseudoboth.90 <- pseudo.ltr(af, 'entryage_55', 'exitage', 'event', tau=90)

start_time <- Sys.time()
af$pseudo.af <- pseudo.ltr(af, 'entryage_55', 'exitage', 'event', 1, tau=95)
end_time <- Sys.time()
end_time - start_time


af$pseudo.death <- pseudo.ltr(af, 'entryage_55', 'exitage', 'event', 2)

# Export for future use
write.csv(af, "af_pseudo.csv", row.names=FALSE)


