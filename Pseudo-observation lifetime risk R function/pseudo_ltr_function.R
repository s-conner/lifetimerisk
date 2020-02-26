#############################################################################################
# Function to derive the pseudo-observations of the lifetime risk
# Requires 'Rcpp' package and 'pseudo_ltr_cpp' file of Rcpp functions 
# Function arguments include:
  # dataset: dataset name
  # entry: entry time
  # times: exit time
  # cause: status at exit time (0 for censor, 1+ for events)
  # eoi: value of cause variable which denotes event of interest. 
  #      i.e., 1 for the disease or 2 for semi-competing risk of death
  # tau: age to calculate lifetime risk. default is the last event time for event of interest

# Please note that you can alternatively calculate the lifetime risk using the ETM package
#############################################################################################



library(Rcpp)
sourceCpp('//restricted//projectnb//conner-thesis//lifetimerisk//sim1nov4//pseudo_ltr_cpp.cpp')


pseudo.lr <- function(dataset, entry, times, cause, eoi, tau=NULL){
  data <- data.frame(entry=dataset[[entry]], times=dataset[[times]], cause=dataset[[cause]])
  jk.lr <- rep(NA, nrow(data))
  pseudo <- rep(NA, nrow(data))
  n <- nrow(data)
  data$jid <- 1:n
  
  if(max(data$times)<tau){print('Event times do not reach tau.')
  }else{
    
    if(is.null(tau)){
      tau <- max(data$times[data$cause==eoi])
    }
    
    event.ages <- data$times[data$cause!=0]
    event.ages <- unique(event.ages[order(event.ages)])
    
    # Obtain lifetime risk in overall sample
    overall.cif <- cif2(start=data$entry, stop=data$times, event=data$cause, eventAges=event.ages, failcode=eoi)
    overall.lr <- overall.cif[which(event.ages==max(event.ages[event.ages<=tau]))]
    
    # Jackknife
    for(i in 1:n){
      samp <- data[data$jid != i, ]
      event.ages <- samp$times[samp$cause!=0]
      event.ages <- unique(event.ages[order(event.ages)])
      
      jk.cif <- cif2(start=samp$entry, stop=samp$times, event=samp$cause, eventAges=event.ages, failcode=eoi)
      jk.lr[i] <- jk.cif[which(event.ages==max(event.ages[event.ages<=tau]))]
      pseudo[i] <- n*overall.lr - (n-1)*jk.lr[i]
    }
    return(pseudo)
  }
}