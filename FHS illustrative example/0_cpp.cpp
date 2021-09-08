#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
NumericVector cif(NumericVector start, NumericVector stop, NumericVector event, NumericVector eventAges, int failcode) {
  
  int n = eventAges.size();
  NumericVector atrisk(n);
  NumericVector eventsall(n);
  NumericVector eoi(n);
  NumericVector hazall(n);
  NumericVector hazeoi(n);
  NumericVector surv(n+1);
  NumericVector ci(n+1);
  NumericVector cif(n);
  double x;
  double y;
  
  // numbers at risk
  for(int i = 0; i < n; ++i){
    atrisk[i]=0;
    for(int j = 0; j < stop.size(); ++j){
      atrisk[i] = atrisk[i] + int(start[j]<eventAges[i] && stop[j]>=eventAges[i]);
    }
  }
  
  // any events of all causes
  for(int i = 0; i < n; ++i){
    eventsall[i]=0;
    for(int j = 0; j < stop.size(); ++j){
      eventsall[i] = eventsall[i] + int(event[j]!=0 && stop[j]==eventAges[i]) ;
    }
  }
  
  // events of interest
  for(int i = 0; i < n; ++i){
    eoi[i]=0;
    for(int j = 0; j < stop.size(); ++j){
       eoi[i] = eoi[i] + int(event[j]==failcode && stop[j]==eventAges[i]);
    }
  }
  
  // hazards (all cause and event of interest)
  hazall = eventsall/atrisk;
  hazeoi = eoi/atrisk;
  
  // All-cause survival
  surv[0]=1;
  for(int i = 1; i < n+1; ++i){
    x=0.;
    for(int j = 0; j < i; ++j){
      x = x + hazall[j] * surv[j];
    }
    surv[i] = 1. - x;
  }
  
  // CIF
  ci[0]=0;
  for(int i = 1; i < n+1; ++i){
    y=0.;
    for(int j = 0; j < i; ++j){
      y = y + hazeoi[j] * surv[j];
    }
    ci[i] = y;
  }
  
  for(int i = 0; i < n; ++i){
    cif[i] = ci[i+1];
  }
  
  return tail(cif, 1);
  
}

