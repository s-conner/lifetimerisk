local taskID : env SGE_TASK_ID

cd "//restricted//projectnb//conner-thesis//lifetimerisk//flex//nph_csh//"

tempfile simres
postfile simres iter ltr1 ltr1_lci ltr1_uci ltr0 ltr0_lci ltr0_uci ltrdiff ltrdiff_lci ltrdiff_uci using loopres`taskID', replace

forvalues i=1/2000 {
  import delimited "//restricted//projectnb//conner-thesis//lifetimerisk//sim1nov4//datacsh//data`taskID'_`i'.csv", clear
  
  stset times, failure(cause==1, 2) id(id) enter(entry) exit(time 40)
  gen cause2 = cond(_d==0,0,cause)
  
  stcrprep, events(cause2) keep(arm) 
  
  stset tstop [iw=weight_c], failure(cause==1) enter(tstart) 
  capture noisily stpm2 arm if failcode == 1, scale(hazard) df(3) tvc(arm) dftvc(1) failconvlininit iterate(100)
  
  if e(converged)==1 { 
    gen tau=40
    gen iter=`i'
    predict ltr1, failure at(arm 1) timevar(tau) ci 
    predict ltr0, failure at(arm 0) timevar(tau) ci 
    predict ltrdiff, sdiff1(arm 1) sdiff2(arm 0) timevar(tau) ci
    keep if _n==1
    keep iter ltr1 ltr1_lci ltr1_uci ltr0 ltr0_lci ltr0_uci ltrdiff ltrdiff_lci ltrdiff_uci
    
    post simres (iter) (ltr1) (ltr1_lci) (ltr1_uci) (ltr0) (ltr0_lci) (ltr0_uci) (ltrdiff) (ltrdiff_lci) (ltrdiff_uci) 
    } 
}

postclose simres
use "//restricted//projectnb//conner-thesis//lifetimerisk//flex//nph_csh//loopres`taskID'.dta", clear
outsheet using "//restricted//projectnb//conner-thesis//lifetimerisk//flex//nph_csh//loopres`taskID'.csv", comma




