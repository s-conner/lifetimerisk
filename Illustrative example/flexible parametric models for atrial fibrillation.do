cd "//restricted//projectnb//conner-thesis//lifetimerisk//fhs//"



/* AF model*/ 
import delimited "//restricted//projectnb//conner-thesis//lifetimerisk//fhs//af_withpseudodeath_v2.csv", clear
  


stset endage_55, failure(event_55==1, 2) id(fid) enter(startage_55) exit(time 95)

gen cause2 = cond(_d==0,0,event_55)
  

stcrprep, events(cause2) keep(startage_55 male sbp_55 dbp_55 hrx_55 currsmk_55 alco_elevated bmi_55 hx_diab_55 pchf_55 pmi_55) 

stset tstop [iw=weight_c], failure(cause==1) enter(tstart) 




stpm2 startage_55 male sbp_55 dbp_55 hrx_55 currsmk_55 alco_elevated bmi_55 hx_diab_55 pchf_55 pmi_55 if failcode == 1, scale(hazard) df(4) tvc(male currsmk_55 hrx_55) dftvc(1) failconvlininit iterate(100)



gen tau=95

predict ltr, failure timevar(tau) ci 


outsheet using "//restricted//projectnb//conner-thesis//lifetimerisk//fhs//rp_predictedltr.csv", comma


/* Death without AF model */
import delimited "//restricted//projectnb//conner-thesis//lifetimerisk//fhs//af_withpseudodeath_v2.csv", clear
  


stset endage_55, failure(event_55==1, 2) id(fid) enter(startage_55) exit(time 95)
gen cause2 = cond(_d==0,0,event_55)
stcrprep, events(cause2) keep(startage_55 male sbp_55 dbp_55 hrx_55 currsmk_55 alco_elevated bmi_55 hx_diab_55 pchf_55 pmi_55) 
stset tstop [iw=weight_c], failure(cause==2) enter(tstart) 

stpm2 startage_55 male sbp_55 dbp_55 hrx_55 currsmk_55 alco_elevated bmi_55 hx_diab_55 pchf_55 pmi_55 if failcode == 2, scale(hazard) df(3) tvc(male sbp_55 hrx_55 hx_diab_55) dftvc(1) failconvlininit 