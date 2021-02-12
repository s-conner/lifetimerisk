# Modeling the lifetime risk 

This project contains R code to model the lifetime risk using pseudo-observations, and to predict lifetime risk with Fine-Gray and flexible parametric models of the subdistribtuion hazard with time-varying effects. We also share code to replicate our simulation study and illustrative example in the Framingham Heart Study (without data).

## Pseudo-observation lifetime risk R functions
pseudo_ltr_function.R
```
This R code derives the lifetime risk at a given timepoint and derives pseudo-observations of the lifetime 
risk with the jackknife approach. Requires the Rcpp package to run pseudo_ltr_cpp.cpp code.

One could alternatively use the ETM package to derive the lifetime risk, and derive pseudo-observations 
from this quantity - results will be identical if the ETM package is used correctly.
```

pseudo_ltr_cpp.cpp
```
C++ functions used to derive the lifetime risk. We use Rcpp to speed up calculations. 
```

## Simulations
generatedata_sim_csh.R
```
R code used to generate cause-specific hazard simulation settings
```

generatedata_sim_sdh.R
```
R code used to generate subdistribution hazard simulation settings
```

pseudo_sim.R
```
R code for modeling the pseudo-observations of the lifetime risk
```



## Illustrative example

pseudo models for atrial fibrillation.R
```
R code for modeling the pseudo-observations of the lifetime risk of atrial fibrillation and death without
atrial fibrillation in the Framingham Heart Study
```

finegray models for atrial fibrillation.sas
```
R code for fitting Fine-Gray models for the cumulative incidence functions of atrial fibrillation and 
death without atrial fibrillation in the Framingham Heart Study
```

flexible parametric models for atrial fibrillation.do
```
R code for fitting fleixble parametric models for the cumulative incidence functions of atrial  
fibrillation and death without atrial fibrillation in the Framingham Heart Study
```
