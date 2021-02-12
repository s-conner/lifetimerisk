# Modeling the lifetime risk 

This project contains R code to model the lifetime risk using pseudo-observations, and to predict lifetime risk with Fine-Gray and flexible parametric models of the subdistribtuion hazard with time-varying effects. We also share code to replicate our simulation study and illustrative example in the Framingham Heart Study (without data).

## Pseudo-observation lifetime risk R functions
pseudo_ltr.R
```
This R code derives the lifetime risk at a given timepoint and derives pseudo-observations of the lifetime 
risk with the jackknife approach. Requires the Rcpp package to run pseudo_cpp.cpp code.

One could alternatively use the ETM package to derive the lifetime risk, and code their own jackknife and pseudo-observations 
from this quantity - pseudo-observations will be identical.
```

pseudo_cpp.cpp
```
C++ functions used to derive the lifetime risk. We use Rcpp to speed up calculations. 
```

## Simulations
0_generatedata_sim_csh.R and 0_generatedata_sim_sdh.R
```
R code used to generate cause-specific hazard data and subdistribution hazard data
```

1_pseudo_csh_55.R and 1_pseudo_sdh_55.R
```
R code used to analyze data with pseudo-observation approach
```

2_fg_csh_55.R and 2_fg_sdh_55.R
```
R code used to fit Fine-Gray model assuming proportional subdistribution hazards to IPCLW data and predict lifetime risk
```

3_fg_logt_csh_55.R and 3_fg_logt_sdh_55.R
```
R code used to fit Fine-Gray model with time-covariate interactions to allow non-proportional subdistribution hazards with IPCLW data, and predict lifetime risk
```

4_flex_csh_55.R and 4_flex_sdh_55.R
```
R code used to fit flexible parametric models to IPCLW data and predict lifetime risk
```

5_compile_results.R, 6_nestedloops_figs.R, 6_nestedloop_fxn.R
```
Compile all simulation results and present as nested loop plots (Rucker et al)
```





