
%let path=%str(/restricted/projectnb/conner-thesis/lifetimerisk/sim1nov4/);
%let i=%scan(&sysparm,1);

%macro runfg;
/* can split this into chunks to increase run time */
%do j=1 %to 2000;
proc import datafile="&path.datacsh/data&i._&j..csv" out=data dbms=csv replace; getnames=yes; run;

data data; set data; if times>40 then do; times=40; cause=0; end; run;

/* Fit Fine-Gray Model, export linear predictor for each arm*/
data pred; arm=0; output; arm=1; output; run;

ods exclude all;
ods output ParameterEstimates=model(keep=Estimate StdErr ProbChiSq);
proc phreg data=data;
model (entry,times)*cause(0) = arm / eventcode=1;
baseline covariates=pred out=predout cif=_all_/ rowid=arm normalsample=1000;
run;

/* Extract LTR and 95% CI per group*/
proc sort data=predout; by arm times; run;
data ltrlong0 ltrlong1; set predout; by arm; if last.arm and arm=0 then output ltrlong0; else if last.arm and arm=1 then output ltrlong1; run;

proc sql; create table ltrwide as 
select a.cif as ltr0, a.stderrcif as se0, a.lowercif as ltr_cil0, a.uppercif as ltr_ciu0, b.cif as ltr1, b.stderrcif as se1, b.lowercif as ltr_cil1, b.uppercif as ltr_ciu1, b.cif-a.cif as ltr_diff
from ltrlong0 a left join ltrlong1 b on 1=1;
quit;


/* Bootstrap CI for difference in LTR */
proc sort data=data; by arm; run;
proc surveyselect data=data out=bootsample noprint method=urs samprate=1 reps=1000; 
strata arm;  
run;
proc sort data=bootsample; by replicate; run;

proc phreg data=bootsample noprint;
by replicate;
freq numberhits;
model (entry,times)*cause(0) = arm / eventcode=1;
baseline covariates=pred out=bootout cif=cif/ rowid=arm ;
run;

proc sort data=bootout; by replicate arm times; run;

data bootlong0 bootlong1; 
set bootout; 
by replicate arm; 
if last.arm and arm=0 then output bootlong0; 
else if last.arm and arm=1 then output bootlong1; 
run;

proc sql; create table bootlong as select a.replicate, a.cif as ltr0, b.cif as ltr1, b.cif-a.cif as ltr_diff 
from bootlong0 a 
left join bootlong1 b on a.replicate=b.replicate; 
quit;
proc means data=bootlong noprint; var ltr_diff; output out=bootse std=ltr_diff_se; run;

data bootse; set bootse; drop _TYPE_ _FREQ_; run;
data ltrall; merge ltrwide bootse; run;
data ltrall; set ltrall; cil=ltr_diff - 1.96*ltr_diff_se; ciu=ltr_diff + 1.96*ltr_diff_se;  run;

/* Export model coefficients and LTR */
proc append base=all_model_res data=model; run;
proc append base=all_ltr_res data=ltrall; run;

%end;
proc export data=all_model_res dbms=csv outfile="&path.fgresa/model&i..csv" replace; run;
proc export data=all_ltr_res dbms=csv outfile="&path.fgresa/ltr&i..csv" replace; run;
%mend;

%runfg;
