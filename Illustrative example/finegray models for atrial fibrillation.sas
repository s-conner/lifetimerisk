/*******************************************************************************************
Code for illustrative example analyses in the Framingham Heart Study

Fits Fine-Gray models for AF and death without AF, and predicts
individuals'  lifetime risk at age 95 (for assessing calibration).
******************************************************************************************/


proc import datafile="Y:\\19SCR-lifetimerisk_varselection\\FHS example\\af_withpseudodeath_v2.csv" out=fhs dbms=csv replace; run;
libname fhs "Y:\\19SCR-lifetimerisk_varselection\\FHS example\\";

proc means data=fhs; var startage_55 sbp_55 dbp_55 bmi_55; run;

****** AF ******;
ods output ParameterEstimates=Param;
proc phreg data=fhs;
title 'AF, PH all covariates';
model (startage_55, endage_55) * event_55(0) = startage_55 male sbp_55 dbp_55 hrx_55 currsmk_55 alco_elevated bmi_55 hx_diab_55 pchf_55 pmi_55/ eventcode=1;
baseline covariates=fhs out=predout cif=cif/ rowid=fid;
hazardratio startage_55 / units=1.0335572;
hazardratio sbp_55 / units=17.08932;
hazardratio dbp_55 / units=9.848818;
hazardratio bmi_55 / units=5.293189;
run;

data fhs_predictedltr; set predout; by fid; if last.fid; run;
proc rank data=fhs_predictedltr groups=10 out=fhs.fhs_predictedltr; var cif; ranks decile; run;

data param; set param;
HR=exp(estimate);
CIL=exp(estimate-1.96*stderr);
CIU=exp(estimate+1.96*stderr);
estci=cat(round(HR,.01), " (", round(CIL,.01), ", ", round(CIU,.01), ")");
run;
proc print; var parameter estci probchisq; title 'AF model'; run;


****** Death without AF ******;
ods output ParameterEstimates=Param2;
proc phreg data=fhs;
title 'Death without AF, PH all covariates';
model (startage_55, endage_55) * event_55(0) = startage_55 male sbp_55 dbp_55 hrx_55 currsmk_55 alco_elevated bmi_55 hx_diab_55 pchf_55  pmi_55/ eventcode=2;
hazardratio startage_55 / units=1.0335572;
hazardratio sbp_55 / units=17.08932;
hazardratio dbp_55 / units=9.848818;
hazardratio bmi_55 / units=5.293189;
run;

data param2; set param2;
HR=exp(estimate);
CIL=exp(estimate-1.96*stderr);
CIU=exp(estimate+1.96*stderr);
estci=cat(round(HR,.01), " (", round(CIL,.01), ", ", round(CIU,.01), ")");
run;

proc print; var parameter estci probchisq; title 'Death without AF model';  run;
