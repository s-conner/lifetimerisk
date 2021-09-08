### Flexible parametric models for FHS - Lifetime risk of AF
### Choose # of knots and include interactions with log(time), determine using AIC
### Updated analyses for revision 7/15/2021

library(survival)
library(mstate)
library(rstpm2) 
options(scipen=999)

# ----- Data prep -----

af <- read.csv('af_pseudo.csv')
vars <- c('male', 'entry_hgt', 'entry_wgt', 'entry_sbp', 'entry_dbp', 'entry_hrx', 
          'entry_smk', 'entry_alc_elev', 'entry_diab', 'entry_pchf', 'entry_pmi')
means <- apply(af[, vars], 2, mean)
sds <- apply(af[, vars], 2, sd)

# Create LTRC weighted dataset
af.fg <- crprep(Tstop="exitage", status="event", data=af, trans=1, cens=0, Tstart="entryage_55", id="fid", keep=vars, shorten=FALSE)
af.fg$weight_ltrc <- af.fg$weight.cens*af.fg$weight.trunc
af.fg$ev1 <- ifelse(af.fg$status==1,1,0)


# ----- Determine # knots in FP model -----

fp.1 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
              entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
            weight=weight_ltrc, data=af.fg, df=1)
fp.2 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
              entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
            weight=weight_ltrc, data=af.fg, df=2)
fp.3 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
              entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
            weight=weight_ltrc, data=af.fg, df=3)
fp.4 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
              entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
            weight=weight_ltrc, data=af.fg, df=4)
fp.5 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
              entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
            weight=weight_ltrc, data=af.fg, df=5)

AIC(fp.1)
AIC(fp.2)
AIC(fp.3)
AIC(fp.4) # AIC 9732.231
AIC(fp.5) # lowest AIC 9731.536

BIC(fp.1)
BIC(fp.2)
BIC(fp.3)
BIC(fp.4) 
BIC(fp.5)


AIC(fp.1) - AIC(fp.2)
AIC(fp.2) - AIC(fp.3)
AIC(fp.3) - AIC(fp.4)
AIC(fp.4) - AIC(fp.5) 

full.aic0 <- AIC(fp.4)


# ----- Round 1: 4 knots, add 1 log(time) interaction at a time ----- 

fp.int1 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(male=1))

fp.int2 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_hgt=1))

fp.int3 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1))

fp.int4 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_sbp=1))

fp.int5 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_dbp=1))

fp.int6 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_hrx=1))

fp.int7 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_smk=1))

fp.int8 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_alc_elev=1))

fp.int9 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_diab=1))

fp.int10 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                    entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                  weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_pchf=1))

fp.int11 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                    entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                  weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_pmi=1))

round1.aic <- data.frame(var=1:11, aic=c(AIC(fp.int1)
                            , AIC(fp.int2)
                            , AIC(fp.int3) 
                            , AIC(fp.int4)
                            , AIC(fp.int5)
                            , AIC(fp.int6)
                            , AIC(fp.int7)
                            , AIC(fp.int8)
                            , AIC(fp.int9)
                            , AIC(fp.int10)
                            , AIC(fp.int11)
))

round1.aic$delta <- full.aic0 - round1.aic$aic
round1.aic[order(round1.aic[,2]),]
full.aic1 <- min(round1.aic$aic)

# Lowest AIC 9690.055
summary(fp.int3) 



# ----- Round 2: 4 knots, entry_wgt + add 1 log(time) interaction at a time ----- 

fp.2.int1 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, male=1))

fp.2.int2 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_hgt=1))

fp.2.int4 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_sbp=1))

fp.2.int5 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_dbp=1))

fp.2.int6 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_hrx=1))

fp.2.int7 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_smk=1))

fp.2.int8 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_alc_elev=1))

fp.2.int9 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1))

fp.2.int10 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                    entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                  weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_pchf=1))

fp.2.int11 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                    entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                  weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_pmi=1))

round2.aic <- data.frame(var=vars[c(1,2,4:11)], aic=c(AIC(fp.2.int1)
                            , AIC(fp.2.int2)
                            , AIC(fp.2.int4)
                            , AIC(fp.2.int5)
                            , AIC(fp.2.int6)
                            , AIC(fp.2.int7)
                            , AIC(fp.2.int8)
                            , AIC(fp.2.int9)
                            , AIC(fp.2.int10)
                            , AIC(fp.2.int11)
))

round2.aic$delta <- full.aic1 - round2.aic$aic
round2.aic[order(round2.aic[,2]),]
sort2 <- round2.aic[order(round2.aic[,2]),]
sort2[sort2$aic>0,]
full.aic2 <- round2.aic$aic[8]

# Delta 29.5 for diab


# ----- Round 3: 4 knots, entry_wgt + entry_diab=1, add 1 log(time) interaction at a time ----- 

fp.3.int1 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, male=1))

fp.3.int2 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))

fp.3.int6 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hrx=1))

fp.3.int8 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_alc_elev=1))

fp.3.int10 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                      entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                    weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_pchf=1))

fp.3.int11 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                      entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                    weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_pmi=1))

round3.aic <- data.frame(var=vars[c(1,2,6,8,10,11)], 
                    aic=c(AIC(fp.3.int1)
                      , AIC(fp.3.int2)
                      , AIC(fp.3.int6)
                      , AIC(fp.3.int8)
                      , AIC(fp.3.int10)
                      , AIC(fp.3.int11)
))


round3.aic$delta <- full.aic2 - round3.aic$aic
sort3 <- round3.aic[order(round3.aic[,2]),]
sort3[sort3$aic>0,]
full.aic3 <- round3.aic$aic[2]

# Delta 21.6 for height


# ----- Round 4: 4 knots, entry_wgt + entry_diab + entry_hgt, add 1 log(time) interaction at a time ----- 


fp.4.int1 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1, male=1))

fp.4.int6 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1, entry_hrx=1))

fp.4.int8 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1, entry_alc_elev=1))

fp.4.int10 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                      entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                    weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1, entry_pchf=1))

fp.4.int11 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                      entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                    weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1, entry_pmi=1))

round4.aic <- data.frame(var=vars[c(1,6,8,10,11)], 
                    aic=c(AIC(fp.4.int1)
                      , AIC(fp.4.int6)
                      , AIC(fp.4.int8)
                      , AIC(fp.4.int10)
                      , AIC(fp.4.int11)
                    ))

round4.aic$delta <- full.aic3 - round4.aic$aic
sort4 <- round4.aic[order(round4.aic[,2]),]
sort4[sort4$aic>0,]
full.aic4 <- round4.aic$aic[2]
  
# delta 3.3, STOP at weight, diab, height



# ----- Round 1 of adding additional df for time-varying terms -----

fp.1b.1 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                      entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                    weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=2, entry_diab=1, entry_hgt=1))
fp.1b.2 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=2, entry_hgt=1))
fp.1b.3 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=2))

round1b.aic <- data.frame(aic=c(AIC(fp.1b.1), AIC(fp.1b.2), AIC(fp.1b.3)))

round1b.aic$delta <- full.aic3 - round1b.aic$aic
sort1b <- round1b.aic[order(round1b.aic[1]),]
sort1b[sort1b$aic>0,]

# None improve the model fit



# ----- Round 1 of adding sex interaction terms -----

fp.1c.1 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_hgt, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.1c.2 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_wgt,  
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.1c.3 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_sbp,  
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.1c.4 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_dbp,  
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.1c.5 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_hrx,  
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.1c.6 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_smk,  
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.1c.7 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_alc_elev,  
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.1c.8 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_diab,  
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.1c.9 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_pchf,  
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.1c.10 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_pmi,  
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))


round1c.aic <- data.frame(vars=vars[c(2:11)], 
                          aic=c(AIC(fp.1c.1), AIC(fp.1c.2), AIC(fp.1c.3),
                                AIC(fp.1c.4), AIC(fp.1c.5), AIC(fp.1c.6),
                                AIC(fp.1c.7), AIC(fp.1c.8), AIC(fp.1c.9), AIC(fp.1c.10)))

round1c.aic$delta <- full.aic3 - round1c.aic$aic
sort1c <- round1c.aic[order(round1c.aic[2]),]
sort1c[sort1c$aic>0,]

summary(fp.1c.5) # did not converge
summary(fp.1c.3) # delta 216 seems huge

full.aic1c <- AIC(fp.1c.3)


# ----- Round 2 of adding sex interaction terms -----

fp.2c.1 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_sbp + male:entry_hgt, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.2c.2 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_sbp + male:entry_wgt,  
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.2c.4 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_sbp + male:entry_dbp,  
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.2c.5 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_sbp + male:entry_hrx,  
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.2c.6 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_sbp + male:entry_smk,  
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.2c.7 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_sbp + male:entry_alc_elev,  
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.2c.8 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_sbp + male:entry_diab,  
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.2c.9 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_sbp + male:entry_pchf,  
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.2c.10 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                    entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_sbp + male:entry_pmi,  
                  weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))

round2c.aic <- data.frame(vars=vars[c(2,3,5:11)], 
                          aic=c(AIC(fp.2c.1), AIC(fp.2c.2), 
                                AIC(fp.2c.4), AIC(fp.2c.5), AIC(fp.2c.6),
                                AIC(fp.2c.7), AIC(fp.2c.8), AIC(fp.2c.9), AIC(fp.2c.10)))

round2c.aic$delta <- full.aic1c - round2c.aic$aic
sort2c <- round2c.aic[order(round2c.aic[2]),]
sort2c[sort2c$aic>0,]

summary(fp.2c.9) # Converged and looks okay, delta is very high 1528.9954
full.aic2c <- AIC(fp.2c.9)




# ----- Round 3 of adding sex interaction terms -----

fp.3c.1 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_sbp + male:entry_pchf + male:entry_hgt, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.3c.2 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_sbp + male:entry_pchf + male:entry_wgt,  
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.3c.4 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_sbp + male:entry_pchf + male:entry_dbp,  
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.3c.5 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_sbp + male:entry_pchf + male:entry_hrx,  
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.3c.6 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_sbp + male:entry_pchf + male:entry_smk,  
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1)) 
fp.3c.7 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_sbp + male:entry_pchf + male:entry_alc_elev,  
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.3c.8 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_sbp + male:entry_pchf + male:entry_diab,  
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))
fp.3c.10 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                    entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_sbp + male:entry_pchf + male:entry_pmi,  
                  weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1))

round3c.aic <- data.frame(vars=vars[c(2,3,5:8,10,11)], 
                          aic=c(AIC(fp.3c.1), AIC(fp.3c.2), 
                                AIC(fp.3c.4), AIC(fp.3c.5), AIC(fp.3c.6),
                                AIC(fp.3c.7), AIC(fp.3c.8), AIC(fp.3c.10)))

round3c.aic$delta <- full.aic2c - round3c.aic$aic
sort3c <- round3c.aic[order(round3c.aic[2]),]
sort3c[sort3c$aic>0,]

# All AICS increased, stop here and proceed w/ Round 2 model!


