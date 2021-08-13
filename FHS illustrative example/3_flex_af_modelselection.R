### Flexible parametric models for FHS - Lifetime risk of AF
### Choose # of knots and include interactions with log(time), determine using AIC
### Updated analyses 1/12/2021

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


# ----- Round 1: 5 knots, add 1 log(time) interaction at a time ----- 

fp.int1 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                  entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                weight=weight_ltrc, data=af.fg, df=5, tvc=list(male=1))

fp.int2 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=5, tvc=list(entry_hgt=1))

fp.int3 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=5, tvc=list(entry_wgt=1))

fp.int4 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=5, tvc=list(entry_sbp=1))

fp.int5 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=5, tvc=list(entry_dbp=1))

fp.int6 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=5, tvc=list(entry_hrx=1))

fp.int7 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=5, tvc=list(entry_smk=1))

fp.int8 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=5, tvc=list(entry_alc_elev=1))

fp.int9 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=5, tvc=list(entry_diab=1))

fp.int10 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=5, tvc=list(entry_pchf=1))

fp.int11 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=5, tvc=list(entry_pmi=1))

round1.aic <- cbind(1:11, c(AIC(fp.int1)
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

round1.aic[order(round1.aic[,2]),]

# Lowest AIC, also didn't converge
summary(fp.int9)

# Next lowest AIC 9264.574, converged!
summary(fp.int3) 
AIC(fp.int3) 



# ----- Round 2: 5 knots, entry_wgt + add 1 log(time) interaction at a time ----- 


fp.2.int1 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=5, tvc=list(entry_wgt=1, male=1))

fp.2.int2 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=5, tvc=list(entry_wgt=1, entry_hgt=1))

fp.2.int6 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=5, tvc=list(entry_wgt=1, entry_hrx=1))

fp.2.int8 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=5, tvc=list(entry_wgt=1, entry_alc_elev=1))

fp.2.int10 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                    entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                  weight=weight_ltrc, data=af.fg, df=5, tvc=list(entry_wgt=1, entry_pchf=1))

fp.2.int11 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                    entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                  weight=weight_ltrc, data=af.fg, df=5, tvc=list(entry_wgt=1, entry_pmi=1))

round2.aic <- cbind(c(1,2,6,8,10,11), 
                    c(AIC(fp.2.int1)
                      , AIC(fp.2.int2)
                      , AIC(fp.2.int6)
                      , AIC(fp.2.int8)
                      , AIC(fp.2.int10)
                      , AIC(fp.2.int11)
))

round2.aic[order(round2.aic[,2]),]

# Lowest AIC didn't converge
summary(fp.2.int11)

# Next lowest AIC 9457.297, is higher than previous round - proceed w/ Round 1 model
summary(fp.2.int2)












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

round1.aic <- cbind(1:11, c(AIC(fp.int1)
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

round1.aic[order(round1.aic[,2]),]

# Lowest AIC 9690.055
summary(fp.int3) 


# ----- Round 1b: 4 knots, add 1 df=2 interaction at a time ----- 

fp.int1 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(male=2))

fp.int2 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_hgt=2))

fp.int3 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=2))

fp.int4 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_sbp=2))

fp.int5 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_dbp=2))

fp.int6 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_hrx=2))

fp.int7 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_smk=2))

fp.int8 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_alc_elev=2))

fp.int9 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_diab=2))

fp.int10 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                    entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                  weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_pchf=2))

fp.int11 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                    entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                  weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_pmi=2))

round1.aic <- cbind(1:11, c(AIC(fp.int1)
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

round1.aic[order(round1.aic[,2]),]

# Lowest AIC df=2 9691.407 is higher than df=1 9690.055
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

round2.aic <- cbind(c(1,2,4:11), c(AIC(fp.2.int1)
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

round2.aic[order(round2.aic[,2]),]

# Lowest AIC 9660.525  
summary(fp.2.int9)


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

round3.aic <- cbind(c(1,2,6,8,10,11), 
                    c(AIC(fp.3.int1)
                      , AIC(fp.3.int2)
                      , AIC(fp.3.int6)
                      , AIC(fp.3.int8)
                      , AIC(fp.3.int10)
                      , AIC(fp.3.int11)
))

round3.aic[order(round3.aic[,2]),]

# Lowest AIC 9638.896
summary(fp.3.int2)


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

round4.aic <- cbind(c(1,6,8,10,11), 
                    c(AIC(fp.4.int1)
                      , AIC(fp.4.int6)
                      , AIC(fp.4.int8)
                      , AIC(fp.4.int10)
                      , AIC(fp.4.int11)
                    ))

round4.aic[order(round4.aic[,2]),]

# vars 8, 1 don't converge
# var 6 aic 9635.606
# var 11  9655.976 is higher

summary(fp.4.int6)
AIC(fp.4.int6)

# ----- Round 5: 4 knots, entry_wgt + entry_diab + entry_hgt + entry_hrx, add 1 log(time) interaction at a time ----- 

fp.5.int10 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                      entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                    weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1, entry_hrx=1, entry_pchf=1))

fp.5.int11 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                      entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                    weight=weight_ltrc, data=af.fg, df=4, tvc=list(entry_wgt=1, entry_diab=1, entry_hgt=1, entry_hrx=1, entry_pmi=1))

round5.aic <- cbind(c(10,11), 
                    c(AIC(fp.5.int10), AIC(fp.5.int11)
                    ))

round5.aic[order(round5.aic[,2]),]

# AIC incresaed, stop & proceed w/ round 4 model
