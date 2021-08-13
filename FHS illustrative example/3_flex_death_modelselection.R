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
af.fg <- crprep(Tstop="exitage", status="event", data=af, trans=2, cens=0, Tstart="entryage_55", id="fid", keep=vars, shorten=FALSE)
af.fg$weight_ltrc <- af.fg$weight.cens*af.fg$weight.trunc
af.fg$ev2 <- ifelse(af.fg$status==2,1,0)


# ----- Determine # knots in FP model -----

fp.1 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
              entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
            weight=weight_ltrc, data=af.fg, df=1)
fp.2 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
              entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
            weight=weight_ltrc, data=af.fg, df=2)
fp.3 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
              entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
            weight=weight_ltrc, data=af.fg, df=3)
fp.4 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
              entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
            weight=weight_ltrc, data=af.fg, df=4)
fp.5 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
              entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
            weight=weight_ltrc, data=af.fg, df=5)

AIC(fp.1)
AIC(fp.2)
AIC(fp.3) #lowest
AIC(fp.4) 
AIC(fp.5) 


# ----- Round 1: 5 knots, add 1 log(time) interaction at a time ----- 

fp.int1 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                  entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1))

fp.int2 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(entry_hgt=1))

fp.int3 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(entry_wgt=1))

fp.int4 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(entry_sbp=1))

fp.int5 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(entry_dbp=1))

fp.int6 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(entry_hrx=1))

fp.int7 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(entry_smk=1))

fp.int8 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(entry_alc_elev=1))

fp.int9 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(entry_diab=1))

fp.int10 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(entry_pchf=1))

fp.int11 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(entry_pmi=1))

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

# Lowest AIC,  didn't converge
summary(fp.int5)
summary(fp.int4) 

# Next lowest AIC 12886.009
summary(fp.int1) 



# ----- Round 2: 5 knots, male + add 1 log(time) interaction at a time ----- 

fp.2.int2 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hgt=1))

fp.2.int3 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_wgt=1))

fp.2.int6 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1))

fp.2.int7 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_smk=1))

fp.2.int8 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_alc_elev=1))

fp.2.int9 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_diab=1))

fp.2.int10 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                    entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                  weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_pchf=1))

fp.2.int11 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                    entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                  weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_pmi=1))

round2.aic <- cbind(c(2,3,6:11), c(AIC(fp.2.int2)
                            , AIC(fp.2.int3) 
                            , AIC(fp.2.int6)
                            , AIC(fp.2.int7)
                            , AIC(fp.2.int8)
                            , AIC(fp.2.int9)
                            , AIC(fp.2.int10)
                            , AIC(fp.2.int11)
))

round2.aic[order(round2.aic[,2]),]

# lowest AIC 12878.11
summary(fp.2.int6)


# ----- Round 3: 5 knots, male + entry_hrx + add 1 log(time) interaction at a time ----- 

fp.3.int2 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_hgt=1))

fp.3.int3 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_wgt=1))

fp.3.int7 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_smk=1))

fp.3.int8 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_alc_elev=1))

fp.3.int9 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1))

fp.3.int10 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                      entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                    weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_pchf=1))

fp.3.int11 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                      entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                    weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_pmi=1))

round3.aic <- cbind(c(2,3,7:11), c(AIC(fp.3.int2)
                                   , AIC(fp.3.int3) 
                                   , AIC(fp.3.int7)
                                   , AIC(fp.3.int8)
                                   , AIC(fp.3.int9)
                                   , AIC(fp.3.int10)
                                   , AIC(fp.3.int11)
))

round3.aic[order(round3.aic[,2]),]

# lowest AIC 12860.41
summary(fp.3.int9)


# ----- Round 4: 5 knots, male + entry_hrx + entry_diab=1, add 1 log(time) interaction at a time ----- 

fp.4.int2 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_hgt=1))

fp.4.int3 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_wgt=1))

fp.4.int7 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1))

fp.4.int8 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_alc_elev=1))

fp.4.int10 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                      entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                    weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_pchf=1))

fp.4.int11 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                      entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                    weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_pmi=1))

round4.aic <- cbind(c(2,3,7,8,10,11), c(AIC(fp.4.int2)
                                   , AIC(fp.4.int3) 
                                   , AIC(fp.4.int7)
                                   , AIC(fp.4.int8)
                                   , AIC(fp.4.int10)
                                   , AIC(fp.4.int11)
))

round4.aic[order(round4.aic[,2]),]

# lowest AIC 12850.54
summary(fp.4.int7)



# ----- Round 5: 5 knots, male + entry_hrx + entry_diab + entry_smk, add 1 log(time) interaction at a time ----- 

fp.5.int2 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_hgt=1))

fp.5.int3 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_wgt=1))

fp.5.int8 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_alc_elev=1))

fp.5.int10 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                      entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                    weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pchf=1))

fp.5.int11 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                      entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                    weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pmi=1))

round5.aic <- cbind(c(2,3,8,10,11), c(AIC(fp.5.int2)
                                        , AIC(fp.5.int3) 
                                        , AIC(fp.5.int8)
                                        , AIC(fp.5.int10)
                                        , AIC(fp.5.int11)
))

round5.aic[order(round5.aic[,2]),]

# lowest AIC 12840.71
summary(fp.5.int10)




# ----- Round 6: 5 knots, male + entry_hrx + entry_diab + entry_smk + entry_pchf, add 1 log(time) interaction at a time ----- 

fp.6.int2 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pchf=1, entry_hgt=1))

fp.6.int3 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pchf=1, entry_wgt=1))

fp.6.int8 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pchf=1, entry_alc_elev=1))

fp.6.int11 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                      entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                    weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pchf=1, entry_pmi=1))

round6.aic <- cbind(c(2,3,8,11), c(AIC(fp.6.int2)
                                      , AIC(fp.6.int3) 
                                      , AIC(fp.6.int8)
                                      , AIC(fp.6.int11)
))

round6.aic[order(round6.aic[,2]),]

# lowest AIC 12838.70
summary(fp.6.int3)



# ----- Round 7: 5 knots, male + entry_hrx + entry_diab + entry_smk + entry_pchf + entry_wgt, add 1 log(time) interaction at a time ----- 

fp.7.int2 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=3, 
                   tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pchf=1, entry_wgt=1, entry_hgt=1))

fp.7.int8 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=3, 
                   tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pchf=1, entry_wgt=1, entry_alc_elev=1))

fp.7.int11 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                      entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                    weight=weight_ltrc, data=af.fg, df=3, 
                    tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pchf=1, entry_wgt=1, entry_pmi=1))

round7.aic <- cbind(c(2,8,11), c(AIC(fp.7.int2)
                                   , AIC(fp.7.int8)
                                   , AIC(fp.7.int11)
))

round7.aic[order(round7.aic[,2]),]

# lowest AIC 12839.03 increased, proceed w/ Round 6 model
summary(fp.6.int3)


