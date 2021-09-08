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

full.aic0 <- AIC(fp.3)


# ----- Round 1: 3 knots, add 1 log(time) interaction at a time ----- 

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

# Lowest AIC,  didn't converge
summary(fp.int5)
summary(fp.int4) 

# Next lowest AIC that converged 12886.009, delta=33
summary(fp.int1) 
full.aic1 <- AIC(fp.int1)


# ----- Round 2: 3 knots, male + add 1 log(time) interaction at a time ----- 

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

round2.aic <- data.frame(var=vars[c(2,3,6:11)], aic=c(AIC(fp.2.int2)
                                                      , AIC(fp.2.int3)
                                                      , AIC(fp.2.int6)
                                                      , AIC(fp.2.int7)
                                                      , AIC(fp.2.int8)
                                                      , AIC(fp.2.int9)
                                                      , AIC(fp.2.int10)
                                                      , AIC(fp.2.int11)
))

round2.aic$delta <- full.aic1 - round2.aic$aic
round2.aic[order(round2.aic[,2]),]

# lowest AIC 12878.11, delta=7.9
summary(fp.2.int6)
full.aic2 <- min(round2.aic$aic)



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

round3.aic <- data.frame(var=vars[c(2,3,7:11)], aic=c(AIC(fp.3.int2)
                                                      , AIC(fp.3.int3)
                                                      , AIC(fp.3.int7)
                                                      , AIC(fp.3.int8)
                                                      , AIC(fp.3.int9)
                                                      , AIC(fp.3.int10)
                                                      , AIC(fp.3.int11)
))

round3.aic$delta <- full.aic2 - round3.aic$aic
round3.aic[order(round3.aic[,2]),]

# lowest AIC 12860.41, delta=17
summary(fp.3.int9)
full.aic3 <- AIC(fp.3.int9)


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

round4.aic <- data.frame(var=vars[c(2,3,7,8,10,11)], aic=c(AIC(fp.4.int2)
                                   , AIC(fp.4.int3) 
                                   , AIC(fp.4.int7)
                                   , AIC(fp.4.int8)
                                   , AIC(fp.4.int10)
                                   , AIC(fp.4.int11)
))
round4.aic$delta <- full.aic3 - round4.aic$aic
round4.aic[order(round4.aic[,2]),]

# lowest AIC 12850.54, delta=9.9
summary(fp.4.int7)
full.aic4 <- AIC(fp.4.int7)


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

round5.aic <- data.frame(var=vars[c(2,3,8,10,11)], aic=c(AIC(fp.5.int2)
                                        , AIC(fp.5.int3) 
                                        , AIC(fp.5.int8)
                                        , AIC(fp.5.int10)
                                        , AIC(fp.5.int11)
))
round5.aic$delta <- full.aic4 - round5.aic$aic
round5.aic[order(round5.aic[,2]),]

# lowest AIC 12840.71, delta=9.8
summary(fp.5.int10)

full.aic5 <- AIC(fp.5.int10)
full.aic5 <- 12840.71


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

# lowest AIC 12838.70, delta=2, proceed w/ Round 5 model
summary(fp.6.int3)



# ----- Consider higher order terms -----


fp.1b.1 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=2, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pchf=1))
fp.1b.2 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=2, entry_diab=1, entry_smk=1, entry_pchf=1))
fp.1b.3 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=2, entry_smk=1, entry_pchf=1))
fp.1b.4 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=2, entry_pchf=1))
fp.1b.5 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pchf=2))
round1b.aic <- data.frame(aic=c(AIC(fp.1b.1), AIC(fp.1b.2), AIC(fp.1b.3), AIC(fp.1b.4), AIC(fp.1b.5)))

round1b.aic$delta <- full.aic5 - round1b.aic$aic
sort1b <- round1b.aic[order(round1b.aic[1]),]
sort1b[sort1b$aic>0,]

# None significantly improved AIC ()


# ----- Round 1 of adding sex interaction terms -----

fp.1c.1 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_hgt, 
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pchf=1))
fp.1c.2 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_wgt,  
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pchf=1))
fp.1c.3 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_sbp,  
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pchf=1))
fp.1c.4 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_dbp,  
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pchf=1))
fp.1c.5 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_hrx,  
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pchf=1))
fp.1c.6 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_smk,  
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pchf=1))
fp.1c.7 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_alc_elev,  
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pchf=1))
fp.1c.8 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_diab,  
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pchf=1))
fp.1c.9 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_pchf,  
                 weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pchf=1))
fp.1c.10 <- stpm2(Surv(Tstart, Tstop, ev2) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                    entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_pmi,  
                  weight=weight_ltrc, data=af.fg, df=3, tvc=list(male=1, entry_hrx=1, entry_diab=1, entry_smk=1, entry_pchf=1))


round1c.aic <- data.frame(vars=vars[c(2:11)], 
                          aic=c(AIC(fp.1c.1), AIC(fp.1c.2), AIC(fp.1c.3),
                                AIC(fp.1c.4), AIC(fp.1c.5), AIC(fp.1c.6),
                                AIC(fp.1c.7), AIC(fp.1c.8), AIC(fp.1c.9), AIC(fp.1c.10)))

round1c.aic$delta <- full.aic5 - round1c.aic$aic
sort1c <- round1c.aic[order(round1c.aic[2]),]
sort1c[sort1c$aic>0,]

# Delta 3.5 for HF, do not add interaction!

