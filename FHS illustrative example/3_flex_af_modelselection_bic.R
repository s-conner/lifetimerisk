### Flexible parametric models for FHS - Lifetime risk of AF
### Choose # of knots and include interactions with log(time), determine using BIC
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

BIC(fp.1)
BIC(fp.2) # BIC decrease by 16.6
BIC(fp.3) # BIC increase
BIC(fp.4) # BIC decrease by 3.6
BIC(fp.5) # BIC increase

full.BIC0 <- BIC(fp.2)


# ----- Round 1: 2 knots, add 1 log(time) interaction at a time ----- 

fp.int1 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(male=1))

fp.int2 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_hgt=1))

fp.int3 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1))

fp.int4 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_sbp=1))

fp.int5 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_dbp=1))

fp.int6 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_hrx=1))

fp.int7 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_smk=1))

fp.int8 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_alc_elev=1))

fp.int9 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_diab=1))

fp.int10 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                    entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                  weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_pchf=1))

fp.int11 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                    entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                  weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_pmi=1))

round1.BIC <- data.frame(var=1:11, BIC=c(BIC(fp.int1)
                            , BIC(fp.int2)
                            , BIC(fp.int3) 
                            , BIC(fp.int4)
                            , BIC(fp.int5)
                            , BIC(fp.int6)
                            , BIC(fp.int7)
                            , BIC(fp.int8)
                            , BIC(fp.int9)
                            , BIC(fp.int10)
                            , BIC(fp.int11)
))

round1.BIC$delta <- full.BIC0 - round1.BIC$BIC
round1.BIC[order(round1.BIC[,2]),]
full.BIC1 <- BIC(fp.int3) 

# Lowest BIC
summary(fp.int3) 



# ----- Round 2: 2 knots, entry_wgt + add 1 log(time) interaction at a time ----- 

fp.2.int1 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, male=1))

fp.2.int2 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_hgt=1))

fp.2.int4 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_sbp=1))

fp.2.int5 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_dbp=1))

fp.2.int6 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_hrx=1))

fp.2.int7 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_smk=1))

fp.2.int8 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_alc_elev=1))

fp.2.int9 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_diab=1))

fp.2.int10 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                    entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                  weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_pchf=1))

fp.2.int11 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                    entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                  weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_pmi=1))

round2.BIC <- data.frame(var=vars[c(1,2,4:11)], BIC=c(BIC(fp.2.int1)
                            , BIC(fp.2.int2)
                            , BIC(fp.2.int4)
                            , BIC(fp.2.int5)
                            , BIC(fp.2.int6)
                            , BIC(fp.2.int7)
                            , BIC(fp.2.int8)
                            , BIC(fp.2.int9)
                            , BIC(fp.2.int10)
                            , BIC(fp.2.int11)
))

round2.BIC$delta <- full.BIC1 - round2.BIC$BIC
round2.BIC[order(round2.BIC[,2]),]
full.BIC2 <- BIC(fp.2.int11) # BIC=

summary(fp.2.int1) # Delta 2681 for male sex but did not converge
summary(fp.2.int2) # Delta 1017 for height but did not converge
summary(fp.2.int11) # Delta 14 for MI



# ----- Round 3: 2 knots, entry_wgt + entry_pmi, add 1 log(time) interaction at a time ----- 

fp.3.int1 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_pmi=1, male=1))

fp.3.int2 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_pmi=1, entry_hgt=1))

fp.3.int4 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_pmi=1, entry_sbp=1))

fp.3.int5 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_pmi=1, entry_dbp=1))

fp.3.int6 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_pmi=1, entry_hrx=1))

fp.3.int7 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_pmi=1, entry_smk=1))

fp.3.int8 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_pmi=1, entry_alc_elev=1))

fp.3.int9 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                     entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                   weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_pmi=1, entry_diab=1))

fp.3.int10 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                      entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                    weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_pmi=1, entry_pchf=1))


round3.BIC <- data.frame(var=vars[c(1,2,4:10)], 
                    BIC=c(BIC(fp.3.int1)
                      , BIC(fp.3.int2)
                      , BIC(fp.3.int4)
                      , BIC(fp.3.int5)
                      , BIC(fp.3.int6)
                      , BIC(fp.3.int7)
                      , BIC(fp.3.int8)
                      , BIC(fp.3.int9)
                      , BIC(fp.3.int10)
))

round3.BIC$delta <- full.BIC2 - round3.BIC$BIC
sort3 <- round3.BIC[order(round3.BIC[,2]),]
sort3[sort3$BIC>0,]
#full.BIC3 <- round3.BIC$BIC[2]

summary(fp.3.int2) # Delta 988 for height, did not converge. 
# Other deltas small or increase, stop here



# ----- Round 1 of adding additional df for time-varying terms -----

fp.1b.1 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                      entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                    weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=2, entry_pmi=1))
fp.1b.2 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi, 
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_pmi=2))

round1b.BIC <- data.frame(BIC=c(BIC(fp.1b.1), BIC(fp.1b.2)))

round1b.BIC$delta <- full.BIC2 - round1b.BIC$BIC
sort1b <- round1b.BIC[order(round1b.BIC[1]),]
sort1b[sort1b$BIC>0,]

summary(fp.1b.1) # Did not converge

# Do not add additional knots to time-varying terms


# ----- Round 1 of adding sex interaction terms -----

fp.1c.1 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_hgt, 
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_pmi=1))
fp.1c.2 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_wgt,  
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_pmi=1))
fp.1c.3 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_sbp,  
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_pmi=1))
fp.1c.4 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_dbp,  
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_pmi=1))
fp.1c.5 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_hrx,  
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_pmi=1))
fp.1c.6 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_smk,  
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_pmi=1))
fp.1c.7 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_alc_elev,  
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_pmi=1))
fp.1c.8 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_diab,  
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_pmi=1))
fp.1c.9 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_pchf,  
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_pmi=1))
fp.1c.10 <- stpm2(Surv(Tstart, Tstop, ev1) ~ male + entry_hgt + entry_wgt + entry_sbp + entry_dbp + entry_hrx + 
                   entry_smk + entry_alc_elev + entry_diab + entry_pchf + entry_pmi + male:entry_pmi,  
                 weight=weight_ltrc, data=af.fg, df=2, tvc=list(entry_wgt=1, entry_pmi=1))


round1c.BIC <- data.frame(vars=vars[c(2, 4:9, 11)], 
                          BIC=c(BIC(fp.1c.1), BIC(fp.1c.3),
                                BIC(fp.1c.4), BIC(fp.1c.5), BIC(fp.1c.6),
                                BIC(fp.1c.7), BIC(fp.1c.8), BIC(fp.1c.10)))

round1c.BIC$delta <- full.BIC2 - round1c.BIC$BIC
sort1c <- round1c.BIC[order(round1c.BIC[2]),]
sort1c[sort1c$BIC>0,] 

# All deltas increased, stop here


