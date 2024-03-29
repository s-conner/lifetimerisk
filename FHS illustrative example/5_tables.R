# Create tables of all methods for manuscript
# 1/25/2021


# ----- AF lifetime risks ----- 

pseudo.ltr <- read.csv('results//pseudologit_af_pred.csv')
fg.ltr <- read.csv('results//finegray_af_diffltr_ci_revision.csv')
flex.ltr <- read.csv('results//flex_af_diffltr_revision.csv')
flex.ltr.bic <- read.csv('results//flex_af_diffltr_revision_bic.csv')

pseudo.ltr <- pseudo.ltr*100
pseudo.fmt <- cbind(ltrdiff=paste0(sprintf('%.2f', pseudo.ltr[,1]), " (", sprintf('%.2f', pseudo.ltr[,2]), ", ", sprintf('%.2f', pseudo.ltr[,3]), ")"))

fg.ltr[,2:4] <- fg.ltr[,2:4]*100
fg.fmt <- cbind(fg.ltr$var, ltrdiff=paste0(sprintf('%.2f', fg.ltr$ltrdiff), " (", sprintf('%.2f', fg.ltr$cil), ", ", sprintf('%.2f', fg.ltr$ciu), ")"))

allpred <- cbind(fg.ltr[,1], pseudo.fmt[2:12,], fg.fmt[,2], flex.ltr[2], flex.ltr.bic[2])
colnames(allpred) <- c('Risk factor', 'Pseudo-observation', 'Fine-Gray', 'Flexible parametric AIC', 'Flexible parametric BIC')
write.csv(allpred, 'results//all_af_diffltr_revision_bic.csv', row.names=FALSE)


# ----- AF models ----- 

pseudo.mod <- read.csv('results//pseudologit_af_mod.csv')
fg.mod <- read.csv('results//finegray_af_mod_revision.csv')
flex.mod <- read.csv('results//flex_af_mod_revision.csv')

# Re order
fg.mod <- fg.mod[c(1:3,12,4:7,14,8:11,13),]
flex.mod <-  flex.mod[c(1,12,13,2,16,3,14,4:9,15,10,11),]

# Pad with NAs, since models have different parameters
blank <- data.frame(matrix('',nrow=5,ncol=2))
colnames(blank) <- colnames(pseudo.mod)
pseudo.mod <- rbind(pseudo.mod, blank)

blank <- data.frame(matrix('',nrow=2,ncol=2))
colnames(blank) <- colnames(fg.mod)
fg.mod <- rbind(fg.mod, blank)
fg.mod <- fg.mod[c(1,15,16,2:14),]

allmod <- cbind(pseudo.mod, fg.mod, flex.mod)

write.csv(allmod[c(1,2,15,3,12,4:7,14,8,9,16,10,11,13),], 'results//all_af_mod.csv', row.names=FALSE)


# ----- Death models ----- 

pseudo.mod.death <- read.csv('results//pseudologit_death_mod.csv')
fg.mod.death <- read.csv('results//finegray_death_mod.csv')
flex.mod.death <- read.csv('results//flex_death_mod.csv')

# Model coefficients
allmod.d <- cbind(flex.mod.death[,1], 
                c(pseudo.mod.death$ltrdiff, NA, NA, NA, NA, NA, NA),
                c(fg.mod.death[c(1,13,2,3),2], NA, fg.mod.death[c(4,5,6,15,7,14,8,9),2], NA, fg.mod.death[c(10:12),2]),
                c(flex.mod.death[c(1,12,2,3,17,4,5,6,13,7,15,8,9,14,10,16,11),2]))

write.csv(allmod.d[,], 'results//all_death_mod.csv', row.names=FALSE)


# ----- Table 1 ----- 

af <- read.csv('af_pseudo.csv')
af$af <- ifelse(af$event==1,1,0)
af$death <- ifelse(af$event==2,1,0)

# median follow-up
median(af$exitage-af$entryage_55)

vars <- c('af', 'death',
          'male', 'entry_hgt', 'entry_wgt', 'entry_sbp', 'entry_dbp', 'entry_hrx', 
          'entry_smk', 'entry_alc_elev', 'entry_diab', 'entry_pchf', 'entry_pmi')

means <- apply(af[, vars], 2, mean)
sds <- apply(af[, vars], 2, sd)
freqs <- apply(af[, vars], 2, sum)
props <- freqs/nrow(af)

means.f <- sprintf('%.2f', means) 
sds.f <- sprintf('%.2f', sds) 
props.f <- paste0(sprintf('%.2f', props*100), '%')

comb1 <- c(freqs[1:3], means.f[4:7], freqs[8:13])
comb2 <- c(props.f[1:3], sds.f[4:7], props.f[8:13])

tab1 <- cbind(vars, paste0(comb1, " (", comb2, ")"))

write.csv(tab1, 'results//table1.csv', row.names=FALSE)


