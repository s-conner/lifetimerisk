# Example of IPCLW data in subset of aidssi data
# Sarah Conner - 8/12/2021

library(mstate)

data(aidssi)
dat0 <- aidssi
set.seed(1)
dat0$entry <- runif(nrow(dat0), 0, 3)
dat <- dat0[dat0$patnr %in% c(5,6,10,14,18,48),]
dat <- dat[dat$entry < dat$time,]

# Unique entry times, event 1 times, and censoring times
round(sort(dat$entry),1)
round(sort(dat$time[dat$status==1]),1)
round(sort(dat$time[dat$status==0]),1)

dat.ipclw <- crprep(Tstop="time", status="status", data=dat, trans=1, cens=0, Tstart="entry", id="patnr", shorten=FALSE)
dat.ipclw$weight <- dat.ipclw$weight.cens*dat.ipclw$weight.trunc
dat.ipclw$ev1 <- ifelse(dat.ipclw$status==1,1,0)

out <- cbind(id=c(1,2,3,4,5,6,6,6,6),round(dat.ipclw,1))
out <- out[,c('id', 'Tstart', 'Tstop', 'status', 'weight.cens', 'weight.trunc', 'weight')]
write.csv(out, '//restricted//projectnb//conner-thesis//lifetimerisk//ipclw_aidssi_example.csv', row.names=FALSE)



