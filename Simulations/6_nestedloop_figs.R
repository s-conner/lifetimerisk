# Nested loop plots of simulation results
# Updated 1/22/2021

library(scales)
source("4_nestedloop_fxn.r")

# Load dataset 'res' with simulation results
myres <- read.csv("resultsv2.csv")
nldata <- nestedloop(myres, varnames=c("setting", "sampsize", "probcens", "ptrunc"),
                     varlabels=c("Setting", "Sample size",
                                 "Censoring", "Left truncation"))

pd2 <- nldata
pd2$probcens <- factor(pd2$probcens, levels=c(0.25, 0.75), labels=c("25%", "75%"))
pd2$ptrunc <- factor(pd2$ptrunc, levels=c(0.5, 0.8), labels=c("50%", "80%"))

nl.lines <- function(x, ...) lines(c(x, tail(x, 1)), ...)

#####  ----- BIAS  ----- #####

pdf("plots//Biasv2.pdf", width=12, height=6)
layout(matrix(c(1, 4, 2, 4, 3, 4), ncol=3), heights=c(3, .5))
par(mai=c(0.25, .4, .75, 0.25), oma=c(0,2,0,0))
ymax <- .075
refmin <- .04
cex.ref <- 1
        
# LTR 1
plot(NULL, 
     xlim=c(1, 64), ylim=c(-.07, ymax),
     ylab="", xlab='', main='Group 1', cex.lab=2, cex.main=2, cex.axis=1.5,
     las=1, xaxt="n")
abline(h=0, lty=2, col='black')
title(ylab="Bias", outer=TRUE, line=.25, cex.lab=2, cex.axis=1.5)
lines(pd2, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd2, which="r", ymin.refline=refmin, ymax.refline=ymax, cex.ref=cex.ref)
nl.lines(pd2$bias.ltr1.pseudo, col=alpha("deepskyblue", .6), type="s", lwd=3)      
nl.lines(pd2$bias.ltr1.fg, col=alpha("green3", .6), type="s", lwd=3, lty=1)  
nl.lines(pd2$bias.ltr1.fglogt, col=alpha("red", .6), type="s", lwd=3, lty=1)  
nl.lines(pd2$bias.ltr1.flex, col=alpha("slateblue1", .6), type="s", lwd=3, lty=1)  

# LTR 0
plot(NULL,
     xlim=c(1, 64), ylim=c(-.07, ymax),
     ylab="", xlab='', main='Group 0', cex.lab=2, cex.main=2, cex.axis=1.5,
     las=1, xaxt="n")
abline(h=0, lty=2, col='black')
lines(pd2, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd2, which="r", ymin.refline=refmin, ymax.refline=ymax, cex.ref=cex.ref)
nl.lines(pd2$bias.ltr0.pseudo, col=alpha("deepskyblue", .6), type="s", lwd=3)      
nl.lines(pd2$bias.ltr0.fg, col=alpha("green3", .6), type="s", lwd=3, lty=1)  
nl.lines(pd2$bias.ltr0.fglogt, col=alpha("red", .6), type="s", lwd=3, lty=1) 
nl.lines(pd2$bias.ltr0.flex, col=alpha("slateblue1", .6), type="s", lwd=3, lty=1)  

# Difference
plot(NULL, 
     xlim=c(1, 64), ylim=c(-.07, ymax),
     ylab="", xlab='', main='Difference (Group 1 - Group 0)', cex.lab=2, cex.main=2, cex.axis=1.5,
     las=1, xaxt="n")
abline(h=0, lty=2, col='black')
lines(pd2, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd2, which="r", ymin.refline=refmin, ymax.refline=ymax, cex.ref=cex.ref)
nl.lines(pd2$bias.ltrdiff.pseudo, col=alpha("deepskyblue", .6), type="s", lwd=3)      
nl.lines(pd2$bias.ltrdiff.fg, col=alpha("green3", .6), type="s", lwd=3, lty=1)   
nl.lines(pd2$bias.ltrdiff.fglogt, col=alpha("red", .6), type="s", lwd=3, lty=1)   
nl.lines(pd2$bias.ltrdiff.flex, col=alpha("slateblue1", .6), type="s", lwd=3, lty=1)    

# Legend
par(mai=c(0,0,0,0))
plot.new()
legend(x="center", ncol=2, lty=1, lwd=3, cex=1.5, pt.cex=1.5,
       col=c(alpha('deepskyblue', .6), alpha('slateblue1', .6), alpha('green3', .6), alpha('red'), .6), 
       legend=c('Pseudo-observation', 'Flexible parametric', 'Fine-Gray', 'Fine-Gray, log(time) interaction'))

dev.off()




##### ----- RMSE ----- #####

pdf("plots//rmsev2.pdf", width=12, height=6)
layout(matrix(c(1, 4, 2, 4, 3, 4), ncol=3), heights=c(4, .5))
par(mai=c(0.25, .4, .75, 0.25), oma=c(0,3.2,0,0))
ymax <- .11
ymin <- .085
cex.ref <- 1

# LTR 1
plot(NULL,
     xlim=c(1, 64), ylim=c(0, ymax),
     ylab="", xlab='', main='Group 1', cex.lab=2, cex.main=2, cex.axis=1.5,
     las=1, xaxt="n")
title(ylab="Root mean squared error", outer=TRUE, line=1.75, cex.lab=2, cex.axis=1.5)
abline(h=0, lty=2, col='black')
lines(pd2, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd2, which="r", ymin.refline=ymin, ymax.refline=ymax, cex.ref=cex.ref)
nl.lines(pd2$rmse.ltr1.pseudo, col=alpha("deepskyblue", .6), type="s", lwd=3)      
nl.lines(pd2$rmse.ltr1.fg, col=alpha("green3", .6), type="s", lwd=3, lty=1)  
nl.lines(pd2$rmse.ltr1.fglogt, col=alpha("red", .6), type="s", lwd=3, lty=1)  
nl.lines(pd2$rmse.ltr1.flex, col=alpha("slateblue1", .6), type="s", lwd=3, lty=1)  

# LTR 0
plot(NULL,
     xlim=c(1, 64), ylim=c(0, ymax),
     ylab="", xlab='', main='Group 0', cex.lab=2, cex.main=2, cex.axis=1.5,
     las=1, xaxt="n")
#title(ylab="Mean squared error", line=5, cex.lab=2, cex.axis=1.5)
abline(h=0, lty=2, col='black')
lines(pd2, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd2, which="r", ymin.refline=ymin, ymax.refline=ymax, cex.ref=cex.ref)
nl.lines(pd2$rmse.ltr0.pseudo, col=alpha("deepskyblue", .6), type="s", lwd=3)      
nl.lines(pd2$rmse.ltr0.fg, col=alpha("green3", .6), type="s", lwd=3, lty=1)  
nl.lines(pd2$rmse.ltr0.fglogt, col=alpha("red", .6), type="s", lwd=3, lty=1)  
nl.lines(pd2$rmse.ltr0.flex, col=alpha("slateblue1", .6), type="s", lwd=3, lty=1)  

# Difference
plot(NULL,
     xlim=c(1, 64), ylim=c(0, ymax),
     ylab="", xlab='', main='Difference (Group 1 - Group 0)', cex.lab=2, cex.main=2, cex.axis=1.5,
     las=1, xaxt="n")
#title(ylab="Mean squared error", line=5, cex.lab=2, cex.axis=1.5)
abline(h=0, lty=2, col='black')
lines(pd2, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd2, which="r", ymin.refline=ymin, ymax.refline=ymax, cex.ref=cex.ref)
nl.lines(pd2$rmse.ltrdiff.pseudo, col=alpha("deepskyblue", .6), type="s", lwd=3)      
nl.lines(pd2$rmse.ltrdiff.fg, col=alpha("green3", .6), type="s", lwd=3, lty=1)    
nl.lines(pd2$rmse.ltrdiff.fglogt, col=alpha("red", .6), type="s", lwd=3, lty=1)    
nl.lines(pd2$rmse.ltrdiff.flex, col=alpha("slateblue1", .6), type="s", lwd=3, lty=1)    

# Legend
par(mai=c(0,0,0,0))
plot.new()
legend(x="center", ncol=2, lty=1, lwd=3, cex=1.5, pt.cex=1.5,
       col=c(alpha('deepskyblue', .6), alpha('slateblue1', .6), alpha('green3', .6), alpha('red'), .6), 
       legend=c('Pseudo-observation', 'Flexible parametric', 'Fine-Gray', 'Fine-Gray, log(time) interaction'))
dev.off()



##### -----  COVERAGE -----  #####

pdf("plots//Coveragev2.pdf", width=12, height=6)
layout(matrix(c(1, 4, 2, 4, 3, 4), ncol=3), heights=c(4, .5))
par(mai=c(0.25, .4, .75, 0.25), oma=c(0,2,0,0))
ymax <- 1
ymin <- .975
cex.ref <- 1

# LTR 1
plot(NULL,
     xlim=c(1, 64), ylim=c(0.9, ymax),
     ylab="", xlab='', main='Group 1', cex.lab=2, cex.main=2, cex.axis=1.5,
     las=1, xaxt="n", yaxt='n')
title(ylab="Coverage", outer=TRUE, line=.25, cex.lab=2, cex.axis=1.5)
abline(h=.95, lty=2, col='black')
axis(side=2, at=c(0.8, 0.9, 0.95, 1), cex.axis=1.5)
lines(pd2, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd2, which="r", ymin.refline=ymin, ymax.refline=ymax, cex.ref=.9)
nl.lines(pd2$cov.ltr1.pseudo, col=alpha("deepskyblue", .6), type="s", lwd=3)      
nl.lines(pd2$cov.ltr1.flex, col=alpha("slateblue1", .6), type="s", lwd=3, lty=1) 

# LTR 0
plot(NULL,
     xlim=c(1, 64), ylim=c(0.9, ymax),
     ylab="", xlab='', main='Group 0', cex.lab=2, cex.main=2, cex.axis=1.5,
     las=1, xaxt="n", yaxt='n')
abline(h=.95, lty=2, col='black')
axis(side=2, at=c(0.8, 0.9, 0.95, 1), cex.axis=1.5)
lines(pd2, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd2, which="r", ymin.refline=ymin, ymax.refline=ymax, cex.ref=.9)
nl.lines(pd2$cov.ltr0.pseudo, col=alpha("deepskyblue", .6), type="s", lwd=3)      
nl.lines(pd2$cov.ltr0.flex, col=alpha("slateblue1", .6), type="s", lwd=3, lty=1)  

# Difference
plot(NULL,
     xlim=c(1, 64), ylim=c(0.9, ymax),
     ylab="", xlab='', main='Difference (Group 1 - Group 0)', cex.lab=2, cex.main=2, cex.axis=1.5,
     las=1, xaxt="n", yaxt='n')
abline(h=.95, lty=2, col='black')
axis(side=2, at=c(0.8, 0.9, 0.95, 1), cex.axis=1.5)
lines(pd2, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd2, which="r", ymin.refline=ymin, ymax.refline=ymax, cex.ref=.9)
nl.lines(pd2$cov.ltrdiff.pseudo, col=alpha("deepskyblue", .6), type="s", lwd=3)      
nl.lines(pd2$cov.ltrdiff.flex, col=alpha("slateblue1", .6), type="s", lwd=3, lty=1)    

# Legend
par(mai=c(0,0,0,0))
plot.new()
legend(x="center", ncol=2, lty=c(1,5), lwd=3, cex=1.5, pt.cex=1,
       col=c(alpha('deepskyblue', .6), alpha('slateblue1', .6)), legend=c('Pseudo-observation', 'Flexible parametric'))
dev.off()



##### Type 1 error & power #####
type1 <- pd2[pd2$setting %in% c(1,5),]
power <- pd2[pd2$setting %in% c(2,3,4,6,7,8),]
power$setting <- factor(power$setting, levels=c(2,3,4,6,7,8), labels=c('2', '3', '4', '6', '7', '8'))

pdf("plots//PowerType1v2.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(NULL,
     xlim=c(1, 16), ylim=c(0, .1),
     ylab="", xlab='', cex.lab=1.5, cex.main=1.5, cex.axis=1, main='Type I error',
     las=1, xaxt="n")
abline(h=.05, lty=2, col='black')
title(ylab="Type I error", outer=TRUE, line=0, cex.lab=1.2, cex.axis=1)
lines(type1, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(type1, which="r", ymin.refline=.07, ymax.refline=.1, cex.ref=.75)
nl.lines(type1$power.pseudo, col=alpha("deepskyblue", .6), type="s", lwd=3)      

plot(NULL,
     xlim=c(0, 48), ylim=c(0, 1.45),
     ylab="", xlab='', cex.lab=1.5, cex.main=1.5, main='Power',
     las=1, yaxt="n", xaxt="n")
axis(side=2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis=1)
title(ylab="Power", outer=TRUE, line=0, cex.lab=1.2, cex.axis=1)
lines(power, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(power, which="r", ymin.refline=1.02, ymax.refline=1.45, cex.ref=.75)
nl.lines(power$power.pseudo, col=alpha("deepskyblue", .6), type="s", lwd=3)      
dev.off()




##### RELATIVE BIAS #####
pdf("plots//RelativeBiasv2.pdf", width=12, height=6)
layout(matrix(c(1, 4, 2, 4, 3, 4), ncol=3), heights=c(4, .5))
par(mai=c(0.25, .4, .75, 0.25),oma=c(0,2,0,0))
ymax <- .45
cex.ref <- 1
refmin <- .2
        
# LTR 1
plot(NULL, 
     xlim=c(1, 48), ylim=c(-.8, ymax),
     ylab="", xlab='', main='Group 1', cex.lab=2, cex.main=2, cex.axis=1.5,
     las=1, xaxt="n")
title(ylab="Relative bias", outer=TRUE, line=.25, cex.lab=2, cex.axis=1.5)
abline(h=0, lty=2, col='black')
lines(power, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(power, which="r", ymin.refline=refmin, ymax.refline=ymax, cex.ref=cex.ref)
nl.lines(power$relbias.ltr1.pseudo, col=alpha("deepskyblue", .6), type="s", lwd=3)      
nl.lines(power$relbias.ltr1.fg, col=alpha("red", .6), type="s", lwd=3, lty=1)  
nl.lines(power$relbias.ltr1.fglogt, col=alpha("green3", .6), type="s", lwd=3, lty=1)  
nl.lines(power$relbias.ltr1.flex, col=alpha("slateblue1", .6), type="s", lwd=3, lty=1) 

# LTR 0
plot(NULL,
     xlim=c(1, 48), ylim=c(-.8, ymax),
     ylab="", xlab='', main='Group 0', cex.lab=2, cex.main=2, cex.axis=1.5,
     las=1, xaxt="n")
abline(h=0, lty=2, col='black')
lines(power, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(power, which="r", ymin.refline=refmin, ymax.refline=ymax, cex.ref=cex.ref)
nl.lines(power$relbias.ltr0.pseudo, col=alpha("deepskyblue", .6), type="s", lwd=3)      
nl.lines(power$relbias.ltr0.fg, col=alpha("red", .6), type="s", lwd=3, lty=1)  
nl.lines(power$relbias.ltr0.fglogt, col=alpha("green3", .6), type="s", lwd=3, lty=1)  
nl.lines(power$relbias.ltr0.flex,  col=alpha("slateblue1", .6), type="s", lwd=3, lty=1)  

# Difference
plot(NULL, 
     xlim=c(1, 48), ylim=c(-.8, ymax),
     ylab="", xlab='', main='Difference (Group 1 - Group 0)', cex.lab=2, cex.main=2, cex.axis=1.5,
     las=1, xaxt="n")
abline(h=0, lty=2, col='black')
lines(power, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(power, which="r", ymin.refline=refmin, ymax.refline=ymax, cex.ref=cex.ref)
nl.lines(power$relbias.ltrdiff.pseudo, col=alpha("deepskyblue",.6), type="s", lwd=3)      
nl.lines(power$relbias.ltrdiff.fg, col=alpha("green3",.6), type="s", lwd=3, lty=1)    
nl.lines(power$relbias.ltrdiff.fglogt, col=alpha("red",.6), type="s", lwd=3, lty=1)    
nl.lines(power$relbias.ltrdiff.flex, col=alpha("slateblue1",.6), type="s", lwd=3, lty=1)    

# Legend
par(mai=c(0,0,0,0))
plot.new()
legend(x="center", ncol=2, lty=1, lwd=3, cex=1.5, pt.cex=1.5,
       col=c(alpha('deepskyblue', .6), alpha('slateblue1', .6), alpha('green3', .6), alpha('red'), .6), 
       legend=c('Pseudo-observation', 'Flexible parametric', 'Fine-Gray', 'Fine-Gray, log(time) interaction'))
dev.off()









##### -----  RELATIVE SE -----  #####

pdf("plots//RelSE.pdf", width=12, height=6)
layout(matrix(c(1, 4, 2, 4, 3, 4), ncol=3), heights=c(4, .5))
par(mai=c(0.25, .4, .75, 0.25), oma=c(0,2,0,0))
ymax <- 1.2
ymin <- 1.12
cex.ref <- 1

# LTR 1
plot(NULL,
     xlim=c(1, 64), ylim=c(0.9, ymax),
     ylab="", xlab='', main='Group 1', cex.lab=2, cex.main=2, cex.axis=1.5,
     las=1, xaxt="n", yaxt='n')
title(ylab="Relative standard errors", outer=TRUE, line=.25, cex.lab=2, cex.axis=1.5)
abline(h=1, lty=2, col='black')
axis(side=2, at=c(0.8, 0.9, 0.95, 1), cex.axis=1.5)
lines(pd2, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd2, which="r", ymin.refline=ymin, ymax.refline=ymax, cex.ref=.9)
nl.lines(pd2$relse.ltr1.pseudo, col=alpha("deepskyblue", .6), type="s", lwd=3)      
nl.lines(pd2$relse.ltr1.flex, col=alpha("slateblue1", .6), type="s", lwd=3, lty=1) 

# LTR 0
plot(NULL,
     xlim=c(1, 64), ylim=c(0.9, ymax),
     ylab="", xlab='', main='Group 0', cex.lab=2, cex.main=2, cex.axis=1.5,
     las=1, xaxt="n", yaxt='n')
abline(h=1, lty=2, col='black')
axis(side=2, at=c(0.8, 0.9, 0.95, 1), cex.axis=1.5)
lines(pd2, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd2, which="r", ymin.refline=ymin, ymax.refline=ymax, cex.ref=.9)
nl.lines(pd2$relse.ltr0.pseudo, col=alpha("deepskyblue", .6), type="s", lwd=3)      
nl.lines(pd2$relse.ltr0.flex, col=alpha("slateblue1", .6), type="s", lwd=3, lty=1)  

# Difference
plot(NULL,
     xlim=c(1, 64), ylim=c(0.9, ymax),
     ylab="", xlab='', main='Difference (Group 1 - Group 0)', cex.lab=2, cex.main=2, cex.axis=1.5,
     las=1, xaxt="n", yaxt='n')
abline(h=1, lty=2, col='black')
axis(side=2, at=c(0.8, 0.9, 0.95, 1), cex.axis=1.5)
lines(pd2, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd2, which="r", ymin.refline=ymin, ymax.refline=ymax, cex.ref=.9)
nl.lines(pd2$relse.ltrdiff.pseudo, col=alpha("deepskyblue", .6), type="s", lwd=3)      
nl.lines(pd2$relse.ltrdiff.flex, col=alpha("slateblue1", .6), type="s", lwd=3, lty=1)    

# Legend
par(mai=c(0,0,0,0))
plot.new()
legend(x="center", ncol=2, lty=c(1,5), lwd=3, cex=1.5, pt.cex=1,
       col=c(alpha('deepskyblue', .6), alpha('slateblue1', .6)), legend=c('Pseudo-observation', 'Flexible parametric'))
dev.off()













