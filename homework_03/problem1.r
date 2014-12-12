##################################
# Problem 1
##################################
d = read.csv('harvest.csv')
N = 67
M = 190
n = nrow(d)

#(muhat = N/M*mean(d$harvest))
#(se_muhat = 1/M*sqrt(N*(N-n)*var(d$harvest)/n))

# part (a)
(murhat = sum(d$harvest)/sum(d$size))
sr2 = sum((d$harvest - murhat*d$size)^2)/(n-1)
(se_murhat = sqrt( (1-n/N)/mean(d$size)^2 * sr2/n ))
(murhat_ci = murhat + c(-1,1)*qnorm(.975)*se_murhat)

# part (b)
(tauhat = N*mean(d$harvest))
(se_tauhat = sqrt(N*(N-n)*var(d$harvest)/n))
(tauhat_ci = tauhat + c(-1,1)*qnorm(.975)*se_tauhat)

# part (c)
plot(d$size,d$harvest,xlab="Household Size",ylab="Household Harvest (in pounds)")
out = lm(d$harvest ~ d$size)
abline(out,col="red")
#dev.copy(pdf,'harvest_ratio.pdf') 
#dev.off() 

(taurhat = M*murhat)
(se_taurhat = M*se_murhat)
(se2_taurhat = sqrt(N^2*(1-n/N) * sr2/n ))

(taureghat = N*out$coefficients[1] + M*out$coefficients[2])
sreg2 = sum(out$residuals^2)/(n-2)
(se_taureghat = N*sqrt((1-n/N) * sreg2/n))
