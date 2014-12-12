rm(list=ls(all=T))
# Homework 2 Problem 3
n = 4
y = c(2,2,5,10)
x = c(1.2,1.2,0.2,0.5)

# Linearization
t = y/x
v = 1/x
(mu_y = mean(t)/mean(v))
(se_mu_y = sqrt(1/ mean(v)^2 * var(t - mean(t)/mean(v)*v)/n))
# Note that this agrees with the other formula
(sqrt(mean(t)^2/mean(v)^4*var(v)/n + 1/mean(v)^2 * var(t)/n - 2 * mean(t)/mean(v)^3 * cor(t,v)*sqrt(var(t)*var(v))/n))

(mu_y.ci = mu_y + c(-1,1)*qt(.975,n-1) * se_mu_y)

# Bootstrapping # uncomment to run, it takes a while
library(boot)
N = 10^5
grab = function(data,i){ mean(data$t[i])/mean(data$v[i]) }
(straps = boot(data.frame(t,v), grab, N))
plot(straps)
dev.copy(png,'pollution_bootstraps.png') 
dev.off()
(boot.ci(straps))

# Estimate average lake size
(muhat_x  = 1/mean(v))
(se_muhat_x = sqrt(var(v)/(n*mean(v)^4) ))

(muhat_x.ci = muhat_x + c(-1,1)*qt(.975,n-1) * se_muhat_x)

N = 10^5
grab = function(data,i){ 1/mean(data$v[i]) }
(straps = boot(data.frame(t,v), grab, N))
plot(straps)
dev.copy(png,'mean_lake_bootstraps.png') 
dev.off()

(boot.ci(straps))
