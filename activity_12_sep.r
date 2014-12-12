my_sample = c(12,70,3,21,90,6,18,63,2,24)
srs = c(50,28,2,4,18,3,49,8,18,25)

tau = mean(my_sample)*100 
s2 = sd(my_sample)^2
tau + c(-1,1)*qt(.975,9)*sqrt(100*99*s2/9) 

library(boot) 
bsamp = boot(data=my_sample,statistic=function(dat,indices){ mean(dat[indices]) },R=1000)
boot.ci(bsamp, type="bca", index=1) 

tau2 = mean(srs)*100 
ss2 = sd(srs)^2
tau2 + c(-1,1)*qt(.975,9)*sqrt(100*99*ss2/9) 
