################################################################################
# Problem 2
################################################################################
rm(list=ls())
source('systematic.R')
flows = scan('fraser.txt')
N = length(flows)
ts.plot(flows[1:(5*12)])
dev.copy(pdf,'fraser_flow.pdf')
dev.off()
k = c(3,6,9,11,12,13,23,24)
systematic_samples = sapply(k, function(k){ systematic(flows,k) })
(round(systematic_samples,2))
(n = floor(length(flows)/k))
srs_flows = function(n){
  samp = sample(flows,n)
  c('y.bar'=mean(samp), 'se.y.bar='=sqrt((1-n/N)*var(flows)/n))
}
srs = sapply(n,srs_flows)
round(srs,2)

#par(mfrow=c(3,2))
#for (i in 1:6){
#  plot(flows[seq(i,N,6)])
#  abline(h=mean(flows[seq(i,N,6)]),col='red') 
#}

