## Homework 2 Problem 6
#rm(list=ls(all=T))

## part (a) ##
census = read.csv('Census.csv') 
plot(census$HousingUnits, census$TotalPop) 
title('Housing Units vs Total Population')
#dev.copy(pdf,'house_pop_compare.pdf')
#dev.off()
plot(log(census$HousingUnits), log(census$TotalPop))
title('Log Housing Units vs Log Total Population') 
#dev.copy(pdf,'log_compare.pdf')
#dev.off()

## part (b) ##
n = 100
N = nrow(census)
mu_x = mean(census$HousingUnits)
tauhat.none = function(y,N) {
  mean(y)*N
}

tauhat.ratio = function(y,x,N,mu_x) {
  N*mean(y)/mean(x) * mu_x
}

tauhat.regression = function(y,x,N,mu_x) {
  b = sum((x - mean(x))*(y-mean(y)))/sum((x-mean(x))^2) 
  N*(mean(y) + b*(mu_x - mean(x)))
}

test.estimators = function() { 
  samp = census[sample(1:N,n),] # extract the relevant columns from the data.frame 
  x = samp$HousingUnits
  y = samp$TotalPop
  list(
    t.none = tauhat.none(y,N),
    t.ratio = tauhat.ratio(y,x,N,mu_x),
    t.regression = tauhat.regression(y,x,N,mu_x)
  )
}
out = replicate(100,test.estimators())

# Plot histograms of each estimator
hist(as.numeric(out['t.none',]),main=expression(paste('Histogram of 100 standard estimates of ',tau)),xlab=expression(hat(tau)["none"]),cex.lab=1.6,cex.main=1.6)
abline(v=sum(census$TotalPop),col='red')
#dev.copy(pdf,'none_hist.pdf')
#dev.off()

hist(as.numeric(out['t.ratio',]),main=expression(paste('Histogram of 100 ratio estimates of ',tau)),xlab=expression(hat(tau)["ratio"]),cex.lab=1.6,cex.main=1.6)
abline(v=sum(census$TotalPop),col='red')
#dev.copy(pdf,'ratio_hist.pdf')
#dev.off()

hist(as.numeric(out['t.regression',]),main=expression(paste('Histogram of 100 regression estimates of ',tau)),xlab=expression(hat(tau)["reg"]),cex.lab=1.6,cex.main=1.6)
abline(v=sum(census$TotalPop),col='red')
#dev.copy(pdf,'regression_hist.pdf')
#dev.off()


## part (c) ##
set.seed(45519)  # Set seed for problem
samp = census[sample(1:N,n),] # extract the relevant columns from the data.frame 
x = samp$HousingUnits
y = samp$TotalPop

hatse.tauhat.none = function(y,N) {
  n = length(y)
  sqrt( N*(N-n)*var(y)/n )
}
hatse.tauhat.ratio = function(y,x,N,mu_x) {
  n = length(y)
  sr2 = (1/(n-1)) * sum((y - mean(y)/mean(x)*x)^2)
  N*sqrt((1-n/N)*sr2/n)
}
hatse.tauhat.regression = function(y,x,N,mu_x) {
  n = length(y)
  reg = lsfit(x,y) # I'll take a cue from Dave and use lsfit to get residuals
  N*sqrt((1-n/N)/(n*(n-2)) * sum(reg$residuals^2))
}
cat(" Specific sample results \n")
(data.frame(
  'none'=c(tauhat.none(y,N), hatse.tauhat.none(y,N)),
  'ratio'=c( tauhat.ratio(y,x,N,mu_x), hatse.tauhat.ratio(y,x,N,mu_x) ),
  'regression'=c(tauhat.regression(y,x,N,mu_x), hatse.tauhat.regression(y,x,N,mu_x) )
))

## part (d) ##
cat(" Population based standard errors \n")
se.muhat.none  = sqrt((1-n/N)*var(census$TotalPop)/n)

R = mean(census$TotalPop)/mean(census$HousingUnits)
SR2 = sum((census$TotalPop - R*census$HousingUnits)^2)/(N-1)
se.muhat.ratio = sqrt((1-n/N)*SR2/n)

reg = lsfit(census$HousingUnits,census$TotalPop)
MSE = sum(reg$residuals^2)/(N-1)
se.muhat.regression = sqrt((1-n/N)*MSE/n)
cat(" Muhat \n=============\n")
(ans.d = data.frame(
'none'  = c(se.muhat.none, hatse.tauhat.none(y,N)/N), 
'ratio' = c(se.muhat.ratio, hatse.tauhat.ratio(y,x,N,mu_x)/N ), 
'regression' = c(se.muhat.regression, hatse.tauhat.regression(y,x,N,mu_x)/N)
))

cat(" Tauhat \n=============\n")
(ans.d*N)

test.se.estimators = function() { 
  samp = census[sample(1:N,n),] # extract the relevant columns from the data.frame 
  x = samp$HousingUnits
  y = samp$TotalPop
  list(
    t.none = hatse.tauhat.none(y,N),
    t.ratio = hatse.tauhat.ratio(y,x,N,mu_x),
    t.regression = hatse.tauhat.regression(y,x,N,mu_x)
  )
}
out = replicate(100,test.se.estimators())
# Plot histograms of each estimator
#dev.new()
hist(as.numeric(out['t.none',]),main=expression(paste('Histogram of 100 standard estimates of ',hat("SE")(hat(tau)["none"]))),xlab=expression(hat(tau)["none"]),cex.lab=1.6,cex.main=1.6)
abline(v=se.muhat.none*N,col='red')
#dev.copy(pdf,'none_hist_se.pdf')
#dev.off()

#dev.new()
hist(as.numeric(out['t.ratio',]),main=expression(paste('Histogram of 100 ratio estimates of ',hat("SE")(hat(tau)["ratio"]))),xlab=expression(hat(tau)["ratio"]),cex.lab=1.6,cex.main=1.6)
abline(v=se.muhat.ratio*N,col='red')
#dev.copy(pdf,'ratio_hist_se.pdf')
#dev.off()

#dev.new()
hist(as.numeric(out['t.regression',]),main=expression(paste('Histogram of 100 regression estimates of ',hat("SE")(hat(tau)["reg"]))),xlab=expression(hat(tau)["reg"]),cex.lab=1.6,cex.main=1.6)
abline(v=se.muhat.regression*N,col='red')
#dev.copy(pdf,'regression_hist_se.pdf')
#dev.off()

## part e ##
z = qnorm(.975)
d = 5e6
required_n = function(sig2) { 
  (d^2/(N^2*z^2*sig2) + 1/N)^(-1)
}
sigs = c(
  none       = var(census$TotalPop),
  ratio      = SR2,
  regression = MSE
)
(required_n(sigs))

## part f ##  # This takes forever. Uncomment to run
#out = replicate(1e5,test.estimators())  
#hist(as.numeric(out['t.ratio',]),main=expression(paste('Histogram of 100000 ratio estimates of ',tau)),xlab=expression(hat(tau)["ratio"]),cex.lab=1.6,cex.main=1.6)
#abline(v=sum(census$TotalPop),col='red')
##dev.copy(pdf,'ratio_hist_big.pdf')
##dev.off()
#
#hist(as.numeric(out['t.regression',]),main=expression(paste('Histogram of 100000 regression estimates of ',tau)),xlab=expression(hat(tau)["reg"]),cex.lab=1.6,cex.main=1.6)
#abline(v=sum(census$TotalPop),col='red')
##dev.copy(pdf,'regression_hist_big.pdf')
##dev.off()
#
#t.ratio.samps = as.numeric(out['t.ratio',])
#t.regression.samps = as.numeric(out['t.regression',])
#
#(sim_results = data.frame(
#  ratio	= c(mean(t.ratio.samps) - sum(census$TotalPop), sd(t.ratio.samps)),
#  regression = c(mean(t.regression.samps) - sum(census$TotalPop), sd(t.regression.samps))
#  )
#)

x = census$HousingUnits
y = census$TotalPop
mu_x = mean(census$HousingUnits)
mu_y = mean(census$TotalPop)
var_barx = (1-n/N)*var(census$HousingUnits)/n
cov_barxbary = (1-n/N)/n * sum((x-mu_x)*(y-mu_y))/(N-1)

(lin.bias = N*mu_y/mu_x^2 *var_barx - N/mu_x *cov_barxbary)

## part g (ughhh) ##
p = census$HousingUnits/sum(census$HousingUnits)
(sd_pps = sqrt(sum(p * (census$TotalPop/p - sum(census$TotalPop))^2)/n))
sprintf("%.4e",sd_pps)
