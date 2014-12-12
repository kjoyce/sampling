# Homework 2, problem #1

## Part (a) ##
N = 13828
n = 250
owned = 158

(ybar = owned/n)
s2 = n/(n-1) * ybar * (1-ybar) 

(se.ybar = sqrt( (1-n/N) * s2/(n-1) ) )
(ybar.ci = ybar + c(-1,1)*qt(.995,n-1)*se.ybar) # This is the normal approximation confidence interval

# The following codes produce a hypergeometric based confidence interval.
tauhat = round(N*ybar) # We actually form a confidence interval for the total number of owned households, tauhat
x = seq(1,n,1)
# This loop indicates that at these levels of N and n, various tau's near tauhat produce a hypergeometric density that is roughly symmetric.
for (i in seq(-500,500,100)) {
  plot(x,dhyper(x,tauhat+i,N-tauhat-i,n)) 
  title(sprintf('tau=%d',tauhat+i))
  print('Strike Anykey')
  readline()
}
alpha1 = alpha2 = .005 # Hence we choose alpha1 = alpha2

taus = seq(1,N,1) # Make all possible integer taus
left_tails = phyper(owned,taus,N-taus,n) # form all cumulative densities for possible tau values
right_tails = 1-left_tails + dhyper(owned,taus,N-taus,n) # form right tails from left ones
idx = which( left_tails > alpha1 & right_tails > alpha2 ) # find taus that satisfy Thompson 5.2
(ybar.hyper_ci = c(min(taus[idx]), max(taus[idx]))/N) # The hypergeometric based CI


## Part (b) ##
d = 0.03
(n0 = qnorm(.995)^2 * ybar*(1-ybar) / d^2)
(n_necessary = (N^-1 + n0^-1)^-1)
(n_max = (N^-1 + 4*d^2/qnorm(.995)^2)^-1)

