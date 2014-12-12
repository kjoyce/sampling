################################################################################
# Problem 6
################################################################################
rm(list=ls())
# part (a)
N = 100
mooses = c( 0, 0, 0, 0,
	    1, 2, 6,11,
	    5, 1, 9, 3,
	    1, 0,10, 4,
	    7,22, 0, 0
	  )
n = length(mooses)
p = .89
(muhat = mean(mooses)/p)
(var.srs = (1/p^2)*((N-n)/N) * var(mooses)/n)
(var.det = (1/p^2)*((1-p)/N) * mean(mooses))
(se.muhat = sqrt(var.srs + var.det))

# part (b)
se.p = c(1,2,3,4,5,6,7,8,10) / 100 
var.phat = 1/p^2 * (mean(mooses)/p)^2 * se.p^2
se2.muhat = sqrt(var.srs + var.det + var.phat)
t(round(data.frame(
   'se.p'	= se.p,
   'se2.muhat' = se2.muhat
   ),3))

