################################################################################
# Problem 4
################################################################################
rm(list=ls())
N = 200
y = c( 0, 0, 0, 0, 0,
       0, 0, 0, 1, 1, 
       1, 1, 2, 2, 3, 
       3, 4, 5, 9, 10 )
n = length(y)

# part (a)
( tau.srs = N * mean(y) )
( se.tau.srs = N*sqrt( (1 - n/N) * var(y)/n ) )

# part (b)
x = y < 1 # x is true when if there is a decayed tooth
nh = table(x) # False comes first
Nh = c(200-60,60)
s2h = tapply(y,x,var) # estimate of sig^2_h
( tau.post = sum(tapply(y,x,mean) * Nh) ) 
var.term1 = N^2 * (1 - n/N)/n * sum( Nh/N * s2h ) 
var.term2 = N^2 * 1/n^2 * (N-n)/(N-1) * sum( (N-Nh)/N * s2h ) 
(se.tau.post.both.terms = sqrt(var.term1 + var.term2))
(se.tau.post.first.term = sqrt(var.term1))
