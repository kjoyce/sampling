# Homework 2 Problem 2
## Part (b) ##
N = 320
n = 4
y = c(2,2,5,10)
p = c(1.2,1.2,0.2,0.5)/80

(mu_p = 1/N * 1/n * sum(y/p))

## Part (c) ##
(var_mu_p = 1/(N^2*n*(n-1)) * sum( (y/p - N*mu_p)^2 )) 
