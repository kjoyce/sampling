y = c(2,26,3, 1,
      1, 1,14,15, 1,1,
      3,19, 1,4, 1,25, 2
     )
net = c( rep(1,4),
	 rep(2,6),
	 rep(3,7))

n = tapply(y,net,sum)
a = c(2/5,3/5,4/5)

idx = list()
idx[[1]] = c(1,2)
idx[[2]] = c(1,2,3)
idx[[3]] = c(3)
idx[[4]] = c(3)
idx[[5]] = c(2,3)

nu.i = c(8,15,13)
tau.pi = sapply(idx,function(i)sum((n/a)[i])) 
nu     = sapply(idx,function(i)sum(nu.i[i]))
(mu.tau.pi = mean(tau.pi))
(sigma.tau.pi = sqrt(sum( (tau.pi - mu.tau.pi)^2 )/5) )

(mu.nu = mean(nu))
(sigma.nu = sqrt(sum( (nu - mu.nu)^2 )/5) )

N = 100
sigma2 = var(c(y,rep(0,N-length(y))))
(sd.tau.srs = sqrt(N*(N- mu.nu)*sigma2/mu.nu))


