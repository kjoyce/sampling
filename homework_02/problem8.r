# Small population example N = 3 n = 2

tauhat.p.sd = function(y,p) { 
  P = outer(p,p)
  taus = outer(y/p,y/p,'+')/2
  Etau = sum( P * taus )
  (VARtau = sum( P*(taus - Etau)^2 ))
  sqrt(abs(VARtau))
}
 
tauhat.pi.sd = function(y,p) { 
  ppi = 1-(1-p)^2
  P = outer(p,p)
  taus = outer(y/ppi,y/ppi,'+')
  diag(taus) = diag(taus)/2 # only add diag once
  Etau = sum(P*taus)
  (VARtau = sum( P*(taus - Etau)^2 ))
  sqrt(abs(VARtau))
}

diff_fun = function(param) {
# The arguments force the p's to sum to 1 
  y = param[1:3]
  p = c(param[4:5],1-sum(param[4:5]))
  (tauhat.pi.sd(y,p) - tauhat.p.sd(y,p))/tauhat.pi.sd(y,p)
}

# Minimize the difference of the estimators with box contraints (not quite right)
param = c(rep(5,3),rep(.2,2)) # Guess for initial params
param = c(1,2,3,.5,.1) # Another guess for initial params
(optim.out = optim(param,diff_fun,lower=rep(1e-4,5),upper=c(rep(10,3),rep(1,2))))
# The constraints are more complicated that box-type, and I just used this to get a ball park of what happens

# In much like the ``circus statistician'' we find that if one observation has lots of weight and is very likely, then the tau_p estimator is terrible, but tau_p is not so bad
y = c(.01,.01,50000)
p = c(.00005,.00005,.9999)
(ans = c(
  tauhat.pi.sd = tauhat.pi.sd(y,p),
  tauhat.p.sd = tauhat.p.sd(y,p)
))

# Build table
I = outer(c(1,2,3),c(1,1,1))
P = outer(p,p)
Y = outer(y,rep(1,3)) 
P = outer(p,p)
ppi = 1-(1-p)^2
tau.pis = outer(y/ppi,y/ppi,'+')
diag(tau.pis) = diag(tau.pis)/2 # only add diag once
tau.ps = outer(y/p,y/p,'+')/2

tab = cbind(as.vector(t(I)),as.vector(I), as.vector(P), as.vector(Y), as.vector(t(Y)), as.vector(tau.ps), as.vector(tau.pis))

tab = rbind(tab,c(0,0,0,0,0,sum(tau.ps*P),sum(tau.pis*P)))
(tab = rbind(tab,c(0,0,0,0,0,sqrt(sum((tau.ps-sum(y))^2*P)),sqrt(sum((tau.pis-sum(y))^2*P)))))

