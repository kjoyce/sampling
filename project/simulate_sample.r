rm(list=ls())
tau = 500  # Population total
n = 4 # Number of transects 
systematic_length = 1/n # Number in systematic sample
sight_radius = .03 # distance from transects that gaurantees counting
buf = 1.3 # Buffer parameter for detectability function
pow = .65 # Power on 1-distance for detectability function

x = runif(tau) # Generate random coordinates
y = runif(tau)

transect_sample = function(x,y) {
#transects = runif(n) # Sample from the bottom of the square  # Randomly placed transects
  transects = (runif(1) + (1:n-1)*systematic_length) %% 1 # A systematic sample
  distance = sapply(transects,function(t)abs(t - x)) # Calculate the distance of every point to each transect  
  idx1 = distance < sight_radius # Include all points within sight_radius
  idx2 = (distance < buf*sight_radius) & # Detectability is random within buffer zone
	 (matrix(runif(tau*n),tau,n) < (1-distance/(buf*sight_radius))^pow) # This implements Bernouli trials
  idx = idx1 | idx2

  xx = matrix(rep(x,n),tau,n,byrow=F) # Replicate x to match dimensions of distance matrix
  yy = matrix(rep(y,n),tau,n,byrow=F) # Replicate y to match dimensions of distance matrix
  tt = matrix(rep(1:n,tau),tau,n,byrow=F) # Replicate transects to match dimensions of distance matrix

  samp = data.frame(
      x = as.vector(xx[idx]),
      y = as.vector(yy[idx]),
      transect = as.vector(tt[idx]), 
      transect_x = transects[as.vector(tt[idx])],
      distance = as.vector(distance[idx])
  )
}
samp = transect_sample(x,y)

# Plot the sample from the simulation
plot(x,y,main=paste('Simulated Systematic Sample n=',n,' transects'),cex.main=1.4)
transects = unique(samp$transect_x)
out = sapply(transects,function(x)abline(v=x,col='red'))
out = sapply(transects,function(x)abline(v=x-sight_radius,lty='dashed'))
out = sapply(transects,function(x)abline(v=x+sight_radius,lty='dashed'))
out = sapply(1:n,function(i)points(samp$x,samp$y,pch='x'))
#dev.copy(pdf,'simulated_sample.pdf',width=10.2,height=10)
#dev.off() 

dev.new()
source('estimator_functions.r')
compare_three_densities(samp$distance)
#dev.copy(pdf,'density_curve_fits.pdf')
#dev.off()

estimate_total = function(samp,f0) {
  yi = tapply(samp$distance,samp$transect,length)
  Li = tapply(samp$distance,samp$transect,function(x)1)
  Di = yi*f0/(2*Li)
  return( mean(Di) )  
}

f.hn = halfnormal.density(samp$distance)
f.beta = me_beta.density(samp$distance)
f.ker = kernel.density(samp$distance)

(tauhats = sapply(c('half-normal'=f.hn(0),'half-beta'=f.beta(0),'kernel'=f.ker[1]), function(f0)estimate_total(samp,f0)))

bootstrapper = function(x,y,pb,i) {
  setTxtProgressBar(pb,i)
  samp = transect_sample(x,y)
  f.hn = halfnormal.density(samp$distance)
  f.beta = me_beta.density(samp$distance)
  f.ker = kernel.density(samp$distance)
  sapply(c('half-normal'=f.hn(0),'half-beta'=f.beta(0),'kernel'=f.ker[1]), function(f0)estimate_total(samp,f0))
}

#N = 5e5 # This took a few hours with my PC
#pb = txtProgressBar(min=1,max=N,style=3)
#bootstraps = sapply(1:N,function(i)bootstrapper(x,y,pb,i))
#close(pb)
#save.image('500k_bootstrapped_simulations.RData')
load('500k_bootstrapped_simulations.RData')

dev.new()
par(mfrow = c(1, 3), mgp = c(2, 0.5, 0),mar = c(4, 3, 2, 0) + 0.1)
xl = c(min(bootstraps),max(bootstraps))
breaks = min(sapply(1:3,function(i)nclass.scott(bootstraps[i,])))
hist(bootstraps['half-normal',],
     xlab=expression(widehat(tau)[hn]), 
     xlim=xl,
     sub=bquote(paste(
	       plain(E),  widehat(tau)[hn]==.(round(mean(bootstraps['half-normal',]),2))," ",
	       plain(SD), widehat(tau)[hn]==.(round(sd(  bootstraps['half-normal',]),2))
	      )),
     main='Half Normal Bootstraps',
     cex.main = 1.3,
     cex.sub  = 1,
     cex.lab  = 1,
     breaks=breaks,
     freq=F
     )
abline(v=tau,col='red')
hist(bootstraps['kernel',],
     xlab=expression(widehat(tau)[k]), 
     xlim=xl,
     sub=bquote(paste(
	       plain(E),  widehat(tau)[k]==.(round(mean(bootstraps['kernel',]),2))," ",
	       plain(SD), widehat(tau)[k]==.(round(sd(  bootstraps['kernel',]),2))
	      )),
     main='Normal Kernel Bootstraps',
     cex.main = 1.3,
     cex.sub  = 1,
     cex.lab  = 1,
     breaks=breaks,
     freq=F
     )
abline(v=tau,col='red')
hist(bootstraps['half-beta',], 
     xlab=expression(widehat(tau)[hb]), 
     xlim=xl,
     sub=bquote(paste(
	       plain(E),  widehat(tau)[hb]==.(round(mean(bootstraps['half-beta',]),2))," ",
	       plain(SD), widehat(tau)[hb]==.(round(sd(  bootstraps['half-beta',]),2))
	      )),
     main='ME Beta Bootstraps',
     cex.main = 1.3,
     cex.sub  = 1,
     cex.lab  = 1,
     breaks=breaks,
     freq=F
     )
abline(v=tau,col='red')
#dev.copy(pdf,'500K_bootstrap_comparison.pdf',width=15,height=5)
#dev.off() 
