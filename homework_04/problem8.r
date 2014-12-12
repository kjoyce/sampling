################################################################################
# Problem 8
################################################################################
rm(list=ls()) 
seedlings = data.frame(
'plot' = c(rep(1,5),rep(2,6),rep(3,6),rep(4,5),rep(5,5)),
'height' = c( 12 , 11 , 11 , 10 , 13 ,    
              10 ,  9 ,  7 ,  8 ,  8 , 10,
               6 ,  5 ,  7 ,  5 ,  6 ,  4,
               7 ,  8 ,  6 ,  7 ,  6 ,   
              10 , 11 , 13 , 12 , 12 ))   


N = 25			    # Total number of primary units
Nh = c(52,56,60,46,49)	    # Known secondary unit totals
nh = table(seedlings$plot)  # Secondary sample sizes
n = length(nh)		    # Primary sample size

# This function does a two stage resample
# where the first stage resamples without 
# replacement from bootstrap population
# made from replicating the sample (N/n) times
# to correct for a small sampling population
twostage.resample = function(y,factors,N,Nh) { 
  factors = factor(factors)		    # Make sure this is a factor level
  #nh = table(factors)			    # Strata sample sizes
  #n = length(nh)			    # Number of primary units sampled
  stageone.pop = rep(1:n,N/n)	            # Construct indices for bootsrp population 
  stageone.idx = sample(stageone.pop,N/n)   # Sample WITHOUT replacement
  ragged_table = tapply(y,factors,	    # Create a list of strata for indexing
		   function(yh){ yh })	    #   using the first stage samples indices
  stagetwo = function(yh) {		    # Function for sampling second stage
    sample(yh,replace=T)		    # Here we "could" implement a finite pop
  }					    # correction, but instead do the regular bootstrp
  yy = lapply(ragged_table[stageone.idx]    # Sample the second stage
	      , stagetwo)		    # Using the stagetwo function
  repidx = rep(stageone.idx,nh[stageone.idx])# Make a vector of indices for making resampled factor and Nh
  ffactors = levels(factors)[repidx]        # Make a vector of the factors sampled at stage one
  yy = unlist(yy,use.names=F)		    # Change yy from a ragged list back to a vector

  if(length(yy[1]) > length(ffactors) ) {
    browser()
  }
  samp = data.frame(
    'resampidx' = rep(1:n,nh[stageone.idx]),# Use this for the "new" strata index
    'factors'   = ffactors,		    # The factor (plot#) that was resampled
    'y'	        = yy,			    # The resampled value
    'Nh'        = Nh[repidx]		    # The true strata size matched up with resample
  )
  if(length(samp$y) != length(samp$resampidx)) {
    browser()
  }

  return(samp)
}

# This function is meant to work with
# the naming convention in twostage.resample
estimate.ybar = function(resampidx,Nh,y) {
  NNh = tapply(Nh,resampidx,mean)     # Using the fact that the mean of repeats is what is repeated
  ybarh = tapply(y,resampidx,mean)	     # Strata wise means
  sum(NNh*ybarh)/sum(NNh)		     # Ratio estimator
}

boot = function() {
  samp = twostage.resample(seedlings$height,seedlings$plot,N,Nh)
  estimate.ybar(samp$resampidx,samp$Nh,samp$y)
}
straps = replicate(10000,boot())

mur = estimate.ybar(seedlings$plot,rep(Nh,nh),seedlings$height)
hist(straps, main="Bootstrap Density of Ratio Estimator for Mean", xlab="mu_r")
abline(v = mean(straps), lty='solid', col='red')
abline(v = mur, lty='dashed', col='red')
dev.copy(pdf,'bootstrap.pdf')
dev.off()
(muboot = mean(straps))
(seboot = sd(straps))
sorted_straps = sort(straps)
(boot.ci = sorted_straps[c(50,9750)])
