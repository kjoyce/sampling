# density = beta_max_entropy(d)
# This function takes a list of positive distances
# and uses quantil maximum entropy to fit a
# half symmetric beta disribution.  The return value
# density is the density function
library('entropy')
beta_max_entropy_density = function(d,bins = NA) {
# These are helper functions to beta_max_entropy
  scale_fun = function(x) { 
    x = x / (2*max(x))
    x = c(-x,x) + .5
  } 
  beta_cdf = function(x,a) { pbeta(x,a,a) } # Fit data to a symmetric beta distribution 
  scale_beta_pdf = function(x,a,max_x) { # This scales beta to [0,max(x)]
    dbeta( x/(2*max_x) + .5,a,a) / max_x
  }
  if(is.na(bins)) bins = nclass.Sturges(x) # Use the Sturges algorithm for deciding the number of bins for the entropy estimator

  x = scale_fun(d)
  cost = function(a) {1/entropy(discretize(beta_cdf(x,a),bins)) } # Define 1/entropy based cost function
  opt_param = optimize(f=cost,interval=c(0,10))
  return( function(x) scale_beta_pdf(x,opt_param$minimum,max(d)) )
}

