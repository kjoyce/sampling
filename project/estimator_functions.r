library('entropy')

# This function takes a list of positive distances
# and uses quantil maximum entropy to fit a
# half symmetric beta disribution.  The return value
# density is the density function
me_beta.density = function(d,breaks = NA) {
  # These are helper functions to beta_max_entropy
  beta_cdf = function(x,a) { pbeta(x,a,a) } # Fit data to a symmetric beta distribution 
  scale_beta_pdf = function(s,a,max_s) { # This scales beta to [0,max(x)]
    dbeta( s/(2*max_s) + .5,a,a) / max_s
  }
  x = d / (2*max(d))
  x = c(-x,x) + .5
  if(is.na(breaks)) breaks = nclass.Sturges(x) # Use the Sturges algorithm for deciding the number of breaks for the entropy estimator
  cost = function(a) {1/entropy(discretize(beta_cdf(x,a),breaks)) } # Define 1/entropy based cost function
  opt_param = optimize(f=cost,interval=c(0,6)) # values over 5 sometimes lead to Inf values
  return( function(s) scale_beta_pdf(s,opt_param$minimum,max(d)) )
}

# These functions were culled from DistanceFunctions.R
halfnormal.density = function(x) {
  what.hn <- sqrt((pi/2) * mean(x^2))	# effective half-width for half-normal
  return( function(xx) exp(( - pi * xx^2)/(4 * what.hn^2))/what.hn )# half-normal fit
}

# h is the window width for the normal kernel method.
# If h is not given then formula (17.12) is used to calculate a window
# width from the data: h=0.9ay^(-1/5) where a = min(s,IQR/1.34).
# kernel.density returns a discretization of the density with n points (default 200)
kernel.density = function(x,h=NA,n=200) {
  y <- length(x)
  xx <- seq(0, max(x), length = n)
  a <- min(sd(x), IQR(x)/1.34)	
  if(is.na(h))
	  h <- (0.9 * a)/y^0.2	# default bandwidth
  f.ker <- rep(0, n)
  for(i in 1:n) f.ker[i] <- f.ker[i] <- (1/(y * h)) * (1/sqrt(2 * pi)) * sum(exp(-0.5 *
		 ((xx[i] - x)/h)^2) + exp(-0.5 * ((xx[i] + x)/h)^2))
  return( f.ker )
}

# This plots densities overlayed on a histogram of the data
# and is also culled from DistanceFunctions.r
compare_three_densities = function(x,h=NA,n=200,breaks=NA) {
  if(is.na(breaks)) { breaks = nclass.Sturges(x) }  
  xx = seq(0,max(x),length=n)
  f.hn = halfnormal.density(x)
  f.beta = me_beta.density(x)
  f.ker = kernel.density(x,h=h,n=n)
  if(is.na(h)){
    a <- min(sd(x), IQR(x)/1.34)	
    h <- (0.9 * a)/length(x)^0.2	# default bandwidth
  }

  oldpar <- par(mfrow = c(2, 2), mgp = c(2, 0.5, 0),mar = c(4, 3, 2, 0) + 0.1,oma=c(0,0,1.1,0))
  ymax <- max(c(f.beta(xx)[!is.infinite(f.beta(xx))], f.hn(xx), f.ker, hist(x,plot = F,breaks = breaks)$density))
  hist(x + 0.0001 * max(x), col = 0, breaks = breaks, freq = F,
      ylim = c(0, ymax), main = paste("Raw data: y = ", length(x),
       " points"), xlab = "Distance")
  hist(x + 0.0001 * max(x), col = 0, breaks = breaks, freq = F,
       ylim = c(0, ymax), main = "Half-normal fit", sub = paste(
       "f(0)=", round(f.hn(0), 6), "    w=", round(1/f.hn(0),
       3)), xlab = "Distance")
  lines(xx, f.hn(xx),col='red')
  hist(x + 0.0001 * max(x), col = 0, breaks = breaks, freq = F,
       ylim = c(0, ymax), main = "Normal kernel fit", sub =
       paste("f(0)=", round(f.ker[1], 6), "    w=", round(
       1/f.ker[1], 3), "    h=", round(h, 4)), xlab = "Distance")
  lines(xx, f.ker,col='red')
  hist(x + 0.0001 * max(x), col = 0, breaks = breaks, freq = F,
       ylim = c(0, ymax), main = "Half-beta fit", sub =
       paste("f(0)=", round(f.beta(0), 6), "    w=", round(
       1/f.beta(0), 3)), xlab = "Distance")
  lines(xx, f.beta(xx),col='red')
}
