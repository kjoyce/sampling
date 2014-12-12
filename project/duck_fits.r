source('estimator_functions.r')
duck = read.csv('ducknests.csv')
crab = read.csv('crabbieMCDS.csv') 
compare_three_densities(duck$distance)
title('Duck nests distances (Anderson Pospahala 1970)',outer=T)
#dev.copy(pdf,'duck_nest_fits.pdf',width=7,height=7)
#dev.off()

scale_fun = function(x) { 
  x = x / (2*max(x))
  x = c(-x,x) + .5
}
x = scale_fun(crab$distance)

beta_cdf = function(x,a) { pbeta(x,a,a) } # Fit data to a symmetric beta distribution 
scale_beta_pdf = function(x,a,max_x) { # This scales beta to [0,max(x)]
  dbeta( x/(2*max_x) + .5,a,a) / max_x
}

cost = function(a) {1/entropy(discretize(beta_cdf(x,a),bins)) } # Define entropy based cost function
opt_param = optimize(f=cost,interval=c(0,4))

# Example of quantile function
f = function(x)pbeta(x,opt_param$minimum,opt_param$minimum)
plot(f,
     ylab = expression(P(X <= x)),
     xlab = expression(x),
     main = "Cumulative Distribuion Function",
     axes = F
     ) 
axis(side=1, at=c(0,.6,1), labels=c(0,expression(x[i]),1))
axis(side=2, at=c(0,f(.6),1), labels=c(0,expression(F(x[i])),1))
idx = seq(1,length(x),4)
rug(x[idx])
rug(pbeta(x[idx],opt_param$minimum,opt_param$minimum),side=2)
segments(c(.6,.6), c(0,f(.6)), c(.6,-.1), c(f(.6),f(.6)))
dev.copy(pdf,'cdf_example.pdf',width=5,height=5)
dev.off()
