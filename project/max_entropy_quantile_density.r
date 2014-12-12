crab = read.csv('crabbieMRDS.csv') 
library('entropy')

# Scale and duplicate data to (0,1) centered at .5
scale_fun = function(x) { 
  x = x / (2*max(x))
  x = c(-x,x) + .5
}
x = scale_fun(crab$distance)

beta_cdf = function(x,a) { pbeta(x,a,a) } # Fit data to a symmetric beta distribution 
scale_beta_pdf = function(x,a,max_x) { # This scales beta to [0,max(x)]
  dbeta( x/(2*max_x) + .5,a,a) / max_x
}

bins = nclass.Sturges(x) # Use the Sturges algorithm for deciding the number of bins for the entropy estimator
cost = function(a) {1/entropy(discretize(beta_cdf(x,a),bins)) } # Define entropy based cost function

# Use nonlinear optimization to find optimal symmetric beta distribution 
# over alpha in [0,4]
opt_parm = optimize(f=cost,interval=c(0,4))

plot(function(x)scale_beta_pdf(x,opt_parm$minimum,max(crab$distance)),from=0,to=max(crab$distance),col='red')
hist(crab$distance,freq=F,add=T)


#crab2 = read.csv('crabbieMCDS.csv')

duck = read.csv('ducknests.csv')
x = scale_fun(duck$distance)
bins = nclass.Sturges(x) # Use the Sturges algorithm for deciding the number of bins for the entropy estimator
cost = function(a) {1/entropy(discretize(beta_cdf(x,a),bins)) } # Define entropy based cost function
opt_parm = optimize(f=cost,interval=c(0,4))

dev.new()
plot(function(x)scale_beta_pdf(x,opt_parm$minimum,max(duck$distance)),from=0,to=max(duck$distance),col='red')
hist(duck$distance,freq=F,add=T)

minke = read.csv('minke.csv')
x = scale_fun(minke$distance[!is.na(minke$distance)])
bins = nclass.Sturges(x) # Use the Sturges algorithm for deciding the number of bins for the entropy estimator
cost = function(a) {1/entropy(discretize(beta_cdf(x,a),bins)) } # Define entropy based cost function
opt_parm = optimize(f=cost,interval=c(0,4))

dev.new()
M = max(minke$distance[!is.na(minke$distance)])
plot(function(x)scale_beta_pdf(x,opt_parm$minimum,M),from=0,to=M,col='red')
hist(minke$distance,freq=F,add=T)


