n = 5
y = c(2,8,4,8,3)

p = 1/625*c(5,28,12,14,13) # selection probabilities
pp = 1-(1-p)**n # inclusion probabilities

# HH estimator
tau.p = mean(y/p)
var.tau.p = 1/n * var(y/p)
se.tau.p = sqrt(var.tau.p)

# HT estimator
tau.pi = n*mean(y/pp)

# Form the joint inclusion probabilities
J = outer(1:5,1:5,function(i,j){ pp[i] + pp[j] - (1 - (1-p[i] - p[j])**n) }) 
diag(J) = pp

# Form product of inclusion probabilities
H = outer(pp,pp)

# The variance estimate is also a quadratic form
var.tau.pi = t(y) %*% ((J - H)/(J*H)) %*% y
se.tau.pi = sqrt(tau.pi)


