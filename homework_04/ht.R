ht <- function(y, pi, p)
{
# Computes Horvitz-Thompson estimate of a population total and
#   estimated variance and SE; Sec. 6.2 of Thompson.
# y is the vector of values of the response variable
# pi is the vector of inclusion probabilities
# p is the matrix of joint inclusion probabilities
# only the values of p[i,j] for j > i need be filled in
# p will be a square matrix with number of rows equal to the
# length of the vector pi
	if(length(y) != length(pi)) stop("y and pi must be same length")
	tauhat <- sum(y/pi)
	# cat("tau.hat = ", tauhat, "\n")
	if(nrow(p) != ncol(p))
		stop("p must be square")
	if(nrow(p) != length(pi))
		stop("no. of rows of p must be same as length of pi")
	c1 <- sum((1/pi^2 - 1/pi) * y^2)
	c2 <- 0
	v <- length(pi)
	for(i in 1:(v - 1)) {
		for(j in (i + 1):v) {
			c2 <- c2 + 2 * ((1/(pi[i] * pi[j]) - 1/p[i, j]) * y[i] * y[j])
		}
	}
	var.tauhat <- c1 + c2
	# cat("Var(tauhat)= ", var.tauhat,"  SE=",sqrt(var.tauhat),"\n")
      list(tauhat=tauhat,SE=sqrt(var.tauhat))
}
