systematic <- function(x, k)
{
# Computes exp. value and std. dev. of ratio estimator of pop. mean of x
# from systematic sample with interval k and random starting point
	N <- length(x)
	if(k > N)
		stop("k must be less than length of x")
	if(k < 1)
		stop("k must be >= 1")
	mn <- rep(0, k)
	for(a in 1:k)
		mn[a] <- mean(x[seq(a, N, k)])
	#cat("N=", N, "  Pop. mean = ", mean(x), "\n")
	#cat("Systematic with k=", k, ":\n  exp. value=", mean(mn), "  std. dev.=", 
	#	sqrt(var(mn)), "\n")
	return(
	  c('evalue'=mean(mn),'stddev'=sqrt(var(mn)),'n'=length(seq(a, N, k)))
	)
}
