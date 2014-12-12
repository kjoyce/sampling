probin <- function(n, b, wk, wkh)
{
# Computes probabilities of inclusion and joint inclusion probabilities
#   for line intercept sampling, SRS of n transects
#  (equations (6) & (8) on pp. 246-247).
# n = the number of transects
# b = the length of the baseline
# wk = the vector of object widths
# wkh = the matrix of object overlaps; only entries above the diagonal are used
	if(!is.vector(wk)) stop("wk must be a vector")
	if(!is.matrix(wkh))
		stop("wkh must be a matrix")
	K <- ncol(wkh)
	if(K != nrow(wkh))
		stop("wkh must be square")
	if(K != length(wk))
		stop("length of wk must be same as dimension of wkh")
	w <- wkh
	w[1, 1] <- wk[1]
	for(i in 2:K) {
		w[i, i] <- wk[i]
		for(j in 1:(i - 1)) {
			w[i, j] <- w[j, i]
		}
	}
	pk <- wk/b
	pik <- 1 - (1 - pk)^n
	j <- length(wk)
	A <- matrix(rep(1, j),ncol=1) %*% t(wk)
	uno <- matrix(1, nrow = j, ncol = j)
	P <- matrix(rep(1, j),ncol=1) %*% t(pik)
	pikh <- P + t(P) - 1 + (uno - (1/b) * (A + t(A) - w))^n
	list(pik = pik, pikh = pikh)
}

