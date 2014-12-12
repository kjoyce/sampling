distance.fit <- function(x, breaks = NA, h = NA)
{
# Computes estimate of f(x) (pdf of observed distances for transect data)
# by three methods: "exp","hnorm","kernel" (Sec. 17.4 and 17.5 of Thompson)
# and shows graphs of fits for all three to histogram of the raw data.
# "exp" is exponential, "hnorm" is half-normal, and "kernel" is kernel
# estimate with normal kernel.
# x is vector of observed distances,
# breaks is number of bars in histogram (may be omitted),
# h is the window width for the normal kernel method.
# If h is not given then formula (17.12) is used to calculate a window
# width from the data: h=0.9ay^(-1/5) where a = min(s,IQR/1.34).
#
# This function is modified to also implement beta max entropy entropy estimator
        n <- 200        # no. of points at which to plot curves
        oldpar <- par(mfrow = c(2, 2), mgp = c(2, 0.5, 0),mar = c(4, 3, 2, 0) + 0.1)
        y <- length(x)
        xx <- seq(0, max(x), length = n)
        what.exp <- mean(x)	# effective half-width for exponential
        f.exp <- exp( - xx/what.exp)/what.exp	# exponential fit
        what.hn <- sqrt((pi/2) * mean(x^2))	# effective half-width for half-normal
        f.hn <- exp(( - pi * xx^2)/(4 * what.hn^2))/what.hn	# half-normal fit
        a <- min(sd(x), IQR(x)/1.34)	
        if(is.na(h))
                h <- (0.9 * a)/y^0.2	# default bandwidth
        f.ker <- rep(0, n)
        for(i in 1:n) f.ker[i] <- f.ker[i] <- (1/(y * h)) * (1/sqrt(2 * pi)) * sum(exp(-0.5 *
                       ((xx[i] - x)/h)^2) + exp(-0.5 * ((xx[i] + x)/h)^2))
        
	  what.ker <- 1/f.ker[1]
	######	
	f.beta_fun = beta_max_entropy_density(x)
	f.beta = f.beta_fun(xx)
	what.beta = 1/f.beta_fun(0)
	######
        if(is.na(breaks)) { breaks = nclass.Sturges(x) }
	ymax <- max(c(f.exp, f.hn, f.ker, hist(x,plot = F,breaks = breaks)$density))
	#hist(x + 0.0001 * max(x), col = 0, breaks = breaks, freq = F,
	#    ylim = c(0, ymax), main = paste("Raw data: y = ", y,
	#     " points"), xlab = "Distance")
	hist(x + 0.0001 * max(x), col = 0, breaks = breaks, freq = F,
	     ylim = c(0, ymax), main = "Exponential fit", sub = paste(
	     "f(0)=", round(1/what.exp, 6), "    w=", round(what.exp,
	     3)), xlab = "Distance")
	lines(xx, f.exp)
	hist(x + 0.0001 * max(x), col = 0, breaks = breaks, freq = F,
	     ylim = c(0, ymax), main = "Half-normal fit", sub = paste(
	     "f(0)=", round(1/what.hn, 6), "    w=", round(what.hn,
	     3)), xlab = "Distance")
	lines(xx, f.hn)
	hist(x + 0.0001 * max(x), col = 0, breaks = breaks, freq = F,
	     ylim = c(0, ymax), main = "Normal kernel fit", sub =
	     paste("f(0)=", round(1/what.ker, 6), "    w=", round(
	     what.ker, 3), "    h=", round(h, 4)), xlab = "Distance")
	lines(xx, f.ker)
	hist(x + 0.0001 * max(x), col = 0, breaks = breaks, freq = F,
	     ylim = c(0, ymax), main = "Half-beta fit", sub =
	     paste("f(0)=", round(1/what.beta, 6), "    w=", round(
	     what.beta, 3), "    h=", round(h, 4)), xlab = "Distance")
	lines(xx, f.beta)

        cat("No. of distances=", y)
        cat("\n  Exponential: f(0)_hat=", round(1/what.exp, 8), " w_hat=",
            round(what.exp, 8))
        cat("\n  Half-normal: f(0)_hat=", round(1/what.hn, 8), " w_hat=",
            round(what.hn, 8))
        cat("\n  Kernel: f(0)_hat=", round(1/what.ker, 8), " w_hat=",
            round(what.ker, 8), "  Window width: h=", round(h, 8))
        cat("\n")
        par(oldpar)
        invisible()
	return(c(
	    'exponential' = 1/what.exp,
	    'half-normal' = 1/what.hn, 
	    'kernel'      = 1/what.ker,
	    'half-beta'   = 1/what.beta
	))
}

distance.boot <- function(x, l, m, method, h = NA)
{
# Compute ratio estimator of density for transect data (Sec. 17.7 of Thompson) along
# with SE from formula on p. 214 and bootstrap SE.
# x is a list of n vectors of distances at which objects were
# observed on n random transects and l is a vector of
# lengths of the transects.
# m is the number of bootstrap replications.  method is the
# method used to estimate the effective half-width:
# "exp" means exponential detectability function,
# "hnorm" means half-normal and "kernel" means kernel method with
# normal kernel. h is the window width for the normal kernel method;
# if h is not given then eq. (12), p. 208 of Thompson is used.
        n <- length(x)
        if(n < 2)
                stop("x must be a list of length >= 2")
        if(length(l) != n)
                stop("length of l must be number of transects")
        if(is.na(match(method, c("exp", "hnorm", "kernel"))))
                stop("method not valid")
        yi <- unlist(lapply(x, length)) # vector of no. of detections
# First compute density estimate from original data
        ux <- unlist(x)
        if(method == "exp")
                f0 <- 1/mean(ux)
        else if(method == "hnorm")
                f0 <- 1/sqrt((pi/2) * mean(ux^2))
        else {
                if(is.na(h))
                  h <- (0.9 * min(sd(ux), IQR(ux)/1.34))/
                       length(ux)^0.2
                f0 <- (2/(h * sqrt(2 * pi))) * mean(exp(-0.5 * (ux/h)^2))
        }
        Dhat <- (length(ux) * f0)/(2 * sum(l))
        se.Dhat <- sqrt(sum(((yi * f0)/2 - Dhat * l)^2)/(n * (n - 1) *
                   mean(l)^2))
        cat("Ratio estimator: Dhat=", round(Dhat, 10),
            "     SE (Sec. 17.7, Var_1)=", round(se.Dhat, 10),
            "\n") # Bootstrap
        dhat <- rep(0, m)
        for(i in 1:m) {
                z <- sample(n, size = n, replace = T)
                xz <- unlist(x[z])
                if(method == "exp")
                        f0 <- 1/mean(xz)
                else if(method == "hnorm")
                        f0 <- 1/sqrt((pi/2) * mean(xz^2))
                else {
                        if(is.na(h))
                          h <- (0.9 * min(sd(xz), IQR(xz)/1.34))/i
                               length(xz)^0.2
                        f0 <- (2/(h * sqrt(2 * pi))) * mean(exp(-0.5 *
                              (xz/h)^2))
                }
                dhat[i] <- (f0 * length(xz))/(2 * sum(l[z]))
        }
        cat("Bootstrap SE=", round(sqrt(var(dhat)), 10),
            ",  Estimated bias=",
            mean(dhat) - Dhat, "\n")
        dhat
}
