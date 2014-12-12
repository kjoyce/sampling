# Line-intercept sampling
# Desk example; random sample of 4 transects.
# Desks are 150 by 44 cm. and there are 16 of them. Baseline is 777 cm. Length of transects is 692 cm.
n <- 4  # number of transects  
transect <- rep(1:4,c(4,3,3,4))
w <- c(146,147,152,150,130.5,153,120,130,124,120,78,155.1,146.8,150.4) # widths of desks
k <- length(w)  # number of objects
y <- rep(1,k)  # response vector (here, to estimate no. of desks)
len <- c(56,47,49,44,69,29,43.5,68,33,76,97.5,45.5,55.2,44.5) # this is the length of the part of the transect intersecting the desk
desksize <- 150*44
A <- 777*692  # area of room

# Separate transects estimate of number of desks
p <- w/777  # probability of an object being intercepted by a single random transect
resp <- y/p
(v <- tapply(resp,transect,sum))
c(mean(v),sqrt(var(v)/length(v)))  # Estimate and SE of number of desks
c(desksize*mean(v)/A,(desksize/A)*sqrt(var(v)/length(v))) # estimate and SE of proportion of room covered by desks
16*desksize/A	# true proportion

# Horvitz-Thompson using only distinct desks
kappa <- 12  # number of distinct desks

# next section creates matrix of overlap widths
leftw <- c(59,126,61,19,33.5,142.5,19,14,18,74,39.5,45.4) # right widths
rightw <- c(87,21,91,131,97,10.5,101,110,60,81.1,107.3,105)  # left widths
pos <- c(rep(273,4),rep(432,3),454,rep(623,4))  # position of transects along walls
aleft <- pos-leftw  # absolute positions of left ends of desks
aright <- pos+rightw # absolute positions of right ends of desks
wk <- aright-aleft  # widths of 12 desks
wkh <- matrix(0,nrow=kappa,ncol=kappa) # matrix of overlaps 
for(i in 1:(kappa-1)){
  for(j in (i+1):kappa) wkh[i,j] <- max(min(aright[i],aright[j])-max(aleft[i],aleft[j]),0)}

pr <- probin(4,777,wk,wkh) # create inclusion and joint inclusion probabilities
pr
yy <- rep(1,kappa)  # to estimate number of desks
(h <- ht(yy,pr$pik,pr$pikh))  # H-T estimate of total

16*desksize/A   # true proportion of room covered by desks
h$tauhat*desksize/A  # H-T esimate of proportion
h$SE*desksize/A  # SE