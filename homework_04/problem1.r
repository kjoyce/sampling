################################################################################
# Problem 1
################################################################################
rm(list=ls())
N = 80

out1 = data.frame(row.names=c('n','tauhat','se.tauhat'))
# SRS: Zor
y <- c(47,41,60,58,60,41,112,67,91,50,75,25,77,49,72,98,35,35,47,52)
n = length(y)
tauhat = N * mean(y)
se.tauhat = N*sqrt((1-n/N)*var(y)/n)
out1['zor'] = c(n,tauhat,se.tauhat)

# SRS: Vernon
y <- c(38,47,57,44,41,31,17,62,35,60,35,50,60,47,24,49,73,35,91,75,56,98)
n = length(y)
tauhat = N * mean(y)
se.tauhat = N*sqrt((1-n/N)*var(y)/n)
out1['vernon'] = c(n,tauhat,se.tauhat)
(round(t(out1),2))

ratio_estimate = function(x,y) { 
  n <- length(y) 
  tau.x <- sum(x)
  x = x[1:n]
  r <- sum(y)/sum(x)
  tau.hat.r <- r*tau.x
  sr2 <- (1/(n-1))*sum((y-r*x)^2)
  var.tau.hat.r <- N^2*(1-n/N)*sr2/n
  se.tau.hat.r = sqrt(var.tau.hat.r)
  return(
    c(
      'tau.hat.r'=tau.hat.r,
      'se.tau.hat.r'=se.tau.hat.r
    )
  )
}

double_sample_estimate = function(x,y) { 
  n <- length(y) 
  x1 = x[1:n]
  np <- length(x)
  r <- sum(y)/sum(x1)
  tau.hat.x <- (N/np)*sum(x)
  tau.hat.r <- r*tau.hat.x
  sr2 <- (1/(n-1))*sum((y-r*x1)^2)
  var.tau.hat.r <- (N*(N-np)*var(y))/np + N^2*((np-n)/np)*sr2/n
  se.tau.hat.r = sqrt(var.tau.hat.r)
  return(
    c(
      'tau.hat.r'=tau.hat.r,
      'se.tau.hat.r'=se.tau.hat.r
    )
  )
}

out = data.frame(row.names=c('n','np','tau.hat.r.','se.tau.hat.r','all.tau.hat.r','all.se.tau.hat.r'))
# Double sampling, every 2nd: Eli and Guedem
x <- c(36,18,52,34,61,11,12,53,28,21,17,31,29,48,25,77,34,28,27,18,11,20,21,27,17,15,27,41,
       32,18,23,11,28,35,16,39,26,46,31,33,67,31)
y <- c(44,35,62,56,102,24,24,67,47,39,35,48,49,47,46,93,54,49,35,25,19)
out['eli_guedem'] = c(length(y),length(x),double_sample_estimate(x,y),NaN,NaN)

# Double sampling, every 2nd: Kevin and Lia
x <- c(60,20,30,45,20,30,23,35,67,27,35,51,10,38,37,30,36,29,19,18,30,25,30,35,50,32,50,71,
       18,38,37,32,51,59,28,41,38,29,41,33)
y <- c(72,31,44,49,25,45,35,52,102,35,52,75,17,58,49,60,50,47,25,17)
x.all <- c(x,35,28,21,32,28,31,28,22,18,14,54,45,31,21,30,43,45,36,32,68,19,40,29,28,25,61,
             50,35,18,20,30,38,40,45,30,59,40,30,20,45)
out['kevin_lia'] = c(length(y),length(x),double_sample_estimate(x,y),ratio_estimate(x.all,y))

# Double sampling, every 4th: Caitlin and Mike
x <- c(58,14,44,36,35,50,28,32,24,74,38,40,34,34,39,41,30,32,28,54,24,47,11,42,72,33,57,30,
       44,41,25,46,22,80,32,35,54,61,27,42,46,27,31,25,43,46,40,22,29,41,30,38,24,53,39,33,
       45,13,67,15,34,37,42,46,53,28,15,53)
y <- c(83,17,58,49,35,60,32,46,28,91,67,47,37,41,56,50,46)
x.all <- c(x,24,80,19,24,43,37,63,48,36,63,47,53)
out['caitlin_mike'] = c(length(y),length(x),double_sample_estimate(x,y),ratio_estimate(x.all,y))

# Double sampling, every 4th: Ann and Jacob
x <- c(45,50,35,35,45,30,80,60,15,80,20,50,70,45,25,70,90,90,20,35,15,35,45,20,40,20,60,30,60,
       40,20,45,65,20,45,10,30,50,55,50,35,30,15,40,50,30,20,35,65,40,30,40,40,45,40,70,30,35,
       35,90,40,25,35,80,45)
y <- c(35,52,35,39,54,28,93,60,24,83,30,50,77,67,35,72,102)
x.all <- c(x,35,60,40,30,40,50,20,15,13,25,70,50,12,45,80)
out['ann_jacob'] = c(length(y),length(x),double_sample_estimate(x,y),ratio_estimate(x.all,y))

# Double sampling, every 5th: Grant and Jon
x <- c(55,75,25,65,35,55,75,80,160,70,20,85,30,35,25,45,20,60,20,30,40,30,25,25,15,30,40,25,60,
       35,45,30,50,25,75,30,45,120,80,40,40,50,80,45,55,60,65,50,25,60,45)
y <- c(56,102,28,75,41,50,83,91,112,60,17)
x.all <- c(x,50,60,20,40,50,35,60,70,35,25,100,50,80,35,35,65,65,90,60,20,35,60,70,15,45,50,60,65,80)
out['grant_jon'] =  c(length(y),length(x),double_sample_estimate(x,y),ratio_estimate(x.all,y))

(round(t(out),2))

# part (e)
# Use Caitlin and Mike
# Double sampling, every 4th: Caitlin and Mike
x <- c(58,14,44,36,35,50,28,32,24,74,38,40,34,34,39,41,30,32,28,54,24,47,11,42,72,33,57,30,
       44,41,25,46,22,80,32,35,54,61,27,42,46,27,31,25,43,46,40,22,29,41,30,38,24,53,39,33,
       45,13,67,15,34,37,42,46,53,28,15,53)
y <- c(83,17,58,49,35,60,32,46,28,91,67,47,37,41,56,50,46)
cost_ratios = c(1/5,1/3,1/2)
n = length(y)
x = x[1:n]
r = sum(y)/sum(x)
sr2 <- (1/(n-1))*sum((y-r*x)^2)
s2 <- var(y)
(sample_ratio = sqrt(cost_ratios * sr2/(s2-sr2)))
