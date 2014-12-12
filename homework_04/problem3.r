################################################################################
# Problem 3
################################################################################
rm(list=ls())
N = 400
nhp = c(30,90)
nh = c(30,20)
ybar = c(20/30,4/20)
w = nhp/sum(nhp)
(ybar.d = sum(w*ybar))
sh2.over.nh = ybar*(1-ybar)/(nh-1)
(
  se.ybar.d = sqrt(
    (N-1)/N * sum( ((nhp-1)/(sum(nhp)-1) - (nh-1)/(N-1)) * (w*sh2.over.nh) )
    + (N - sum(nhp))/(N*(sum(nhp)-1)) * sum( nhp/sum(nhp) * (ybar - ybar.d)^2 )
  )
)

n = 50
(p = 24/n) 
(se.p = sqrt((1-n/N) * p*(1-p)/(n-1)) )
