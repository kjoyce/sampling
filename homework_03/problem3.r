#####################################
# Problem 3
#####################################
N = 4000
Nh = c(1800,700,1500)
Rh = c(80,120,200)
ch = c(20,25,35)

sigh = (Rh/6) # Using P(|X|< 3 sig) = .997
d = 5

total_n = function(wh) { # Use the general equation and just input wh
  sum( (Nh*sigh)^2/wh ) / ( (N*d/qnorm(.975))^2 + sum(Nh*sigh^2) )
}

equal_wh = rep(1/3,3)
equal_n = total_n(equal_wh)
equal_nh = equal_n * equal_wh

prop_wh = Nh/N
prop_n = total_n(prop_wh)
prop_nh = prop_n * prop_wh

equalcost_wh = Nh * sigh / sum( Nh * sigh )
equalcost_n = total_n(equalcost_wh)
equalcost_nh = equalcost_n * equalcost_wh

unequalcost_wh = Nh*sigh/sqrt(ch) / sum( Nh*sigh * sqrt(ch) )
cstar = total_n(unequalcost_wh)
unequalcost_nh = cstar*unequalcost_wh
unequalcost_n = sum(unequalcost_nh)

(t(data.frame(
'equal'      = c( round(equal_nh)      , sum(round(equal_nh))     , sum(round(equal_nh     )*ch)),
'prop'       = c( round(prop_nh)       , sum(round(prop_nh))      , sum(round(prop_nh      )*ch )),
'equalcost'  = c( round(equalcost_nh)  , sum(round(equalcost_nh)) , sum(round(equalcost_nh )*ch)),
'unequalcost'= c( round(unequalcost_nh), sum(round(unequalcost_n)), sum(round(unequalcost_nh)*ch))
)))
