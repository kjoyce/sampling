################################################################################
# Problem 5
################################################################################
clovers = read.csv('SweetCloverTransects.csv')
width = 108*12 # We take inches to be the units
y = tapply(clovers$width,clovers$transect,function(w) { sum(1/(w/width)) })
( v_st=mean(y) )
( se.v_st=sqrt(var(y)/length(y)) )
( ci.v_st = v_st + c(-1,1)*qnorm(.975)*se.v_st )

# part (b)
d = 100
(n_suff = (qnorm(.975)/d)^2*var(y) )

