w = c(); p = c(); var.p = c();
w[1] = .5
p[1] = .3
var.p[1] = w[1]^2*(p[1]*(1-p[1]))/100

w[2] = .3
p[2] = .42
var.p[2] = .08^2

w[3] = .2
p[3] = .1
var.p[3] = .04^2

(phat = sum(p*w))
(sdhat = sqrt(sum(var.p)))
