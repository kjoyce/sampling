#####################################
# Problem 4
#####################################
seedlings = data.frame(
'plot' = c(rep(1,5),rep(2,6),rep(3,6),rep(4,5),rep(5,5)),
'height' = c( 12 , 11 , 11 , 10 , 13 ,    
              10 ,  9 ,  7 ,  8 ,  8 , 10,
               6 ,  5 ,  7 ,  5 ,  6 ,  4,
               7 ,  8 ,  6 ,  7 ,  6 ,   
              10 , 11 , 13 , 12 , 12 ))   

Mi = c(52,56,60,46,49)
mi = table(seedlings$plot)
hatyi = Mi*tapply(seedlings$height,seedlings$plot,mean)
(mur = sum(hatyi) / sum(Mi))

n = 5
mbar = mean(Mi)
N = 25
si2 = tapply(seedlings$height,seedlings$plot,var)

(varbetween = (1 - n/N) * 1/(mbar^2 * n) * 1/(n -1) * sum((hatyi - Mi*mur)^2) )
(varwithin = 1/(N*mbar^2*n) *sum(Mi*(Mi - mi)*si2/mi))
(se_mur = sqrt(varbetween + varwithin))

(mur_ci = mur + c(-1,1)*se_mur*qt(.975,n-1))
