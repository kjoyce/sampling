######################################
# Problem 2
######################################
otters = read.csv('otters.csv') 

nh = table(otters$Habitat) 
Nh = c('1'=89,'2'=61,'3'=40,'4'=47)
ybarh = tapply(otters$Holts,otters$Habitat,mean)  
sh2 = tapply(otters$Holts,otters$Habitat,var)

(tauhat = sum(Nh*ybarh))
(se_tauhat = sqrt(sum(Nh^2*(1-nh/Nh)*sh2/nh)))

ah = Nh*(Nh-nh)/nh
satterthwaite_df = sum(ah*sh2)^2/sum((ah*sh2)^2/nh-1)
(tauhat_ci = tauhat + c(-1,1)*se_tauhat*qnorm(.975))
(tauhat_sat_ci = tauhat + c(-1,1)*se_tauhat*qt(.975,satterthwaite_df))
