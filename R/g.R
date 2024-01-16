####################################################
## Define g() function ##
g <- function(v,n1,n2,sig){
eu = -digamma(1) # Euler constant
vtmp = v/sig +eu
cdf =  pgamma(vtmp,n1,n2)

if (min(cdf)>.999999999999) cdf = .999999999999 * cdf 
if (max(cdf)< 0.00000000001) cdf = 0.00000000001 

return(qnorm(cdf)) 
}
####################################################

