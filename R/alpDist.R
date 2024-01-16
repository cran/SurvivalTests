####################################################
alpDist <- function(n, m, nsim){

alpA <- betA <- rep(0,nsim)

for (j in 1:nsim) {
Wdat = sort(rexp(n)) # This is for sampling from Exp(1)

# Handle censored values
Wdat[Wdat>Wdat[m]] <- Wdat[m]
    
Zi <- log(Wdat)
alp0 <- pi* var(Zi)^(-1/2) /sqrt(6) # Approx MLE
alpout <- optimize(f, ys=Wdat, m=m, lower=alp0/3 , upper=2*alp0) # Alpha hat for Weibil(1,1)
alph <- alpout$minimum
    
# Prepare samaples to return
alpA[j] <- alph
betA[j] <- sum(Wdat^alph / m)^(1/alph) 
}

return(cbind(alpA,betA))
}
####################################################

