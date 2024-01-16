####################################################
## Function to compute p-value
getP <- function(Tobs, sighat, s2, nn, alpi, Z, nM){

sqn <- sqrt(nn)
NN <- nrow(Z)

Tvec <- rep(0, nM)
ind <- sample(c(1:nM),NN)

## Compute p-value and  Count hits ##
for (i in 1:nM) {	
Zvec <- Z[ind[i],]	
Alpvec <- alpi[ind[i],]
Wt <- nn * Alpvec^2 /s2 
factor <- sighat/sqn  # This is due to XB=sig *Z/sqrt(n)
mubar <- sum(factor*Wt*Zvec)/sum(Wt)
Tvec[i] <- sum(Wt*(factor*Zvec - mubar)^2) 
}

pval = mean(Tvec>Tobs)
return(pval)
}# End of getP function	
####################################################

