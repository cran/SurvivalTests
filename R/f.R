####################################################
## Function to find MLE of alpha ##
f <- function(x, ys, m) { # ys data is supposed to be Weibull
lys <- log(ys)
f <- 1/x - sum(ys^x * lys)/sum(ys^x) + sum(lys[1:m])/m # This is to be solved, but numerically unstable
return(f^2) # This is to be minimized 
}
####################################################


