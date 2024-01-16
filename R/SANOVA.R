####################################################
SANOVA <- function(formula, data, nM = 5000, seed = 123, alpha = 0.05, na.rm = TRUE, verbose = TRUE){

set.seed(seed)
data <- model.frame(formula, data)
dp <- as.character(formula)

METHOD <- "Generalized Test for Survival ANOVA"

coln <- colnames(data)
DNAME <- paste(coln[1], "and", coln[2], "| status =", coln[3])

if (na.rm) {
        completeObs <- complete.cases(data)
        data <- data[completeObs, ]
    }

Time = data[,1]
Group = data[,2]
Mind =  data[,3]


if (!is.factor(Group)) stop("The group variable must be a factor.")
if (!(is.numeric(Mind)|is.integer(Mind))) stop("The status variable must be a numeric or an integer.")
if (!is.numeric(Time)) stop("The time variable must be a numeric.")


if(!all((names(table(data$Mind))) == c("0","1"))) stop("The status variable must be a numeric or an integer including 0 and 1. 0 must be censored, 1 must be non-censored.")

Mind <- 1 - Mind

data <- data.frame(Time,Group,Mind) 	


n0 <- table(data$Group)
grps <- names(n0)
n <- as.numeric(n0)
k <- length(n)
ym <- yv <- ydel <- gym <- gvec <- XB <- alphat <- M <- rep(0,k)

k <- length(n0)
eu <- -digamma(1) # Euler constant
    
for (i in 1:k) {
Dati <-  data[data$Group==grps[i], ]  
n[i] <- nrow(Dati)  
M[i] <- sum(1-Dati$Mind) # Mind=1 if censores    
Yi <- Dati$Time   # $Time 
lnYi <- log(Yi) # Testing location invariance
Dati$logTime <- lnYi
alp0 <- pi* var(lnYi)^(-1/2)/sqrt(6) # Approx MLE
ot <- optimize(f, ys=Yi, m=M[i] , lower=alp0/3 , upper=2*alp0) 
alphat[i] <- ot$minimum # estimate is minumum 
ym[i] <- mean(Dati$logTime) # Sample means for k groups
yv[i] <- var(Dati$logTime) # Sample variances for k groups    
}

sd <- sqrt(yv)

sighat <- sigma <- 1/alphat
nn <- M
s20 <- sighat^2 #Use this for mean correction only
wt0 <- nn/s20 # weight is simple weighted mean
sqn <- sqrt(nn)
factor <- sighat/sqn

scale <- ym/gamma(1+1/alphat)
muhat <- sum(wt0*ym)/sum(wt0)

eu1 <- 1-eu # mean of standard Gumbel
sd1 <- sqrt(1/nn)


ydel0 <- rep(0,k)
for (i in 1:k) {# Standardize data
    ydel0[i] =  ym[i] - muhat + eu1 
}

# Standardize sd
ydel0 <- ydel0 * sd1/sd(ydel0)  # corrected variance was very important

# Now apply g transformation
for (i in 1:k) gvec[i] = g(sighat[i]*ydel0[i],nn[i],nn[i],sighat[i]) 

wgvec <- factor*gvec

alpi <- matrix(0, nM, k) 
Z <- matrix(0, nM, k)
 
for (i in 1:k) {
paras <- alpDist(n[i],M[i], nM) # alpDist use Exp(1) 
alpi[,i] <- paras[,1] # Use only alpha 	
}

SQ <- s2 <- nn * (alpi * sigma)^2 ## like ns^2/chi^2
sd <- sqrt(s2)
	
alpi <- matrix(0, nM, k)
Z <- matrix(0, nM, k) 

for (i in 1:k) {
Z[,i] <- rnorm(nM) 
paras <- alpDist(n[i],M[i], nM) # Exp(1) data
alpi[,i] <- paras[,1] # Use only alpha 
}

SQ <- nn * (alpi * sigma)^2 ## like ns^2/chi^2
sd <- sqrt(SQ)

xb <- wgvec 
s2 <- colMeans(SQ) 
sd <- sqrt(s2)
  
wt <- nn/s2
xbs <- sum(wt*xb)/sum(wt) # common sample mean
Tobs <- sum(wt*(xb-xbs)^2) # Observed statistic

## Compute p-value and  Count hits ##

pval <- getP(Tobs, sighat, s2, nn, alpi, Z, nM)

  if (verbose) {
    cat("\n", "",METHOD, paste("(alpha = ",alpha,")",sep = ""), "\n", 
        sep = " ")
    line<-paste(rep("-",nchar(METHOD)+28),sep = "")
    cat(line, 
        "\n", sep = "")
    cat("  data :", DNAME, "\n\n", sep = " ")
    
    
    cat("  p.value    :", pval, "\n\n", sep = " ")
    cat(if (pval > alpha) {
      "  Result     : Difference is not statistically significant."
    }
    else {
      "  Result     : Difference is statistically significant."
    }, "\n")
    cat(line, 
        "\n\n", sep = "")
  }

data$Mind <- 1 - data$Mind 	
colnames(data) <- coln

result <- list()
result$p.value <- pval
result$alpha <- alpha
result$method <- METHOD
result$data <- data
result$formula <- formula
result$seed <- seed
attr(result, "class") <- "survtests"
invisible(result)
}
####################################################


