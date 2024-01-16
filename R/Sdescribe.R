####################################################
Sdescribe <- function(formula, data, level = 0.95, nM = 5000, na.rm = TRUE, verbose = TRUE){

data <- model.frame(formula, data)
dp <- as.character(formula)


if ((level<=0)|(level>=1)) stop("Confidence level must be between 0 and 1")


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

y.n <- tapply(Time, Group, length)
n.event <- tapply(Mind, Group, sum)


Mind <- 1 - Mind

data <- data.frame(Time,Group,Mind) 	


n0 <- table(data$Group)
grps <- names(n0)
n <- as.numeric(n0)
k <- length(n)
ym <- yv <- ydel <- gym <- gvec <- XB <- alphat <- M <- rep(0,k)

    alphatA=rep(0,k)  ##@  #For Kris Interval
    thetahatA=rep(0,k)  ##@  #For Kris Interval

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

        alp0A <- pi* var(lnYi)^(-1/2) /sqrt(6) # Approx MLE ##@
        otA = optimize(f, ys=Yi, m=M[i] , lower=alp0A/6 , upper=6*alp0A) ##@
        alphatA[i] = otA$minimum # estimate is minumum ##@  #MLE for alpha


ym[i] <- mean(Dati$logTime) # Sample means for k groups
yv[i] <- var(Dati$logTime) # Sample variances for k groups    
}

	for (j in 1:k) {	##@
	YiA=Yi^alphatA[j] ##@
	thetahatA[j]= sum(YiA[1:M[j]]/M[j])^(1/alphatA[j]) ##@  #MLE for theta
	} ##@

	#print("MLEs of alpha, theta(scale), and means")
	#print(alphatA)  ##@ 
	#print(thetahatA) ##@
	#print( thetahatA*gamma( 1 + 1/alphatA) )  ##@
    means_groups <- thetahatA*gamma( 1 + 1/alphatA) 

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

NN=nM
	#####################  ##@
	##@ For Kris calculations	
	alpiA = matrix(0,NN,k)  ##@
	thetaiA = matrix(0,NN,k) ##@
	Gmu = matrix(0,NN,k) ##@
	for (i in 1:k) { ##@
	  paras = alpDist(n[i],M[i], NN) # Exp(1) data ##@
	  alpiA[,i] = paras[,1]  ##@
	  thetaiA[,i]= paras[,2] ##@
	} ##@
	for (i in 1:k) { ##@
	alpiA[,i] = alphatA[i]/alpiA[,i]  ##@
	thetaiA[,i]= thetahatA[i]*(1/thetaiA[,i])^(1/alpiA[,i]) ##@
	} ##@	
	
	#print("For parameters") ##@
	#print("90% CIs for Alpha and Theta - just to check") ##@
	for (i in 1:k) { ##@
	alpiA[,i] = sort(alpiA[,i]) ##@
	#print( c( alpiA[NN*((1-level)/2),i], alpiA[NN*(1-(1-level)/2),i] )) ##@
	
	thetaiA[,i] = sort (thetaiA[,i]) ##@
	#print( c( thetaiA[NN*((1-level)/2),i], thetaiA[NN*(1-(1-level)/2),i] )) ##@
	}
	
	##@ Kris CI for means
	Gmu= thetaiA * gamma(1+1/alpiA) ##@
	for (i in 1:k) { ##@
	  Gmu[,i]=sort(Gmu[,i]) ##@
	} ##@
	
	#print(".........")
	#print("90% CIs for k Means - Need to make Bonferoni adjustment if necessary")
	##@   Kris is reporting 95% interval for his data
	confIntervals <- NULL
	for (i in 1:k) { ##@
	#  print( c( Gmu[NN*.025,i], Gmu[NN*.975,i] ) ) ##@
	confIntervals <- rbind(confIntervals, c( Gmu[NN*(1-level)/2,i], Gmu[NN*(1-(1-level)/2),i] ))
	} ##@

desc_stat <- cbind(y.n,n.event,means_groups,confIntervals)



colnames(desc_stat)<-c("n", "nE","Mean", paste((1-level)/2*100, "%",sep = ""), paste((1-(1-level)/2)*100, "%",sep = ""))
rownames(desc_stat)<-levels(data$Group)

attr(desc_stat, "class") <- "Sdescribe"
return(desc_stat)

}
####################################################




