Spaircomp <- function(x, ...) UseMethod("Spaircomp")

Spaircomp.default <- function(x, ...) Spaircomp.survtests(x, ...)


Spaircomp.survtests<- function(x, adjust.method = c("bonferroni", "holm", "hochberg", "hommel", "BH", 
                                             "BY", "fdr", "none"), verbose = TRUE, ...){
  
  if(x$p.value>x$alpha) stop(paste("Pairwise comparisons could not be performed since difference is not statistically significant (alpha = ",x$alpha,").",sep = ""))
  
  data<-x$data
  

  
  Time = data[,1]
  Group = data[,2]
  Mind =  data[,3]
  
  if (length(levels(Group))==2) stop("The group variable must have more than two levels to make pairwise comparison.")
  

  alpha <- x$alpha
  id <- levels(Group)
  comb <- t(combn((id), 2))
  comb2 <- dim(comb)[1]
  
  pval <- NULL
  for (i in 1:comb2){
    data_sub <- data[(Group == comb[i,1])|(Group == comb[i,2]),]
	data_sub[,2] <- factor(data_sub[,2])
	set.seed(x$seed)
    if (x$method == "Generalized Test for Survival ANOVA") pval <- c(pval, SANOVA(x$formula,data_sub,verbose=F)$p.value)
	
  }
  
  adjust.method <- match.arg(adjust.method)
  
  if (adjust.method == "bonferroni") method.name = paste("Bonferroni Correction (alpha = ",alpha,")",sep = "")
  if (adjust.method == "holm") method.name = paste("Holm Correction (alpha = ",alpha,")",sep = "")
  if (adjust.method == "hochberg") method.name = paste("Hochberg Correction (alpha = ",alpha,")",sep = "")
  if (adjust.method == "hommel") method.name = paste("Hommel Correction (alpha = ",alpha,")",sep = "")
  if ((adjust.method == "fdr")|(adjust.method == "BH")) method.name = paste("Benjamini-Hochberg Correction (alpha = ",alpha,")",sep = "")
  if (adjust.method == "BY") method.name = paste("Benjamini-Yekutieli Correction (alpha = ",alpha,")",sep = "")
  if (adjust.method == "none") method.name = paste("No Correction (alpha = ",alpha,")",sep = "")
  
  
  padj <- p.adjust(pval, method = adjust.method)
  
  
  
  store = data.frame(matrix(NA, nrow = comb2, ncol = 4))
  
  store$X1 = comb[,1]
  store$X2 = comb[,2]
  store$X3 = padj
  store$X4 = ifelse(store$X3 <= alpha, "Reject", "Not reject")
  colnames(store) = c("Level (a)", "Level (b)", "p.value", "  No difference")
  
  if (verbose) {
    cat("\n", "",method.name, "\n", sep = " ")
    cat("---------------------------------------------------------------", "\n", sep = " ")
    print(store)
    cat("---------------------------------------------------------------", "\n\n", sep = " ")
  }
  
  invisible(store)
  
}

