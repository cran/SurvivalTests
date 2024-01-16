


weibull.test <- function(formula, data, alpha = 0.05, na.rm = TRUE, verbose = TRUE){ 


  data <- model.frame(formula, data)
  dp <- as.character(formula)
  DNAME <- paste(dp[[2L]], "and", dp[[3L]])

if (any(colnames(data)==dp[[3L]])==FALSE) stop("The name of group variable does not match the variable names in the data. The group variable must be one factor.")
if (any(colnames(data)==dp[[2L]])==FALSE) stop("The name of response variable does not match the variable names in the data.")

y = data[[dp[[2L]]]]
group = data[[dp[[3L]]]]


if (!(is.factor(group)|is.character(group))) stop("The group variable must be a factor or a character.") 
if (is.character(group)) group <- as.factor(group)
if (!is.numeric(y)) stop("The response must be a numeric variable.") 

  dname1<-dp[[2L]]
  dname2<-dp[[3L]]

    group = as.factor(group)
    k = length(levels(group))
    if (length(y) != length(group)) {
        stop("The lengths of y and group must be equal")
    }



if (na.rm) {
        completeObs <- complete.cases(y, group)
        y <- y[completeObs]
        group <- group[completeObs]
    }



 kk <- tapply(y, group, wp.test)
 stor1 <- as.numeric(unlist(lapply(1:k, function(i) kk[[i]]$statistic)))
 stor2 <- as.numeric(unlist(lapply(1:k, function(i) kk[[i]]$p)))
 method.name = "Weibullness Test from Weibull Plot"




   store = data.frame(matrix(NA, nrow = k, ncol = 4))
    store$X1 = levels(group)
    store$X2 = stor1
    store$X3 = stor2
    store$X4 = ifelse(store$X3 > alpha, "Not reject", "Reject")
    colnames(store) = c("Level", "Statistic", "p.value", "  Weibullness")


 if (verbose) {
        cat("\n", "",method.name, paste("(alpha = ",alpha,")",sep = ""), "\n", sep = " ")
        cat("---------------------------------------------------------------", 
            "\n", sep = " ")
        cat("  data :", dname1, "and", dname2, "\n\n", sep = " ")
          
        print(store)
          cat("---------------------------------------------------------------", 
            "\n\n", sep = " ")
    }

    
invisible(store)
}

