\name{weibull.test}
\alias{weibull.test}
\title{
Weibullness Test
}
\description{
\code{weibull.test} performs weibullness test from weibull plot in each group.
}

\usage{weibull.test(formula, data, alpha = 0.05, na.rm = TRUE, verbose = TRUE)}

\arguments{
  \item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs} gives the survival times and \code{rhs} the corresponding groups.}
  \item{data}{a tibble or data frame containing the variables in \code{formula}.} 
  \item{alpha}{the level of significance to assess normality. Default is set to alpha = 0.05.}
  \item{na.rm}{a logical value indicating whether NA values should be stripped before the computation proceeds.}
  \item{verbose}{a logical for printing output to R console.}

}


\value{
A data frame gives the test results for the weibullness of groups.
}


\author{
Osman Dag
}

\seealso{
\code{\link[weibullness]{wp.test}}  
}



\examples{

library(survival)
lung$status <- ifelse(lung$status == 2, 1, 0)
lung$age <- arules::discretize(lung$age, breaks = 3, labels = c("Low","Medium","High"))

library(SurvivalTests)

weibull.test(time~age, lung)


}

\keyword{functions}
