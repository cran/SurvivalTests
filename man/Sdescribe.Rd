\name{Sdescribe}
\alias{Sdescribe}
\title{
Descriptive Statistics
}
\description{
\code{Sdescribe} produces number of observations (n), number of events (nE), and the mean survival with confidence limits for the mean.
}

\usage{Sdescribe(formula, data, level = 0.95, nM = 5000, na.rm = TRUE, 
	  verbose = TRUE)}

\arguments{

  \item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs} gives the survival times and \code{rhs} the corresponding groups + the status variable (1: non-censored, 0: censored).}
  \item{data}{a data frame containing the variables in \code{formula}.} 
  \item{level}{the confidence level.}
  \item{nM}{a number of bootstrap samples.}  
  \item{na.rm}{a logical value indicating whether NA values should be stripped before the computation proceeds.}
  \item{verbose}{a logical for printing output to R console.}


}


\value{
Returns a data.frame of output with class "Sdescribe".
}



\author{
Malwane Ananda, Osman Dag
}



\examples{

library(survival)
lung$status <- ifelse(lung$status == 2, 1, 0)
lung$age <- arules::discretize(lung$age, breaks = 3, labels = c("Low","Medium","High"))

library(SurvivalTests)
\donttest{
Sdescribe(time~age+status, lung)
}
}

\keyword{functions}
