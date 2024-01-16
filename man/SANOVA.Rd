\name{SANOVA}
\alias{SANOVA}
\title{
Generalized Test for Survival ANOVA
}
\description{
\code{SANOVA} performs generalized test for survival ANOVA.
}

\usage{SANOVA(formula, data, nM = 5000, seed = 123, alpha = 0.05, na.rm = TRUE, verbose = TRUE)}

\arguments{
  \item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs} gives the survival times and \code{rhs} the corresponding groups + the status variable (1: non-censored, 0: censored).}
  \item{data}{a data frame containing the variables in \code{formula}.} 
  \item{nM}{a number of bootstrap samples.}  
  \item{seed}{a seed number for the reproducibility of results. Default is set to 123.}  
  \item{alpha}{the level of significance to assess the statistical difference. Default is set to alpha = 0.05.} 
  \item{na.rm}{a logical value indicating whether NA values should be stripped before the computation proceeds.}
  \item{verbose}{a logical for printing output to R console.}

}

\details{
SANOVA performs parametric survival ANOVA when the underlying data are distributed as Weibull or Gumbel. SANOVA tests are based on generalized p-value approach (cf. Tsui and Weerahandi (1989) and Weerahandi (2004)) extended for situations of affine invariance.

}


\value{
A list with class "survtests" containing the following components:

\item{p.value}{the p-value of the test.}
\item{alpha}{the level of significance to assess the statistical difference.} 
\item{method}{the character string "Generalized Test for Survival ANOVA".}
\item{data}{a data frame containing the variables in which NA values (if exist) are removed.}
\item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs} gives the sample values and \code{rhs} the corresponding groups + the status variable (1: non-censored, 0: censored).}
\item{seed}{a seed number for the reproducibility of results.}

}

\references{

Tsui K. and Weerahandi, S. (1989). Generalized P-values in Significance Testing of Hypotheses in the Presence of Nuisance Parameters. \emph{Journal of the American Statistical Association}, \bold{84}, 602-607. 

Weerahandi, S. (2004). Generalized Inference in Repeated Measures: Exact Methods in MANOVA and Mixed Models, Wiley. 

}


\author{
Sam Weerahandi, Malwane Ananda, Osman Dag
}



\examples{

library(survival)
lung$status <- ifelse(lung$status == 2, 1, 0)
lung$age <- arules::discretize(lung$age, breaks = 3, labels = c("Low","Medium","High"))

library(SurvivalTests)
\donttest{
SANOVA(time~age+status, lung, alpha = 0.05)
}
}

\keyword{functions}
