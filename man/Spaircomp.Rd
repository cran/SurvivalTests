\name{Spaircomp}
\alias{Spaircomp}
\alias{Spaircomp.default}
\alias{Spaircomp.survtests}
\title{
Pairwise Comparisons of Survival Tests for One-Way Layout
}
\description{
\code{Spaircomp} is a generic function for pairwise comparisons of survival tests in one-way layout by adjusting p-values.
}

\usage{
\method{Spaircomp}{survtests}(x, adjust.method = c("bonferroni", "holm", "hochberg", "hommel", "BH", 
  "BY", "fdr", "none"), verbose = TRUE, ...)}

\arguments{
  \item{x}{a \code{survtests} object.}
  \item{adjust.method}{Method for adjusting p values (see \code{\link[stats]{p.adjust}}). Default is set to "bonferroni".} 
  \item{verbose}{a logical for printing output to R console.}
  \item{...}{Additional arguments affecting multiple comparisons of groups in one-way independent designs.}

}

\value{
Returns a data.frame of output.
}

\author{
Osman Dag
}


\examples{

library(survival)
lung$status <- ifelse(lung$status == 2, 1, 0)
lung$age <- arules::discretize(lung$age, breaks = 3, labels = c("Low","Medium","High"))

library(SurvivalTests)
\donttest{
out = SANOVA(time~age+status, lung, alpha = 0.05)

Spaircomp(out, adjust.method = "bonferroni")


# to illustrate pairwise comparison, alpha is set to 0.5
out2 = SANOVA(time~age+status, lung, alpha = 0.5)

Spaircomp(out2, adjust.method = "bonferroni")

Spaircomp(out2, adjust.method = "hochberg")

Spaircomp(out2, adjust.method = "hommel")

Spaircomp(out2, adjust.method = "holm")
}
}

\keyword{functions}
