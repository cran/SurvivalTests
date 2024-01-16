\name{plot.Sdescribe}
\alias{plot.Sdescribe}
\title{
Error Bars with Means and Confidence Limits
}
\description{
\code{plot.Sdescribe} produce error bars with means and confidence limits of the given grouped values.
}

\usage{
\method{plot}{Sdescribe}(x, ylim = NULL, xlab = NULL, ylab = NULL, title = NULL, width = NULL, ...)}


\arguments{
  \item{x}{a Sdescribe object.}
  \item{ylim}{a limit for y axes.} 
  \item{xlab}{a label for the x axis.} 
  \item{ylab}{a label for the y axis.} 
  \item{title}{a main title for the plot.} 
  \item{width}{the little lines at the tops and bottoms of the error bars (defaults to 0.80).} 
  \item{...}{additional arguments.} 

}


\value{No return value, called for side effects}


\author{
Osman Dag
}



\examples{

library(survival)
lung$status <- ifelse(lung$status == 2, 1, 0)
lung$age <- arules::discretize(lung$age, breaks = 3, labels = c("Low","Medium","High"))


library(SurvivalTests)
\donttest{
out <- Sdescribe(time~age+status, lung)

plot(out, ylim = c(0,500))
plot(out)
}

}

\keyword{functions}
