% Generated by roxygen2 (4.1.0.9001): do not edit by hand
% Please edit documentation in R/protractedbirthdeath.R
\name{plotsim}
\alias{plotsim}
\title{A plot of repeats with means and sds.}
\usage{
plotsim(x)
}
\arguments{
\item{x}{an output from summaryrepsim}
}
\value{
a ggplot
}
\description{
A plot of repeats with means and sds.
}
\examples{
\dontrun{
plotsim(summaryrepsim(c(0.2,0.2,0.2,0.1,0.1),15,15))
}
}
\seealso{
\code{\link{repsim2}} which produces the inputs \code{\link{plotsim}}
  which plots this functions output
}

