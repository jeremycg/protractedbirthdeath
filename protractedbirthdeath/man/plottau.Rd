% Generated by roxygen2 (4.1.0.9001): do not edit by hand
% Please edit documentation in R/protractedbirthdeath.R
\name{plottau}
\alias{plottau}
\title{A plot of tau (time to speciation)}
\usage{
plottau(holding, vars, max)
}
\arguments{
\item{holding}{an output of tauloop}

\item{vars}{variable to fix}

\item{max}{the value to clip tau on the plot at}
}
\value{
a countour ggplot of tau
}
\description{
A plot of tau (time to speciation)
}
\examples{
\dontrun{
z<-tauloop(0.5,"incipext")
plottau(z,"incipext",5)
}
}
\seealso{
\code{\link{tauloop}} which produces the inputs
  which this function plots
}

