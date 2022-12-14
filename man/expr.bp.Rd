% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QC_plot_fun.R
\name{expr.bp}
\alias{expr.bp}
\title{Plot Expression Value Boxplots for Samples}
\usage{
expr.bp(Y, n = 10)
}
\arguments{
\item{Y}{a matrix of expression values with row names corresponding to genes and column names to sample names}

\item{n}{the number of samples used to plot the boxplots. Default is 10}
}
\value{
boxplots for the expression values with each representing those of a sample
}
\description{
Plot Expression Value Boxplots for Samples
}
\examples{
data(example)
expr.bp(Y = example)
expr.bp(Y = example, n = 5)

}
