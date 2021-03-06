% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-helpers.R
\name{sampl_sort}
\alias{sampl_sort}
\title{function defining all of the sort orders}
\usage{
sampl_sort(sort_idx, vals_e, vals_j, n)
}
\arguments{
\item{sort_idx}{integer value, see documentation for \code{splicegrahm}}

\item{vals_e}{matrix of exon coverages}

\item{vals_j}{matrix of junction coverages}

\item{n}{number of samples}
}
\description{
function defining all of the sort orders
}
\details{
sort_idx can take values of either:
\itemize{
\item{\code{1}}: sort on first exon
\item{\code{2}}: sort on PC 2
\item{\code{3}}: sort on mean exon coverage
\item{\code{4}}: sort on mean exon log-coverage
\item{\code{5}}: rev sort on mean exon coverage
\item{\code{6}}: rev sorton mean exon log-coverage
}
}
\author{
Patrick Kimes
}
\keyword{internal}

