% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-helpers.R
\name{find_annotations}
\alias{find_annotations}
\title{construct annotation track for concomp object}
\usage{
find_annotations(obj, txlist, txdb, orgdb, eps)
}
\arguments{
\item{obj}{\code{concomp} object}

\item{txlist}{see \code{splicegrahm} documentation}

\item{txdb}{see \code{splicegrahm} documentation}

\item{orgdb}{see \code{splicegrahm} documentation}

\item{eps}{see \code{splicegrahm} documentation}
}
\value{
list containing a \code{GRangesList}, \code{tx_match}, and a
\code{ggplot} object, \code{annot_track}
}
\description{
construct annotation track for concomp object
}
\author{
Patrick Kimes
}
\keyword{internal}

