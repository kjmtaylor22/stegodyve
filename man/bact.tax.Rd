% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bact_tax.R
\name{bact.tax}
\alias{bact.tax}
\title{Bacterial Taxonomy}
\usage{
bact.tax(taxonomy, database = NULL)
}
\arguments{
\item{taxonomy}{Taxonomy file output by \code{biom_as_csv} function.}

\item{database}{Character string; either 'silva' or 'green'.}
}
\value{
Parsed taxonomy table; also output to file (.RD).
}
\description{
Parses and formats data from BIOM table.
}
\author{
Kara J.M. Taylor
}
