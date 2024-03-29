% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pooling_mapping.R
\name{pool_bin_stat}
\alias{pool_bin_stat}
\title{Calculate mean or median per feature}
\usage{
pool_bin_stat(bin_tibble, operator = "mean")
}
\arguments{
\item{bin_tibble}{a tibble of coverage values (numeric) that need to be pooled,
column names should contain \verb{Sequence Id}, \verb{Bin Id}, \code{Coverage}}

\item{operator}{a character vector either of "mean" or "median" value; identifies pooling strategy.}
}
\value{
A tibble with first column "Feature" that contains bin IDs, and the rest of the columns represent samples with bins' pooled values.
}
\description{
pool_bin_stat calculates per bin mean or median scaffold value (coverage or mapped read numbers).
}
