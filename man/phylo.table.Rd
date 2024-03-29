% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_qSIP_MAGs.R
\name{phylo.table}
\alias{phylo.table}
\title{Master phyloseq object using the MAG phyloseq objects}
\usage{
phylo.table(mag, taxa, samples)
}
\arguments{
\item{mag}{phyloseq-styled MAG table}

\item{taxa}{phyloseq-styled taxa table}

\item{samples}{sample information table}
}
\value{
phyloseq object for MAGs
}
\description{
Creates a phyloseq-style object using processed phyloseq objects for otu
table (here, MAG table), taxa table, and sample table
}
