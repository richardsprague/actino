% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/actinoJSON.R
\name{full_taxa}
\alias{full_taxa}
\title{given a row taxaRow in df, return a character vector showing the full set of tax_ranks}
\usage{
full_taxa(df, taxaRow)
}
\arguments{
\item{df}{dataframe that contains all the taxa}

\item{taxaRow}{a list of one row in the dataframe (the one you want to convert)}
}
\description{
Needed by functions that convert tax_names to QIIME (or other) format.
}
