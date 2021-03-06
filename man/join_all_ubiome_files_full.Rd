% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/actinoJSON.R
\name{join_all_ubiome_files_full}
\alias{join_all_ubiome_files_full}
\title{returns a big dataframe joining all tax_names from a list of JSON files}
\usage{
join_all_ubiome_files_full(flist, count.normalized = FALSE, site = "gut")
}
\arguments{
\item{flist}{list of file pathnames, each pointing to a uBiome JSON file}

\item{count.normalized}{use the count_norm field instead of raw read numbers (default is FALSE)}

\item{site}{character vector}
}
\value{
data frame where first col is all taxa names, other cols are samples
}
\description{
returns a big dataframe joining all tax_names from a list of JSON files
}
