% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Networks_analysis.R
\name{vector_smartmerge}
\alias{vector_smartmerge}
\title{Merge dataframes by shared rows
Dataframes to be merged need to have rownames with the TFs or other entities being merged on.}
\usage{
vector_smartmerge(df1, df2, margin = "row")
}
\arguments{
\item{df1}{Dataframe 1, essentially a named vector in dataframe format}

\item{df2}{Dataframe 2, essentially a named vector in dataframe format with}

\item{margin}{The direction of merging dataframes only implemented for rows}
}
\description{
Merge dataframes by shared rows
Dataframes to be merged need to have rownames with the TFs or other entities being merged on.
}
