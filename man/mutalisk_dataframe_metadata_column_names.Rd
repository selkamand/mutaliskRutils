% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ParseMutaliskResults.R
\name{mutalisk_dataframe_metadata_column_names}
\alias{mutalisk_dataframe_metadata_column_names}
\title{Mutalisk Dataframe}
\usage{
mutalisk_dataframe_metadata_column_names(mutalisk_dataframe)
}
\arguments{
\item{mutalisk_dataframe}{a dataframe that can be produced using mutalisk_best_signature_directory_to_dataframe. Can also just make it yourself, if you want to visualise non-mutalisk data. Dataframe just needs 3 columns:
\enumerate{
\item \strong{SampleID}: a sample identifier.
\item \strong{Signatures}: an identifier for a particular signature.
\item \strong{Contributions}: the percentage contribution of the signature to the patients genetic profile (0.1 = 10\%).
}}
}
\value{
a character vector containing names of metadata columns. If no metadata columns have been added, returns a zero length character vector. (character)
}
\description{
Get a vector of metadata columns from a mutalisk_dataframe.
}
