% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/6_zzz.R
\name{extract_ms1_database}
\alias{extract_ms1_database}
\title{Extract MS1 Information from a Database}
\usage{
extract_ms1_database(database)
}
\arguments{
\item{database}{A `databaseClass` object containing MS1 information in the `spectra.info` slot.}
}
\value{
A data frame with MS1 information, where the `mz` and `RT` columns are numeric and rows with missing (`NA`) m/z values are removed. If the input is not a `databaseClass` object, the function returns `NULL`.
}
\description{
This function extracts MS1 data (m/z and retention time) from a `databaseClass` object and returns it as a filtered data frame. It ensures that the m/z values are numeric and removes rows where the m/z values are missing (`NA`).
}
\details{
The function extracts the `spectra.info` slot from the `databaseClass` object, converts the `mz` and `RT` columns to numeric, and filters out rows with missing m/z values. This ensures that only valid MS1 data is returned.
}
\examples{
\dontrun{
# Load a sample database
my_database <- load_database("path/to/database")

# Extract MS1 data from the database
ms1_data <- extract_ms1_database(my_database)
head(ms1_data)
}


}
