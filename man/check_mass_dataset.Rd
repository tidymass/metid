% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/52_check_parameters_functions.R
\name{check_mass_dataset}
\alias{check_mass_dataset}
\title{Validate the mass_dataset Object for Metabolite Annotation}
\usage{
check_mass_dataset(object, based_on = c("ms1", "rt", "ms2"))
}
\arguments{
\item{object}{A `mass_dataset` object containing MS1, RT, and/or MS2 data.}

\item{based_on}{Character vector. Specifies the criteria to base the validation on. Can include `"ms1"`, `"rt"`, and `"ms2"`. Default is `c("ms1", "rt", "ms2")`.}
}
\value{
The function does not return a value but throws an error if any validation checks fail.
}
\description{
This function checks if the provided `mass_dataset` object is valid for metabolite annotation based on the specified criteria (`ms1`, `rt`, and/or `ms2`). It ensures that the object contains the necessary data for the chosen annotation method.
}
\details{
The function checks the following:
* Ensures that `object` is a valid `mass_dataset` object.
* If `"ms2"` is included in `based_on`, it checks that the `mass_dataset` object contains MS2 data (`object@ms2_data`). If no MS2 data is present, the function throws an error.

If any of these conditions are not met, the function will stop with an appropriate error message.
}
\examples{
\dontrun{
# Load a sample mass dataset
my_data <- load_mass_dataset("path/to/dataset")

# Validate the dataset for MS1 and MS2 information
check_mass_dataset(object = my_data, based_on = c("ms1", "ms2"))
}


}
