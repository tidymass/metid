% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/50_annotate_metabolites_functions.R
\name{annotate_metabolites}
\alias{annotate_metabolites}
\title{Annotate Metabolites}
\usage{
annotate_metabolites(
  object,
  database,
  based_on = c("ms1", "rt", "ms2"),
  polarity = c("positive", "negative"),
  column = c("rp", "hilic"),
  adduct.table = NULL,
  ce = "all",
  ms1.match.ppm = 25,
  ms2.match.ppm = 30,
  mz.ppm.thr = 400,
  ms2.match.tol = 0.5,
  fraction.weight = 0.3,
  dp.forward.weight = 0.6,
  dp.reverse.weight = 0.1,
  rt.match.tol = 30,
  ms1.match.weight = 0.25,
  rt.match.weight = 0.25,
  ms2.match.weight = 0.5,
  total.score.tol = 0.5,
  candidate.num = 3,
  remove_fragment_intensity_cutoff = 0,
  return_format = c("mass_dataset", "data.frame"),
  threads = 3
)
}
\arguments{
\item{object}{A `mass_dataset` object containing MS1, RT, and/or MS2 data.}

\item{database}{A `databaseClass` object used for metabolite annotation.}

\item{based_on}{Character vector. Specifies the matching criteria to be used for annotation. Can include `"ms1"`, `"rt"`, and/or `"ms2"`. Default is `c("ms1", "rt", "ms2")`.}

\item{polarity}{Character. Ionization mode, either `"positive"` or `"negative"`. Default is `"positive"`.}

\item{column}{Character. The chromatographic column type, either `"hilic"` or `"rp"` (reversed-phase). Default is `"hilic"`.}

\item{adduct.table}{A data frame specifying the adduct table for metabolite annotation. If `NULL`, a default adduct table is loaded based on polarity and column type.}

\item{ce}{Character. Collision energy used in MS2. Default is `"all"`.}

\item{ms1.match.ppm}{Numeric. Mass tolerance in parts per million (ppm) for MS1 peak matching. Default is 25.}

\item{ms2.match.ppm}{Numeric. Mass tolerance in ppm for MS2 peak matching. Default is 30.}

\item{mz.ppm.thr}{Numeric. m/z threshold for ppm calculation. Default is 400.}

\item{ms2.match.tol}{Numeric. Retention time tolerance for MS2 fragment matching. Default is 0.5.}

\item{fraction.weight}{Numeric. Weight for the fraction of matched fragments in MS2 spectra. Default is 0.3.}

\item{dp.forward.weight}{Numeric. Weight for the forward dot product score in MS2 matching. Default is 0.6.}

\item{dp.reverse.weight}{Numeric. Weight for the reverse dot product score in MS2 matching. Default is 0.1.}

\item{rt.match.tol}{Numeric. Retention time matching tolerance in seconds. Default is 30.}

\item{ms1.match.weight}{Numeric. Weight for MS1 matching score in the overall annotation score. Default is 0.25.}

\item{rt.match.weight}{Numeric. Weight for retention time matching score in the overall annotation score. Default is 0.25.}

\item{ms2.match.weight}{Numeric. Weight for MS2 matching score in the overall annotation score. Default is 0.5.}

\item{total.score.tol}{Numeric. Tolerance for the total matching score. Default is 0.5.}

\item{candidate.num}{Numeric. Maximum number of candidate annotations to retain per metabolite. Default is 3.}

\item{remove_fragment_intensity_cutoff}{Numeric. Cutoff to remove low-intensity MS2 fragments. Default is 0.}

\item{return_format}{Character. Specifies the format of the output. Can be `"mass_dataset"` or `"data.frame"`. Default is `"mass_dataset"`.}

\item{threads}{Numeric. Number of threads to use for parallel processing. Default is 3.}
}
\value{
A modified `mass_dataset` object with annotated metabolites added to the `annotation_table` slot.
}
\description{
This function annotates metabolites in the provided object based on MS1, retention time (RT), and/or MS2 spectra data using a specified database. It allows for customization of matching parameters such as m/z match tolerance, retention time tolerance, and MS2 matching criteria.
}
\details{
This function performs metabolite annotation using a combination of MS1, retention time, and MS2 data (if available) from the provided object. The function allows users to customize the matching process, including setting tolerances for MS1 and MS2 matching, adjusting the weights of different scoring components, and selecting a specific chromatographic column and adduct table.

If `ms2` is included in the `based_on` argument, the function extracts both MS1 and MS2 information for annotation. The final annotations are filtered based on the specified score thresholds and only the top `candidate.num` annotations are retained for each metabolite.
}
\examples{
\dontrun{
# Load a sample dataset and database
my_data <- load_mass_dataset("path/to/data")
my_database <- load_database("path/to/database")

# Annotate metabolites using MS1 and MS2 data
annotated_data <- annotate_metabolites(
  object = my_data,
  database = my_database,
  based_on = c("ms1", "ms2"),
  polarity = "positive",
  column = "rp",
  ms1.match.ppm = 20,
  ms2.match.ppm = 25,
  candidate.num = 5,
  threads = 4
)
}


}
