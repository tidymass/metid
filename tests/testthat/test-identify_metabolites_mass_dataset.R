library(massdataset)
library(magrittr)
library(dplyr)
ms1_data =
  readr::read_csv(file.path(
    system.file("ms1_peak", package = "metid"),
    "ms1.peak.table.csv"
  ))

ms1_data = data.frame(ms1_data, sample1 = 1, sample2 = 2)

expression_data = ms1_data %>%
  dplyr::select(-c(name:rt))

variable_info =
  ms1_data %>%
  dplyr::select(name:rt) %>%
  dplyr::rename(variable_id = name)

sample_info =
  data.frame(
    sample_id = colnames(expression_data),
    injection.order = c(1, 2),
    class = c("Subject", "Subject"),
    group = c("Subject", "Subject")
  )
rownames(expression_data) = variable_info$variable_id

object = create_mass_dataset(
  expression_data = expression_data,
  sample_info = sample_info,
  variable_info = variable_info
)

object

data("hmdb_ms1_database0.0.3", package = "metid")

data("snyder_database_rplc0.0.3", package = "metid")

database = snyder_database_rplc0.0.3

object1 <-
  annotate_metabolites_mass_dataset(object = object[1:10,],
                                    database = snyder_database_rplc0.0.3)


test_that("identify_metabolites_mass_dataset", {
  testthat::expect_s4_class(object = object1, "mass_dataset")
  # testthat::expect_true(massdataset::check_mass_dataset_class(object = object1))
})
