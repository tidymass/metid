##------------------------------------------------------------------------------
#' @title Construct public MS2 database from MoNA with msp format.
#' @description Construct MS2 spectra database according to mzXML data and compound information table (csv format).
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param file The file name of MassBank or MoNA database (mgf format).
#' @param only.remain.ms2 Only remain the metabolites with MS2 spectra?
#' @param path Work directory.
#' @param version The version of you database. Default is 0.0.1.
#' @param source The source of your database.
#' @param link Website link of the source.
#' @param creater Creater name. For example, Xiaotao Shen.
#' @param email email address.
#' @param rt Do the metabolites have RT information or not?. If not, set it as FALSE.
#' @param threads The number of threads
#' @importFrom magrittr %>%
#' @importFrom crayon red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @importFrom readr cols
#' @importFrom pbapply pblapply
#' @return A databaseClass object.
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/metid.html}
#' @export

# # ####MoNA library
# sxtTools::setwd_project()
# setwd("other_files/all_ms2_database/mona/2021_6_10/")
# mona_database = construct_mona_database(file = "MoNA-export-LC-MS-MS_Spectra.msp")
# 
# mz <- 167.035
# rt <- 450
# ms2 <- data.frame(
#   mz = c(
#     91.0187,
#     95.0138,
#     95.0503,
#     99.0088,
#     108.0217,
#     121.0297,
#     123.045,
#     152.0113,
#     167.0348
#   ),
#   intensity = c(75.6,
#                 81.8,
#                 61.7,
#                 69.9,
#                 4877.1,
#                 46.5,
#                 8682.4,
#                 25133.4,
#                 44848.5),
#   stringsAsFactors = FALSE
# )
# 
# save(mona_database, file = "mona_database", compress = "xz")
# 
# annotation_result <-
#   metid::identify_single_peak(
#     ms1.mz = mz,
#     ms1.rt = rt,
#     ms2 = ms2,
#     ms1.match.ppm = 15,
#     rt.match.tol = 30,
#     ms2.match.tol = 0.5,
#     database = "mona_database",
#     path = ".",
#     polarity = "negative"
#   )
# 
# which_has_identification(annotation_result)
# get_identification_table(annotation_result, type = "new")
# ms2plot(object = annotation_result,
#         database = mona_database,
#         which.peak = "mz167.035rt450")

construct_mona_database = function(file,
                                   only.remain.ms2 = TRUE,
                                   path = ".",
                                   version = "0.0.1",
                                   source = "MoNA",
                                   link = "https://mona.fiehnlab.ucdavis.edu/",
                                   creater = "Xiaotao Shen",
                                   email = "shenxt1990@163.com",
                                   rt = FALSE,
                                   threads = 5) {
  mona_database = read_msp_mona(file = file)
  
  all_metabolite_names =
    purrr::map(mona_database, function(x) {
      rownames(x$info)
    }) %>%
    unlist() %>%
    unique()
  
  metabolite_info =
    mona_database %>%
    purrr::map(function(x) {
      x = as.data.frame(x$info)
      new_x = x[, 1]
      # new_x = as.data.frame(matrix(data = new_x, nrow = 1))
      names(new_x) = rownames(x)
      new_x = new_x[all_metabolite_names]
      names(new_x) = all_metabolite_names
      new_x
    })  %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  colnames(metabolite_info) = all_metabolite_names
  
  ###remove the metabolites without MS2 spectra
  if (only.remain.ms2) {
    remain_idx =
      which(metabolite_info$Spectrum_type == "MS2")
    metabolite_info =
      metabolite_info[remain_idx,]
    mona_database = mona_database[remain_idx]
  }
  
  metabolite_info =
    metabolite_info %>%
    dplyr::select(
      Compound.name = Name,
      mz = ExactMass,
      Formula,
      MoNA.ID = `DB#`,
      dplyr::everything()
    )
  
  metabolite_info =
    metabolite_info %>%
    dplyr::mutate(
      Lab.ID = paste("MoNA", seq_len(nrow(metabolite_info)), sep = "_"),
      RT = NA,
      CAS.ID = NA,
      HMDB.ID = NA,
      KEGG.ID = NA,
      mz.pos = NA,
      mz.neg = NA,
      Submitter = "MoNA",
      Family = NA,
      Sub.pathway = NA,
      Note = NA
    ) %>%
    dplyr::select(
      Lab.ID,
      Compound.name,
      mz,
      RT,
      CAS.ID,
      HMDB.ID,
      KEGG.ID,
      Formula,
      mz.pos,
      mz.neg,
      Submitter,
      Family,
      Sub.pathway,
      Note,
      dplyr::everything()
    )
  
  metabolite_info$Collision_energy[is.na(metabolite_info$Collision_energy)] = "not_available"
  metabolite_info$Collision_energy[metabolite_info$Collision_energy == ""] = "not_available"
  
  #####create metid database format
  positive_idx = which(metabolite_info$Ion_mode == "P")
  negative_idx = which(metabolite_info$Ion_mode == "N")
  
  Spectra.positive = mona_database[positive_idx]
  Spectra.negative = mona_database[negative_idx]
  
  names(Spectra.positive) = metabolite_info$Lab.ID[positive_idx]
  names(Spectra.negative) = metabolite_info$Lab.ID[negative_idx]
  
  Spectra.positive =
    purrr::map2(
      .x = Spectra.positive,
      .y = metabolite_info$Collision_energy[positive_idx],
      .f = function(x, y) {
        x = x$spec
        x = list(x)
        names(x) = y
        x
      }
    )
  
  Spectra.negative =
    purrr::map2(
      .x = Spectra.negative,
      .y = metabolite_info$Collision_energy[negative_idx],
      .f = function(x, y) {
        x = x$spec
        x = list(x)
        names(x) = y
        x
      }
    )
  
  database.info <- list(
    "Version" = version,
    "Source" = source,
    "Link" = link,
    "Creater" = creater,
    "Email" = email,
    "RT" = rt
  )
  
  spectra.info <- as.data.frame(metabolite_info)
  rm(list = "metabolite_info")
  
  Spectra <- list("Spectra.positive" = Spectra.positive,
                  "Spectra.negative" = Spectra.negative)
  
  database <- new(
    Class = "databaseClass",
    database.info = database.info,
    spectra.info = spectra.info,
    spectra.data = Spectra
  )
  
  database@database.info$RT <-
    ifelse(all(is.na(database@spectra.info$RT)), FALSE, TRUE)
  message(crayon::bgRed("All done!\n"))
  return(database)
}
