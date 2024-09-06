#' Construct a Mass Spectra Database
#'
#' This function constructs a spectral database from provided metabolite information and MS2 spectra in both positive and negative modes. The database is created as a `databaseClass` object, including metadata, spectra information, and spectral data for metabolites.
#'
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @param path A character string specifying the path to the folder containing the metabolite information and MS2 data files. Defaults to the current working directory (`"."`).
#' @param version A character string specifying the version of the database. Defaults to `"0.0.1"`.
#' @param metabolite.info.name A character string specifying the file name of the metabolite information (CSV format). Defaults to `"metabolite.info.csv"`.
#' @param source A character string specifying the source of the database. Defaults to `"Shen Lab"`.
#' @param link A character string specifying the link to the source of the database. Defaults to `"http://snyderlab.stanford.edu/"`.
#' @param creater A character string specifying the creator of the database. Defaults to `"Xiaotao Shen"`.
#' @param email A character string specifying the email of the creator. Defaults to `"xiaotao.shen@outlook.com"`.
#' @param rt A logical value indicating whether the retention time (RT) information is available. Defaults to `TRUE`.
#' @param mz.tol A numeric value specifying the tolerance for matching the mass-to-charge ratio (m/z). Defaults to `15`.
#' @param rt.tol A numeric value specifying the tolerance for matching the retention time (RT). Defaults to `30`.
#' @param threads A numeric value specifying the number of threads to use for parallel processing. Defaults to `3`.
#'
#' @return A `databaseClass` object containing the constructed database, including metadata, spectra information, and MS2 spectral data for both positive and negative modes.
#'
#' @details
#' The function reads the metabolite information and MS2 spectra from the specified path. It checks for the required files (`metabolite.info.csv`, `POS`, and `NEG` folders) and loads the MS2 data from mzXML or MGF files. The spectra are matched to the metabolites based on mass-to-charge ratio (m/z) and retention time (RT) tolerances.
#'
#' - **Positive Mode**: The function reads and processes spectra files in the `POS` folder.
#' - **Negative Mode**: The function reads and processes spectra files in the `NEG` folder.
#'
#' The function uses matching algorithms to associate the metabolites with their corresponding MS2 spectra, based on the provided tolerances.
#'
#' @examples
#' \dontrun{
#' # Construct a database from the current working directory
#' db <- construct_database(path = ".", version = "1.0.0", metabolite.info.name = "metabolite_info.csv")
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom crayon red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @importFrom readr cols
#' @importFrom pbapply pblapply
#' @export

construct_database <-
  function(path = ".",
           version = "0.0.1",
           metabolite.info.name = "metabolite.info.csv",
           source = "Shen Lab",
           link = "http://snyderlab.stanford.edu/",
           creater = "Xiaotao Shen",
           email = "xiaotao.shen@outlook.com",
           rt = TRUE,
           mz.tol = 15,
           rt.tol = 30,
           threads = 3) {
    ##check data first
    
    file <- dir(path)
    if (all(file != metabolite.info.name)) {
      message(crayon::red("No", metabolite.info.name, "in your", path, "\n"))
      return(NULL)
    }
    
    if (all(file != "POS")) {
      message(crayon::red("No POS file in your", path, "\n"))
    } else{
      file_pos <- dir(file.path(path, "POS"))
      if (length(file_pos) == 0) {
        message(crayon::red("No mzXML/mgf files in POS folder."))
      } else{
        if (sum(stringr::str_detect(file_pos, "mzXML|mgf")) == 0) {
          message(crayon::red("No mzXML/mgf files in POS folder."))
        }
      }
    }
    
    if (all(file != "NEG")) {
      message(crayon::red("No NEG file in your", path, "\n"))
    } else{
      file_neg <- dir(file.path(path, "NEG"))
      if (length(file_neg) == 0) {
        message(crayon::red("No mzXML/mgf files in NEG folder."))
      } else{
        if (sum(stringr::str_detect(file_neg, "mzXML|mgf")) == 0) {
          message(crayon::red("No mzXML/mgf files in NEG folder."))
        }
      }
    }
    
    ##read metabolite information
    message(crayon::green("Reading metabolite information..."))
    metabolite.info <-
      readTable(file = file.path(path, metabolite.info.name))
    
    message(crayon::green("Reading positive MS2 data..."))
    
    file.pos <-
      dir(file.path(path, 'POS'), full.names = TRUE)
    
    if (length(file.pos) > 0) {
      ####metabolite.info has to have mz.pos column
      if (all(colnames(metabolite.info) != "mz.pos")) {
        stop("You need to provide mz.pos column in metabolite.info")
      } else{
        metabolite.info$mz.pos <-
          as.numeric(metabolite.info$mz.pos)
        
        if (all(is.na(metabolite.info$mz.pos))) {
          stop("All are NAs in the mz.pos column.")
        }
        
        if (any(is.na(metabolite.info$mz.pos))) {
          warning("NA in the mz.pos column.")
        }
      }
      
      if (stringr::str_detect(file.pos, "mzXML")[1]) {
        ms2.data.pos <-
          masstools::read_mzxml(file = file.pos, threads = threads)
      }
      
      if (stringr::str_detect(file.pos, "mgf")[1]) {
        ms2.data.pos <-
          masstools::read_mgf4database(file = file.pos)
      }
      
      ms1.info.pos <- lapply(ms2.data.pos, function(x) {
        x[[1]]
      })
      
      ms1.info.pos <- do.call(rbind, ms1.info.pos) %>%
        as.data.frame()
      
      ms1.info.pos$file <- basename(ms1.info.pos$file)
      
      ms2.info.pos <- lapply(ms2.data.pos, function(x) {
        x[[2]]
      })
      
      rm(list = "ms2.data.pos")
      
      message(crayon::red("OK."))
    } else{
      ms1.info.pos <- NULL
    }
    
    message(crayon::green("Reading negative MS2 data..."))
    
    file.neg <-
      dir(file.path(path, 'NEG'), full.names = TRUE)
    
    if (length(file.neg) > 0) {
      if (all(colnames(metabolite.info) != "mz.neg")) {
        stop("You need to provide mz.neg column in metabolite.info")
      } else{
        metabolite.info$mz.neg <-
          as.numeric(metabolite.info$mz.neg)
        
        if (all(is.na(metabolite.info$mz.neg))) {
          stop("All are NAs in the mz.neg column.")
        }
        
        if (any(is.na(metabolite.info$mz.pos))) {
          warning("NA in the mz.neg column.")
        }
      }
      
      if (stringr::str_detect(file.neg, "mzXML")[1]) {
        ms2.data.neg <-
          masstools::read_mzxml(file = file.neg, threads = threads)
      }
      
      if (stringr::str_detect(file.neg, "mgf")[1]) {
        ms2.data.neg <-
          masstools::read_mgf4database(file = file.neg)
      }
      
      ms1.info.neg <- lapply(ms2.data.neg, function(x) {
        x[[1]]
      })
      
      ms1.info.neg <- do.call(rbind, ms1.info.neg) %>%
        as.data.frame()
      
      ms1.info.neg$file <- basename(ms1.info.neg$file)
      
      ms2.info.neg <- lapply(ms2.data.neg, function(x) {
        x[[2]]
      })
      
      rm(list = "ms2.data.neg")
      message(crayon::red("OK."))
    } else{
      ms1.info.neg <- NULL
    }
    
    ###---------------------------------------------------------------------------
    message(crayon::green("Matching metabolites with MS2 spectra (positive)..."))
    
    if (!is.null(ms1.info.pos)) {
      match.result.pos <-
        masstools::mz_rt_match(
          data1 = as.data.frame(metabolite.info[, c("mz.pos", "RT")]),
          data2 = ms1.info.pos[, c(2, 3)],
          mz.tol = mz.tol,
          rt.tol = rt.tol,
          rt.error.type = "abs"
        )
      
      match.result.pos <-
        data.frame(match.result.pos,
                   "file" = ms1.info.pos$file[match.result.pos[, 2]],
                   stringsAsFactors = FALSE)
      
      if (nrow(match.result.pos) == 0) {
        warning("No metabolites matched MS2 spectra.")
      }
      
      unique.idx1 <- unique(match.result.pos[, 1])
      
      spectra.pos <-
        pbapply::pblapply(unique.idx1, function(idx) {
          temp.match.result.pos <-
            match.result.pos[which(match.result.pos$Index1 == idx), , drop = FALSE]
          if (nrow(temp.match.result.pos) == 0) {
            return(NULL)
          }
          temp.submitter <- metabolite.info$Submitter[idx]
          if (!is.na(temp.submitter) &
              length(grep(temp.submitter, temp.match.result.pos[, 9])) > 0) {
            temp.match.result.pos <-
              temp.match.result.pos[grep(temp.submitter, temp.match.result.pos[, 9]), ]
          }
          
          if (nrow(temp.match.result.pos) == 0) {
            return(NULL)
          }
          
          if (nrow(temp.match.result.pos) == 1) {
            temp.ms2.pos <- ms2.info.pos[temp.match.result.pos[1, 2]]
            names(temp.ms2.pos) <-
              stringr::str_extract(string = temp.match.result.pos[1, 9], pattern = "NCE[0-9]{1,3}")
            return(temp.ms2.pos)
          }
          
          unique.file.name <-
            unique(temp.match.result.pos$file)
          
          temp.ms2.pos <-
            lapply(unique.file.name, function(temp.name) {
              temp.x <-
                temp.match.result.pos[which(temp.match.result.pos$file == temp.name), , drop = FALSE]
              temp.idx <-
                which.max(unlist(lapply(ms2.info.pos[temp.x[, 2]], function(y) {
                  sum(y[, 2])
                })))
              ms2.info.pos[[temp.x[temp.idx, 2]]]
            })
          
          names(temp.ms2.pos) <-
            stringr::str_extract(string = unique.file.name, pattern = "NCE[0-9]{1,3}")
          temp.ms2.pos
        })
      
      names(spectra.pos) <-
        metabolite.info$Lab.ID[unique.idx1]
      
      if (length(spectra.pos) == 0) {
        spectra.pos <- list()
      } else{
        spectra.pos <-
          spectra.pos[which(!unlist(lapply(spectra.pos, is.null)))]
      }
      message(crayon::red("OK."))
    } else{
      spectra.pos <- NULL
    }
    
    ###---------------------------------------------------------------------------
    message(crayon::green("Matching metabolites with MS2 spectra (negative)..."))
    if (!is.null(ms1.info.neg)) {
      match.result.neg <-
        masstools::mz_rt_match(
          data1 = as.data.frame(metabolite.info[, c("mz.neg", "RT")]),
          data2 = ms1.info.neg[, c(2, 3)],
          mz.tol = mz.tol,
          rt.tol = rt.tol,
          rt.error.type = "abs"
        )
      
      match.result.neg <- data.frame(match.result.neg,
                                     "file" = ms1.info.neg$file[match.result.neg[, 2]],
                                     stringsAsFactors = FALSE)
      
      if (nrow(match.result.neg) == 0) {
        warning("No metabolites matched MS2 spectra.")
      }
      
      unique.idx1 <- unique(match.result.neg[, 1])
      
      spectra.neg <-
        pbapply::pblapply(unique.idx1, function(idx) {
          temp.match.result.neg <-
            match.result.neg[which(match.result.neg$Index1 == idx), , drop = FALSE]
          if (nrow(temp.match.result.neg) == 0) {
            return(NULL)
          }
          
          temp.submitter <- metabolite.info$Submitter[idx]
          if (!is.na(temp.submitter) &
              length(grep(temp.submitter, temp.match.result.neg[, 9])) > 0) {
            temp.match.result.neg <-
              temp.match.result.neg[grep(temp.submitter, temp.match.result.neg[, 9]), ]
          }
          
          if (nrow(temp.match.result.neg) == 0) {
            return(NULL)
          }
          
          if (nrow(temp.match.result.neg) == 1) {
            temp.ms2.neg <- ms2.info.neg[temp.match.result.neg[1, 2]]
            names(temp.ms2.neg) <-
              stringr::str_extract(string = temp.match.result.neg[1, 9], pattern = "NCE[0-9]{1,3}")
            return(temp.ms2.neg)
          }
          
          unique.file.name <-
            unique(temp.match.result.neg$file)
          
          temp.ms2.neg <-
            lapply(unique.file.name, function(temp.name) {
              temp.x <-
                temp.match.result.neg[which(temp.match.result.neg$file == temp.name), , drop = FALSE]
              temp.idx <-
                which.max(unlist(lapply(ms2.info.neg[temp.x[, 2]], function(y) {
                  sum(y[, 2])
                })))
              ms2.info.neg[[temp.x[temp.idx, 2]]]
            })
          
          names(temp.ms2.neg) <-
            stringr::str_extract(string = unique.file.name, pattern = "NCE[0-9]{1,3}")
          temp.ms2.neg
        })
      
      names(spectra.neg) <-
        metabolite.info$Lab.ID[unique.idx1]
      
      if (length(spectra.neg) == 0) {
        spectra.neg <- list()
      } else{
        spectra.neg <-
          spectra.neg[which(!unlist(lapply(spectra.neg, is.null)))]
      }
      
      message(crayon::red("OK."))
    } else{
      spectra.neg <- NULL
    }
    
    if (is.null(spectra.pos) & is.null(spectra.neg)) {
      Spectra <- list()
    } else{
      Spectra <- list("Spectra.positive" = spectra.pos,
                      "Spectra.negative" = spectra.neg)
    }
    
    database.info <- list(
      "Version" = version,
      "Source" = source,
      "Link" = link,
      "Creater" = creater,
      "Email" = email,
      "RT" = rt
    )
    
    spectra.info <- as.data.frame(metabolite.info)
    rm(list = "metabolite.info")
    
    msDatabase0.0.1 <- new(
      Class = "databaseClass",
      database.info = database.info,
      spectra.info = spectra.info,
      spectra.data = Spectra
    )
    
    msDatabase0.0.1@database.info$RT <-
      ifelse(all(is.na(msDatabase0.0.1@spectra.info$RT)), FALSE, TRUE)
    message(crayon::bgRed("All done!"))
    return(msDatabase0.0.1)
  }









##------------------------------------------------------------------------------

# # ####MassBank library
# setwd(masstools::get_project_wd())
# setwd("other_files/all_ms2_database/massbank/2021_3_5/")
# massbank_database = construct_massbank_database(file = "MassBank_NIST.msp")
#
# mz <- 167.035
# rt <- 450
# ms2 <- data.frame(
#   mz = c(91.0187,
#          95.0138,
#          95.0503,
#          99.0088,
#          108.0217,
#          121.0297,
#          123.045,
#          152.0113,
#          167.0348
#   ),
#   intensity = c(75.6,
#                 81.8,
#                 61.7,
#                 69.9,
#                 4877.1,
#                 46.5,
#                 8682.4,
#                 25133.4,
#                 44848.5
#   ),
#   stringsAsFactors = FALSE
# )
#
# save(massbank_database, file = "massbank_database", compress = "xz")
#
# annotation_result <-
#   metid::identify_single_peak(ms1.mz = mz,
#                               ms1.rt = rt,
#                               ms2 = ms2,
#                               ms1.match.ppm = 15,
#                               rt.match.tol = 30,
#                               ms2.match.tol = 0.5,
#                               database = "massbank_database",
#                               path = ".", polarity = "negative")
#
# which_has_identification(annotation_result)
# get_identification_table(annotation_result, type = "new")
# ms2plot(object = annotation_result, database = massbank_database, which.peak = "mz167.035rt450")


#' Construct a Spectral Database from MassBank Data
#'
#' This function constructs a spectral database from MassBank data, specifically formatted as a `databaseClass` object. The database contains metabolite information and MS2 spectral data for both positive and negative ionization modes. 
#'
#' @param file A character string specifying the path to the MassBank data file (MSP format).
#' @param only.remain.ms2 A logical value indicating whether to retain only the metabolites with MS2 spectra. Defaults to `TRUE`.
#' @param path A character string specifying the path where the database will be saved. Defaults to the current working directory (`"."`).
#' @param version A character string specifying the version of the database. Defaults to `"0.0.1"`.
#' @param source A character string specifying the source of the database. Defaults to `"MassBank"`.
#' @param link A character string specifying the URL of the source database. Defaults to `"https://massbank.eu/MassBank/"`.
#' @param creater A character string specifying the creator of the database. Defaults to `"Xiaotao Shen"`.
#' @param email A character string specifying the email of the creator. Defaults to `"xiaotao.shen@outlook.com"`.
#' @param rt A logical value indicating whether retention time (RT) information is available. Defaults to `FALSE`.
#' @param threads An integer specifying the number of threads to use for parallel processing. Defaults to `5`.
#'
#' @return A `databaseClass` object containing the MassBank database, with metabolite information and MS2 spectral data.
#'
#' @details
#' The function reads MassBank data in MSP format and constructs a `databaseClass` object. It processes both positive and negative ionization modes, and optionally filters out metabolites that do not have MS2 spectra if `only.remain.ms2` is set to `TRUE`.
#'
#' The MassBank data is organized into metabolite information and corresponding MS2 spectra. The spectra are stored in the `Spectra.positive` and `Spectra.negative` slots based on their ionization modes.
#'
#' @examples
#' \dontrun{
#' # Construct a database from a MassBank file
#' massbank_db <- construct_massbank_database(file = "MassBank.msp", only.remain.ms2 = TRUE)
#' }
#'
#' @importFrom purrr map map2
#' @importFrom dplyr select mutate everything
#' @importFrom crayon bgRed
#' @seealso \code{\link{databaseClass}}, \code{\link{read_msp_mona}}
#' @export

construct_massbank_database <-
  function(file,
           only.remain.ms2 = TRUE,
           path = ".",
           version = "0.0.1",
           source = "MassBank",
           link = "https://massbank.eu/MassBank/",
           creater = "Xiaotao Shen",
           email = "xiaotao.shen@outlook.com",
           rt = FALSE,
           threads = 5) {
    massbank_database = read_msp_mona(file = file)
    
    all_metabolite_names =
      purrr::map(massbank_database, function(x) {
        rownames(x$info)
      }) %>%
      unlist() %>%
      unique()
    
    metabolite_info =
      massbank_database %>%
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
        metabolite_info[remain_idx, ]
      massbank_database = massbank_database[remain_idx]
    }
    
    metabolite_info =
      metabolite_info %>%
      dplyr::select(
        Compound.name = Name,
        mz = ExactMass,
        Formula,
        MassBank.ID = `DB#`,
        dplyr::everything()
      )
    
    metabolite_info =
      metabolite_info %>%
      dplyr::mutate(
        Lab.ID = paste("MassBank", seq_len(nrow(metabolite_info)), sep = "_"),
        RT = NA,
        CAS.ID = NA,
        HMDB.ID = NA,
        KEGG.ID = NA,
        mz.pos = NA,
        mz.neg = NA,
        Submitter = "MassBank",
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
    positive_idx = which(metabolite_info$Ion_mode == "POSITIVE")
    negative_idx = which(metabolite_info$Ion_mode == "NEGATIVE")
    
    Spectra.positive = massbank_database[positive_idx]
    Spectra.negative = massbank_database[negative_idx]
    
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
    message(crayon::bgRed("All done!"))
    return(database)
  }












#' Construct a Spectral Database from MoNA Data
#'
#' This function constructs a spectral database from MoNA (MassBank of North America) data, specifically formatted as a `databaseClass` object. The database contains metabolite information and MS2 spectral data for both positive and negative ionization modes.
#'
#' @param file A character string specifying the path to the MoNA data file (MSP format).
#' @param only.remain.ms2 A logical value indicating whether to retain only the metabolites with MS2 spectra. Defaults to `TRUE`.
#' @param path A character string specifying the path where the database will be saved. Defaults to the current working directory (`"."`).
#' @param version A character string specifying the version of the database. Defaults to `"0.0.1"`.
#' @param source A character string specifying the source of the database. Defaults to `"MoNA"`.
#' @param link A character string specifying the URL of the source database. Defaults to `"https://mona.fiehnlab.ucdavis.edu/"`.
#' @param creater A character string specifying the creator of the database. Defaults to `"Xiaotao Shen"`.
#' @param email A character string specifying the email of the creator. Defaults to `"xiaotao.shen@outlook.com"`.
#' @param rt A logical value indicating whether retention time (RT) information is available. Defaults to `FALSE`.
#' @param threads An integer specifying the number of threads to use for parallel processing. Defaults to `5`.
#'
#' @return A `databaseClass` object containing the MoNA database, with metabolite information and MS2 spectral data.
#'
#' @details
#' The function reads MoNA data in MSP format and constructs a `databaseClass` object. It processes both positive and negative ionization modes and optionally filters out metabolites that do not have MS2 spectra if `only.remain.ms2` is set to `TRUE`.
#'
#' The MoNA data is organized into metabolite information and corresponding MS2 spectra. The spectra are stored in the `Spectra.positive` and `Spectra.negative` slots based on their ionization modes.
#'
#' @examples
#' \dontrun{
#' # Construct a database from a MoNA file
#' mona_db <- construct_mona_database(file = "MoNA.msp", only.remain.ms2 = TRUE)
#' }
#'
#' @importFrom purrr map map2
#' @importFrom dplyr select mutate everything
#' @importFrom crayon bgRed
#' @seealso \code{\link{databaseClass}}, \code{\link{read_msp_mona}}
#' @export

# # # ####MoNA library
# setwd(masstools::get_project_wd())
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

construct_mona_database <-
  function(file,
           only.remain.ms2 = TRUE,
           path = ".",
           version = "0.0.1",
           source = "MoNA",
           link = "https://mona.fiehnlab.ucdavis.edu/",
           creater = "Xiaotao Shen",
           email = "xiaotao.shen@outlook.com",
           rt = FALSE,
           threads = 5) {
    if (missing(file)) {
      stop("please prove file.")
    }
    if (all(dir(path) != file)) {
      stop(file, " is not in ", path)
    }
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
        metabolite_info[remain_idx, ]
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
