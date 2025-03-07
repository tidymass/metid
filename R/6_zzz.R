.onAttach <- function(libname, pkgname) {
  # needed <- core[!is_attached(core)]
  # if (length(needed) == 0)
  #   return()
  #
  crayon::num_colors(TRUE)
  metid_attach()
  #
  # if (!"package:conflicted" %in% search()) {
  #   x <- metid_conflicts()
  #   msg(metid_conflict_message(x), startup = TRUE)
  # }
  packageStartupMessage(paste0("metid ", metid_version, " (", update_date, ')'))
}

is_attached <- function(x) {
  paste0("package:", x) %in% search()
  
}


####database for annotation
#' Extract the Name of the Database
#'
#' This function extracts and formats the name of a database from a `databaseClass` object by combining the source and version information.
#'
#' @param database A `databaseClass` object containing database information.
#'
#' @return A character string representing the name of the database, formatted as `Source_Version`. Returns `NULL` if the input is not a `databaseClass` object.
#'
#' @details
#' The function extracts the `Source` and `Version` fields from the `database.info` slot of the `databaseClass` object and combines them using an underscore (`_`). This formatted string represents the name of the database.
#'
#' @examples
#' \dontrun{
#' # Load a sample database
#' my_database <- load_database("path/to/database")
#'
#' # Extract the name of the database
#' db_name <- extract_database_name(my_database)
#' print(db_name)  # Output might be "HMDB_4.0"
#' }
#'
#'
#' @export


extract_database_name <-
  function(database) {
    ##database is a databaseClass object
    if (is(database, "databaseClass")) {
      paste(database@database.info$Source,
            database@database.info$Version,
            sep = "_")
    }
  }


#' Extract MS1 Information from a Database
#'
#' This function extracts MS1 data (m/z and retention time) from a `databaseClass` object and returns it as a filtered data frame. It ensures that the m/z values are numeric and removes rows where the m/z values are missing (`NA`).
#'
#' @param database A `databaseClass` object containing MS1 information in the `spectra.info` slot.
#'
#' @return A data frame with MS1 information, where the `mz` and `RT` columns are numeric and rows with missing (`NA`) m/z values are removed. If the input is not a `databaseClass` object, the function returns `NULL`.
#'
#' @details
#' The function extracts the `spectra.info` slot from the `databaseClass` object, converts the `mz` and `RT` columns to numeric, and filters out rows with missing m/z values. This ensures that only valid MS1 data is returned.
#'
#' @examples
#' \dontrun{
#' # Load a sample database
#' my_database <- load_database("path/to/database")
#'
#' # Extract MS1 data from the database
#' ms1_data <- extract_ms1_database(my_database)
#' head(ms1_data)
#' }
#'
#'
#' @export

extract_ms1_database <-
  function(database) {
    ###database is a databaseClass object
    if (is(database, "databaseClass")) {
      database@spectra.info %>%
        dplyr::mutate(mz = as.numeric(mz)) %>%
        dplyr::mutate(RT = as.numeric(RT)) %>%
        dplyr::filter(!is.na(mz))
    }
  }

#' Extract MS2 Spectra from a Database
#'
#' This function extracts MS2 spectra from a `databaseClass` object based on the specified polarity and collision energy (CE) values. It supports both positive and negative ionization modes, and allows users to filter by specific CE values or use all available CE values in the database.
#'
#' @param database A `databaseClass` object containing MS2 spectra data.
#' @param polarity Character. The ionization mode, either `"positive"` or `"negative"`. Default is `"positive"`.
#' @param ce Character or vector. The collision energy (CE) values to extract MS2 spectra for. Set to `"all"` to use all CE values in the database. Default is `"all"`.
#'
#' @return A list containing the extracted MS2 spectra for the specified polarity and CE values. If no matching spectra are found, the function returns `NULL`.
#'
#' @details
#' The function checks the polarity and collision energy values against the available data in the provided `databaseClass` object. It extracts MS2 spectra corresponding to the specified polarity and filters them based on the given CE values. If the `ce` argument is set to `"all"`, all CE values in the database are used. If none of the requested CE values are available, the function returns `NULL`.
#'
#' @examples
#' \dontrun{
#' # Load a database
#' my_database <- load_database("path/to/database")
#'
#' # Extract all MS2 spectra in positive ion mode
#' ms2_spectra <- extract_ms2_database(database = my_database,
#' polarity = "positive", ce = "all")
#'
#' # Extract MS2 spectra for specific CE values in negative ion mode
#' ms2_spectra <- extract_ms2_database(database = my_database,
#' polarity = "negative", ce = c("10", "20"))
#' }
#'
#'
#' @export


extract_ms2_database <-
  function(database,
           polarity = c("positive", "negative"),
           ce = "all") {
    ###database is a databaseClass object
    if (!is(database, "databaseClass")) {
      stop("database should be a databaseClass object")
    }
    
    polarity <- match.arg(polarity)
    
    spectra_pos_number <-
      database@spectra.data[['Spectra.positive']] %>%
      length()
    
    spectra_neg_number <-
      database@spectra.data[['Spectra.negative']] %>%
      length()
    
    spectra_number <-
      ifelse(polarity == "positive",
             spectra_pos_number,
             spectra_neg_number)
    
    if (spectra_number == 0) {
      return(NULL)
    }
    
    ##CE values
    ## Check if CE values (collision energies) match with the database
    ce.list.pos <- unique(unlist(lapply(
      database@spectra.data$Spectra.positive, names
    )))
    ce.list.neg <- unique(unlist(lapply(
      database@spectra.data$Spectra.negative, names
    )))
    
    if (polarity == "positive") {
      ce.list <- ce.list.pos
    } else{
      ce.list <- ce.list.neg
    }
    
    ## Stop execution if provided CE values are not in the database
    if (all(!ce %in% ce.list) & !("all" %in% ce)) {
      message(
        "The CE values you can chose are: ",
        paste(ce.list, collapse = ", "),
        "\n",
        "or just set it as ",
        "'all'"
      )
      stop("All ce values you set are not in the database. Please check it.\n")
    }
    
    if (polarity == "positive") {
      spectra.data <- database@spectra.data$Spectra.positive
    } else{
      spectra.data <- database@spectra.data$Spectra.negative
    }
    
    ##get the MS2 spectra within the CE values
    if (any(ce == "all")) {
      message("\n")
      message("Use all CE values.\n")
      ce <-
        unique(unlist(lapply(spectra.data, function(x) {
          names(x)
        })))
    } else{
      spectra.data <-
        lapply(spectra.data, function(x) {
          x <-
            x[which(names(x) %in% ce)]
          if (length(x) == 0)
            return(NULL)
          return(x)
        })
    }
    
    ##remove some metabolites which have no spectra
    spectra.data <-
      spectra.data[which(!unlist(lapply(spectra.data, is.null)))]
    
    if (length(spectra.data) == 0) {
      return(NULL)
      stop("No spectra with CE: ",
           paste(ce, collapse = ", "),
           " in you database.\n")
    }
    return(spectra.data)
  }


#####mass_dataset for annotation

#' Extract MS1 Information from a mass_dataset Object
#'
#' This function extracts MS1-related data from a `mass_dataset` object, including MS2 spectrum ID, m/z, retention time (RT), file ID, and variable ID.
#'
#' @param object A `mass_dataset` object containing MS2 data.
#'
#' @return A data frame containing the following columns:
#' \describe{
#'   \item{ms2_spectrum_id}{A combined ID made up of the MS2 data ID and spectrum ID.}
#'   \item{mz}{The m/z values for the corresponding MS2 spectra.}
#'   \item{rt}{The retention time (RT) values for the corresponding MS2 spectra.}
#'   \item{ms2_files_id}{The file ID associated with each MS2 spectrum.}
#'   \item{variable_id}{The variable ID for each MS2 spectrum.}
#' }
#'
#' @details
#' The function extracts MS1-related information from the MS2 data stored within a `mass_dataset` object. It combines the MS2 data ID and spectrum ID to create a unique `ms2_spectrum_id`, and returns the corresponding m/z, RT, file ID, and variable ID in a data frame.
#'
#' @examples
#' \dontrun{
#' # Load a sample mass dataset
#' my_dataset <- load_mass_dataset("path/to/dataset")
#'
#' # Extract MS1 information
#' ms1_info <- extract_ms1_info(my_dataset)
#' head(ms1_info)
#' }
#'
#'
#' @export

extract_ms1_info <-
  function(object) {
    if (!is(object, "mass_dataset")) {
      stop("object should be a massdataset object")
    } else{
      ms1.info <-
        purrr::map2(
          .x = names(object@ms2_data),
          .y = object@ms2_data,
          .f = function(temp_ms2_data_id, temp_ms2_data) {
            data.frame(
              ms2_spectrum_id = paste(temp_ms2_data_id, temp_ms2_data@ms2_spectrum_id, sep = "_"),
              mz = temp_ms2_data@ms2_mz,
              rt = temp_ms2_data@ms2_rt,
              # ms2_files_id = temp_ms2_data@ms2_file,
              ms2_files_id = temp_ms2_data_id,
              variable_id = temp_ms2_data@variable_id
            )
          }
        ) %>%
        do.call(rbind, .) %>%
        as.data.frame()
    }
    
    ms1.info <-
      ms1.info %>%
      dplyr::filter(variable_id %in% object@variable_info$variable_id)
    return(ms1.info)
    
  }


#' Extract MS2 Information from a mass_dataset Object
#'
#' This function extracts MS2 spectral data from a `mass_dataset` object, returning the MS2 spectra with combined file and spectrum IDs.
#'
#' @param object A `mass_dataset` object containing MS2 data.
#'
#' @return A named list where each element is an MS2 spectrum (a matrix of m/z and intensity values) with names representing the combined file and spectrum IDs.
#'
#' @details
#' The function extracts the MS2 spectra from a `mass_dataset` object, combines the file and spectrum IDs to create unique identifiers, and returns a list of MS2 spectra. Each MS2 spectrum contains m/z and intensity data, and the list is named according to the file and spectrum IDs.
#'
#' @examples
#' \dontrun{
#' # Load a sample mass dataset
#' my_dataset <- load_mass_dataset("path/to/dataset")
#'
#' # Extract MS2 information
#' ms2_info <- extract_ms2_info(my_dataset)
#' head(names(ms2_info))
#' }
#'
#'
#' @export

extract_ms2_info <-
  function(object) {
    if (!is(object, "mass_dataset")) {
      stop("object should be a massdataset object")
    } else{
      ms2.info <-
        purrr::map2(
          .x = names(object@ms2_data),
          .y = object@ms2_data,
          .f = function(temp_ms2_data_id, temp_ms2_data) {
            temp_ms2_spectra <-
              temp_ms2_data@ms2_spectra
            names(temp_ms2_spectra) <-
              paste(temp_ms2_data_id,
                    names(temp_ms2_data@ms2_spectra),
                    sep = "_")
            temp_ms2_spectra
          }
        )
      
      names(ms2.info) <- NULL
      
      ms2.info <-
        ms2.info %>%
        do.call(c, .)
      
      return(ms2.info)
      
    }
  }



#' Calculate m/z Match Score
#'
#' This function calculates the m/z match score based on the m/z error and the MS1 match tolerance in parts per million (ppm). The score is calculated using a Gaussian-like function.
#'
#' @param mz.error Numeric vector. The m/z error values between experimental and reference m/z values.
#' @param ms1.match.ppm Numeric. The MS1 match tolerance in parts per million (ppm). Default is 25.
#'
#' @return A numeric vector representing the m/z match score for each m/z error. The score is calculated using the formula: \eqn{exp(-0.5 * (mz.error / ms1.match.ppm)^2)}.
#'
#' @details
#' The m/z match score is computed using a Gaussian-like function where the error is normalized by the MS1 match tolerance (`ms1.match.ppm`). Lower m/z errors result in higher match scores.
#'
#' @examples
#' \dontrun{
#' # Example m/z errors
#' mz_errors <- c(10, 15, 5, 20)
#'
#' # Calculate m/z match scores with a default tolerance of 25 ppm
#' scores <- calculate_mz_match_score(mz.error = mz_errors, ms1.match.ppm = 25)
#' print(scores)
#' }
#'
#'
#' @export


calculate_mz_match_score <-
  function(mz.error, ms1.match.ppm = 25) {
    exp(-0.5 * (mz.error / (ms1.match.ppm))^2)
  }

#' Calculate Retention Time (RT) Match Score
#'
#' This function calculates the retention time (RT) match score based on the RT error and the RT match tolerance. The score is calculated using a Gaussian-like function.
#'
#' @param RT.error Numeric vector. The RT error values between experimental and reference RT values.
#' @param rt.match.tol Numeric. The RT match tolerance in seconds. Default is 30.
#'
#' @return A numeric vector representing the RT match score for each RT error. The score is calculated using the formula: \eqn{exp(-0.5 * (RT.error / rt.match.tol)^2)}.
#'
#' @details
#' The RT match score is computed using a Gaussian-like function where the RT error is normalized by the RT match tolerance (`rt.match.tol`). Lower RT errors result in higher match scores.
#'
#' @examples
#' \dontrun{
#' # Example RT errors
#' rt_errors <- c(5, 10, 2, 25)
#'
#' # Calculate RT match scores with a default tolerance of 30 seconds
#' scores <- calculate_rt_match_score(RT.error = rt_errors, rt.match.tol = 30)
#' print(scores)
#' }
#'
#'
#' @export

calculate_rt_match_score <-
  function(RT.error, rt.match.tol = 30) {
    exp(-0.5 * (RT.error / (rt.match.tol))^2)
  }

#' Calculate the Total Annotation Score
#'
#' This function calculates the total annotation score based on the matching scores for m/z, retention time (RT), and spectral similarity (SS). The total score is calculated as a weighted sum of these matching scores, with weights determined by the `based_on` argument.
#'
#' @param mz.match.score Numeric vector. The m/z match scores.
#' @param RT.match.score Numeric vector. The retention time (RT) match scores.
#' @param SS Numeric vector. The MS2 spectral similarity (SS) scores.
#' @param ms1.match.weight Numeric. The weight assigned to the MS1 (m/z) match score. Default is 0.25.
#' @param rt.match.weight Numeric. The weight assigned to the RT match score. Default is 0.25.
#' @param ms2.match.weight Numeric. The weight assigned to the MS2 (spectral similarity) match score. Default is 0.5.
#' @param based_on Character vector. Specifies the criteria to base the total score on. Can include `"ms1"`, `"rt"`, and/or `"ms2"`. Default is `c("ms1", "rt", "ms2")`.
#'
#' @return A numeric vector representing the total score for each match, computed as the weighted sum of m/z, RT, and MS2 scores.
#'
#' @details
#' The function calculates the total score by applying weights to the individual match scores (`mz.match.score`, `RT.match.score`, `SS`). If the `based_on` argument specifies only one or two matching criteria (e.g., `"ms1"` or `"ms1"`, `"rt"`), the weights are adjusted accordingly.
#'
#' The weights for `ms1.match.weight`, `rt.match.weight`, and `ms2.match.weight` must sum to 1. If any of the matching scores (m/z, RT, SS) contain `NA` values, those scores are treated as 0 in the total score calculation.
#'
#' @examples
#' \dontrun{
#' # Example matching scores
#' mz_scores <- c(0.9, 0.8, 0.7)
#' rt_scores <- c(0.8, 0.7, 0.9)
#' ss_scores <- c(0.7, 0.9, 0.8)
#'
#' # Calculate the total score using the default weights
#' total_scores <- calculate_total_score(
#'   mz.match.score = mz_scores,
#'   RT.match.score = rt_scores,
#'   SS = ss_scores,
#'   based_on = c("ms1", "rt", "ms2")
#' )
#' print(total_scores)
#' }
#'
#'
#' @export

calculate_total_score <-
  function(mz.match.score,
           RT.match.score,
           SS,
           ms1.match.weight = 0.25,
           rt.match.weight = 0.25,
           ms2.match.weight = 0.5,
           based_on = c("ms1", "rt", "ms2")) {
    check_parameters4calculate_total_score(
      mz.match.score = mz.match.score,
      RT.match.score = RT.match.score,
      SS = SS,
      ms1.match.weight = ms1.match.weight,
      rt.match.weight = rt.match.weight,
      ms2.match.weight = ms2.match.weight
    )
    
    if (length(based_on) == 1) {
      if (based_on == "ms1") {
        ms1.match.weight <- 1
        rt.match.weight <- 0
        ms2.match.weight <- 0
      }
      
      if (based_on == "rt") {
        ms1.match.weight <- 0
        rt.match.weight <- 1
        ms2.match.weight <- 0
      }
      
      if (based_on == "ms2") {
        ms1.match.weight <- 0
        rt.match.weight <- 0
        ms2.match.weight <- 1
      }
    }
    
    if (length(based_on) == 2) {
      if ("ms1" %in% based_on & "rt" %in% based_on) {
        ms1.match.weight <- ms1.match.weight / (ms1.match.weight + rt.match.weight)
        rt.match.weight <- rt.match.weight / (ms1.match.weight + rt.match.weight)
        ms2.match.weight <- 0
      }
      
      if ("ms1" %in% based_on & "ms2" %in% based_on) {
        ms1.match.weight <- ms1.match.weight / (ms1.match.weight + ms2.match.weight)
        rt.match.weight <- 0
        ms2.match.weight <- ms2.match.weight / (ms1.match.weight + ms2.match.weight)
      }
      
      if ("rt" %in% based_on & "ms2" %in% based_on) {
        ms1.match.weight <- 0
        rt.match.weight <- rt.match.weight / (rt.match.weight + ms2.match.weight)
        ms2.match.weight <- ms2.match.weight / (rt.match.weight + ms2.match.weight)
      }
    }
    
    mz_score <- mz.match.score * ms1.match.weight
    rt_score <- RT.match.score * rt.match.weight
    SS_score <- SS * ms2.match.weight
    
    mz_score[is.na(mz_score)] <- 0
    rt_score[is.na(rt_score)] <- 0
    SS_score[is.na(SS_score)] <- 0
    
    total.score <-
      mz_score + rt_score + SS_score
    
    return(total.score)
  }

#' Calculate Confidence Level for Metabolite Annotations
#'
#' This function calculates the confidence level for each annotation based on the presence or absence of m/z error, retention time (RT) error, and spectral similarity (SS).
#'
#' @param annotation_result A data frame containing the annotation results, including the columns `mz.error`, `RT.error`, and `SS`.
#'
#' @details
#' The function assigns a confidence level to each annotation according to the following rules:
#' \describe{
#'   \item{1}{All of `mz.error`, `RT.error`, and `SS` are available (not NA).}
#'   \item{2}{Any two of the three (`mz.error`, `RT.error`, or `SS`) are available.}
#'   \item{3}{Only one or none of `mz.error`, `RT.error`, or `SS` is available.}
#' }
#'
#' @return A numeric vector of confidence levels (1, 2, or 3) for each annotation in the input data.
#'
#' @examples
#' \dontrun{
#' # Assuming `my_annotation_result` is the annotation data frame:
#' confidence_levels <- calculate_confidence_level(my_annotation_result)
#' }
#'
#' @seealso \code{\link{annotate_metabolites}}, \code{\link{plot_ms2_matching}}
#' @export

calculate_confidence_level <-
  function(annotation_result) {
    annotation_result %>%
      dplyr::select(mz.error, RT.error, SS) %>%
      t() %>%
      as.data.frame() %>%
      purrr::map(function(x) {
        if (all(!is.na(x))) {
          return(1)
        }
        if (!is.na(x[3]) & !is.na(x[1])) {
          return(2)
        }
        if (!is.na(x[3]) & !is.na(x[2])) {
          return(2)
        }
        if (!is.na(x[1]) & !is.na(x[2])) {
          return(2)
        }
        return(3)
      }) %>%
      unlist()
  }


#' Load the Adduct Table Based on Polarity and Column Type
#'
#' This function loads the appropriate adduct table based on the specified ionization polarity and chromatographic column type. The adduct tables contain information about the different adducts used in mass spectrometry.
#'
#' @param polarity Character. The ionization mode, either `"positive"` or `"negative"`. Default is `"positive"`.
#' @param column Character. The chromatographic column type, either `"rp"` (reversed-phase) or `"hilic"` (hydrophilic interaction chromatography). Default is `"rp"`.
#'
#' @return A data frame representing the adduct table for the specified `polarity` and `column`.
#'
#' @details
#' The function loads the appropriate adduct table from the environment based on the ionization mode (`polarity`) and column type (`column`). It supports four combinations:
#' * Positive polarity with HILIC (`"hilic.pos"`)
#' * Positive polarity with RP (`"rp.pos"`)
#' * Negative polarity with HILIC (`"hilic.neg"`)
#' * Negative polarity with RP (`"rp.neg"`)
#'
#' The corresponding adduct table is returned as a data frame.
#'
#' @examples
#' \dontrun{
#' # Load adduct table for positive polarity and reversed-phase column
#' adduct_table <- load_adduct_table(polarity = "positive", column = "rp")
#' head(adduct_table)
#' }
#'
#'
#' @export

load_adduct_table <-
  function(polarity = c("positive", "negative"),
           column = c("rp", "hilic")) {
    polarity <- match.arg(polarity)
    column <- match.arg(column)
    
    if (polarity == "positive" & column == "hilic") {
      data("hilic.pos", envir = environment())
      adduct.table <- hilic.pos
    }
    
    if (polarity == "positive" & column == "rp") {
      data("rp.pos", envir = environment())
      adduct.table <- rp.pos
    }
    
    if (polarity == "negative" & column == "hilic") {
      data("hilic.neg", envir = environment())
      adduct.table <- hilic.neg
    }
    
    if (polarity == "negative" & column == "rp") {
      data("rp.neg", envir = environment())
      adduct.table <- rp.neg
    }
    
    return(adduct.table)
  }


#' Remove Impossible Annotations Based on Adducts
#'
#' This function removes impossible metabolite annotations based on the chemical formula and detected adducts. It checks for impossible adducts such as multiple losses of water (e.g., `-H2O`, `-2H2O`) that cannot be supported by the molecular formula's hydrogen and oxygen content.
#'
#' @param match_result A data frame containing matched metabolites, including the columns `Formula` and `Adduct`.
#'
#' @return A data frame with impossible annotations removed. The `Formula` column is dropped from the final output.
#'
#' @details
#' The function analyzes the `Formula` and `Adduct` columns to identify impossible matches. For adducts that suggest a loss of water (e.g., `-H2O`, `-2H2O`), the function compares the number of hydrogens and oxygens in the molecular formula to see if the loss of water is feasible. If the molecular formula cannot support the loss, the annotation is removed.
#'
#' Specifically, the function:
#' * Extracts the number of hydrogens (`H`) and oxygens (`O`) from the molecular formula.
#' * Calculates whether the number of hydrogens and oxygens can support the indicated loss of water in the adduct.
#' * Removes annotations where the molecular formula has insufficient hydrogen or oxygen atoms for the adduct.
#'
#' @examples
#' \dontrun{
#' # Example data frame with matched metabolites
#' match_result <- data.frame(
#'   Formula = c("C6H12O6", "C5H10O5", "C7H14O7"),
#'   Adduct = c("(M-H2O)+", "(M-2H2O)+", "(M+H)+"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Remove impossible annotations
#' cleaned_results <- remove_impossible_annotations(match_result)
#' print(cleaned_results)
#' }
#'
#' @export


remove_impossible_annotations <-
  function(match_result) {
    ####remove some impossible adducts
    adduct_check <-
      match_result %>%
      dplyr::select(Formula, Adduct) %>%
      dplyr::filter(!is.na(Formula)) %>%
      dplyr::distinct(Formula, Adduct) %>%
      dplyr::filter(stringr::str_detect(Adduct, "(-H2O)|(-2H2O)")) %>%
      dplyr::mutate(minus_h2o_number = stringr::str_extract(Adduct, "(-H2O)|(-2H2O)")) %>%
      dplyr::mutate(minus_h2o_number = stringr::str_replace(minus_h2o_number, "(H2O)", "")) %>%
      dplyr::mutate(minus_h2o_number = stringr::str_replace(minus_h2o_number, "-", "")) %>%
      dplyr::mutate(minus_h2o_number = as.numeric(minus_h2o_number))
    
    adduct_check$minus_h2o_number[is.na(adduct_check$minus_h2o_number)] <-
      1
    
    extract_element_count <-
      function(formula, element = "H") {
        match <- regmatches(formula, regexec(paste0(element, "(\\d*)"), formula))
        
        # Extract matched value, if missing assume '1', otherwise convert to numeric
        if (length(match[[1]]) > 1 && match[[1]][2] != "") {
          return(as.numeric(match[[1]][2]))
        } else if (grepl("O", formula)) {
          return(1)
        } else {
          return(0)
        }
      }
    
    adduct_check$h_number <-
      sapply(adduct_check$Formula, extract_element_count, element = "H")
    adduct_check$o_number <-
      sapply(adduct_check$Formula, extract_element_count, element = "O")
    # adduct_check <-
    #   adduct_check %>%
    #   dplyr::mutate(
    #     h_number = stringr::str_extract(Formula, "H[0-9]{1,2}") %>%
    #       stringr::str_replace("H", "") %>%
    #       as.numeric()
    #   ) %>%
    #   dplyr::mutate(
    #     o_number = stringr::str_extract(Formula, "O[0-9]{0,2}") %>%
    #     stringr::str_replace("O", "") %>%
    #     as.numeric())
    #
    adduct_check$h_number[is.na(adduct_check$h_number)] <- 0
    adduct_check$o_number[is.na(adduct_check$o_number)] <- 0
    
    adduct_check <-
      adduct_check %>%
      dplyr::mutate(h_diff = h_number - minus_h2o_number * 2,
                    o_diff = o_number - minus_h2o_number) %>%
      dplyr::filter(h_diff < 0 | o_diff < 0)
    
    if (nrow(adduct_check) > 0) {
      match_result <-
        match_result %>%
        dplyr::anti_join(adduct_check, by = c("Formula", "Adduct"))
    }
    match_result <-
      match_result %>%
      dplyr::select(-Formula) %>%
      as.data.frame()
    
    return(match_result)
  }