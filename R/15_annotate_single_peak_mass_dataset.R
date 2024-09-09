#' Annotate a Single Peak in a mass_dataset Object
#'
#' This function annotates a single peak in a `mass_dataset` object using MS1 and optionally MS2 data. It allows for matching with a reference database based on m/z, retention time (RT), and MS2 fragmentation patterns.
#'
#' @param object A `mass_dataset` object containing the peak data to annotate.
#' @param variable_id The ID of the peak to annotate. Either `variable_id` or `variable_index` must be provided.
#' @param variable_index The index of the peak to annotate. Either `variable_id` or `variable_index` must be provided.
#' @param based_on_rt Logical, if `TRUE` (default), retention time will be used for matching.
#' @param based_on_ms2 Logical, if `TRUE` (default), MS2 spectra will be used for matching if available.
#' @param add_to_annotation_table Logical, if `TRUE`, the annotation results will be added to the annotation table in the object. Defaults to `FALSE`.
#' @param ms1.match.ppm Numeric, mass accuracy threshold for MS1 matching in parts per million (ppm). Defaults to `25`.
#' @param ms2.match.ppm Numeric, mass accuracy threshold for MS2 matching in ppm. Defaults to `30`.
#' @param mz.ppm.thr Numeric, m/z threshold in ppm for matching MS1 and MS2. Defaults to `400`.
#' @param ms2.match.tol Numeric, tolerance for MS2 fragment ion matching. Defaults to `0.5`.
#' @param fraction.weight Numeric, weight for the MS2 fragmentation score. Defaults to `0.3`.
#' @param dp.forward.weight Numeric, weight for the forward dot product in MS2 matching. Defaults to `0.6`.
#' @param dp.reverse.weight Numeric, weight for the reverse dot product in MS2 matching. Defaults to `0.1`.
#' @param rt.match.tol Numeric, retention time matching tolerance in seconds. Defaults to `30`.
#' @param polarity Character, ionization mode, either `"positive"` or `"negative"`. Defaults to `"positive"`.
#' @param ce Character, collision energy for MS2 matching. Defaults to `"all"`.
#' @param column Character, chromatographic column type, either `"rp"` (reverse phase) or `"hilic"`. Defaults to `"rp"`.
#' @param ms1.match.weight Numeric, weight of MS1 matching in total score calculation. Defaults to `0.25`.
#' @param rt.match.weight Numeric, weight of RT matching in total score calculation. Defaults to `0.25`.
#' @param ms2.match.weight Numeric, weight of MS2 matching in total score calculation. Defaults to `0.5`.
#' @param total.score.tol Numeric, threshold for the total score. Defaults to `0.5`.
#' @param candidate.num Numeric, the number of top candidates to retain for each peak. Defaults to `3`.
#' @param database A `databaseClass` object containing the reference spectral database for annotation.
#' @param threads Numeric, the number of threads to use for parallel processing. Defaults to `3`.
#'
#' @return Either the `mass_dataset` object with updated annotation table or a data frame containing the annotation results for the specified peak.
#'
#' @details
#' This function performs peak annotation using MS1 data and optionally MS2 data for a single peak in a `mass_dataset` object. The matching process is based on m/z, retention time, and MS2 spectra comparison with a reference database. The results can either be returned as a data frame or added to the annotation table in the `mass_dataset` object.
#'
#' @examples
#' \dontrun{
#' # Annotate a single peak using MS1 and MS2 data
#' annotation <- annotate_single_peak_mass_dataset(
#'   object = mass_object,
#'   variable_id = "P001",
#'   database = reference_database,
#'   based_on_rt = TRUE,
#'   based_on_ms2 = TRUE
#' )
#' }
#'
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @export

annotate_single_peak_mass_dataset <-
  function(object,
           variable_id,
           variable_index,
           based_on_rt = TRUE,
           based_on_ms2 = TRUE,
           add_to_annotation_table = FALSE,
           ms1.match.ppm = 25,
           ms2.match.ppm = 30,
           mz.ppm.thr = 400,
           ms2.match.tol = 0.5,
           fraction.weight = 0.3,
           dp.forward.weight = 0.6,
           dp.reverse.weight = 0.1,
           rt.match.tol = 30,
           polarity = c("positive", "negative"),
           ce = "all",
           column = c("rp", "hilic"),
           ms1.match.weight = 0.25,
           rt.match.weight = 0.25,
           ms2.match.weight = 0.5,
           total.score.tol = 0.5,
           candidate.num = 3,
           database,
           threads = 3) {
    massdataset::check_object_class(object = object, class = "mass_dataset")
    
    ###check parameters
    if (!is.numeric(ms1.match.ppm)) {
      stop("ms1.match.ppm should be numeric.\n")
    } else{
      if (ms1.match.ppm <= 0 | ms1.match.ppm >= 500) {
        stop("ms1.match.ppm should > 0 and < 500\n")
      }
    }
    
    if (!is.numeric(ms2.match.ppm)) {
      stop("ms2.match.ppm should be numeric.\n")
    } else{
      if (ms2.match.ppm <= 0 | ms2.match.ppm >= 500) {
        stop("ms2.match.ppm should > 0 and < 500\n")
      }
    }
    
    if (!is.numeric(candidate.num)) {
      stop("candidate.num should be numeric.\n")
    } else{
      if (candidate.num <= 0) {
        stop("candidate.num should > 0.\n")
      }
    }
    
    if (!based_on_rt) {
      rt.match.tol = 1000000
    }
    
    if (is.na(rt.match.tol)) {
      rt.match.tol = 1000000
    }
    
    ###Check data
    if (missing(database)) {
      stop("No database is provided.\n")
    }
    
    ##parameter specification
    polarity <- match.arg(polarity)
    column <- match.arg(column)
    
    if (!is(database, "databaseClass")) {
      stop("database should be databaseClass object.\n")
    }
    
    database.name = paste(database@database.info$Source,
                          database@database.info$Version,
                          sep = "_")
    
    ######check variable_id and variable_index
    if (missing(variable_id) & missing(variable_index)) {
      stop("provide variable_id or variable_index.\n")
    }
    
    if (!missing(variable_id)) {
      purrr::walk(
        variable_id,
        .f = function(temp_variable_id) {
          if (!temp_variable_id %in% object@variable_info$variable_id) {
            stop(paste(temp_variable_id, "is not in variable_info.\n"))
          }
        }
      )
      variable_index = match(variable_id, object@variable_info$variable_id)
    } else{
      purrr::walk(
        variable_index,
        .f = function(temp_variable_index) {
          if (temp_variable_index <= 0 |
              temp_variable_index > nrow(object@variable_info)) {
            stop(
              "variable_index ",
              temp_variable_index,
              " should be range from 1 to ",
              nrow(object@variable_info)
            )
          }
        }
      )
    }
    
    temp_object = object[variable_index, ]
    
    ######NO MS2 in object
    if (length(object@ms2_data) == 0 | !based_on_ms2) {
      message(crayon::yellow("No MS2 data in object, so only use mz and/or RT for matching."))
      annotation_result =
        mzIdentify_mass_dataset(
          object = temp_object,
          rt.match.tol = rt.match.tol,
          ms1.match.ppm = ms1.match.ppm,
          polarity = polarity,
          column = column,
          candidate.num = candidate.num,
          database = database,
          threads = threads
        )
      annotation_result$SS = NA
    }
    
    ######MS2 in object
    if (length(object@ms2_data) > 0 & based_on_ms2) {
      spectra_pos_number =
        database@spectra.data[['Spectra.positive']] %>%
        length()
      
      spectra_neg_number =
        database@spectra.data[['Spectra.negative']] %>%
        length()
      
      if (polarity == "positive") {
        spectra_number = spectra_pos_number
      } else{
        spectra_number = spectra_neg_number
      }
      
      ######NO MS2 in database
      if (spectra_number == 0) {
        message(crayon::yellow(
          "No MS2 data in database, so only use mz and/or RT for matching."
        ))
        annotation_result =
          mzIdentify_mass_dataset(
            object = temp_object,
            rt.match.tol = rt.match.tol,
            ms1.match.ppm = ms1.match.ppm,
            polarity = polarity,
            column = column,
            candidate.num = candidate.num,
            database = database,
            threads = threads
          )
        annotation_result$SS = NA
      } else{
        ######MS2 in database
        annotation_result =
          metIdentify_mass_dataset(
            object = temp_object,
            ms1.match.ppm = ms1.match.ppm,
            ms2.match.ppm = ms2.match.ppm,
            mz.ppm.thr = mz.ppm.thr,
            ms2.match.tol = ms2.match.tol,
            fraction.weight = fraction.weight,
            dp.forward.weight = dp.forward.weight,
            dp.reverse.weight = dp.reverse.weight,
            rt.match.tol = rt.match.tol,
            polarity = polarity,
            ce = ce,
            column = column,
            ms1.match.weight = ms1.match.weight,
            rt.match.weight = rt.match.weight,
            ms2.match.weight = ms2.match.weight,
            total.score.tol = total.score.tol,
            candidate.num = candidate.num,
            database = database,
            threads = threads
          )
      }
    }
    
    Level =
      annotation_result %>%
      dplyr::select(RT.error, SS) %>%
      t() %>%
      as.data.frame() %>%
      purrr::map(function(x) {
        x = as.numeric(x)
        if (sum(is.na(x)) == 2) {
          return(3)
        }
        
        if (sum(is.na(x)) == 1) {
          return(2)
        }
        
        if (sum(is.na(x)) == 0) {
          return(1)
        }
      }) %>%
      unlist()
    
    annotation_result$Level = Level
    
    # browser()
    
    annotation_result =
      annotation_result %>%
      dplyr::arrange(variable_id, Level, dplyr::desc(Total.score))
    
    annotation_result =
      annotation_result %>%
      dplyr::filter(variable_id %in% temp_object@variable_info$variable_id)
    
    
    annotation_result =
      annotation_result[, c(
        "variable_id",
        "ms2_files_id",
        "ms2_spectrum_id",
        "Compound.name",
        "CAS.ID",
        "HMDB.ID",
        "KEGG.ID",
        "Lab.ID",
        "Adduct",
        "mz.error",
        "mz.match.score",
        "RT.error",
        "RT.match.score",
        "CE",
        "SS",
        "Total.score",
        "Database",
        "Level"
      )]
    
    
    if (add_to_annotation_table) {
      if (nrow(object@annotation_table) == 0) {
        object@annotation_table = as.data.frame(annotation_result)
      } else{
        object@annotation_table =
          rbind(object@annotation_table, annotation_result) %>%
          dplyr::arrange(variable_id, Level, dplyr::desc(Total.score))
        
        ###only remain top annotations
        object@annotation_table =
          object@annotation_table %>%
          dplyr::group_by(variable_id) %>%
          dplyr::slice_head(n = candidate.num) %>%
          dplyr::ungroup() %>%
          dplyr::distinct(.keep_all = TRUE) %>%
          as.data.frame()
      }
      
      ###processing information
      process_info = object@process_info
      
      parameter <- new(
        Class = "tidymass_parameter",
        pacakge_name = "metid",
        function_name = "identify_metabolites_mass_dataset()",
        parameter = list(
          ms1.match.ppm = ms1.match.ppm,
          ms2.match.ppm = ms2.match.ppm,
          mz.ppm.thr = mz.ppm.thr,
          ms2.match.tol = ms2.match.tol,
          fraction.weight = fraction.weight,
          dp.forward.weight = dp.forward.weight,
          dp.reverse.weight = dp.reverse.weight,
          rt.match.tol = rt.match.tol,
          polarity = polarity,
          ce = ce,
          column = column,
          ms1.match.weight = ms1.match.weight,
          rt.match.weight = rt.match.weight,
          ms2.match.weight = ms2.match.weight,
          total.score.tol = total.score.tol,
          candidate.num = candidate.num,
          database = database.name,
          threads = threads
        ),
        time = Sys.time()
      )
      
      if (all(names(process_info) != "identify_metabolites_mass_dataset")) {
        process_info$identify_metabolites_mass_dataset = parameter
      } else{
        process_info$identify_metabolites_mass_dataset = c(process_info$identify_metabolites_mass_dataset,
                                                           parameter)
      }
      
      object@process_info = process_info
      return(object)
    } else{
      return(annotation_result)
    }
  }
