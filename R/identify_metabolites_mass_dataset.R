#' @title Identify metabolites based on MS1 or MS/MS database
#' @description Identify metabolites based on MS1 or MS/MS database.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param object A mass_dataset class obejct.
#' @param ms1.match.ppm Precursor match ppm tolerance.
#' @param ms2.match.ppm Fragment ion match ppm tolerance.
#' @param mz.ppm.thr Accurate mass tolerance for m/z error calculation.
#' @param ms2.match.tol MS2 match (MS2 similarity) tolerance.
#' @param fraction.weight The weight for matched fragments.
#' @param dp.forward.weight Forward dot product weight.
#' @param dp.reverse.weight Reverse dot product weight.
#' @param remove_fragment_intensity_cutoff remove_fragment_intensity_cutoff
#' @param rt.match.tol RT match tolerance.
#' @param polarity The polarity of data, "positive"or "negative".
#' @param ce Collision energy. Please confirm the CE values in your database.
#' Default is "all".
#' @param column "hilic" (HILIC column) or "rp" (reverse phase).
#' @param ms1.match.weight The weight of MS1 match for total score calculation.
#' @param rt.match.weight The weight of RT match for total score calculation.
#' @param ms2.match.weight The weight of MS2 match for total score calculation.
#' @param total.score.tol Total score tolerance. The total score are referring
#' to MS-DIAL.
#' @param candidate.num The number of candidate.
#' @param database MS2 database name or MS database.
#' @param threads Number of threads
#' @return A metIdentifyClass object.
#' @importFrom crayon yellow green red bgRed
#' @importFrom magrittr %>%
#' @export
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/metid.html}
#' @examples
#' library(massdataset)
#' library(magrittr)
#' library(dplyr)
#' ms1_data =
#'   readr::read_csv(file.path(
#'     system.file("ms1_peak", package = "metid"),
#'     "ms1.peak.table.csv"
#'   ))
#'
#' ms1_data = data.frame(ms1_data, sample1 = 1, sample2 = 2)
#'
#' expression_data = ms1_data %>%
#'   dplyr::select(-c(name:rt))
#'
#' variable_info =
#'   ms1_data %>%
#'   dplyr::select(name:rt) %>%
#'   dplyr::rename(variable_id = name)
#'
#' sample_info =
#'   data.frame(
#'     sample_id = colnames(expression_data),
#'     injection.order = c(1, 2),
#'     class = c("Subject", "Subject"),
#'     group = c("Subject", "Subject")
#'   )
#' rownames(expression_data) = variable_info$variable_id
#'
#' object = create_mass_dataset(
#'   expression_data = expression_data,
#'   sample_info = sample_info,
#'   variable_info = variable_info
#' )
#'
#' object
#'
#' data("hmdb_ms1_database0.0.3", package = "metid")
#'
#' object1 =
#'   annotate_metabolites_mass_dataset(object = object,
#'                                     database = hmdb_ms1_database0.0.3)
#'
#' data("snyder_database_rplc0.0.3", package = "metid")
#'
#' database = snyder_database_rplc0.0.3
#'
#' object2 =
#'   annotate_metabolites_mass_dataset(object = object1,
#'                                     database = snyder_database_rplc0.0.3)
#' head(object2@annotation_table)
#' extract_variable_info(object = object)

annotate_metabolites_mass_dataset <-
  function(object,
           ms1.match.ppm = 25,
           ms2.match.ppm = 30,
           mz.ppm.thr = 400,
           ms2.match.tol = 0.5,
           fraction.weight = 0.3,
           dp.forward.weight = 0.6,
           dp.reverse.weight = 0.1,
           remove_fragment_intensity_cutoff = 0,
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
    
    if (class(database) != "databaseClass") {
      stop("database should be databaseClass object.\n")
    }
    
    database.name = paste(database@database.info$Source,
                          database@database.info$Version,
                          sep = "_")
    
    ######NO MS2 in object
    if (length(object@ms2_data) == 0) {
      message(crayon::yellow("No MS2 data in object, so only use mz and/or RT for matching."))
      annotation_result =
        mzIdentify_mass_dataset(
          object = object,
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
    if (length(object@ms2_data) > 0) {
      spectra_pos_number <-
        database@spectra.data[['Spectra.positive']] %>%
        length()
      
      spectra_neg_number <-
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
        annotation_result <-
          mzIdentify_mass_dataset(
            object = object,
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
        annotation_result <-
          metIdentify_mass_dataset(
            object = object,
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
            threads = threads,
            remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff
          )
      }
    }
    
    ###processing information
    process_info <- object@process_info
    
    parameter <- new(
      Class = "tidymass_parameter",
      pacakge_name = "metid",
      function_name = "annotate_metabolites_mass_dataset()",
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
    
    if (all(names(process_info) != "annotate_metabolites_mass_dataset")) {
      process_info$annotate_metabolites_mass_dataset <- parameter
    } else{
      process_info$annotate_metabolites_mass_dataset <-
        c(process_info$annotate_metabolites_mass_dataset,
          parameter)
    }
    
    object@process_info <- process_info
    
    if (nrow(annotation_result) == 0) {
      return(object)
    }
    
    Level <-
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
    
    annotation_result$Level <- Level
    
    annotation_result <-
      annotation_result %>%
      dplyr::arrange(variable_id, Level, dplyr::desc(Total.score))
    
    annotation_result =
      annotation_result %>%
      dplyr::filter(variable_id %in% object@variable_info$variable_id)
    
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
    
    if (nrow(object@annotation_table) == 0) {
      object@annotation_table <- annotation_result
    } else{
      object@annotation_table <-
        rbind(object@annotation_table,
              annotation_result) %>%
        dplyr::arrange(variable_id, Level, dplyr::desc(Total.score))
      
      ###only remain top annotations
      object@annotation_table =
        object@annotation_table %>%
        dplyr::group_by(variable_id) %>%
        dplyr::slice_head(n = candidate.num) %>%
        dplyr::ungroup() %>%
        dplyr::distinct(.keep_all = TRUE)
    }
    return(object)
  }
