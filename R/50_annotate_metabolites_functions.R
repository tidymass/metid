#' Annotate Metabolites
#'
#' This function annotates metabolites in the provided object based on MS1, retention time (RT), and/or MS2 spectra data using a specified database. It allows for customization of matching parameters such as m/z match tolerance, retention time tolerance, and MS2 matching criteria.
#'
#' @param object A `mass_dataset` object containing MS1, RT, and/or MS2 data.
#' @param database A `databaseClass` object used for metabolite annotation.
#' @param based_on Character vector. Specifies the matching criteria to be used for annotation. Can include `"ms1"`, `"rt"`, and/or `"ms2"`. Default is `c("ms1", "rt", "ms2")`.
#' @param polarity Character. Ionization mode, either `"positive"` or `"negative"`. Default is `"positive"`.
#' @param column Character. The chromatographic column type, either `"hilic"` or `"rp"` (reversed-phase). Default is `"hilic"`.
#' @param adduct.table A data frame specifying the adduct table for metabolite annotation. If `NULL`, a default adduct table is loaded based on polarity and column type.
#' @param ce Character. Collision energy used in MS2. Default is `"all"`.
#' @param ms1.match.ppm Numeric. Mass tolerance in parts per million (ppm) for MS1 peak matching. Default is 25.
#' @param ms2.match.ppm Numeric. Mass tolerance in ppm for MS2 peak matching. Default is 30.
#' @param mz.ppm.thr Numeric. m/z threshold for ppm calculation. Default is 400.
#' @param ms2.match.tol Numeric. Retention time tolerance for MS2 fragment matching. Default is 0.5.
#' @param fraction.weight Numeric. Weight for the fraction of matched fragments in MS2 spectra. Default is 0.3.
#' @param dp.forward.weight Numeric. Weight for the forward dot product score in MS2 matching. Default is 0.6.
#' @param dp.reverse.weight Numeric. Weight for the reverse dot product score in MS2 matching. Default is 0.1.
#' @param rt.match.tol Numeric. Retention time matching tolerance in seconds. Default is 30.
#' @param ms1.match.weight Numeric. Weight for MS1 matching score in the overall annotation score. Default is 0.25.
#' @param rt.match.weight Numeric. Weight for retention time matching score in the overall annotation score. Default is 0.25.
#' @param ms2.match.weight Numeric. Weight for MS2 matching score in the overall annotation score. Default is 0.5.
#' @param total.score.tol Numeric. Tolerance for the total matching score. Default is 0.5.
#' @param candidate.num Numeric. Maximum number of candidate annotations to retain per metabolite. Default is 3.
#' @param remove_fragment_intensity_cutoff Numeric. Cutoff to remove low-intensity MS2 fragments. Default is 0.
#' @param return_format Character. Specifies the format of the output. Can be `"mass_dataset"` or `"data.frame"`. Default is `"mass_dataset"`.
#' @param threads Numeric. Number of threads to use for parallel processing. Default is 3.
#'
#' @return A modified `mass_dataset` object with annotated metabolites added to the `annotation_table` slot.
#'
#' @details
#' This function performs metabolite annotation using a combination of MS1, retention time, and MS2 data (if available) from the provided object. The function allows users to customize the matching process, including setting tolerances for MS1 and MS2 matching, adjusting the weights of different scoring components, and selecting a specific chromatographic column and adduct table.
#'
#' If `ms2` is included in the `based_on` argument, the function extracts both MS1 and MS2 information for annotation. The final annotations are filtered based on the specified score thresholds and only the top `candidate.num` annotations are retained for each metabolite.
#'
#' @examples
#' \dontrun{
#' # Load a sample dataset and database
#' my_data <- load_mass_dataset("path/to/data")
#' my_database <- load_database("path/to/database")
#'
#' # Annotate metabolites using MS1 and MS2 data
#' annotated_data <- annotate_metabolites(
#'   object = my_data,
#'   database = my_database,
#'   based_on = c("ms1", "ms2"),
#'   polarity = "positive",
#'   column = "rp",
#'   ms1.match.ppm = 20,
#'   ms2.match.ppm = 25,
#'   candidate.num = 5,
#'   threads = 4
#' )
#' }
#'
#'
#' @export

annotate_metabolites <-
  function(object,
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
           threads = 3) {
    ###Check data
    if (missing(object)) {
      stop("No object is provided.\n")
    }
    # browser()
    if (missing(database)) {
      stop("No database is provided.\n")
    }
    
    if (!is(database, "databaseClass")) {
      stop("The database is not a databaseClass object.\n")
    }
    
    ## Match argument values for polarity and column type
    polarity <- match.arg(polarity)
    column <- match.arg(column)
    return_format <- match.arg(return_format)
    based_on <- match.arg(based_on, several.ok = TRUE)
    
    ###adduct table
    if (is.null(adduct.table)) {
      message("No adduct table is provided. Use the default adduct table.\n")
      #------------------------------------------------------------------
      ## Load the appropriate adduct table based on polarity and column type
      adduct.table <-
        load_adduct_table(polarity = polarity, column = column)
      
    }
    
    check_adduct_table(adduct.table)
    
    ###check parameters
    message("Checking parameters...\n")
    check_parameters4annotate_metabolites(
      object = object,
      database = database,
      based_on = based_on,
      polarity = polarity,
      column = column,
      adduct.table = adduct.table,
      ce = ce,
      ms1.match.ppm = ms1.match.ppm,
      ms2.match.ppm = ms2.match.ppm,
      mz.ppm.thr = mz.ppm.thr,
      ms2.match.tol = ms2.match.tol,
      fraction.weight = fraction.weight,
      dp.forward.weight = dp.forward.weight,
      dp.reverse.weight = dp.reverse.weight,
      rt.match.tol = rt.match.tol,
      ms1.match.weight = ms1.match.weight,
      rt.match.weight = rt.match.weight,
      ms2.match.weight = ms2.match.weight,
      total.score.tol = total.score.tol,
      candidate.num = candidate.num,
      remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff,
      threads = threads
    )
    
    message("All parameters are correct.\n")
    
    if ("ms2" %in% based_on) {
      ms1.info <-
        extract_ms1_info(object = object)
      
      ms2.info <-
        extract_ms2_info(object = object)
      
      if (!is.null(ms2.info)) {
        if (length(ms2.info) > 0) {
          ms2.info <-
            ms2.info[which(names(ms2.info) %in% ms1.info$ms2_spectrum_id)]
        }
      }
      
    } else{
      ms1.info <-
        extract_variable_info(object = object)
      ms2.info <- NULL
    }
    
    ##### Annotate metabolites for each MS2 data set in the object
    annotation_result <-
      annotate_peaks_mz_rt_ms2(
        ms1.info = ms1.info,
        ms2.info = ms2.info,
        database = database,
        based_on = based_on,
        polarity = polarity,
        column = column,
        adduct.table = adduct.table,
        ce = ce,
        ms1.match.ppm = ms1.match.ppm,
        ms2.match.ppm = ms2.match.ppm,
        mz.ppm.thr = mz.ppm.thr,
        ms2.match.tol = ms2.match.tol,
        fraction.weight = fraction.weight,
        dp.forward.weight = dp.forward.weight,
        dp.reverse.weight = dp.reverse.weight,
        rt.match.tol = rt.match.tol,
        ms1.match.weight = ms1.match.weight,
        rt.match.weight = rt.match.weight,
        ms2.match.weight = ms2.match.weight,
        total.score.tol = total.score.tol,
        candidate.num = candidate.num,
        remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff,
        threads = threads
      )
    
    ##change the ms2_spectrum_id
    if (!is.null(annotation_result)) {
      if (nrow(annotation_result) > 0) {
        annotation_result$ms2_spectrum_id <-
          stringr::str_replace(
            annotation_result$ms2_spectrum_id,
            paste0(annotation_result$ms2_files_id, "_"),
            ""
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
        database = extract_database_name(database),
        based_on = based_on,
        polarity = polarity,
        column = column,
        adduct.table = adduct.table,
        ce = ce,
        ms1.match.ppm = ms1.match.ppm,
        ms2.match.ppm = ms2.match.ppm,
        mz.ppm.thr = mz.ppm.thr,
        ms2.match.tol = ms2.match.tol,
        fraction.weight = fraction.weight,
        dp.forward.weight = dp.forward.weight,
        dp.reverse.weight = dp.reverse.weight,
        rt.match.tol = rt.match.tol,
        ms1.match.weight = ms1.match.weight,
        rt.match.weight = rt.match.weight,
        ms2.match.weight = ms2.match.weight,
        total.score.tol = total.score.tol,
        candidate.num = candidate.num,
        remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff,
        threads = threads
      ),
      time = Sys.time()
    )
    
    
    if (all(names(process_info) != "annotate_metabolites")) {
      process_info$annotate_metabolites <- parameter
    } else{
      process_info$annotate_metabolites <-
        c(process_info$annotate_metabolites, parameter)
    }
    
    object@process_info <- process_info
    
    annotation_result <-
      annotation_result %>%
      dplyr::filter(variable_id %in% object@variable_info$variable_id)
    
    if (is.null(annotation_result)) {
      return(object)
    }
    
    if (nrow(annotation_result) == 0) {
      return(object)
    }
    
    annotation_result <-
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
        rbind(object@annotation_table, annotation_result) %>%
        dplyr::arrange(variable_id, Level, dplyr::desc(Total.score))
      
      ###only remain top annotations
      object@annotation_table <-
        object@annotation_table %>%
        dplyr::group_by(variable_id) %>%
        dplyr::slice_head(n = candidate.num) %>%
        dplyr::ungroup() %>%
        dplyr::distinct(.keep_all = TRUE)
    }
    
    message("All done.")
    if (return_format == "mass_dataset") {
      return(object)
    } else{
      return(annotation_result)
    }
    return(object)
  }
