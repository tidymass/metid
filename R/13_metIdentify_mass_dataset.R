#' Metabolite Identification in a mass_dataset Object Using MS1 and MS2 Data
#'
#' This function identifies potential metabolites in a `mass_dataset` object by matching MS1 and MS2 data with a reference spectral database. The function uses both MS1 (m/z) and MS2 (fragment ions) matching for more accurate identification.
#'
#' @param object A `mass_dataset` object that contains MS1 and MS2 data.
#' @param ms1.match.ppm A numeric value specifying the mass accuracy threshold for MS1 matching in parts per million (ppm). Defaults to `25`.
#' @param ms2.match.ppm A numeric value specifying the mass accuracy threshold for MS2 matching in ppm. Defaults to `30`.
#' @param mz.ppm.thr A numeric value specifying the m/z threshold in ppm for matching MS1 and MS2. Defaults to `400`.
#' @param ms2.match.tol A numeric value specifying the tolerance for MS2 fragment ion matching. Defaults to `0.5`.
#' @param fraction.weight A numeric value specifying the weight for the MS2 fragmentation score. Defaults to `0.3`.
#' @param dp.forward.weight A numeric value specifying the weight for the forward dot product in MS2 matching. Defaults to `0.6`.
#' @param dp.reverse.weight A numeric value specifying the weight for the reverse dot product in MS2 matching. Defaults to `0.1`.
#' @param rt.match.tol A numeric value specifying the retention time matching tolerance in seconds. Defaults to `30`.
#' @param polarity A character string specifying the ionization mode. It can be either `"positive"` or `"negative"`. Defaults to `"positive"`.
#' @param ce A character string specifying the collision energy for MS2 matching. Defaults to `"all"`.
#' @param column A character string specifying the chromatographic column type, either `"hilic"` (hydrophilic interaction) or `"rp"` (reverse phase). Defaults to `"hilic"`.
#' @param ms1.match.weight A numeric value specifying the weight of MS1 matching in the total score calculation. Defaults to `0.25`.
#' @param rt.match.weight A numeric value specifying the weight of RT matching in the total score calculation. Defaults to `0.25`.
#' @param ms2.match.weight A numeric value specifying the weight of MS2 matching in the total score calculation. Defaults to `0.5`.
#' @param total.score.tol A numeric value specifying the threshold for the total score. Defaults to `0.5`.
#' @param candidate.num A numeric value specifying the number of top candidates to retain per feature. Defaults to `3`.
#' @param database A `databaseClass` object containing the reference spectral database for annotation.
#' @param threads An integer specifying the number of threads to use for parallel processing. Defaults to `3`.
#' @param remove_fragment_intensity_cutoff A numeric value specifying the intensity cutoff for removing fragments in MS2 matching. Defaults to `0`.
#'
#' @return A data frame containing the metabolite identification results, including m/z error, RT error, MS2 matching scores, and information about the identified compounds.
#'
#' @details
#' This function performs MS1 and MS2-based matching between the experimental data in the `mass_dataset` object and a reference spectral database. The matching process is based on mass-to-charge ratio (m/z), retention time (RT), and MS2 fragmentation patterns. The function supports both positive and negative ionization modes and can work with either HILIC or reverse-phase columns.
#'
#' The matching process can be fine-tuned by adjusting the weights of MS1, MS2, and RT matching, as well as the tolerance parameters for m/z and MS2 matching.
#'
#' @examples
#' \dontrun{
#' # Perform MS1 and MS2-based metabolite identification in a mass_dataset object
#' identification_result <- metIdentify_mass_dataset(
#'   object = mass_object,
#'   ms1.match.ppm = 20,
#'   ms2.match.ppm = 25,
#'   rt.match.tol = 20,
#'   database = reference_database,
#'   threads = 4
#' )
#' }
#'
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @importFrom dplyr filter mutate select left_join bind_rows
#' @importFrom purrr map2
#' @importFrom crayon yellow green bgRed
#' @export

metIdentify_mass_dataset <-
  function(object,
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
           column = c("hilic", "rp"),
           ms1.match.weight = 0.25,
           rt.match.weight = 0.25,
           ms2.match.weight = 0.5,
           total.score.tol = 0.5,
           candidate.num = 3,
           database,
           threads = 3,
           remove_fragment_intensity_cutoff = 0) {
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
    
    #load MS2 database
    database.name = paste(database@database.info$Source,
                          database@database.info$Version,
                          sep = "_")
    
    
    if (!is(database, "databaseClass")) {
      stop("database must be databaseClass object\n")
    }
    
    ce.list.pos <-
      unique(unlist(lapply(
        database@spectra.data$Spectra.positive, names
      )))
    
    ce.list.neg <-
      unique(unlist(lapply(
        database@spectra.data$Spectra.negative, names
      )))
    
    ce.list <-
      ifelse(polarity == "positive", ce.list.pos, ce.list.neg)
    
    if (all(ce %in% ce.list) & ce != "all") {
      stop("All ce values you set are not in database. Please check it.\n")
      ce <- ce[ce %in% ce.list]
    }
    
    rm(list = c("ce.list.pos", "ce.list.neg", "ce.list"))
    
    ##ce values
    if (all(ce != "all")) {
      if (polarity == "positive") {
        ce.list <-
          unique(unlist(
            lapply(database@spectra.data$Spectra.positive, function(x) {
              names(x)
            })
          ))
        if (length(grep("Unknown", ce.list)) > 0) {
          ce <-
            unique(c(ce, grep(
              pattern = "Unknown", ce.list, value = TRUE
            )))
        }
      } else{
        ce.list <-
          unique(unlist(
            lapply(database@spectra.data$Spectra.negative, function(x) {
              names(x)
            })
          ))
        if (length(grep("Unknown", ce.list)) > 0) {
          ce <-
            unique(c(ce, grep(
              pattern = "Unknown", ce.list, value = TRUE
            )))
        }
      }
    }
    
    ##RT in database or not
    if (!database@database.info$RT) {
      message(crayon::yellow(
        "No RT information in database.\nThe weight of RT have been set as 0."
      ))
    }
    
    #------------------------------------------------------------------
    ##load adduct table
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
    
    if (length(object@ms2_data) == 0) {
      stop("No MS2 in you object.\n")
    }
    
    if (lapply(object@ms2_data, function(x) {
      length(x@ms2_spectra)
    }) %>%
    unlist() %>%
    sum() %>%
    `==`(0)) {
      stop("No MS2 in you object.\n")
    }
    
    #####annotation result for each set MS2 data
    annotation_result <-
      purrr::map2(.x = names(object@ms2_data), .y = object@ms2_data, function(temp_ms2_data_id, temp_ms2_data) {
        message(crayon::yellow(temp_ms2_data_id, "file:"))
        message(crayon::green(length(temp_ms2_data@ms2_spectra), "MS2 spectra."))
        
        ms1.info = data.frame(
          name = temp_ms2_data@ms2_spectrum_id,
          mz = temp_ms2_data@ms2_mz,
          rt = temp_ms2_data@ms2_rt,
          file = temp_ms2_data@ms2_file,
          variable_id = temp_ms2_data@variable_id
        )
        
        ms2.info = temp_ms2_data@ms2_spectra
        
        ms2_matchresult <-
          metIdentification(
            ms1.info = ms1.info,
            ms2.info = ms2.info,
            polarity = polarity,
            ce = ce,
            database = database,
            ms1.match.ppm = ms1.match.ppm,
            ms2.match.ppm = ms2.match.ppm,
            mz.ppm.thr = mz.ppm.thr,
            ms2.match.tol = ms2.match.tol,
            rt.match.tol = rt.match.tol,
            column = column,
            ms1.match.weight = ms1.match.weight,
            rt.match.weight = rt.match.weight,
            ms2.match.weight = ms2.match.weight,
            total.score.tol = total.score.tol,
            candidate.num = candidate.num,
            adduct.table = adduct.table,
            threads = threads,
            fraction.weight = fraction.weight,
            dp.forward.weight = dp.forward.weight,
            dp.reverse.weight = dp.reverse.weight,
            remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff
          )
        
        ms2_matchresult =
          purrr::map2(
            .x = names(ms2_matchresult),
            .y = ms2_matchresult,
            .f = function(temp_ms2_id, temp_annotation_result) {
              data.frame(
                ms2_files_id = temp_ms2_data_id,
                ms2_spectrum_id = temp_ms2_id,
                temp_annotation_result
              ) %>%
                dplyr::left_join(ms1.info[, c("name", "variable_id")], by = c("ms2_spectrum_id" = "name")) %>%
                dplyr::select(variable_id, dplyr::everything())
              
            }
          ) %>%
          dplyr::bind_rows()
        ms2_matchresult
      })
    
    annotation_result =
      annotation_result %>%
      dplyr::bind_rows() %>%
      as.data.frame() %>%
      dplyr::mutate(Database = database.name)
    
    message(crayon::bgRed("All done."))
    return(annotation_result)
  }
