#' Annotate Peaks Based on MS1, Retention Time (RT), and MS2
#'
#' This function annotates metabolites by matching MS1, retention time (RT), and MS2 spectra against a provided database. It allows for the use of custom parameters such as mass tolerance, MS2 fragment matching criteria, and retention time tolerance.
#'
#' @param ms1.info A data frame containing MS1 peak information (m/z and RT). If `based_on` includes `"ms1"` or `"rt"`, this argument is required.
#' @param ms2.info A list containing MS2 spectra for each corresponding `ms2_spectrum_id`. If `based_on` includes `"ms2"`, this argument is required.
#' @param database A `databaseClass` object containing the reference database for metabolite annotation.
#' @param based_on Character vector. Specifies which criteria to base the matching on. Can include `"ms1"`, `"rt"`, and `"ms2"`. Default is `c("ms1", "rt", "ms2")`.
#' @param polarity Character. The ionization mode, either `"positive"` or `"negative"`. Default is `"positive"`.
#' @param ce Character. Collision energy used in MS2 spectra. Default is `"all"`.
#' @param column Character. The chromatographic column type, either `"hilic"` or `"rp"` (reversed-phase). Default is `"hilic"`.
#' @param adduct.table A data frame containing the adducts to use in the matching process. If `NULL`, a default table is loaded based on the `polarity` and `column`.
#' @param ms1.match.ppm Numeric. The mass tolerance in parts per million (ppm) for MS1 peak matching. Default is 25.
#' @param mz.ppm.thr Numeric. m/z threshold for ppm calculation. Default is 400.
#' @param rt.match.tol Numeric. Retention time matching tolerance in seconds. Default is 30.
#' @param ms2.match.ppm Numeric. The mass tolerance in ppm for MS2 peak matching. Default is 30.
#' @param ms2.match.tol Numeric. The retention time tolerance for MS2 fragment matching. Default is 0.5.
#' @param fraction.weight Numeric. Weight for the fraction of matched fragments in MS2 spectra. Default is 0.3.
#' @param dp.forward.weight Numeric. Weight for the forward dot product score in MS2 matching. Default is 0.6.
#' @param dp.reverse.weight Numeric. Weight for the reverse dot product score in MS2 matching. Default is 0.1.
#' @param remove_fragment_intensity_cutoff Numeric. Intensity cutoff for removing low-intensity MS2 fragments. Default is 0.
#' @param ms1.match.weight Numeric. Weight for MS1 matching score in the total score calculation. Default is 0.25.
#' @param rt.match.weight Numeric. Weight for RT matching score in the total score calculation. Default is 0.25.
#' @param ms2.match.weight Numeric. Weight for MS2 matching score in the total score calculation. Default is 0.5.
#' @param total.score.tol Numeric. Threshold for the total score. Only results with a score above this value are retained. Default is 0.5.
#' @param candidate.num Numeric. Maximum number of top candidate annotations to retain per metabolite. Default is 3.
#' @param threads Numeric. Number of threads to use for parallel processing. Default is 3.
#'
#' @return A data frame with annotated metabolites, including columns for matched m/z, retention time, MS2 spectra, and the calculated scores for each match.
#'
#' @details
#' The function uses a combination of MS1 peak information (m/z and retention time), MS2 spectra, and a reference database to annotate metabolites. The matching process can be customized by adjusting the mass tolerance, retention time tolerance, and MS2 fragment matching parameters.
#'
#' If `based_on` includes `"ms1"` or `"rt"`, the MS1 information is extracted from the `ms1.info` data frame. If `based_on` includes `"ms2"`, the function uses the provided MS2 spectra in `ms2.info` to perform fragment matching. The function calculates individual scores for m/z, retention time, and MS2 fragment matches, which are then combined into a total score. Annotations with total scores above `total.score.tol` are retained, and only the top `candidate.num` annotations are kept for each metabolite.
#'
#' @examples
#' \dontrun{
#' # Example MS1 and MS2 data
#' ms1_info <- data.frame(
#'   variable_id = c("id1", "id2"),
#'   mz = c(150.08, 180.12),
#'   rt = c(12.5, 14.7)
#' )
#' ms2_info <- list(
#'   id1 = matrix(c(75, 1000, 80, 2000),
#'   ncol = 2, byrow = TRUE,
#'   dimnames = list(NULL, c("mz", "intensity"))),
#'   id2 = matrix(c(85, 3000, 90, 1500),
#'   ncol = 2, byrow = TRUE,
#'   dimnames = list(NULL, c("mz", "intensity")))
#' )
#'
#' # Example database
#' database <- load_database("path/to/database")
#'
#' # Annotate metabolites using MS1, RT, and MS2 data
#' annotations <- annotate_peaks_mz_rt_ms2(
#'   ms1.info = ms1_info,
#'   ms2.info = ms2_info,
#'   database = database,
#'   based_on = c("ms1", "rt", "ms2"),
#'   polarity = "positive",
#'   column = "rp",
#'   ms1.match.ppm = 20,
#'   rt.match.tol = 20,
#'   candidate.num = 5,
#'   threads = 4
#' )
#'
#' print(annotations)
#' }
#'
#'
#' @export


annotate_peaks_mz_rt_ms2 <-
  function(ms1.info = NULL,
           ms2.info = NULL,
           database = NULL,
           based_on = c("ms1", "rt", "ms2"),
           polarity = c("positive", "negative"),
           ce = "all",
           column = c("hilic", "rp"),
           adduct.table = NULL,
           ms1.match.ppm = 25,
           mz.ppm.thr = 400,
           rt.match.tol = 30,
           ms2.match.ppm = 30,
           ms2.match.tol = 0.5,
           fraction.weight = 0.3,
           dp.forward.weight = 0.6,
           dp.reverse.weight = 0.1,
           remove_fragment_intensity_cutoff = 0,
           ms1.match.weight = 0.25,
           rt.match.weight = 0.25,
           ms2.match.weight = 0.5,
           total.score.tol = 0.5,
           candidate.num = 3,
           threads = 3) {
    options(warn = -1)
    
    based_on <- match.arg(based_on, several.ok = TRUE)
    
    if (is.null(database)) {
      stop("No database is provided.\n")
    }
    
    ##check ms1.file and ms2.file
    if (!is(database, "databaseClass")) {
      stop("database should be databaseClass object.\n")
    }
    
    ###if based_on contains ms2, the ms1.info should contains the ms2_spectrum_id
    ###and file
    
    ##parameter specification
    polarity <- match.arg(polarity)
    column <- match.arg(column)
    # ce <- match.arg(ce)
    
    if (is.null(adduct.table)) {
      adduct.table <-
        load_adduct_table(polarity = polarity, column = column)
    }
    
    check_adduct_table(adduct.table)
    
    database.name <-
      extract_database_name(database = database)
    
    ####------------------------------------------------------------------
    ####ms1 and/or rt match
    ####------------------------------------------------------------------
    if ("ms1" %in% based_on | "rt" %in% based_on) {
      ###Check data
      if (is.null(ms1.info)) {
        stop("No ms1.info is provided.\n")
      }
      
      ###check ms.info
      check_ms1_ms2_info(ms1.info = ms1.info,
                         ms2.info = ms2.info,
                         based_on = based_on)
      
      ms1_database <-
        extract_ms1_database(database = database)
      
      ###calculate the mz_rt_matrix for all the compound with all adducts
      mz_matrix <-
        seq_len(nrow(adduct.table)) %>%
        purrr::map(function(i) {
          temp_n <-
            stringr::str_extract(string = adduct.table$adduct[i], pattern = "[0-9]{1}M") %>%
            stringr::str_replace("M", "") %>%
            as.numeric()
          temp_n[is.na(temp_n)] <- 1
          temp <-
            data.frame(as.numeric(adduct.table$mz[i]) + temp_n * as.numeric(ms1_database$mz))
          colnames(temp) <- i
          temp
        }) %>%
        dplyr::bind_cols()
      
      colnames(mz_matrix) <-
        adduct.table$adduct
      
      rownames(mz_matrix) <-
        ms1_database$Lab.ID
      
      mz_matrix <-
        mz_matrix %>%
        tibble::rownames_to_column(var = "Lab.ID") %>%
        tidyr::pivot_longer(cols = -Lab.ID,
                            names_to = "Adduct",
                            values_to = "mz") %>%
        as.data.frame()
      
      mz_rt_matrix <-
        mz_matrix %>%
        dplyr::left_join(ms1_database[, c("Lab.ID", "RT")], by = "Lab.ID")
      
      ###mz match
      # pb <- progress::progress_bar$new(total = nrow(ms1.info))
      future::plan(strategy = future::multisession, workers = threads)
      
      match_result <-
        seq_len(nrow(ms1.info)) %>%
        furrr::future_map(function(i) {
          # pb$tick()
          # cat(i, " ")
          temp <-
            mz_rt_matrix %>%
            dplyr::mutate(variable_id = ms1.info$variable_id[i]) %>%
            dplyr::mutate(
              mz.error = abs(ms1.info$mz[i] - mz) * 10 ^ 6 / ifelse(ms1.info$mz[i] < 400, 400, ms1.info$mz[i])
            ) %>%
            dplyr::mutate(RT.error = abs(ms1.info$rt[i] - RT)) %>%
            dplyr::filter(mz.error < ifelse("ms1" %in% based_on, ms1.match.ppm, max(mz.error) + 100)) %>%
            dplyr::select(variable_id, Lab.ID, Adduct, mz.error, RT.error) %>%
            dplyr::arrange(mz.error, RT.error)
          
          if ("rt" %in% based_on) {
            temp <-
              temp %>%
              dplyr::filter(RT.error < rt.match.tol)
          }
          temp
        }, .progress = TRUE)
      
      match_result <-
        match_result %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(Database = database.name) %>%
        dplyr::select(variable_id, dplyr::everything())
      
      if (nrow(match_result) == 0) {
        message("No matched compounds.")
        return(NULL)
      }
      
      match_result$mz.match.score <-
        calculate_mz_match_score(mz.error = match_result$mz.error, ms1.match.ppm = ms1.match.ppm)
      
      match_result$RT.match.score <-
        calculate_rt_match_score(RT.error = match_result$RT.error, rt.match.tol = rt.match.tol)
      
      match_result <-
        match_result %>%
        dplyr::left_join(ms1_database, by = "Lab.ID")
      
      match_result <-
        match_result %>%
        dplyr::select(
          "variable_id",
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
          "Database",
          dplyr::everything()
        )
      
      ####remove some impossible adducts
      match_result <-
        remove_impossible_annotations(match_result)
      
      match_result <-
        match_result %>%
        dplyr::arrange(variable_id,
                       dplyr::desc(mz.match.score),
                       dplyr::desc(RT.match.score))
      
      if (nrow(match_result) == 0) {
        message("No matched compounds.")
        return(match_result)
      }
    } else{
      match_result <- NULL
    }
    
    ####------------------------------------------------------------------
    ####ms2 and/or ms1/rt match
    ####------------------------------------------------------------------
    if ("ms2" %in% based_on) {
      if (missing(ms2.info)) {
        stop("No ms2.info is provided.\n")
      }
      
      check_ms1_ms2_info(ms1.info = ms1.info,
                         ms2.info = ms2.info,
                         based_on = based_on)
      
      ##extract MS2 data
      spectra.data <-
        extract_ms2_database(database = database,
                             polarity = polarity,
                             ce = ce)
      
      if (is.null(spectra.data)) {
        stop("No MS2 spectra in the database.\n")
      }
      
      if (masstools::get_os() == "windows") {
        bpparam = BiocParallel::SnowParam(workers = threads, progressbar = TRUE)
      } else{
        bpparam = BiocParallel::MulticoreParam(workers = threads, progressbar = TRUE)
      }
      
      ###remove some metabolites from spectra.data according to ms1 and/or rt matching
      if (!is.null(match_result)) {
        match_result <-
          match_result %>%
          dplyr::left_join(ms1.info[, c("variable_id", "ms2_files_id", "ms2_spectrum_id")], by = "variable_id")
      }
      
      ms2.info <-
        ms2.info[which(names(ms2.info) %in% match_result$ms2_spectrum_id)]
      
      if (length(ms2.info) == 0) {
        return(NULL)
      }
      
      ###debug
      # for (i in seq_len(length(ms2.info))) {
      #   cat(i, " ")
      #   match_ms2_temp(
      #     idx = i,
      #     ms2.info = ms2.info,
      #     pre_match_result = match_result,
      #     spectra.data = spectra.data,
      #     ms2.match.ppm = ms2.match.ppm,
      #     mz.ppm.thr = mz.ppm.thr,
      #     ms2.match.tol = ms2.match.tol,
      #     candidate.num = candidate.num,
      #     fraction.weight = fraction.weight,
      #     dp.forward.weight = dp.forward.weight,
      #     dp.reverse.weight = dp.reverse.weight,
      #     remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff
      #   )
      # }
      
      match_result_ms2 <-
        suppressMessages(
          BiocParallel::bplapply(
            seq_len(length(ms2.info)),
            FUN = match_ms2_temp,
            BPPARAM = bpparam,
            pre_match_result = match_result,
            ms2.info = ms2.info,
            spectra.data = spectra.data,
            ms2.match.ppm = ms2.match.ppm,
            mz.ppm.thr = mz.ppm.thr,
            ms2.match.tol = ms2.match.tol,
            candidate.num = candidate.num,
            fraction.weight = fraction.weight,
            dp.forward.weight = dp.forward.weight,
            dp.reverse.weight = dp.reverse.weight,
            remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff
          )
        )
      
      names(match_result_ms2) <- names(ms2.info)
      
      match_result_ms2 <-
        match_result_ms2[which(!unlist(lapply(match_result_ms2, function(x)
          all(is.null(x)))))]
      
      match_result_ms2 <-
        match_result_ms2 %>%
        do.call(rbind, .) %>%
        as.data.frame()
      
      if (nrow(match_result_ms2) == 0) {
        return(NULL)
      }
      
      rownames(match_result_ms2) <- NULL
      
      if ("ms1" %in% based_on | "rt" %in% based_on) {
        match_result <-
          match_result_ms2 %>%
          dplyr::left_join(match_result, by = c("ms2_spectrum_id", "Lab.ID")) %>%
          dplyr::select(
            "variable_id",
            "ms2_spectrum_id",
            "ms2_files_id",
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
            "SS",
            "Database",
            dplyr::everything()
          )
      } else{
        match_result_ms2$Database <- database.name
        match_result <-
          match_result_ms2 %>%
          dplyr::left_join(ms1.info, by = "ms2_spectrum_id") %>%
          dplyr::left_join(ms1_database, by = c("Lab.ID")) %>%
          dplyr::mutate(
            Adduct = NA,
            mz.error = NA,
            mz.match.score = NA,
            RT.error = NA,
            RT.match.score = NA,
            SS = NA
          ) %>%
          dplyr::select(
            "variable_id",
            "ms2_spectrum_id",
            "ms2_files_id",
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
            "SS",
            "Database",
            dplyr::everything()
          )
      }
    } else{
      match_result <-
        match_result %>%
        dplyr::mutate(
          ms2_spectrum_id = NA,
          ms2_files_id = NA,
          CE = NA,
          SS = NA
        )
    }
    
    ###calculate the total score
    Total.score <-
      calculate_total_score(
        mz.match.score = match_result$mz.match.score,
        RT.match.score = match_result$RT.match.score,
        SS = match_result$SS,
        ms1.match.weight = ms1.match.weight,
        rt.match.weight = rt.match.weight,
        ms2.match.weight = ms2.match.weight,
        based_on = based_on
      )
    
    match_result$Total.score <- Total.score
    
    match_result <-
      match_result %>%
      dplyr::filter(Total.score > total.score.tol)
    
    if (nrow(match_result) == 0) {
      message("No matched compounds.")
      return(NULL)
    }
    
    Level <-
      calculate_confidence_level(annotation_result = match_result)
    
    match_result$Level <- Level
    
    match_result <-
      match_result %>%
      dplyr::arrange(
        variable_id,
        Level,
        dplyr::desc(Total.score),
        dplyr::desc(SS),
        dplyr::desc(mz.match.score),
        dplyr::desc(RT.match.score)
      ) %>%
      dplyr::group_by(variable_id) %>%
      dplyr::slice_head(n = candidate.num) %>%
      dplyr::ungroup() %>%
      as.data.frame()
    
    return(match_result)
    
  }


#' Match MS2 Spectra to Database Compounds
#'
#' This function matches MS2 spectra from experimental data against MS2 spectra in a database. It calculates MS2 matching scores for each compound in the database, and selects the best matches based on the provided matching parameters.
#'
#' @param idx Numeric. The index of the current MS2 spectrum in the `ms2.info` list.
#' @param ms2.info A list of MS2 spectra, where each element is a matrix with two columns (`mz` and `intensity`).
#' @param pre_match_result A data frame containing the preliminary match results from MS1 or RT matching, including compound identifiers (`Lab.ID`).
#' @param spectra.data A list of MS2 spectra from the database, where each element is a matrix with two columns (`mz` and `intensity`).
#' @param ms2.match.ppm Numeric. Mass tolerance in parts per million (ppm) for MS2 peak matching. Default is 30.
#' @param mz.ppm.thr Numeric. m/z threshold for ppm calculation. Default is 400.
#' @param ms2.match.tol Numeric. Score threshold for MS2 matches. Only matches with scores above this value are retained. Default is 0.5.
#' @param candidate.num Numeric. The number of top candidate annotations to retain per MS2 spectrum. Default is 3.
#' @param fraction.weight Numeric. Weight for the fraction of matched fragments in the total MS2 score calculation. Default is 0.3.
#' @param dp.forward.weight Numeric. Weight for the forward dot product score in the total MS2 score calculation. Default is 0.6.
#' @param dp.reverse.weight Numeric. Weight for the reverse dot product score in the total MS2 score calculation. Default is 0.1.
#' @param remove_fragment_intensity_cutoff Numeric. Intensity cutoff to remove low-intensity MS2 fragments from matching. Default is 0.
#' @param ... Additional arguments to pass to internal functions.
#'
#' @return A data frame containing the top candidate matches for the given MS2 spectrum. The output includes:
#' \describe{
#'   \item{Lab.ID}{The database compound ID.}
#'   \item{CE}{The collision energy used in the MS2 spectrum.}
#'   \item{SS}{The matching score for the compound.}
#'   \item{ms2_spectrum_id}{The ID of the experimental MS2 spectrum.}
#' }
#'
#' @details
#' The function takes an MS2 spectrum from the experimental data and matches it against MS2 spectra in a reference database. It computes matching scores based on the provided mass tolerance (`ms2.match.ppm`), fragment matching weight parameters (`fraction.weight`, `dp.forward.weight`, `dp.reverse.weight`), and intensity cutoff. The best-matching spectra from the database are retained and returned as candidates.
#'
#' If preliminary MS1/RT matches are provided (`pre_match_result`), the function only considers those compounds for MS2 matching. Otherwise, all compounds in the database are considered.
#'
#' @export


match_ms2_temp <-
  function(idx,
           ms2.info,
           pre_match_result,
           spectra.data,
           ms2.match.ppm = 30,
           mz.ppm.thr = 400,
           ms2.match.tol = 0.5,
           candidate.num = 3,
           fraction.weight = 0.3,
           dp.forward.weight = 0.6,
           dp.reverse.weight = 0.1,
           remove_fragment_intensity_cutoff = 0,
           ...) {
    peak_ms2_spectrum <-
      as.data.frame(ms2.info[[idx]])
    ms2_spectrum_id <- names(ms2.info)[idx]
    rm(list = c("ms2.info"))
    
    if (length(peak_ms2_spectrum) == 0) {
      return(NA)
    }
    
    if (!is.null(pre_match_result)) {
      library_compound_id <-
        pre_match_result$Lab.ID[which(pre_match_result$ms2_spectrum_id == ms2_spectrum_id)]
    } else{
      library_compound_id <- names(spectra.data)
    }
    
    library_compound_id <-
      library_compound_id[which(library_compound_id %in% names(spectra.data))]
    
    if (length(library_compound_id) == 0) {
      return(NULL)
    } else{
      ###MS2 spectra match
      ms2_match_score <-
        seq_len(length(spectra.data[library_compound_id])) %>%
        purrr::map(function(i) {
          temp_spectra_data <- spectra.data[library_compound_id][[i]]
          score <-
            lapply(temp_spectra_data, function(y) {
              y <- as.data.frame(y)
              y$mz <- as.numeric(y$mz)
              y$intensity <- as.numeric(y$intensity)
              calculate_ms2_matching_score(
                experimental.spectrum = peak_ms2_spectrum,
                library.spectrum = y,
                ms2.match.ppm = ms2.match.ppm,
                mz.ppm.thr = mz.ppm.thr,
                fraction.weight = fraction.weight,
                dp.forward.weight = dp.forward.weight,
                dp.reverse.weight = dp.reverse.weight,
                remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff
              )
            })
          score <- score[which.max(unlist(score))]
          score <- unlist(score)
          data.frame(
            Lab.ID = names(spectra.data[library_compound_id])[i],
            "CE" = names(score),
            "SS" = score,
            stringsAsFactors = FALSE
          )
        }) %>%
        dplyr::bind_rows()
      
      rownames(ms2_match_score) <- NULL
      
      ms2_match_score <-
        ms2_match_score %>%
        dplyr::mutate(ms2_spectrum_id = ms2_spectrum_id) %>%
        dplyr::arrange(dplyr::desc(SS)) %>%
        head(candidate.num) %>%
        dplyr::filter(SS > ms2.match.tol)
    }
    
    # rm(list = c("ms2_match_score"))
    if (nrow(ms2_match_score) == 0) {
      return(NULL)
    } else{
      return(ms2_match_score)
    }
  }
