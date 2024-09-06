#' @title Identify peaks based on MS1 database
#' @description Identify peaks based on MS1 database.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param object A mass_dataset class object.
#' @param ms1.match.ppm Precursor match ppm tolerance.
#' @param rt.match.tol RT match tolerance.
#' @param polarity The polarity of data, "positive"or "negative".
#' @param column "hilic" (HILIC column) or "rp" (reverse phase).
#' @param candidate.num The number of candidates.
#' @param database MS1 database name or MS1 database.
#' @param threads Number of threads
#' @return A mzIdentifyClass or metIdentifyClass object.
#' @importFrom magrittr %>%
#' @importFrom dplyr pull filter
#' @export
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/metid.html}


mzIdentify_mass_dataset2 <-
  function(object,
           ms1.match.ppm = 25,
           rt.match.tol = 30,
           polarity = c("positive", "negative"),
           column = c("hilic", "rp"),
           candidate.num = 3,
           database,
           threads = 3) {
    options(warn = -1)
    ###Check data
    if (missing(database)) {
      stop("No database is provided.\n")
    }
    
    ##parameter specification
    polarity <- match.arg(polarity)
    column <- match.arg(column)
    ##check ms1.file and ms2.file
    if (!is(database, "databaseClass")) {
      stop("database should be databaseClass object.\n")
    }
    
    database.name <-
      paste(database@database.info$Source,
            database@database.info$Version,
            sep = "_")
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
    
    ms1.data <-
      object@variable_info %>%
      dplyr::rename(name = variable_id)
    
    if (rt.match.tol > 10000) {
      message(crayon::yellow("You set rt.match.tol as NA, so RT will not be used for matching."))
    } else{
      message(
        crayon::yellow(
          "You set rt.match.tol < 10,000, so if the compounds have RT,  RTs will be used for matching."
        )
      )
    }
    
    temp.fun <-
      function(idx,
               ms1.data,
               ms1.match.ppm = 25,
               rt.match.tol = 30,
               database,
               adduct.table,
               candidate.num = 3) {
        temp_mz <-
          as.numeric(ms1.data$mz[idx])
        temp_rt <-
          as.numeric(ms1.data$rt[idx])
        
        rm(list = c("ms1.data"))
        
        temp_mz_diff1 <- abs(temp_mz - as.numeric(database$mz))
        temp_mz_diff2 <- abs(temp_mz - as.numeric(database$mz) * 2)
        temp_mz_diff3 <- abs(temp_mz - as.numeric(database$mz) * 3)
        
        temp_mz_diff1[is.na(temp_mz_diff1)] <- 100000
        temp_mz_diff2[is.na(temp_mz_diff2)] <- 100000
        temp_mz_diff3[is.na(temp_mz_diff3)] <- 100000
        
        max_mz_diff <-
          max(abs(adduct.table$mz)) + 1
        
        database <-
          database[which(
            temp_mz_diff1 < max_mz_diff |
              temp_mz_diff2 < max_mz_diff |
              temp_mz_diff3 < max_mz_diff
          )
          , , drop = FALSE]
        
        rm(list = c("temp_mz_diff1", "temp_mz_diff2", "temp_mz_diff3"))
        
        if (nrow(database) == 0) {
          return(NA)
        }
        
        # spectra_mz <-
        #   purrr::map(as.data.frame(t(adduct.table)),
        #              function(x) {
        #                temp_n <-
        #                  stringr::str_extract(string = as.character(x[1]),
        #                                       pattern = "[0-9]{1}M")
        #                temp_n <-
        #                  as.numeric(stringr::str_replace(
        #                    string = temp_n,
        #                    pattern = "M",
        #                    replacement = ""
        #                  ))
        #                temp_n[is.na(temp_n)] <- 1
        #                as.numeric(x[2]) + temp_n * as.numeric(database$mz)
        #              }) %>%
        #   do.call(cbind, .)
        
        spectra_mz <-
          seq_len(nrow(adduct.table)) %>%
          purrr::map(function(i) {
            temp_n <-
              stringr::str_extract(string = adduct.table$adduct[i],
                                   pattern = "[0-9]{1}M") %>%
              stringr::str_replace("M", "") %>%
              as.numeric()
            temp_n[is.na(temp_n)] <- 1
            as.numeric(adduct.table$mz[i]) + temp_n * as.numeric(database$mz)
          }) %>%
          do.call(cbind, .) %>%
          as.data.frame()
        
        colnames(spectra_mz) <- adduct.table$adduct
        rownames(spectra_mz) <- database$Lab.ID
        
        ###mz match
        temp <-
          abs(spectra_mz - temp_mz) * 10 ^ 6 / ifelse(temp_mz < 400, 400, temp_mz)
        
        temp_idx <-
          which(temp < ms1.match.ppm, arr.ind = TRUE) %>%
          as.data.frame()
        
        if (nrow(temp_idx) == 0) {
          return(NA)
        }
        
        match_idx <-
          seq_len(nrow(temp_idx)) %>%
          purrr::map(function(i) {
            data.frame(
              "Lab.ID" = rownames(spectra_mz)[temp_idx$row[i]],
              "Addcut" = colnames(spectra_mz)[temp_idx$col[i]],
              "mz.error" = temp[temp_idx$row[i], temp_idx$col[i]],
              stringsAsFactors = FALSE
            )
          })
        
        rm(list = c("spectra_mz", "adduct.table", "temp", "temp_idx"))
        
        ##remove some none matched
        match_idx <-
          match_idx[which(unlist(lapply(match_idx, function(x) {
            nrow(x)
          })) != 0)]
        
        if (length(match_idx) == 0) {
          return(NA)
        }
        
        match_idx <-
          match_idx %>%
          dplyr::bind_rows() %>%
          dplyr::arrange(mz.error)
        
        # match_idx <- data.frame(rownames(match_idx),
        # match_idx, stringsAsFactors = FALSE)
        colnames(match_idx) <-
          c("Lab.ID", "Adduct", "mz.error")
        
        ##rt match
        RT.error <-
          abs(temp_rt - as.numeric(database$RT)[match(match_idx$Lab.ID, database$Lab.ID)])
        
        match_idx <- data.frame(match_idx, RT.error,
                                stringsAsFactors = FALSE)
        
        match_idx <-
          dplyr::filter(match_idx,
                        is.na(RT.error) | RT.error < rt.match.tol)
        
        if (nrow(match_idx) == 0) {
          return(NA)
        }
        
        if (nrow(match_idx) > candidate.num) {
          match_idx <- match_idx[seq_len(candidate.num), , drop = FALSE]
        }
        
        match_idx <-
          data.frame(match_idx,
                     database[match(match_idx$Lab.ID, database$Lab.ID),
                              c("Compound.name", "CAS.ID", "HMDB.ID", "KEGG.ID"),
                              drop = FALSE],
                     stringsAsFactors = FALSE)
        
        match_idx <-
          match_idx[, c(
            "Compound.name",
            "CAS.ID",
            "HMDB.ID",
            "KEGG.ID",
            "Lab.ID",
            "Adduct",
            "mz.error",
            'RT.error'
          )]
        
        rownames(match_idx) <- NULL
        return(match_idx)
      }
    
    if (masstools::get_os() == "windows") {
      bpparam = BiocParallel::SnowParam(workers = threads,
                                        progressbar = TRUE)
    } else{
      bpparam = BiocParallel::MulticoreParam(workers = threads,
                                             progressbar = TRUE)
    }
    
    if (is(database, "databaseClass")) {
      ms1_database <-
        database@spectra.info %>%
        dplyr::mutate(mz = as.numeric(mz)) %>%
        dplyr::filter(!is.na(mz)) %>%
        dplyr::distinct(Compound.name, .keep_all = TRUE)
    } else{
      ms1_database <-
        database
    }
    
    match_result <-
      BiocParallel::bplapply(
        seq_len(nrow(ms1.data)),
        FUN = temp.fun,
        BPPARAM = bpparam,
        ms1.data = ms1.data,
        ms1.match.ppm = ms1.match.ppm,
        rt.match.tol = rt.match.tol,
        database = ms1_database,
        adduct.table = adduct.table,
        candidate.num = candidate.num
      )
    
    names(match_result) <- ms1.data$name
    
    temp_idx <-
      which(unlist(lapply(match_result, function(x) {
        all(is.na(x))
      })))
    
    if (length(temp_idx) > 0) {
      match_result <- match_result[-temp_idx]
    }
    
    match_result <-
      seq_along(match_result) %>%
      purrr::map(function(i) {
        data.frame(variable_id = names(match_result)[i],
                   match_result[[i]])
      }) %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(Database = database.name) %>%
      dplyr::mutate(ms2_files_id = NA,
                    ms2_spectrum_id = NA) %>%
      dplyr::select(variable_id,
                    ms2_files_id,
                    ms2_spectrum_id,
                    dplyr::everything())
    
    match_result$mz.match.score <-
      exp(-0.5 * (match_result$mz.error / (ms1.match.ppm)) ^ 2)
    
    if (rt.match.tol > 10000) {
      match_result$RT.error <- NA
      match_result$RT.match.score <- NA
      match_result$Total.score <- match_result$mz.match.score
    } else{
      match_result$RT.match.score <-
        exp(-0.5 * (match_result$RT.error / (rt.match.tol)) ^ 2)
      match_result$Total.score <-
        match_result$mz.match.score * 0.5 +
        match_result$RT.match.score * 0.5
      match_result$Total.score[is.na(match_result$Total.score)] <-
        match_result$mz.match.score[is.na(match_result$Total.score)]
    }
    
    match_result$CE <-
      match_result$SS <-
      NA
    
    match_result <-
      match_result %>%
      dplyr::select(
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
        dplyr::everything()
      )
    
    message(crayon::bgRed("All done."))
    return(match_result)
  }









#' @title Identify peaks based on MS1 database
#' @description Identify peaks based on MS1 database.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param object A mass_dataset class object.
#' @param ms1.match.ppm Precursor match ppm tolerance.
#' @param rt.match.tol RT match tolerance.
#' @param polarity The polarity of data, "positive"or "negative".
#' @param column "hilic" (HILIC column) or "rp" (reverse phase).
#' @param candidate.num The number of candidates.
#' @param database MS1 database name or MS1 database.
#' @param threads Number of threads
#' @return A mzIdentifyClass or metIdentifyClass object.
#' @importFrom magrittr %>%
#' @importFrom dplyr pull filter
#' @export
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/metid.html}


mzIdentify_mass_dataset <-
  function(object,
           ms1.match.ppm = 25,
           rt.match.tol = 30,
           polarity = c("positive", "negative"),
           column = c("hilic", "rp"),
           candidate.num = 3,
           database,
           threads = 3) {
    options(warn = -1)
    ###Check data
    if (missing(database)) {
      stop("No database is provided.\n")
    }
    
    ##parameter specification
    polarity <- match.arg(polarity)
    column <- match.arg(column)
    ##check ms1.file and ms2.file
    if (!is(database, "databaseClass")) {
      stop("database should be databaseClass object.\n")
    }
    
    database.name <-
      paste(database@database.info$Source,
            database@database.info$Version,
            sep = "_")
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
    
    ms1.data <-
      object@variable_info %>%
      dplyr::rename(name = variable_id)
    
    if (rt.match.tol > 10000) {
      message(crayon::yellow("You set rt.match.tol as NA, so RT will not be used for matching."))
    } else{
      message(
        crayon::yellow(
          "You set rt.match.tol < 10,000, so if the compounds have RT, RTs will be used for matching."
        )
      )
    }
    
    ms1_database <-
      database@spectra.info %>%
      dplyr::mutate(mz = as.numeric(mz)) %>%
      dplyr::filter(!is.na(mz)) %>%
      dplyr::distinct(Compound.name, .keep_all = TRUE)
    
    ###calculate the mz_matrix for all the compound with all adducts
    mz_matrix <-
      seq_len(nrow(adduct.table)) %>%
      purrr::map(function(i) {
        temp_n <-
          stringr::str_extract(string = adduct.table$adduct[i],
                               pattern = "[0-9]{1}M") %>%
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
    
    mz_rt_matrix <-
      mz_matrix %>%
      tibble::rownames_to_column(var = "Lab.ID") %>%
      tidyr::pivot_longer(cols = -c(Lab.ID),
                          names_to = "Adduct",
                          values_to = "mz") %>%
      dplyr::left_join(ms1_database[, c("Lab.ID", "RT")],
                       by = "Lab.ID")
    
    rm(list = "mz_matrix")
    gc()
    
    ###mz and rt match
    # pb <- progress::progress_bar$new(total = nrow(ms1.data))
    future::plan(strategy = future::multisession, workers = threads)
    
    match_result <-
      seq_len(nrow(ms1.data)) %>%
      furrr::future_map(function(i) {
        # pb$tick()
        # cat(i, " ")
        mz_rt_matrix %>%
          dplyr::mutate(variable_id = ms1.data$name[i]) %>%
          dplyr::mutate(mz.error = abs(ms1.data$mz[i] - mz) * 10 ^ 6 / ifelse(ms1.data$mz[i] < 400, 400, ms1.data$mz[i])) %>%
          dplyr::mutate(RT.error = abs(ms1.data$rt[i] - RT)) %>%
          dplyr::filter(mz.error < ms1.match.ppm) %>%
          dplyr::filter(is.na(RT.error) |
                          RT.error < rt.match.tol) %>%
          dplyr::select(variable_id, Lab.ID, Adduct, mz.error, RT.error) %>%
          dplyr::arrange(mz.error, RT.error) %>%
          head(candidate.num)
      }, .progress = TRUE)
    
    match_result <-
      match_result %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(Database = database.name) %>%
      dplyr::mutate(ms2_files_id = NA,
                    ms2_spectrum_id = NA) %>%
      dplyr::select(variable_id,
                    ms2_files_id,
                    ms2_spectrum_id,
                    dplyr::everything())
    
    match_result$mz.match.score <-
      exp(-0.5 * (match_result$mz.error / (ms1.match.ppm)) ^ 2)
    
    if (rt.match.tol > 10000) {
      match_result$RT.error <- NA
      match_result$RT.match.score <- NA
      match_result$Total.score <- match_result$mz.match.score
    } else{
      match_result$RT.match.score <-
        exp(-0.5 * (match_result$RT.error / (rt.match.tol)) ^ 2)
      match_result$Total.score <-
        match_result$mz.match.score * 0.5 +
        match_result$RT.match.score * 0.5
      match_result$Total.score[is.na(match_result$Total.score)] <-
        match_result$mz.match.score[is.na(match_result$Total.score)]
    }
    
    match_result$CE <-
      match_result$SS <-
      NA
    
    match_result <-
      match_result %>%
      dplyr::left_join(ms1_database[, c("Lab.ID",
                                        "Compound.name",
                                        "CAS.ID",
                                        "HMDB.ID",
                                        "KEGG.ID",
                                        "Formula")],
                       by = "Lab.ID")
    
    match_result <-
      match_result %>%
      dplyr::select(
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
        dplyr::everything()
      )
    
    ####remove some impossible adducts
    adduct_check <-
      match_result %>%
      dplyr::select(Formula, Adduct) %>%
      dplyr::distinct(Formula, Adduct) %>%
      dplyr::filter(stringr::str_detect(Adduct, "(-H2O)|(-2H2O)")) %>%
      dplyr::mutate(minus_h2o_number = stringr::str_extract(Adduct, "(-H2O)|(-2H2O)")) %>%
      dplyr::mutate(minus_h2o_number = stringr::str_replace(minus_h2o_number, "(H2O)", "")) %>%
      dplyr::mutate(minus_h2o_number = stringr::str_replace(minus_h2o_number, "-", "")) %>%
      dplyr::mutate(minus_h2o_number = as.numeric(minus_h2o_number))
    
    adduct_check$minus_h2o_number[is.na(adduct_check$minus_h2o_number)] <-
      1
    
    adduct_check <-
      adduct_check %>%
      dplyr::mutate(
        h_number = stringr::str_extract(Formula, "H[0-9]{1,2}") %>%
          stringr::str_replace("H", "") %>%
          as.numeric()
      ) %>%
      dplyr::mutate(
        o_number = stringr::str_extract(Formula, "O[0-9]{1,2}") %>%
          stringr::str_replace("O", "") %>%
          as.numeric()
      )
    
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
    message()
    message(crayon::bgRed("All done."))
    return(match_result)
  }
