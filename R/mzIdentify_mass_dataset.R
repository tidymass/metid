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


mzIdentify_mass_dataset =
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
    if (class(database) != "databaseClass") {
      stop("database should be databaseClass object.\n")
    }
    
    database.name = paste(database@database.info$Source,
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
      message(
        crayon::yellow(
          "You set rt.match.tol as NA, so RT will not be used for matching."
        )
      )
    } else{
      message(
        crayon::yellow(
          "You set rt.match.tol < 10,000, so if the compounds have RT,  RTs will be used for matching."
        )
      )
    }
    
    temp.fun <- function(idx,
                         ms1.data,
                         ms1.match.ppm = 25,
                         rt.match.tol = 30,
                         database,
                         adduct.table,
                         candidate.num = 3) {
      temp.mz <-
        as.numeric(dplyr::pull(.data = ms1.data[idx,], var = "mz"))
      temp.rt <-
        as.numeric(dplyr::pull(.data = ms1.data[idx,], var = "rt"))
      
      rm(list = c("ms1.data"))
      
      if (class(database) == "databaseClass") {
        database <- database@spectra.info
      }
      
      temp.mz.diff1 <- abs(temp.mz - as.numeric(database$mz))
      temp.mz.diff2 <- abs(temp.mz - as.numeric(database$mz) * 2)
      temp.mz.diff3 <- abs(temp.mz - as.numeric(database$mz) * 3)
      
      temp.mz.diff1[is.na(temp.mz.diff1)] <- 100000
      temp.mz.diff2[is.na(temp.mz.diff2)] <- 100000
      temp.mz.diff3[is.na(temp.mz.diff3)] <- 100000
      
      database <-
        database[which(
          temp.mz.diff1 < max(abs(range(
            adduct.table$mz
          ))) + 1 |
            temp.mz.diff2 < max(abs(range(
              adduct.table$mz
            ))) + 1 |
            temp.mz.diff3 < max(abs(range(
              adduct.table$mz
            ))) + 1
        )
        , , drop = FALSE]
      
      rm(list = c("temp.mz.diff1", "temp.mz.diff2", "temp.mz.diff3"))
      
      if (nrow(database) == 0)
        return(NA)
      
      spectra.mz <- purrr::map(as.data.frame(t(adduct.table)),
                               function(x) {
                                 temp.n <-
                                   stringr::str_extract(string = as.character(x[1]), pattern = "[0-9]{1}M")
                                 temp.n <-
                                   as.numeric(stringr::str_replace(
                                     string = temp.n,
                                     pattern = "M",
                                     replacement = ""
                                   ))
                                 temp.n[is.na(temp.n)] <- 1
                                 as.numeric(x[2]) + temp.n * as.numeric(database$mz)
                               }) %>%
        do.call(cbind, .)
      
      # spectra.mz <- apply(adduct.table, 1, function(x) {
      #   temp.n <-
      #     stringr::str_extract(string = as.character(x[1]), pattern = "[0-9]{1}M")
      #   temp.n <-
      #     as.numeric(stringr::str_replace(
      #       string = temp.n,
      #       pattern = "M",
      #       replacement = ""
      #     ))
      #   temp.n[is.na(temp.n)] <- 1
      #   as.numeric(x[2]) + temp.n * as.numeric(database$mz)
      # })
      
      colnames(spectra.mz) <- adduct.table[, 1]
      rownames(spectra.mz) <- database$Lab.ID
      
      ###mz match
      temp <-
        abs(spectra.mz - temp.mz) * 10 ^ 6 / ifelse(temp.mz < 400, 400, temp.mz)
      temp.idx <-
        which(temp < ms1.match.ppm, arr.ind = TRUE)
      if (nrow(temp.idx) == 0)
        return(NA)
      
      match.idx <- apply(temp.idx, 1, function(x) {
        data.frame(
          "Lab.ID" = rownames(spectra.mz)[x[1]],
          "Addcut" = colnames(spectra.mz)[x[2]],
          "mz.error" = temp[x[1], x[2]],
          stringsAsFactors = FALSE
        )
      })
      
      rm(list = c("spectra.mz", "adduct.table", "temp", "temp.idx"))
      
      ##remove some none matched
      match.idx <-
        match.idx[which(unlist(lapply(match.idx, function(x) {
          nrow(x)
        })) != 0)]
      
      if (length(match.idx) == 0) {
        return(NA)
      }
      
      match.idx <- do.call(rbind, match.idx)
      # match.idx <- data.frame(rownames(match.idx), match.idx, stringsAsFactors = FALSE)
      colnames(match.idx) <-
        c("Lab.ID", "Adduct", "mz.error")
      
      match.idx <-
        match.idx[order(match.idx$mz.error, decreasing = FALSE), , drop = FALSE]
      
      ##rt match
      RT.error <-
        abs(temp.rt - as.numeric(database$RT)[match(match.idx$Lab.ID, database$Lab.ID)])
      
      match.idx <- data.frame(match.idx, RT.error,
                              stringsAsFactors = FALSE)
      
      match.idx <-
        dplyr::filter(match.idx,
                      is.na(RT.error) | RT.error < rt.match.tol)
      
      if (nrow(match.idx) == 0) {
        return(NA)
      }
      
      if (nrow(match.idx) > candidate.num) {
        match.idx <- match.idx[seq_len(candidate.num), , drop = FALSE]
      }
      
      match.idx <- data.frame(match.idx,
                              database[match(match.idx$Lab.ID, database$Lab.ID),
                                       c("Compound.name", "CAS.ID", "HMDB.ID", "KEGG.ID"), drop = FALSE],
                              stringsAsFactors = FALSE)
      
      match.idx <-
        match.idx[, c(
          "Compound.name",
          "CAS.ID",
          "HMDB.ID",
          "KEGG.ID",
          "Lab.ID",
          "Adduct",
          "mz.error",
          'RT.error'
        )]
      
      rownames(match.idx) <- NULL
      return(match.idx)
    }
    
    
    if (masstools::get_os() == "windows") {
      bpparam = BiocParallel::SnowParam(workers = threads,
                                        progressbar = TRUE)
    } else{
      bpparam = BiocParallel::MulticoreParam(workers = threads,
                                             progressbar = TRUE)
    }
    
    match.result <-
      BiocParallel::bplapply(
        seq_len(nrow(ms1.data)),
        FUN = temp.fun,
        BPPARAM = bpparam,
        ms1.data = ms1.data,
        ms1.match.ppm = ms1.match.ppm,
        rt.match.tol = rt.match.tol,
        database = database,
        adduct.table = adduct.table,
        candidate.num = candidate.num
      )
    
    names(match.result) <- ms1.data$name
    
    temp.idx <-
      which(unlist(lapply(match.result, function(x) {
        all(is.na(x))
      })))
    
    if (length(temp.idx) > 0) {
      match.result <- match.result[-temp.idx]
    }
    
    match.result <-
      lapply(match.result, function(x) {
        if (is.null(x)) {
          return(x)
        } else{
          x$mz.match.score <-
            exp(-0.5 * (x$mz.error / (ms1.match.ppm)) ^ 2)
          if (rt.match.tol > 10000) {
            # cat(crayon::yellow("You set rt.match.tol as NA, so RT will not be used for matching.\n"))
            x$RT.error <- NA
            x$RT.match.score <- NA
            x$Total.score <- x$mz.match.score
          } else{
            x$RT.match.score <-
              exp(-0.5 * (x$RT.error / (rt.match.tol)) ^ 2)
            
            x$Total.score <- x$mz.match.score * 0.5 +
              x$RT.match.score * 0.5
            x$Total.score[is.na(x$Total.score)] <-
              x$mz.match.score[is.na(x$Total.score)]
          }
          x$CE <- NA
          x$SS <- NA
          return(x)
        }
      })
    
    
    annotation_result =
      purrr::map2(
        .x = names(match.result),
        .y = match.result,
        .f = function(x, y) {
          data.frame(variable_id = x, y)
        }
      ) %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(Database = database.name) %>%
      dplyr::mutate(ms2_files_id = NA,
                    ms2_spectrum_id = NA) %>%
      dplyr::select(variable_id,
                    ms2_files_id,
                    ms2_spectrum_id,
                    dplyr::everything())
    
    message(crayon::bgRed("All done."))
    return(annotation_result)
  }