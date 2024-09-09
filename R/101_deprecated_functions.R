




#' @title Get MS2 spectra of peaks from databaseClass object
#' @description Get MS2 spectra of peaks from databaseClass object.
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @param lab.id The lab ID of metabolite.
#' @param database Database (databaseClass object).
#' @param polarity positive or negative.
#' @param ce Collision value.
#' @importFrom crayon red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @importFrom readr cols
#' @importFrom pbapply pblapply
#' @return A MS2 spectrum (data.frame).
#' @export
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/metid.html}

getMS2spectrum <-
  function(lab.id,
           database,
           polarity = c("positive", "negative"),
           ce = "30") {
    message("`getMS2spectrum()` is deprecated, use `get_ms2_spectrum()`.")
    polarity <- match.arg(polarity)
    if (!is(database, "databaseClass")) {
      stop("The database must be databaseClass object.\n")
    }
    pol <- ifelse(polarity == "positive", 1, 2)
    temp <-
      database@spectra.data[[pol]][[match(lab.id, names(database@spectra.data[[pol]]))]]
    temp[[match(ce, names(temp))]]
  }





#' @title Generate the mzIdentify parameter list
#' @description Generate the mzIdentify parameter list.
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @param ms1.match.ppm Precursor match ppm tolerance.
#' @param polarity The polarity of data, "positive"or "negative".
#' @param column "hilic" (HILIC column) or "rp" (reverse phase).
#' @param candidate.num The number of candidates.
#' @param database MS1 database name.
#' @param threads Number of threads
#' @return A mzIdentifyClass object.
#' @export
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/metid.html}

mzIdentifyParam <-
  function(ms1.match.ppm = 25,
           polarity = c("positive", "negative"),
           column = c("rp", "hilic"),
           candidate.num = 1,
           database,
           threads = 3) {
    message("`mzIdentifyParam()` is deprecated.")
    
    if (missing(database)) {
      stop("The database name must be provided.\n")
    }
    if (database != "HMDB.metabolite.data") {
      stop("database can only be HMDB.metabolite.data\n")
    }
    polarity <- match.arg(polarity)
    column <- match.arg(column)
    
    param <- list(
      ms1.match.ppm = ms1.match.ppm,
      polarity = polarity,
      column = column,
      candidate.num = candidate.num,
      database = database,
      threads = threads
    )
    list("mzIdentifyParam" = param)
  }


#' @title Show the base information of metid pacakge
#' @description Show the base information of metid pacakge.
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @return A ASCII log of metid
#' @export
#' @examples
#' metid()

metid <-
  function() {
    message("`metid()` is deprecated, please use `metid_logo()`.")
    cat(crayon::green(
      c(
        "                _    _____  ___ ",
        " _ __ ___   ___| |_  \\_   \\/   \\",
        "| '_ ` _ \\ / _ \\ __|  / /\\/ /\\ /",
        "| | | | | |  __/ |_/\\/ /_/ /_// ",
        "|_| |_| |_|\\___|\\__\\____/___,'  ",
        "                                "
      )
      
    ), sep = "\n")
  }






# library(cowsay)
# # https://onlineasciitools.com/convert-text-to-ascii-art
# # writeLines(capture.output(say("Hello"), type = "message"), con = "ascii_art.txt")
# art <- readLines("logo.txt")
# dput(art)
# metid_logo <-
#   c("                _    _____  ___ ", " _ __ ___   ___| |_  \\_   \\/   \\",
#     "| '_ ` _ \\ / _ \\ __|  / /\\/ /\\ /", "| | | | | |  __/ |_/\\/ /_/ /_// ",
#     "|_| |_| |_|\\___|\\__\\____/___,'  ", "                                "
#   )
# cat(metid_logo, sep = "\n")









#' @title Identify peaks based on MS1 database
#' @description Identify peaks based on MS1 database.
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @param ms1.data The name of ms1 peak table (csv format). Column 1 is "name", column 2 is
#' "mz" and column 3 is "rt" (retention time, second).
#' @param ms1.match.ppm Precursor match ppm tolerance.
#' @param rt.match.tol RT match tolerance.
#' @param polarity The polarity of data, "positive"or "negative".
#' @param column "hilic" (HILIC column) or "rp" (reverse phase).
#' @param path Work directory.
#' @param candidate.num The number of candidates.
#' @param database MS1 database name or MS1 database.
#' @param threads Number of threads
#' @param silence.deprecated Silenc the deprecated information or not.
#' @return A mzIdentifyClass or metIdentifyClass object.
#' @importFrom magrittr %>%
#' @importFrom dplyr pull filter
#' @export
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/metid.html}


mzIdentify <-
  function(ms1.data,
           ##csv format
           ms1.match.ppm = 25,
           rt.match.tol = 30,
           polarity = c("positive", "negative"),
           column = c("hilic", "rp"),
           path = ".",
           candidate.num = 3,
           database,
           threads = 3,
           silence.deprecated = FALSE) {
    if (!silence.deprecated) {
      message("`mzIdentify()` is deprecated.")
    }
    options(warn = -1)
    ###Check data
    if (missing(database)) {
      stop("No database is provided.\n")
    }
    
    ##parameter specification
    polarity <- match.arg(polarity)
    column <- match.arg(column)
    ##check ms1.file and ms2.file
    file <- dir(path)
    
    if (!all(ms1.data %in% file)) {
      stop("MS1 data is not in the directory, please check it.\n")
    }
    
    if (!is(database, "databaseClass")) {
      if (!all(database %in% file)) {
        stop("Database is not in this directory, please check it.\n")
      }
    }
    
    #load database
    if (!is(database, "databaseClass")) {
      database.name <- database
      load(file.path(path, database.name))
      database <- get(database.name)
    } else{
      database.name = paste(database@database.info$Source,
                            database@database.info$Version,
                            sep = "_")
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
    
    ms1.data <-
      readr::read_csv(file = file.path(path, ms1.data),
                      col_types = readr::cols())
    
    colnames(ms1.data)[seq_len(3)] <- c("name", "mz", "rt")
    
    
    if (rt.match.tol > 10000) {
      message(crayon::yellow(
        "You set rt.match.tol > 10,000, so RT will not be used for matching."
      ))
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
        as.numeric(dplyr::pull(.data = ms1.data[idx, ], var = "mz"))
      temp.rt <-
        as.numeric(dplyr::pull(.data = ms1.data[idx, ], var = "rt"))
      
      rm(list = c("ms1.data"))
      
      if (is(database, "databaseClass")) {
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
      
      spectra.mz <- purrr::map(as.data.frame(t(adduct.table)), function(x) {
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
      
      match.idx <- data.frame(match.idx, RT.error, stringsAsFactors = FALSE)
      
      match.idx <-
        dplyr::filter(match.idx, is.na(RT.error) |
                        RT.error < rt.match.tol)
      
      if (nrow(match.idx) == 0) {
        return(NA)
      }
      
      if (nrow(match.idx) > candidate.num) {
        match.idx <- match.idx[seq_len(candidate.num), , drop = FALSE]
      }
      
      match.idx <- data.frame(match.idx, database[match(match.idx$Lab.ID, database$Lab.ID), c("Compound.name", "CAS.ID", "HMDB.ID", "KEGG.ID"), drop = FALSE], stringsAsFactors = FALSE)
      
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
      bpparam = BiocParallel::SnowParam(workers = threads, progressbar = TRUE)
    } else{
      bpparam = BiocParallel::MulticoreParam(workers = threads, progressbar = TRUE)
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
    
    if (is(database, "databaseClass")) {
      return.result <- new(
        Class = "metIdentifyClass",
        ms1.data = ms1.data,
        # ms1.info = ms1.info,
        # ms2.info = ms2.info,
        identification.result = match.result,
        # match.result = match.result,
        adduct.table = adduct.table,
        # ms1.ms2.match.mz.tol = ms1.ms2.match.mz.tol,
        # ms1.ms2.match.rt.tol = ms1.ms2.match.rt.tol,
        ms1.match.ppm = ms1.match.ppm,
        # ms2.match.ppm = ms2.match.ppm,
        # ms2.match.tol = ms2.match.tol,
        rt.match.tol = rt.match.tol,
        polarity = polarity,
        # ce = paste(ce, collapse = ";"),
        column = column,
        # ms1.match.weight = ms1.match.weight,
        # rt.match.weight = rt.match.weight,
        # ms2.match.weight = ms2.match.weight,
        path = path,
        # total.score.tol = total.score.tol,
        candidate.num = candidate.num,
        database = database.name,
        threads = threads,
        version = "1.0.0"
      )
      
      return.result@identification.result <-
        lapply(return.result@identification.result, function(x) {
          if (is.null(x)) {
            return(x)
          } else{
            x$mz.match.score <-
              exp(-0.5 * (x$mz.error / (ms1.match.ppm)) ^ 2)
            if (rt.match.tol > 10000) {
              # cat(crayon::yellow("You set rt.match.tol > 10,000, so RT will not be used for matching.\n"))
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
            x$SS <- 0
            return(x)
          }
        })
    } else{
      return.result <- new(
        Class = "mzIdentifyClass",
        ms1.data = ms1.data,
        identification.result = match.result,
        adduct.table = adduct.table,
        ms1.match.ppm = ms1.match.ppm,
        polarity = polarity,
        column = column,
        path = path,
        candidate.num = candidate.num,
        database = database.name,
        threads = threads
      )
    }
    message(crayon::bgRed("All done."))
    return(return.result)
  }



##------------------------------------------------------------------------------
#' @title get_identification_table_all
#' @description Extract the identifications from multiple results of `identify_metabolite_all()`.
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @param ... One or multiple results from `identify_metabolite_all()`.
#' @param candidate.num candidate.num
#' @param level_condition Condition for level assign.
#' @return A identification table (data.frame).
#' @export
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom magrittr %>%
#' @importFrom stringr str_split
#' @importFrom purrr map map2
#' @importFrom crayon green red yellow
#' @importFrom dplyr select mutate everything filter left_join
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/metid.html}

get_identification_table_all <-
  function(...,
           candidate.num = 1,
           level_condition = c(
             "mz&rt&ms2" = 1,
             "mz&rt" = 2,
             "mz&ms2" = 2,
             "mz" = 3
           )) {
    message("This function is deprecated")
    options(warn = -1)
    result = list(...)
    ##rename result
    result =
      purrr::map2(
        .x = result,
        .y = seq_along(result),
        .f = function(x, y) {
          names(x) = paste(y, names(x), sep = "_")
          x
        }
      )
    
    result = unlist(result)
    
    database_level =
      result %>%
      purrr::map(function(y) {
        database = y@database
        mz_match = ifelse(is.na(y@identification.result[[1]]$mz.error[1]), "", "mz")
        rt_match = ifelse(is.na(y@identification.result[[1]]$RT.error[1]), "", "rt")
        ms2_match = ifelse(all(y@identification.result[[1]]$SS == 0), "", "ms2")
        final = paste(mz_match, rt_match, ms2_match, sep = "&")
        final = stringr::str_split(string = final, pattern = "&")[[1]]
        final = paste(final[final != ""], collapse = "&")
        level = level_condition[match(final, names(level_condition))]
        c(database, level)
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    colnames(database_level)[seq_len(2)] = c("database", "level")
    database_level =
      database_level %>%
      tibble::rownames_to_column(var = "name")
    
    database_level$level[is.na(database_level$level)] = 5
    
    ####annotation table level 1
    if (any(database_level$level == 1)) {
      message(crayon::green("Level 1 table..."))
      annotation_table1 <-
        get_identification_table(result[which(database_level$level == 1)], type = "new", candidate.num = candidate.num)
      annotation_table1 <-
        data.frame(annotation_table1,
                   Level = 1,
                   stringsAsFactors = FALSE)
      message(crayon::red("OK"))
    } else{
      message(crayon::yellow("No level 1."))
      annotation_table1 = NULL
    }
    
    ####annotation table level 2
    if (any(database_level$level == 2)) {
      message(crayon::green("Level 2 table..."))
      annotation_table2 <-
        get_identification_table(result[which(database_level$level == 2)], type = "new", candidate.num = candidate.num)
      annotation_table2 <-
        data.frame(annotation_table2,
                   Level = 2,
                   stringsAsFactors = FALSE)
      message(crayon::red("OK"))
    } else{
      message(crayon::yellow("No level 2."))
      annotation_table2 = NULL
    }
    
    ####annotation table level 3
    if (any(database_level$level == 3)) {
      message(crayon::green("Level 3 table..."))
      annotation_table3 <-
        get_identification_table(result[which(database_level$level == 3)], type = "new", candidate.num = candidate.num)
      annotation_table3 <-
        data.frame(annotation_table3,
                   Level = 3,
                   stringsAsFactors = FALSE)
      message(crayon::red("OK"))
    } else{
      message(crayon::yellow("No level 3."))
      annotation_table3 = NULL
    }
    
    ###remove the redundant result
    if (!is.null(annotation_table1)) {
      annotation_table1 =
        annotation_table1 %>%
        dplyr::filter(!is.na(Compound.name))
    }
    
    if (!is.null(annotation_table2)) {
      annotation_table2 =
        annotation_table2 %>%
        dplyr::filter(!is.na(Compound.name))
      if (!is.null(annotation_table1)) {
        annotation_table2 =
          annotation_table2 %>%
          dplyr::filter(!name %in% annotation_table1$name)
      }
    }
    
    if (!is.null(annotation_table3)) {
      annotation_table3 =
        annotation_table3 %>%
        dplyr::filter(!is.na(Compound.name))
      
      if (!is.null(annotation_table1)) {
        annotation_table3 =
          annotation_table3 %>%
          dplyr::filter(!name %in% annotation_table1$name)
      }
      
      if (!is.null(annotation_table2)) {
        annotation_table3 =
          annotation_table3 %>%
          dplyr::filter(!name %in% annotation_table2$name)
      }
    }
    
    annotation_table =
      rbind(annotation_table1, annotation_table2, annotation_table3)
    
    ms1_peak = result[[1]]@ms1.data
    
    diff_name =
      setdiff(ms1_peak$name, annotation_table$name)
    
    if (length(diff_name) > 0) {
      new_matrix =
        ms1_peak[match(diff_name, ms1_peak$name), ]
      new_matrix =
        new_matrix %>%
        dplyr::mutate(
          MS2.spectra.name = NA,
          Candidate.number = NA,
          Compound.name = NA,
          CAS.ID = NA,
          HMDB.ID = NA,
          KEGG.ID = NA,
          Lab.ID = NA,
          Adduct = NA,
          mz.error = NA,
          mz.match.score = NA,
          RT.error = NA,
          RT.match.score = NA,
          CE = NA,
          SS = NA,
          Total.score = NA,
          Database = NA,
          Level = NA
        )
      annotation_table =
        rbind(annotation_table, new_matrix) %>%
        as.data.frame()
    }
    
    ms1_peak_name = colnames(ms1_peak)
    ms1_peak_name = ms1_peak_name[ms1_peak_name != "name"]
    
    annotation_table =
      annotation_table %>%
      dplyr::select(-ms1_peak_name) %>%
      dplyr::left_join(ms1_peak, by = "name") %>%
      dplyr::select(name, ms1_peak_name, dplyr::everything())
    
    if (candidate.num == 1) {
      annotation_table =
        annotation_table[match(result[[1]]@ms1.data$name, annotation_table$name), ]
    }
    
    message(crayon::red("All done."))
    return(tibble::as_tibble(annotation_table))
  }





##------------------------------------------------------------------------------
#' @title Get identification table from a metIdentifyClass object
#' @description Get identification table from a metIdentifyClass object.
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @param ... One or multiple metIdentifyClass objects.
#' @param candidate.num The number of candidates.
#' @param type The type of identification table.
#' @return A identification table (data.frame).
#' @export
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate everything lag filter
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/metid.html}

get_identification_table <-
  function(...,
           candidate.num = 3,
           type = c("old", "new")) {
    message("This function is deprecated")
    options(warn = -1)
    candidate.num <- round(candidate.num)
    if (candidate.num <= 0) {
      candidate.num <- 1
    }
    
    object <- list(...)
    
    ###if object is a one list or not
    if (length(object) == 1 & is.list(object[[1]])) {
      object = object[[1]]
    }
    
    if (any(unique(unlist(lapply(object, class))) != "metIdentifyClass")) {
      stop("Only for metIdentifyClass\n")
    }
    
    type <- match.arg(type)
    
    ##if object number is 1
    if (length(object) == 1) {
      object <- object[[1]]
      database <- object@database
      
      identification.result <- object@identification.result
      
      if (is.null(identification.result[[1]])) {
        return(NULL)
      }
      
      if (nrow(object@match.result) == 0) {
        message(crayon::yellow("The object is identified without MS2 spectra."))
        return(
          getIdentificationTable2(
            object = object,
            candidate.num = candidate.num,
            type = type,
            silence.deprecated = TRUE
          )
        )
      }
      
      ##add database information
      identification.result <-
        lapply(identification.result, function(x) {
          if (nrow(x) > candidate.num) {
            x <- x[seq_len(candidate.num), , drop = FALSE]
          }
          data.frame(x,
                     "Database" = object@database,
                     stringsAsFactors = FALSE)
        })
      
      peak.table <- object@ms1.data
      match.result <- object@match.result
      
      if (type == "old") {
        identification.table <-
          as.data.frame(matrix(nrow = nrow(peak.table), ncol = 3))
        colnames(identification.table) <-
          c("MS2.spectra.name",
            "Candidate.number",
            "Identification")
        identification.table[match.result[, 1], 1] <-
          object@ms1.info$name[match.result[, 2]]
        item <- colnames(identification.result[[1]])
        identification.result <-
          lapply(identification.result, function(x) {
            paste(apply(x, 1, function(y) {
              paste(paste(item, as.character(y), sep = ":"), collapse = ";")
            }), collapse = "{}")
          })
        
        identification.table$Identification[match(names(identification.result),
                                                  identification.table$MS2.spectra.name)] <-
          unlist(identification.result)
        
        identification.table$Candidate.number <-
          lapply(identification.table$Identification, function(x) {
            if (is.na(x))
              return(0)
            return(length(stringr::str_split(
              string = x, pattern = "\\{\\}"
            )[[1]]))
          }) %>%
          unlist()
        identification.table <-
          data.frame(peak.table, identification.table, stringsAsFactors = FALSE)
      } else {
        identification.table <-
          vector(mode = "list", length = nrow(peak.table))
        names(identification.table)[match.result[, 1]] <-
          object@ms1.info$name[match.result[, 2]]
        identification.table[match(names(identification.result),
                                   names(identification.table))] <-
          identification.result
        peak.table <- apply(peak.table, 1, list)
        peak.table <- lapply(peak.table, unlist)
        
        identification.table <-
          mapply(
            FUN = function(x, y) {
              if (all(is.na(y))) {
                temp <-
                  as.data.frame(matrix(c(x, rep(NA, 14)), nrow = 1), stringsAsFactors = FALSE)
                colnames(temp) <-
                  c(
                    names(x),
                    c(
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
                      "Database"
                    )
                  )
                list(temp)
              } else{
                temp <-
                  as.data.frame(matrix(rep(x, nrow(y)), nrow = nrow(y), byrow = TRUE),
                                stringsAsFactors = FALSE)
                if (nrow(temp) > 1) {
                  temp[2:nrow(temp), 2:ncol(temp)] <- ""
                }
                colnames(temp) <- names(x)
                temp <-
                  data.frame(temp, y, stringsAsFactors = FALSE)
                list(temp)
              }
            },
            x = peak.table,
            y = identification.table
          )
        
        identification.table <-
          as.data.frame(do.call(rbind, identification.table))
        
        identification.table$mz <-
          identification.table$mz %>%
          stringr::str_trim(side = 'both') %>%
          as.numeric()
        
        identification.table$rt <-
          identification.table$rt %>%
          stringr::str_trim(side = 'both') %>%
          as.numeric()
        
        identification.table$mz.error <-
          identification.table$mz.error %>%
          stringr::str_trim(side = 'both') %>%
          as.numeric()
        
        identification.table$mz.match.score <-
          identification.table$mz.match.score %>%
          stringr::str_trim(side = 'both') %>%
          as.numeric()
        
        identification.table$RT.error <-
          identification.table$RT.error %>%
          stringr::str_trim(side = 'both') %>%
          as.numeric()
        
        identification.table$RT.match.score <-
          identification.table$RT.match.score %>%
          stringr::str_trim(side = 'both') %>%
          as.numeric()
        
        identification.table$SS <-
          identification.table$SS %>%
          stringr::str_trim(side = 'both') %>%
          as.numeric()
        
        identification.table$Total.score <-
          identification.table$Total.score %>%
          stringr::str_trim(side = 'both') %>%
          as.numeric()
        
        ###add Candidate.number
        Candidate.number = purrr::map2(
          .x = object@identification.result,
          .y = names(object@identification.result),
          .f = function(x, y) {
            data.frame(
              MS2.spectra.name = y,
              Candidate.number = nrow(x),
              stringsAsFactors = FALSE
            )
          }
        ) %>%
          do.call(rbind, .) %>%
          as.data.frame()
        Candidate.number =
          Candidate.number %>%
          dplyr::left_join(object@match.result[, c("MS1.peak.name", "MS2.spectra.name")], by = "MS2.spectra.name") %>%
          dplyr::select(MS1.peak.name, dplyr::everything())
        
        identification.table =
          identification.table %>%
          dplyr::left_join(Candidate.number, by = c("name" = "MS1.peak.name")) %>%
          dplyr::select(
            colnames(object@ms1.data),
            MS2.spectra.name,
            Candidate.number,
            dplyr::everything()
          )
        
        identification.table$MS2.spectra.name[which(is.na(identification.table$Compound.name))] <- NA
        identification.table$Candidate.number[which(is.na(identification.table$Compound.name))] <- NA
        
      }
      rownames(identification.table) <- NULL
      return(tibble::as_tibble(identification.table))
    } else{
      len <- lapply(object, function(x) {
        nrow(x@match.result)
      }) %>% unlist()
      
      if (any(len == 0)) {
        warning(
          crayon::yellow(
            "Some objects are identified without MS2 spectra.\nPlease make sure that you want to compare them together.\n"
          )
        )
      }
      
      ###if object number is not 1
      peak.table <- object[[1]]@ms1.data
      database.name <-
        unlist(lapply(object, function(x)
          x@database))
      
      identification.result <-
        lapply(object, function(x) {
          if (nrow(x@match.result) != 0) {
            iden.result <- x@identification.result
            iden.result <- mapply(function(x, y) {
              if (nrow(x) > candidate.num) {
                x <- x[seq_len(candidate.num), , drop = FALSE]
              }
              x <- data.frame(x,
                              "MS2.spectra.name" = y,
                              stringsAsFactors = FALSE)
              list(x)
            }, x = iden.result, y = names(iden.result))
            match.result <- x@match.result
            # names(iden.result) <-
            #   match.result$MS1.peak.name[match(names(iden.result), match.result$MS2.spectra.name)]
            iden.result <-
              do.call(rbind, iden.result)
            Peak.name <-
              match.result$MS1.peak.name[match(iden.result$MS2.spectra.name,
                                               match.result$MS2.spectra.name)]
            iden.result <-
              data.frame(Peak.name, iden.result, stringsAsFactors = FALSE)
            rownames(iden.result) <- NULL
            return(iden.result)
          } else{
            iden.result <- x@identification.result
            iden.result <- mapply(function(x, y) {
              if (nrow(x) > candidate.num) {
                x <- x[seq_len(candidate.num), , drop = FALSE]
              }
              x <- data.frame(x,
                              "MS2.spectra.name" = y,
                              stringsAsFactors = FALSE)
              list(x)
            }, x = iden.result, y = names(iden.result))
          }
          iden.result <-
            do.call(rbind, iden.result)
          
          iden.result <-
            iden.result %>%
            dplyr::mutate(Peak.name = MS2.spectra.name, MS2.spectra.name = NA) %>%
            dplyr::select(Peak.name, dplyr::everything())
          
        })
      
      identification.result <- mapply(
        function(x, y, z) {
          list(data.frame(
            x,
            "Database" = y,
            "Dataset.index" = z,
            stringsAsFactors = FALSE
          ))
        },
        x = identification.result,
        y = database.name,
        z = seq_along(identification.result)
      )
      
      identification.result <-
        do.call(rbind, identification.result)
      identification.result <-
        as.data.frame(identification.result)
      
      identification.result <-
        plyr::dlply(
          .data = identification.result,
          .variables = plyr::.(Peak.name),
          .fun = function(x) {
            x <- x[order(x$SS, decreasing = TRUE), , drop = FALSE]
            if (nrow(x) > candidate.num) {
              x <- x[seq_len(candidate.num), , drop = FALSE]
            }
            x
          }
        )
      
      # if(type == "old"){
      identification.table <-
        as.data.frame(matrix(nrow = nrow(peak.table), ncol = 3))
      colnames(identification.table) <-
        c("MS2.spectra.name", "Candidate.number", "Identification")
      # identification.table[match.result[,1],1] <- object@ms1.info$name[match.result[,2]]
      item <- colnames(identification.result[[1]])[-1]
      identification.result <-
        lapply(identification.result, function(x) {
          paste(apply(x[, -1], 1, function(y) {
            paste(paste(item, as.character(y), sep = ":"), collapse = ";")
          }), collapse = "{}")
        })
      
      identification.table$Identification[match(names(identification.result), peak.table$name)] <-
        unlist(identification.result)
      
      identification.table$Candidate.number <-
        lapply(identification.table$Identification, function(x) {
          if (is.na(x))
            return(0)
          return(length(stringr::str_split(
            string = x, pattern = "\\{\\}"
          )[[1]]))
        }) %>%
        unlist()
      
      identification.table$MS2.spectra.name <-
        unlist(lapply(identification.table$Identification, function(x) {
          if (is.na(x))
            return(NA)
          temp.name <-
            stringr::str_extract_all(string = x, pattern = "MS2.spectra.name:mz[0-9\\.\\-\\_]{1,30}rt[0-9\\.\\-\\_]{1,30}")[[1]]
          temp.name <-
            stringr::str_replace_all(string = temp.name,
                                     pattern = "MS2.spectra.name:",
                                     replacement = "")
          paste(temp.name, collapse = ";")
        }))
      
      identification.table <-
        data.frame(peak.table, identification.table, stringsAsFactors = FALSE)
      if (type == "new") {
        identification.table <-
          trans_to_new_style(identification.table = identification.table)
        identification.table$mz <-
          identification.table$mz %>%
          stringr::str_trim(side = 'both') %>%
          as.numeric()
        
        identification.table$rt <-
          identification.table$rt %>%
          stringr::str_trim(side = 'both') %>%
          as.numeric()
        
        identification.table$mz.error <-
          identification.table$mz.error %>%
          stringr::str_trim(side = 'both') %>%
          as.numeric()
        
        identification.table$mz.match.score <-
          identification.table$mz.match.score %>%
          stringr::str_trim(side = 'both') %>%
          as.numeric()
        
        identification.table$RT.error <-
          identification.table$RT.error %>%
          stringr::str_trim(side = 'both') %>%
          as.numeric()
        
        identification.table$RT.match.score <-
          identification.table$RT.match.score %>%
          stringr::str_trim(side = 'both') %>%
          as.numeric()
        
        identification.table$SS <-
          identification.table$SS %>%
          stringr::str_trim(side = 'both') %>%
          as.numeric()
        
        identification.table$Total.score <-
          identification.table$Total.score %>%
          stringr::str_trim(side = 'both') %>%
          as.numeric()
      }
      return(tibble::as_tibble(identification.table))
    }
  }

##------------------------------------------------------------------------------
#' @title Transform old style identification table to new style
#' @description Transform old style identification table to new style.
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @param identification.table Identification table from get_identification_table.
#' @return A identification table (data.frame).
#' @export
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/metid.html}

trans_to_new_style <-
  function(identification.table) {
    message("This function is deprecated")
    if (all(colnames(identification.table) != "Identification")) {
      return(identification.table)
    } else{
      identification <- identification.table$Identification
      identification <-
        pbapply::pblapply(identification, function(x) {
          if (is.na(x)) {
            return(
              data.frame(
                "Compound.name" = NA,
                "CAS.ID" = NA,
                "HMDB.ID" = NA,
                "KEGG.ID" = NA,
                "Lab.ID" = NA,
                "Adduct" = NA,
                "mz.error" = NA,
                "mz.match.score" = NA,
                "RT.error" = NA,
                "RT.match.score" = NA,
                "CE" = NA,
                "SS" = NA,
                "Total.score" = NA,
                "Database" = NA,
                stringsAsFactors = FALSE
              )
            )
          } else{
            x <- stringr::str_split(string = x, pattern = "\\{\\}")[[1]][1]
            Compound.name <-
              stringr::str_split(string = x,
                                 pattern = "\\;CAS\\.ID",
                                 n = 2)[[1]][1]
            # require(magrittr)
            Compound.name <- Compound.name %>%
              stringr::str_replace(string = .,
                                   pattern = "Compound.name\\:",
                                   replacement = "")
            
            x <-
              stringr::str_split(string = x, pattern = ";")[[1]]
            
            x <-
              lapply(c(
                "CAS.ID\\:",
                "HMDB.ID\\:",
                "KEGG.ID\\:",
                "Lab.ID\\:",
                "Adduct\\:",
                "mz.error\\:",
                "mz.match.score\\:",
                "RT.error\\:",
                "RT.match.score\\:",
                "CE\\:",
                "SS\\:",
                "Total.score\\:",
                "Database\\:"
              ), function(y) {
                y <- grep(y, x, value = TRUE) %>%
                  stringr::str_replace(
                    string = .,
                    pattern = paste(y, "\\:", sep = ""),
                    replacement = ""
                  )
                if (length(y) == 0) {
                  return(NA)
                }
                y
                
              }) %>%
              unlist()
            x <- stringr::str_replace(
              string = x,
              pattern = c(
                "CAS.ID\\:",
                "HMDB.ID\\:",
                "KEGG.ID\\:",
                "Lab.ID\\:",
                "Adduct\\:",
                "mz.error\\:",
                "mz.match.score\\:",
                "RT.error\\:",
                "RT.match.score\\:",
                "CE\\:",
                "SS\\:",
                "Total.score\\:",
                "Database\\:"
              ),
              replacement = ""
            )
            x <-
              matrix(c(Compound.name, x),
                     nrow = 1,
                     byrow = TRUE)
            colnames(x) <-
              c(
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
                "Database"
              )
            return(as.data.frame(x, stringsAsFactors = FALSE))
          }
        })
      
      identification <- do.call(rbind, identification)
      identification.table <-
        tibble::as_tibble(identification.table)
      identification.table <-
        dplyr::select(.data = identification.table, -(Identification))
      identification.table <-
        cbind(identification.table, identification)
      
      invisible(tibble::as_tibble(identification.table))
      # identification1 <- dplyr::bind_rows(identification)
    }
  }





##------------------------------------------------------------------------------
#' @title Get identification table from a mzIdentifyClass object
#' @description Get identification table from a mzIdentifyClass object.
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @param object mzIdentifyClass object.
#' @param candidate.num The number of candidates.
#' @param type The type of identification table.
#' @param silence.deprecated Silenc the deprecated information or not.
#' @return A identification table (data.frame).
#' @export
#' @importFrom magrittr %>%
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/metid.html}

# object = result_rp_pos25[[8]]

getIdentificationTable2 <-
  function(object,
           candidate.num = 3,
           type = c("old", "new"),
           silence.deprecated = FALSE) {
    if (!silence.deprecated) {
      message("`getIdentificationTable2()` is deprecated, use `get_identification_table()`.")
    }
    
    if (!is(object, "mzIdentifyClass") &
        !is(object, "metIdentifyClass")) {
      stop("Only for mzIdentifyClass or metIdentifyClass\n")
    }
    
    type <- match.arg(type)
    database <- object@database
    
    identification.result <- object@identification.result
    
    if (is.null(identification.result[[1]])) {
      return(NULL)
    }
    
    ##add database information
    identification.result <-
      lapply(identification.result, function(x) {
        if (nrow(x) > candidate.num) {
          x <- x[seq_len(candidate.num), , drop = FALSE]
        }
        data.frame(x,
                   "Database" = object@database,
                   stringsAsFactors = FALSE)
      })
    
    peak.table <- object@ms1.data
    
    if (type == "old") {
      identification.table <-
        as.data.frame(matrix(nrow = nrow(peak.table), ncol = 2))
      colnames(identification.table) <-
        c("Candidate.number", "Identification")
      item <- colnames(identification.result[[1]])
      identification.result <-
        lapply(identification.result, function(x) {
          paste(apply(x, 1, function(y) {
            paste(paste(item, as.character(y), sep = ":"), collapse = ";")
          }), collapse = "{}")
        })
      
      identification.table$Identification[match(names(identification.result), peak.table$name)] <-
        unlist(identification.result)
      
      identification.table$Candidate.number <-
        lapply(identification.table$Identification, function(x) {
          if (is.na(x))
            return(0)
          return(length(stringr::str_split(
            string = x, pattern = "\\{\\}"
          )[[1]]))
        }) %>%
        unlist()
      identification.table <-
        data.frame(peak.table, identification.table, stringsAsFactors = FALSE)
    } else{
      identification.table <-
        vector(mode = "list", length = nrow(peak.table))
      
      names(identification.table) <- object@ms1.data$name
      
      identification.table[match(names(identification.result),
                                 names(identification.table))] <-
        identification.result
      
      peak.table <- apply(peak.table, 1, list)
      peak.table <- lapply(peak.table, unlist)
      
      identification.table <-
        mapply(
          FUN = function(x, y) {
            if (all(is.na(y)) | is.null(y)) {
              temp <-
                as.data.frame(matrix(c(x, rep(NA, 14)), nrow = 1), stringsAsFactors = FALSE)
              
              colnames(temp) <-
                c(
                  names(x),
                  c(
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
                    "Database"
                  )
                )
              list(temp)
            } else{
              temp <-
                as.data.frame(matrix(rep(x, nrow(y)), nrow = nrow(y), byrow = TRUE),
                              stringsAsFactors = FALSE)
              if (nrow(temp) > 1) {
                temp[2:nrow(temp), 2:ncol(temp)] <- ""
              }
              colnames(temp) <- names(x)
              temp <-
                data.frame(temp, y, stringsAsFactors = FALSE)
              list(temp)
            }
          },
          x = peak.table,
          y = identification.table
        )
      
      identification.table <-
        as.data.frame(do.call(rbind, identification.table))
      ###as.numeric for several column
      identification.table$mz =
        identification.table$mz %>%
        stringr::str_trim(side = "both") %>%
        as.numeric()
      
      identification.table$rt =
        identification.table$rt %>%
        stringr::str_trim(side = "both") %>%
        as.numeric()
      
      identification.table$mz.error =
        identification.table$mz.error %>%
        stringr::str_trim(side = "both") %>%
        as.numeric()
      
      identification.table$RT.error =
        identification.table$RT.error %>%
        stringr::str_trim(side = "both") %>%
        as.numeric()
      
      identification.table$mz.match.score =
        identification.table$mz.match.score %>%
        stringr::str_trim(side = "both") %>%
        as.numeric()
      
      identification.table$RT.match.score =
        identification.table$RT.match.score %>%
        stringr::str_trim(side = "both") %>%
        as.numeric()
      
      identification.table$Total.score =
        identification.table$Total.score %>%
        stringr::str_trim(side = "both") %>%
        as.numeric()
      
      identification.table$SS =
        identification.table$SS %>%
        stringr::str_trim(side = "both") %>%
        as.numeric()
      
      ##add MS2.spectra.name and Candidate.number
      identification.table$MS2.spectra.name = NA
      identification.table$Candidate.number = NA
      identification.table =
        identification.table %>%
        dplyr::select(
          colnames(object@ms1.data),
          MS2.spectra.name,
          Candidate.number,
          Compound.name,
          CAS.ID,
          HMDB.ID,
          KEGG.ID,
          Lab.ID,
          Adduct,
          mz.error,
          mz.match.score,
          RT.error,
          RT.match.score,
          CE,
          SS,
          Total.score,
          Database
        )
      
    }
    
    rownames(identification.table) <- NULL
    
    
    return(tibble::as_tibble(identification.table))
  }







#' @title Identify metabolites using multiple databases one time
#' @description Identify metabolites using multiple databases one time.
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @param ms1.data The name of ms1 peak table (csv format). Column 1 is "name", column 2 is
#' "mz" and column 3 is "rt" (second).
#' @param ms2.data MS2 data, must be mgf, msp or mzXML format. For example, ms2.data = c("test.mgf", "test2.msp").
#' @param parameter.list A list contains paramters for each processing.
#' The parameter must get using metIdentifyParam or mzIdentifyParam.
#' @param path Work directory.
#' @return A list containing mzIdentifyClass object.
#' @export
#' @importFrom magrittr %>%
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/multiple_databases.html}

identify_metabolite_all <-
  function(ms1.data, ms2.data, parameter.list, path = ".") {
    message("This function is deprecated")
    dir.create(path = path, showWarnings = FALSE)
    old.path <- path
    path <- file.path(path, "Result")
    dir.create(path = path, showWarnings = FALSE)
    
    threads = parameter.list[[1]]$threads
    ms1.data.name <- ms1.data
    ms2.data.name <- ms2.data
    
    intermediate_path <-
      file.path(path, "intermediate_data")
    dir.create(intermediate_path, showWarnings = FALSE)
    
    file <- dir(intermediate_path)
    
    if (all(c("ms1.info", "ms2.info") %in% file)) {
      message(crayon::yellow("Use old data"))
      load(file.path(intermediate_path, "ms1.info"))
      load(file.path(intermediate_path, "ms2.info"))
    } else{
      ##read MS2 data
      message(crayon::green("Reading MS2 data..."))
      temp.ms2.type <-
        stringr::str_split(string = ms2.data.name, pattern = "\\.")[[1]]
      temp.ms2.type <- temp.ms2.type[length(temp.ms2.type)]
      
      if (temp.ms2.type %in% c("mzXML", "mzML")) {
        ms2.data <-
          masstools::read_mzxml(file = file.path(old.path, ms2.data.name),
                                threads = threads)
      } else{
        ms2.data <- lapply(ms2.data.name, function(temp.ms2.data) {
          temp.ms2.type <- stringr::str_split(string = temp.ms2.data, pattern = "\\.")[[1]]
          temp.ms2.type <-
            temp.ms2.type[length(temp.ms2.type)]
          if (!temp.ms2.type %in% c("mgf", "msp"))
            stop("We only support mgf or msp.\n")
          if (temp.ms2.type == "msp") {
            temp.ms2.data <- readMSP(file = file.path(old.path, temp.ms2.data))
          } else{
            temp.ms2.data <-
              masstools::read_mgf(file = file.path(old.path, temp.ms2.data))
          }
          temp.ms2.data
        })
        
        names(ms2.data) <- ms2.data.name
        ###prepare data for metIdentification function
        message(crayon::green("Preparing MS2 data for identification..."))
        ms2.data <-
          mapply(
            FUN = function(temp.ms2.data, temp.ms2.data.name) {
              temp.ms2.data <- lapply(temp.ms2.data, function(x) {
                info <- x$info
                info <-
                  data.frame(
                    name = paste("mz", info[1], "rt", info[2], sep = ""),
                    "mz" = info[1],
                    "rt" = info[2],
                    "file" = temp.ms2.data.name,
                    stringsAsFactors = FALSE
                  )
                rownames(info) <- NULL
                x$info <- info
                x
              })
              temp.ms2.data
            },
            temp.ms2.data = ms2.data,
            temp.ms2.data.name = ms2.data.name
          )
        
        if (is.matrix(ms2.data)) {
          ms2.data <- ms2.data[, 1]
        } else{
          ms2.data <- do.call(what = c, args = ms2.data)
        }
      }
      
      ms1.info <- lapply(ms2.data, function(x) {
        x[[1]]
      })
      
      ms2.info <- lapply(ms2.data, function(x) {
        x[[2]]
      })
      
      ms1.info <- do.call(what = rbind, args = ms1.info)
      ms1.info <- as.data.frame(ms1.info)
      rownames(ms1.info) <- NULL
      
      duplicated.name <-
        unique(ms1.info$name[duplicated(ms1.info$name)])
      if (length(duplicated.name) > 0) {
        lapply(duplicated.name, function(x) {
          ms1.info$name[which(ms1.info$name == x)] <-
            paste(x, seq_len(sum(ms1.info$name == x)), sep = "_")
        })
      }
      
      names(ms2.info) <- ms1.info$name
      ##save intermediate data
      save(ms1.info,
           file = file.path(intermediate_path, "ms1.info"),
           compress = "xz")
      save(ms2.info,
           file = file.path(intermediate_path, "ms2.info"),
           compress = "xz")
    }
    
    database_class =
      purrr::map(parameter.list, function(x) {
        class(x$database)
      }) %>%
      unlist()
    
    database.name <-
      unlist(lapply(parameter.list, function(x) {
        if (!is(x$database, "databaseClass")) {
          x$database
        } else{
          paste(x$database@database.info$Source,
                x$database@database.info$Version,
                sep = "_")
        }
      }))
    
    ##check databases with same names
    database.name = make.unique(names = database.name, sep = "_")
    
    ##output database information to intermediate_data path
    database_info =
      data.frame(database.name, database_class, parameter = seq_along(database.name))
    
    write.csv(
      database_info,
      file = file.path(intermediate_path, "database_info.csv"),
      row.names = FALSE
    )
    
    if (!all(database.name[database_class != "databaseClass"] %in% dir(old.path))) {
      stop(
        "The database: ",
        paste(database.name[!database.name %in% dir(old.path)], collapse = ", "),
        "\n",
        " you want to use are not in you directory: \n",
        old.path
      )
    }
    
    identification.result <-
      vector(mode = "list", length = length(database.name))
    
    names(identification.result) <- database.name
    
    for (i in seq_along(database.name)) {
      message(crayon::yellow("-------------------------------"))
      message(crayon::yellow('Database ', i, ": ", database.name[i]))
      message(crayon::yellow("-------------------------------"))
      
      new.path <-
        file.path(path, paste(database.name[i], "Result", sep =  '_'))
      
      dir.create(new.path, showWarnings = FALSE)
      
      if (any(dir(new.path) == "result")) {
        load(file.path(new.path, "result"))
        identification.result[[i]] <- result
        rm(list = "result")
        next()
      }
      
      if (is(parameter.list[[i]]$database, "databaseClass")) {
        temp_database =
          parameter.list[[i]]$database
      } else{
        temp_database <-
          load(file.path(old.path, parameter.list[[i]]$database))
        
        temp_database <-
          get(temp_database)
      }
      
      if (length(temp_database@spectra.data) == 0) {
        rm(list = parameter.list[[i]]$database)
        result <- mzIdentify(
          ms1.data = ms1.data.name,
          ms1.match.ppm = parameter.list[[i]]$ms1.match.ppm,
          rt.match.tol = parameter.list[[i]]$rt.match.tol,
          polarity = parameter.list[[i]]$polarity,
          column = parameter.list[[i]]$column,
          path = old.path,
          candidate.num = parameter.list[[i]]$candidate.num,
          database = parameter.list[[i]]$database,
          threads = parameter.list[[i]]$threads,
          silence.deprecated = TRUE
        )
        
      } else{
        # rm(list = parameter.list[[i]]$database)
        result <- identify_metabolites(
          ms1.data = ms1.data.name,
          ms2.data = ms2.data.name,
          ms1.ms2.match.mz.tol = parameter.list[[i]]$ms1.ms2.match.mz.tol,
          ms1.ms2.match.rt.tol = parameter.list[[i]]$ms1.ms2.match.rt.tol,
          ms1.match.ppm = parameter.list[[i]]$ms1.match.ppm,
          ms2.match.ppm = parameter.list[[i]]$ms2.match.ppm,
          mz.ppm.thr = parameter.list[[i]]$mz.ppm.thr,
          ms2.match.tol = parameter.list[[i]]$ms2.match.tol,
          fraction.weight = parameter.list[[i]]$fraction.weight,
          dp.forward.weight = parameter.list[[i]]$dp.forward.weight,
          dp.reverse.weight = parameter.list[[i]]$dp.reverse.weight,
          rt.match.tol = parameter.list[[i]]$rt.match.tol,
          polarity = parameter.list[[i]]$polarity,
          ce = parameter.list[[i]]$ce,
          column = parameter.list[[i]]$column,
          ms1.match.weight = parameter.list[[i]]$ms1.match.weight,
          rt.match.weight = parameter.list[[i]]$rt.match.weight,
          ms2.match.weight = parameter.list[[i]]$ms2.match.weight,
          path = old.path,
          total.score.tol = parameter.list[[i]]$total.score.tol,
          candidate.num = parameter.list[[i]]$candidate.num,
          database = parameter.list[[i]]$database,
          threads = parameter.list[[i]]$threads
        )
      }
      # unlink(x = new.path, recursive = TRUE, force = TRUE)
      identification.result[[i]] <- result
      save(result, file = file.path(new.path, "result"))
      rm(list = "result")
    }
    invisible(identification.result)
  }


#' @title Generate the parameter list for identify_metabolites function
#' @description Generate the parameter list for metIdentify function.
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @param ms1.ms2.match.mz.tol MS1 peak and MS2 spectrum matching m/z tolerance. Default is 25 pm.
#' @param ms1.ms2.match.rt.tol MS1 peak and MS2 spectrum matching RT tolerance. Default is 10 s.
#' @param ms1.match.ppm Precursor match ppm tolerance.
#' @param ms2.match.ppm Fragment ion match ppm tolerance.
#' @param mz.ppm.thr Accurate mass tolerance for m/z error calculation.
#' @param ms2.match.tol MS2 match (MS2 similarity) tolerance.
#' @param fraction.weight The weight for matched fragments.
#' @param dp.forward.weight Forward dot product weight.
#' @param dp.reverse.weight Reverse dot product weight.
#' @param rt.match.tol RT match tolerance.
#' @param polarity The polarity of data, "positive"or "negative".
#' @param ce Collision energy. Please confirm the CE values in your database. Default is "all".
#' @param column "hilic" (HILIC column) or "rp" (reverse phase).
#' @param ms1.match.weight The weight of MS1 match for total score calculation.
#' @param rt.match.weight The weight of RT match for total score calculation.
#' @param ms2.match.weight The weight of MS2 match for total score calculation.
#' @param total.score.tol Total score tolerance. The total score are refering to MS-DIAL.
#' @param candidate.num The number of candidate.
#' @param database MS2 database name or MS2 database.
#' @param threads Number of threads
#' @return A metIdentifyClass object.
#' @export
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/multiple_databases.html}
#' @examples
#'  param1 <-
#' identify_metabolites_params(
#'   ms1.match.ppm = 15,
#'   rt.match.tol = 15,
#'   polarity = "positive",
#'   ce = "all",
#'   column = "rp",
#'   total.score.tol = 0.5,
#'   candidate.num = 3,
#'   threads = 3,
#'   database = "msDatabase_rplc0.0.2"
#' )
#' param1

identify_metabolites_params <-
  function(ms1.ms2.match.mz.tol = 25,
           ms1.ms2.match.rt.tol = 10,
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
           threads = 3) {
    message("This function is deprecated")
    if (missing(database)) {
      stop("The database or database name must be provided.\n")
    }
    
    polarity <- match.arg(polarity)
    column <- match.arg(column)
    param <-
      list(
        ms1.ms2.match.mz.tol = ms1.ms2.match.mz.tol,
        ms1.ms2.match.rt.tol = ms1.ms2.match.rt.tol,
        ms1.ms2.match.mz.tol = ms1.ms2.match.mz.tol,
        ms1.ms2.match.rt.tol = ms1.ms2.match.rt.tol,
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
    list("metIdentifyParam" = param)
  }


















#' @title Identify metabolites based on MS1 or MS/MS database
#' @description Identify metabolites based on MS1 or MS/MS database.
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @param ms1.data The name of ms1 peak table (csv format). Column 1 is "name", Column 2 is
#' "mz" and column is "rt" (second).
#' @param ms2.data MS2 data, must be mgf, msp or mzXML format. For example, ms2.data = c("test.mgf", "test2.msp").
#' @param ms1.ms2.match.mz.tol MS1 peak and MS2 spectrum matching m/z tolerance. Default is 25 pm.
#' @param ms1.ms2.match.rt.tol MS1 peak and MS2 spectrum matching RT tolerance. Default is 10 s.
#' @param ms1.match.ppm Precursor match ppm tolerance.
#' @param ms2.match.ppm Fragment ion match ppm tolerance.
#' @param mz.ppm.thr Accurate mass tolerance for m/z error calculation.
#' @param ms2.match.tol MS2 match (MS2 similarity) tolerance.
#' @param fraction.weight The weight for matched fragments.
#' @param dp.forward.weight Forward dot product weight.
#' @param dp.reverse.weight Reverse dot product weight.
#' @param rt.match.tol RT match tolerance.
#' @param polarity The polarity of data, "positive"or "negative".
#' @param ce Collision energy. Please confirm the CE values in your database. Default is "all".
#' @param column "hilic" (HILIC column) or "rp" (reverse phase).
#' @param ms1.match.weight The weight of MS1 match for total score calculation.
#' @param rt.match.weight The weight of RT match for total score calculation.
#' @param ms2.match.weight The weight of MS2 match for total score calculation.
#' @param path Work directory.
#' @param total.score.tol Total score tolerance. The total score are refering to MS-DIAL.
#' @param candidate.num The number of candidate.
#' @param database MS2 database name or MS database.
#' @param threads Number of threads
#' @return A metIdentifyClass object.
#' @importFrom crayon yellow green red bgRed
#' @importFrom magrittr %>%
#' @export
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/metid.html}

identify_metabolites <-
  function(ms1.data,
           ms2.data = NULL,
           ms1.ms2.match.mz.tol = 25,
           ms1.ms2.match.rt.tol = 10,
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
           path = ".",
           total.score.tol = 0.5,
           candidate.num = 3,
           database,
           threads = 3) {
    ###Check data
    if (missing(database)) {
      stop("No database is provided.\n")
    }
    
    if (missing(ms1.data)) {
      stop("Please provide MS1 data name.\n")
    }
    
    ##parameter specification
    polarity <- match.arg(polarity)
    column <- match.arg(column)
    ##check ms1.file and ms2.file
    file <- dir(path)
    
    if (!all(ms1.data %in% file)) {
      stop("MS1 data is not in the directory, please check it.\n")
    }
    
    if (!is.null(ms2.data)) {
      if (!all(ms2.data %in% file)) {
        stop("Some MS2 data are not in the directory, please check it.\n")
      }
    }
    
    if (!is(database, "databaseClass")) {
      if (!all(database %in% file)) {
        stop("Database is not in this directory, please check it.\n")
      }
    }
    
    if (is.null(ms2.data)) {
      message(crayon::yellow(
        "You don't provide MS2 data, so only use mz and/or RT for matching."
      ))
      mzIdentify(
        ms1.data = ms1.data,
        rt.match.tol = rt.match.tol,
        ms1.match.ppm = ms1.match.ppm,
        polarity = polarity,
        column = column,
        path = path,
        candidate.num = candidate.num,
        database = database,
        threads = threads,
        silence.deprecated = TRUE
      )
    } else{
      metIdentify(
        ms1.data = ms1.data,
        ms2.data = ms2.data,
        ##only msp and mgf and mz(X)ML are supported
        ms1.ms2.match.mz.tol = ms1.ms2.match.mz.tol,
        ms1.ms2.match.rt.tol = ms1.ms2.match.rt.tol,
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
        path = path,
        total.score.tol = total.score.tol,
        candidate.num = candidate.num,
        database = database,
        threads = threads,
        silence.deprecated = TRUE
      )
    }
  }




#' @title Identify metabolites based on MS/MS database.
#' @description Identify metabolites based on MS/MS database.
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @param ms1.data The name of ms1 peak table (csv format). Column 1 is "name", Column 2 is
#' "mz" and column is "rt" (second).
#' @param ms2.data MS2 data, must be mgf, msp or mzXML format. For example, ms2.data = c("test.mgf", "test2.msp").
#' @param ms1.ms2.match.mz.tol MS1 peak and MS2 spectrum matching m/z tolerance. Default is 25 pm.
#' @param ms1.ms2.match.rt.tol MS1 peak and MS2 spectrum matching RT tolerance. Default is 10 s.
#' @param ms1.match.ppm Precursor match ppm tolerance.
#' @param ms2.match.ppm Fragment ion match ppm tolerance.
#' @param mz.ppm.thr Accurate mass tolerance for m/z error calculation.
#' @param ms2.match.tol MS2 match (MS2 similarity) tolerance.
#' @param fraction.weight The weight for matched fragments.
#' @param dp.forward.weight Forward dot product weight.
#' @param dp.reverse.weight Reverse dot product weight.
#' @param rt.match.tol RT match tolerance.
#' @param polarity The polarity of data, "positive"or "negative".
#' @param ce Collision energy. Please confirm the CE values in your database. Default is "all".
#' @param column "hilic" (HILIC column) or "rp" (reverse phase).
#' @param ms1.match.weight The weight of MS1 match for total score calculation.
#' @param rt.match.weight The weight of RT match for total score calculation.
#' @param ms2.match.weight The weight of MS2 match for total score calculation.
#' @param path Work directory.
#' @param total.score.tol Total score tolerance. The total score are refering to MS-DIAL.
#' @param candidate.num The number of candidate.
#' @param database MS2 database name or MS2 database.
#' @param threads Number of threads
#' @param silence.deprecated Silenc the deprecated information or not.
#' @return A metIdentifyClass object.
#' @importFrom crayon yellow green red bgRed
#' @export
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/metid.html}

metIdentify <-
  function(ms1.data,
           ##csv format
           ms2.data = NULL,
           ##only msp and mgf and mz(X)ML are supported
           ms1.ms2.match.mz.tol = 25,
           ms1.ms2.match.rt.tol = 10,
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
           path = ".",
           total.score.tol = 0.5,
           candidate.num = 3,
           database,
           threads = 3,
           silence.deprecated = FALSE) {
    if (!silence.deprecated) {
      message("`metIdentify()` is deprecated, use `identify_metabolites()`.")
    }
    
    ###Check data
    if (missing(database)) {
      stop("No database is provided.\n")
    }
    
    if (missing(ms1.data)) {
      stop("Please provide MS1 data name.\n")
    }
    
    ##parameter specification
    polarity <- match.arg(polarity)
    column <- match.arg(column)
    
    ##check ms1.file and ms2.file
    file <- dir(path)
    intermediate_path <- file.path(path, "intermediate_data")
    dir.create(intermediate_path, showWarnings = FALSE)
    
    if (!all(ms1.data %in% file)) {
      stop("MS1 data is not in the directory, please check it.\n")
    }
    
    if (!is.null(ms2.data)) {
      if (!all(ms2.data %in% file)) {
        stop("Some MS2 data are not in the directory, please check it.\n")
      }
    }
    
    if (!is(database, "databaseClass")) {
      if (!all(database %in% file)) {
        stop("Database is not in this directory, please check it.\n")
      }
    }
    
    #load MS2 database
    if (!is(database, "databaseClass")) {
      database.name <- database
      load(file.path(path, database.name))
      database <- get(database.name)
    } else{
      database.name = paste(database@database.info$Source,
                            database@database.info$Version,
                            sep = "_")
    }
    
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
    
    if (all(c("ms1.info", "ms2.info") %in% dir(intermediate_path))) {
      message(crayon::yellow("Use old data."))
      load(file.path(intermediate_path, "ms1.info"))
      load(file.path(intermediate_path, "ms2.info"))
    } else{
      ##read MS2 data
      # cat(crayon::green("Reading MS2 data...\n"))
      ms2.data.name <- ms2.data
      temp.ms2.type <-
        stringr::str_split(string = ms2.data.name, pattern = "\\.")[[1]]
      temp.ms2.type <- temp.ms2.type[length(temp.ms2.type)]
      
      if (temp.ms2.type %in% c("mzXML", "mzML")) {
        ms2.data <-
          masstools::read_mzxml(file = file.path(path, ms2.data.name),
                                threads = threads)
      } else{
        ms2.data <- lapply(ms2.data.name, function(temp.ms2.data) {
          temp.ms2.type <- stringr::str_split(string = temp.ms2.data, pattern = "\\.")[[1]]
          temp.ms2.type <-
            temp.ms2.type[length(temp.ms2.type)]
          if (!temp.ms2.type %in% c("mgf", "msp"))
            stop("We only support mgf or msp.\n")
          if (temp.ms2.type == "msp") {
            temp.ms2.data <- readMSP(file = file.path(path, temp.ms2.data))
          } else{
            temp.ms2.data <-
              masstools::read_mgf(file = file.path(path, temp.ms2.data))
          }
          temp.ms2.data
        })
        
        names(ms2.data) <- ms2.data.name
        ###prepare data for metidentification function
        message(crayon::green("Preparing MS2 data for identification..."))
        ms2.data <-
          mapply(
            FUN = function(temp.ms2.data, temp.ms2.data.name) {
              temp.ms2.data <- lapply(temp.ms2.data, function(x) {
                info <- x$info
                info <-
                  data.frame(
                    name = paste("mz", info[1], "rt", info[2], sep = ""),
                    "mz" = info[1],
                    "rt" = info[2],
                    "file" = temp.ms2.data.name,
                    stringsAsFactors = FALSE
                  )
                rownames(info) <- NULL
                x$info <- info
                x
              })
              temp.ms2.data
            },
            temp.ms2.data = ms2.data,
            temp.ms2.data.name = ms2.data.name
          )
        
        if (is(ms2.data, "matrix")) {
          ms2.data <- ms2.data[, 1]
        } else{
          ms2.data <- do.call(what = c, args = ms2.data)
        }
      }
      
      ms1.info <- lapply(ms2.data, function(x) {
        x[[1]]
      })
      
      ms2.info <- lapply(ms2.data, function(x) {
        x[[2]]
      })
      
      ms1.info <- do.call(what = rbind, args = ms1.info)
      ms1.info <- as.data.frame(ms1.info)
      rownames(ms1.info) <- NULL
      
      duplicated.name <-
        unique(ms1.info$name[duplicated(ms1.info$name)])
      
      if (length(duplicated.name) > 0) {
        lapply(duplicated.name, function(x) {
          ms1.info$name[which(ms1.info$name == x)] <-
            paste(x, seq_len(sum(ms1.info$name == x)), sep = "_")
        })
      }
      
      names(ms2.info) <- ms1.info$name
      
      ##save intermediate data
      save(ms1.info,
           file = file.path(intermediate_path, "ms1.info"),
           compress = "xz")
      
      save(ms2.info,
           file = file.path(intermediate_path, "ms2.info"),
           compress = "xz")
      message(crayon::red("OK."))
    }
    
    if (!missing(ms1.data)) {
      message(crayon::green("Matching peak table with MS2 spectrum..."))
      ##check ms1 data format
      if (length(grep("csv", ms1.data)) == 0) {
        stop("Only support csv format ms1 data.\n")
      }
      
      ms1.data <-
        readr::read_csv(file = file.path(path, ms1.data),
                        col_types = readr::cols())
      
      ##check for the ms1 data
      if (ncol(ms1.data) < 3) {
        stop(
          "MS1 data should have there columns.
           See here: \n https://tidymass.github.io/metid/articles/metabolite_annotation_using_MS1.html"
        )
      }
      
      if (colnames(ms1.data)[1] != "name" |
          colnames(ms1.data)[2] != "mz" |
          colnames(ms1.data)[3] != "rt") {
        stop("The columns should be name, mz and rt, respectively.\n")
      }
      
      colnames(ms1.data)[seq_len(3)] <- c("name", "mz", "rt")
      match.result <-
        masstools::mz_rt_match(
          data1 = ms1.data[, c(2, 3)],
          data2 = ms1.info[, c(2, 3)],
          mz.tol = ms1.ms2.match.mz.tol,
          rt.tol = ms1.ms2.match.rt.tol,
          rt.error.type = "abs"
        )
      if (is.null(match.result))
        return("No peaks are matched with MS2 spectra.\n")
      if (nrow(match.result) == 0)
        return("No peaks are matched with MS2 spectra.\n")
      message(crayon::green(
        length(unique(match.result[, 1])),
        "out of",
        nrow(ms1.data),
        "peaks have MS2 spectra."
      ))
      
      ###if one peak matches multiple peaks, select the more reliable MS2 spectrum
      message(crayon::green("Selecting the most intense MS2 spectrum for each peak..."))
      temp.idx <- unique(match.result[, 1])
      
      match.result <- lapply(temp.idx, function(idx) {
        idx2 <- match.result[which(match.result[, 1] == idx), 2]
        if (length(idx2) == 1) {
          return(c(idx, idx2))
        } else{
          temp.ms2.info <- ms2.info[idx2]
          return(c(idx, idx2[which.max(unlist(lapply(temp.ms2.info, function(y) {
            y <- y[order(y[, 2], decreasing = TRUE), , drop = FALSE]
            if (nrow(y) > 5)
              y <- y[seq_len(5), ]
            sum(y[, 2])
          })))]))
        }
      })
      
      match.result <- do.call(rbind, match.result)
      match.result <- as.data.frame(match.result)
      colnames(match.result) <- c("Index1", "Index2")
      match.result <- data.frame(match.result,
                                 ms1.data$name[match.result$Index1],
                                 ms1.info$name[match.result$Index2],
                                 stringsAsFactors = FALSE)
      colnames(match.result) <-
        c("Index1.ms1.data",
          "Index.ms2.spectra",
          "MS1.peak.name",
          "MS2.spectra.name")
      ms1.info <-
        ms1.info[unique(match.result[, 2]), , drop = FALSE]
      
      ms2.info <- ms2.info[unique(match.result[, 2])]
      
      match.result$Index.ms2.spectra <-
        match(match.result$MS2.spectra.name, ms1.info$name)
      
      save(
        match.result,
        file = file.path(intermediate_path, "match.result"),
        compress = "xz"
      )
      message(crayon::red("OK."))
    } else{
      stop("Please provide MS1 data name.\n")
    }
    
    ms2Matchresult <-
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
        dp.reverse.weight = dp.reverse.weight
      )
    
    return.result <- new(
      Class = "metIdentifyClass",
      ms1.data = ms1.data,
      ms1.info = ms1.info,
      ms2.info = ms2.info,
      identification.result = ms2Matchresult,
      match.result = match.result,
      adduct.table = adduct.table,
      ms1.ms2.match.mz.tol = ms1.ms2.match.mz.tol,
      ms1.ms2.match.rt.tol = ms1.ms2.match.rt.tol,
      ms1.match.ppm = ms1.match.ppm,
      ms2.match.ppm = ms2.match.ppm,
      ms2.match.tol = ms2.match.tol,
      rt.match.tol = rt.match.tol,
      polarity = polarity,
      ce = paste(ce, collapse = ";"),
      column = column,
      ms1.match.weight = ms1.match.weight,
      rt.match.weight = rt.match.weight,
      ms2.match.weight = ms2.match.weight,
      path = path,
      total.score.tol = total.score.tol,
      candidate.num = candidate.num,
      database = database.name,
      threads = threads,
      version = "1.0.0"
    )
    message(crayon::bgRed("All done."))
    return(return.result)
  }




#' ---------------------------------------------------------------------------
#' @title readMSP_MoNA
#' @description Read MSP data from MoNA.
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @param file The vector of names of ms2 files. MS2 file must be msp. The msp data must from MoNA.
#' @return Return ms2 data. This is a list.
#' @export

readMSP_MoNA <-
  function(file) {
    message("`readMSP_MoNA()` is deprecated, use `read_msp_mona()`.")
    message(crayon::green("Reading MSP data..."))
    msp.data <- readr::read_lines(file)
    n.null <- which(msp.data == '')
    temp.idx1 <- c(1, n.null[-length(n.null)])
    temp.idx2 <- n.null - 1
    
    temp.idx <- data.frame(temp.idx1, temp.idx2, stringsAsFactors = FALSE)
    
    rm(list = c("temp.idx1", "temp.idx2"))
    gc(verbose = FALSE)
    
    temp.idx <- apply(temp.idx, 1, list)
    
    temp.idx <- lapply(temp.idx, unlist)
    temp.idx <- temp.idx[which(unlist(lapply(temp.idx, function(x) {
      x[1] != x[2]
    })))]
    
    pbapply::pboptions(style = 1)
    # fix bug
    message(crayon::yellow("Transforming..."))
    info.spec <- pbapply::pblapply(temp.idx, function(idx) {
      if (idx[1] == idx[2])
        return(NULL)
      temp.msp.data <- msp.data[idx[1]:idx[2]]
      temp.msp.data <- temp.msp.data[temp.msp.data != ""]
      info.idx <- grep("[A-Za-z]", temp.msp.data)
      temp.info <- temp.msp.data[info.idx]
      temp.info <-
        stringr::str_split(string = temp.info,
                           pattern = ":",
                           n = 2)
      temp.info <- do.call(rbind, temp.info)
      temp.info <- data.frame(temp.info, stringsAsFactors = FALSE)
      temp.info[, 2] <-
        stringr::str_trim(temp.info[, 2], side = "both")
      temp.info[, 1] <-
        stringr::str_trim(temp.info[, 1], side = "both")
      colnames(temp.info) <- rownames(temp.info) <- NULL
      #combine synons
      if (length(grep("Synon", temp.info[, 1])) > 1) {
        Synon <-
          stringr::str_c(temp.info[stringr::str_which(string = temp.info[, 1], pattern = "Synon"), 2, drop = TRUE], collapse = ";")
        temp.info[which(temp.info[, 1] == "Synon")[1], 2] <- Synon
        temp.info <- temp.info[!duplicated(temp.info[, 1]), ]
        
      }
      
      rownames(temp.info) <- temp.info[, 1]
      temp.info <- temp.info[, -1, drop = FALSE]
      temp.spec <- temp.msp.data[-info.idx]
      
      if (length(temp.spec) != 0) {
        if (length(grep(" ", temp.spec[1])) == 1) {
          temp.spec <- strsplit(temp.spec, split = ' ')
        }
        
        if (length(grep("\t", temp.spec[1])) == 1) {
          temp.spec <- strsplit(x = temp.spec, split = "\t")
        }
        
        temp.spec <- do.call(rbind, temp.spec)
        temp.spec <- data.frame(temp.spec, stringsAsFactors = FALSE)
        colnames(temp.spec) <- c('mz', 'intensity')
        rownames(temp.spec) <- NULL
        temp.spec$mz <- as.numeric(as.character(temp.spec$mz))
        temp.spec$intensity <- as.numeric(temp.spec$intensity)
        temp.spec <- temp.spec[temp.spec$intensity != 0, ]
      } else{
        temp.spec <- NULL
      }
      
      list('info' = temp.info, 'spec' = temp.spec)
    })
    
    rm(list = c("msp.data", "temp.idx"))
    gc()
    ##remove NULL
    info.spec <- info.spec[!unlist(lapply(info.spec, is.null))]
    
    remove.idx <-
      which(unlist(lapply(info.spec, function(x)
        is.null(x[[2]]))))
    if (length(remove.idx) > 0) {
      info.spec <- info.spec[-remove.idx]
    }
    rm(list = c("remove.idx"))
    gc()
    info.spec <- info.spec
  }




#######------------------------------------------------------------------------
#' @title readMGF
#' @description Read MGF data.
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @param file The vector of names of ms2 files. MS2 file must be mgf.
#' @return Return ms2 data. This is a list.
#' @export

readMGF <-
  function(file) {
    message("`readMGF()` is deprecated, use `read_mgf()`.")
    pbapply::pboptions(style = 1)
    message(crayon::green("Reading MS2 data..."))
    # mgf.data.list <- pbapply::pblapply(file, list_mgf)
    ms2 <- pbapply::pblapply(file, function(mgf.data) {
      mgf.data <- list_mgf(mgf.data)
      # nl.spec <- grep('^\\d', mgf.data)
      nl.spec <-
        lapply(mgf.data, function(x)
          grep('^\\d', x))
      info.mz <-
        lapply(mgf.data, function(x)
          grep('^PEPMASS', x, value = TRUE))
      info.rt <-
        lapply(mgf.data, function(x)
          grep('^RTINSECONDS', x, value = TRUE))
      
      info.mz <- unlist(info.mz)
      #for orbitrap data, the intensity of precursor ion should be removed
      info.mz <-
        unlist(lapply(strsplit(x = info.mz, split = " "), function(x)
          x[1]))
      info.mz <-
        as.numeric(gsub(pattern = "\\w+=", "", info.mz))
      info.rt <- unlist(info.rt)
      info.rt <-
        as.numeric(gsub(pattern = "\\w+=", "", info.rt))
      
      if (length(mgf.data) == 1) {
        spec <- mapply(function(x, y) {
          temp <- do.call(rbind, strsplit(x[y], split = " "))
          list(temp)
        }, x = mgf.data, y = nl.spec)
      } else{
        spec <- mapply(function(x, y) {
          do.call(rbind, strsplit(x[y], split = " "))
        }, x = mgf.data, y = nl.spec)
      }
      
      spec <- lapply(spec, function(x) {
        temp <- cbind(as.numeric(x[, 1]), as.numeric(x[, 2]))
        temp <- matrix(temp, ncol = 2)
        # if(nrow(temp) > 0) temp <- temp[temp[,2] >= max(temp[,2])*0.01,]
        temp <- matrix(temp, ncol = 2)
        colnames(temp) <- c("mz", "intensity")
        temp
      })
      
      ms2 <- mapply(function(x, y, z) {
        info <- c(y, z)
        names(info) <- c("mz", "rt")
        spectrum <- as.matrix(x)
        temp <- list(info, spectrum)
        names(temp) <- c("info", "spec")
        list(temp)
      },
      x = spec,
      y = info.mz,
      z = info.rt)
      
      ms2
      
    })
    
    
    spec.info <- ms2[[1]]
    if (length(ms2) > 1) {
      for (i in 2:length(ms2)) {
        spec.info <- c(spec.info, ms2[[i]])
      }
    }
    
    remove.idx <-
      which(unlist(lapply(spec.info, function(x)
        nrow(x[[2]]))) == 0)
    if (length(remove.idx) != 0)
      spec.info <- spec.info[-remove.idx]
    # ##remove noise
    # cat("\n")
    # cat("Remove noise of MS/MS spectra...\n")
    # spec.info <- pbapply::pblapply(spec.info, function(x){
    #   temp.spec <- x[[2]]
    #   temp.spec <- remove_noise(temp.spec)
    #   x[[2]] <- temp.spec
    #   x
    # })
    
    spec.info <- spec.info
  }



#---------------------------------------------------------------------------
#' @title readMSP
#' @description Read MSP data.
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @param file The vector of names of ms2 files. MS2 file must be msp format.
#' @param silence.deprecated Silenc the deprecated information or not.
#' @return Return ms2 data. This is a list.
#' @export

readMSP <-
  function(file, silence.deprecated = TRUE) {
    if (!silence.deprecated) {
      message("`readMSP()` is deprecated, use `read_msp()`.")
    }
    
    msp.data <- readLines(file)
    if (length(grep("BEGIN IONS", msp.data)) > 0) {
      msp.data <- msp.data[msp.data != ""]
      temp.idx1 <- grep("BEGIN IONS", msp.data)
      temp.idx2 <- grep("END IONS", msp.data)
      if (length(temp.idx2) < length(temp.idx1)) {
        temp.idx2 <- c(temp.idx2, length(msp.data))
      }
      
      temp.idx <-
        purrr::map2(.x = temp.idx1, temp.idx2, function(x, y) {
          c(x + 1, y - 1)
        })
      
      ms2_spec <- purrr::map(
        .x = temp.idx,
        .f = function(x) {
          temp_spec <- msp.data[x[1]:x[2]]
          temp_spec <- temp_spec
          spec_info <-
            temp_spec[stringr::str_detect(temp_spec, "[A-Za-z]")]
          spec <-
            temp_spec[!stringr::str_detect(temp_spec, "[A-Za-z]")]
          spec_info <- stringr::str_split(spec_info, "\\=") %>%
            do.call(rbind, .)
          mz <-
            as.numeric(spec_info[grep("MASS|MZ", spec_info[, 1]), 2])
          rt <-
            as.numeric(spec_info[grep("RT|RETETION", spec_info[, 1]), 2])
          spec_info <- c(mz = mz, rt = rt)
          
          spec <- purrr::map(
            .x = spec,
            .f = function(x) {
              stringr::str_split(x, " ")[[1]] %>% as.numeric()
            }
          ) %>%
            do.call(rbind, .)
          
          spec <-
            spec %>% as.matrix()
          
          colnames(spec) <- c("mz", "intensity")
          
          spec <- list(info = spec_info, spec = spec)
          spec
        }
      )
      return(ms2_spec)
    } else{
      # n.tot <- length(msp.data)
      n.null <- which(msp.data == '')
      
      temp.idx1 <- c(1, n.null[-length(n.null)])
      temp.idx2 <- n.null - 1
      
      temp.idx <- data.frame(temp.idx1, temp.idx2, stringsAsFactors = FALSE)
      temp.idx <- apply(temp.idx, 1, list)
      
      temp.idx <- lapply(temp.idx, unlist)
      
      # n.spec <- which(grepl('^\\d', msp.data))
      # n.info <- seq(n.tot)[-c(n.spec, n.null)]
      
      pbapply::pboptions(style = 1)
      info.spec <- pbapply::pblapply(temp.idx, function(idx) {
        temp.msp.data <- msp.data[idx[1]:idx[2]]
        
        temp.msp.data <- temp.msp.data[temp.msp.data != ""]
        info.idx <- grep("[A-Za-z]", temp.msp.data)
        temp.info <- temp.msp.data[info.idx]
        temp.info <- strsplit(temp.info, split = ":")
        temp.info <- do.call(rbind, temp.info)
        temp.info <- data.frame(temp.info, stringsAsFactors = FALSE)
        temp.info[, 2] <- stringr::str_trim(temp.info[, 2])
        colnames(temp.info) <- rownames(temp.info) <- NULL
        rownames(temp.info) <- temp.info[, 1]
        temp.info <- temp.info[, -1, drop = FALSE]
        
        temp.spec <- temp.msp.data[-info.idx]
        
        if (length(temp.spec) != 0) {
          if (length(grep(" ", temp.spec[1])) == 1) {
            temp.spec <- strsplit(temp.spec, split = ' ')
          }
          
          if (length(grep("\t", temp.spec[1])) == 1) {
            temp.spec <- strsplit(x = temp.spec, split = "\t")
          }
          
          temp.spec <- do.call(rbind, temp.spec)
          temp.spec <- data.frame(temp.spec, stringsAsFactors = FALSE)
          colnames(temp.spec) <- c('mz', 'intensity')
          rownames(temp.spec) <- NULL
          temp.spec$mz <- as.numeric(as.character(temp.spec$mz))
          temp.spec$intensity <- as.numeric(temp.spec$intensity)
          temp.spec <- temp.spec[temp.spec$intensity != 0, ]
        } else{
          temp.spec <- NULL
        }
        
        list('info' = temp.info, 'spec' = temp.spec)
      })
      
      mz.idx <- grep("[Mm][Zz]", rownames(info.spec[[1]][[1]]))
      rt.idx <-
        grep("Time|TIME|time|RT|rt|Rt", rownames(info.spec[[1]][[1]]))
      
      ##fix bug in msp data from metAnalyzer
      if (length(rt.idx) == 0) {
        message(crayon::yellow("The msp data are from MetAnalyzer software."))
        rt.idx <-
          grep("NAME|Name|name", rownames(info.spec[[1]][[1]]))
        ##rt.idx is the name of peak
        info.spec <- lapply(info.spec, function(x) {
          info <- x[[1]]
          mz <- as.numeric(info[mz.idx, 1])
          rt <- as.character(info[rt.idx, 1])
          info <- c(mz, rt)
          names(info) <- c("mz", "rt")
          x[[1]] <- info
          x
        })
      } else{
        info.spec <- lapply(info.spec, function(x) {
          info <- x[[1]]
          mz <- as.numeric(info[mz.idx, 1])
          rt <- as.numeric(info[rt.idx, 1])
          info <- c(mz, rt)
          names(info) <- c("mz", "rt")
          x[[1]] <- info
          x
        })
      }
      
    }
    
    remove.idx <-
      which(unlist(lapply(info.spec, function(x)
        is.null(x[[2]]))))
    if (length(remove.idx) > 0) {
      info.spec <- info.spec[-remove.idx]
    }
    
    info.spec <- info.spec
  }


#' @title readMZXML
#' @description Read mzXML data.
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @param file The vector of names of ms2 files. MS2 file must be mzXML or mzML.
#' @param threads Thread number
#' @return Return ms2 data. This is a list.
#' @export

readMZXML <-
  function(file, threads = 3) {
    message("`readMZXML()` is deprecated, use `read_mzxml()`.")
    # pbapply::pboptions(style = 1)
    message(crayon::green("Reading MS2 data..."))
    # mzxml.data.list <- pbapply::pblapply(file, list_mgf)
    ms2 <-
      MSnbase::readMSData(files = file,
                          msLevel. = 2,
                          mode = "onDisk")
    message(crayon::green("Processing..."))
    
    new.ms2 <- ProtGenerics::spectra(object = ms2)
    rm(list = c("ms2"))
    #
    temp.fun <- function(idx, ms2) {
      temp.ms2 <- ms2[[idx]]
      rm(list = c("ms2"))
      info <-
        data.frame(
          name = paste("mz", temp.ms2@precursorMz, "rt", temp.ms2@rt, sep = ""),
          "mz" = temp.ms2@precursorMz,
          "rt" = temp.ms2@rt,
          "file" = file[temp.ms2@fromFile],
          stringsAsFactors = FALSE
        )
      duplicated.name <-
        unique(info$name[duplicated(info$name)])
      if (length(duplicated.name) > 0) {
        lapply(duplicated.name, function(x) {
          info$name[which(info$name == x)] <-
            paste(x, seq_len(sum(info$name == x)), sep = "_")
        })
      }
      
      rownames(info) <- NULL
      spec <- data.frame(
        "mz" = temp.ms2@mz,
        "intensity" = temp.ms2@intensity,
        stringsAsFactors = FALSE
      )
      list(info = info, spec = spec)
    }
    
    new.ms2 <-
      BiocParallel::bplapply(
        X = seq_along(new.ms2),
        FUN = temp.fun,
        BPPARAM = BiocParallel::MulticoreParam(workers = threads, progressbar = TRUE),
        ms2 = new.ms2
      )
    
    new.ms2 <- new.ms2
  }



#' @title Identify single peak based on database.
#' @description We can use this function to identify single peak, you can just provide m/z or rt, or you can also provide MS2 spectrum for this peak.
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @param ms1.mz m/z value of the peaks
#' @param ms1.rt rt value of the peaks
#' @param ms2 MS2 spectra of the peaks. It must be a two column data frame, the
#' first column is m/z and the second column is the intensity.
#' @param ms1.match.ppm Precursor match ppm tolerance.
#' @param ms2.match.ppm Fragment ion match ppm tolerance.
#' @param mz.ppm.thr Accurate mass tolerance for m/z error calculation.
#' @param ms2.match.tol MS2 match (MS2 similarity) tolerance.
#' @param fraction.weight The weight for matched fragments.
#' @param dp.forward.weight Forward dot product weight.
#' @param dp.reverse.weight Reverse dot product weight.
#' @param rt.match.tol RT match tolerance.
#' @param polarity The polarity of data, "positive"or "negative".
#' @param ce Collision energy. Please confirm the CE values in your database. Default is "all".
#' @param column "hilic" (HILIC column) or "rp" (reverse phase).
#' @param ms1.match.weight The weight of MS1 match for total score calculation.
#' @param rt.match.weight The weight of RT match for total score calculation.
#' @param ms2.match.weight The weight of MS2 match for total score calculation.
#' @param path Work directory.
#' @param total.score.tol Total score tolerance. The total score are refering to MS-DIAL.
#' @param candidate.num The number of candidate.
#' @param database MS2 database name.
#' @param threads Number of threads
#' @return A metIdentifyClass object.
#' @export

identify_single_peak <-
  function(ms1.mz,
           ms1.rt,
           ms2,
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
           path = ".",
           total.score.tol = 0.5,
           candidate.num = 3,
           database,
           threads = 3) {
    message("This function is deprecated")
    ###Check data
    if (missing(database)) {
      stop("No database is provided.\n")
    }
    
    if (missing(ms1.mz) | missing(ms1.rt) | missing(ms2)) {
      stop("ms1.mz, ms1.rt or ms2 is not provided.\n")
    }
    
    ##parameter specification
    polarity <- match.arg(polarity)
    column <- match.arg(column)
    
    #load MS2 database
    database.name <- database
    load(file.path(path, database.name))
    database <- get(database.name)
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
    
    
    name <- paste(paste("mz", ms1.mz, sep = ""), paste("rt", ms1.rt, sep = ""), sep = "")
    
    file <- "test"
    
    ms1.info <- data.frame(name,
                           mz = ms1.mz,
                           rt = ms1.rt,
                           file,
                           stringsAsFactors = FALSE)
    
    if (is.matrix(ms2) | is.data.frame(ms2)) {
      ms2.info <- list(ms2)
    } else{
      ms2.info <- ms2
    }
    
    names(ms2.info) <- ms1.info$name
    
    if (c(length(ms1.mz), length(ms1.rt), length(ms2.info)) %>%
        unique() %>%
        length() != 1) {
      stop("Length of ms1.mz, ms1.rt and ms2 must be same.\n")
    }
    
    ms2Matchresult <-
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
        dp.reverse.weight = dp.reverse.weight
      )
    
    
    ms1.data <- ms1.info %>%
      dplyr::select(-file)
    
    match.result <- data.frame(
      Index1.ms1.data = seq_len(nrow(ms1.info)),
      Index.ms2.spectra = seq_len(nrow(ms1.info)),
      MS1.peak.name = ms1.info$name,
      MS2.spectra.name = ms1.info$name,
      stringsAsFactors = FALSE
    )
    
    return.result <- new(
      Class = "metIdentifyClass",
      ms1.data = ms1.data,
      ms1.info = ms1.info,
      ms2.info = ms2.info,
      identification.result = ms2Matchresult,
      match.result = match.result,
      adduct.table = adduct.table,
      ms1.ms2.match.mz.tol = 0,
      ms1.ms2.match.rt.tol = 0,
      ms1.match.ppm = ms1.match.ppm,
      ms2.match.ppm = ms2.match.ppm,
      ms2.match.tol = ms2.match.tol,
      rt.match.tol = rt.match.tol,
      polarity = polarity,
      ce = paste(ce, collapse = ";"),
      column = column,
      ms1.match.weight = ms1.match.weight,
      rt.match.weight = rt.match.weight,
      ms2.match.weight = ms2.match.weight,
      path = path,
      total.score.tol = total.score.tol,
      candidate.num = candidate.num,
      database = database.name,
      threads = threads,
      version = "1.0.0"
    )
    message(crayon::bgRed("All done."))
    return(return.result)
  }
