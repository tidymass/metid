#' @title Identify peaks based on MS1 database
#' @description Identify peaks based on MS1 database.
#' \lifecycle{deprecated}
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
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
#' https://tidymass.github.io/metid/articles/metid.html


mzIdentify =
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
      cat(crayon::yellow(
        "`mzIdentify()` is deprecated, use `identify_metabolites()`."
      ))
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
    
    if(class(database) != "databaseClass"){
      if (!all(database %in% file)) {
        stop("Database is not in this directory, please check it.\n")
      }  
    }
    
    #load database
    if(class(database) != "databaseClass"){
      database.name <- database
      load(file.path(path, database.name))
      database <- get(database.name) 
    }else{
      database.name = paste(database@database.info$Source, 
                            database@database.info$Version, sep = "_")
    }
    # if(class(database) != "data.frame"){
    #   stop("The database must be HMDB.metabolite.data provided.\n")
    # }
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
    
    colnames(ms1.data)[1:3] <- c("name", "mz", "rt")
    
    
    if (rt.match.tol > 10000) {
      cat(
        crayon::yellow(
          "You set rt.match.tol > 10,000, so RT will not be used for matching.\n"
        )
      )
    } else{
      cat(
        crayon::yellow(
          "You set rt.match.tol < 10,000, so if the compounds have RT,  RTs will be used for matching\n"
        )
      )
    }
    
    # ##debug
    # for(i in 1:nrow(ms1.data)){
    #   cat(i, " ")
    #   result = temp.fun(idx = i, ms1.data = ms1.data,
    #                     ms1.match.ppm = ms1.match.ppm,
    #                     rt.match.tol = rt.match.tol,
    #                     database = database,
    #                     adduct.table = adduct.table, candidate.num = candidate.num)
    # }
    
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
        match.idx <- match.idx[1:candidate.num, , drop = FALSE]
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
    
    
    if(tinytools::get_os() == "windows"){
      bpparam = BiocParallel::SnowParam(workers = threads, 
                                        progressbar = TRUE)
    }else{
      bpparam = BiocParallel::MulticoreParam(workers = threads, 
                                        progressbar = TRUE)
    }
    
    match.result <-
      BiocParallel::bplapply(
        1:nrow(ms1.data),
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
    
    if (class(database) == "databaseClass") {
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
    cat(crayon::bgRed("All done.\n"))
    return(return.result)
  }


#------------------------------------------------------------------------------
#' @title Get parameters from a mzIdentifyClass object
#' @description Get parameters from a mzIdentifyClass object.
#' \lifecycle{deprecated}
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param object A mzIdentifyClass object.
#' @return A data frame contains all the parameters of this metIdentifiyClass object.
#' @export
#' @seealso The example and demo data of this function can be found
#' https://tidymass.github.io/metid/articles/metid.html

getParams2 = function(object){
  cat(crayon::yellow(
    "`getParams2()` is deprecated, use `get_parameters_metid()`."
  ))
  
  if (class(object) != "mzIdentifyClass")
    stop("Only for mzIdentifyClass\n")
  data.frame(
    "Parameter" = c(
      "ms1.match.ppm",
      "polarity",
      "column",
      "path",
      "candidate.num",
      "database",
      "threads"
    ),
    "Meaning" = c(
      "MS1 match tolerance (ppm)",
      "Polarity",
      "Column",
      "Work directory",
      "Candidate number",
      "MS2 database",
      "Thread number"
    ),
    "Value" = c(
      object@ms1.match.ppm,
      object@polarity,
      object@column,
      object@path,
      object@candidate.num,
      object@database,
      object@threads
    ),
    stringsAsFactors = FALSE
  )
}





##------------------------------------------------------------------------------
#' @title Get identification table from a mzIdentifyClass object
#' @description Get identification table from a mzIdentifyClass object.
#' \lifecycle{deprecated}
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param object mzIdentifyClass object.
#' @param candidate.num The number of candidates.
#' @param type The type of identification table.
#' @param silence.deprecated Silenc the deprecated information or not.
#' @return A identification table (data.frame).
#' @export
#' @importFrom magrittr %>%
#' @seealso The example and demo data of this function can be found
#' https://tidymass.github.io/metid/articles/metid.html

# object = result_rp_pos25[[8]]

getIdentificationTable2 = function(object,
                                   candidate.num = 3,
                                   type = c("old", "new"),
                                   silence.deprecated = FALSE) {
  if (!silence.deprecated) {
    cat(
      crayon::yellow(
        "`getIdentificationTable2()` is deprecated, use `get_identification_table()`."
      )
    )
  }
  
  if (class(object) != "mzIdentifyClass" &
      class(object) != "metIdentifyClass") {
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
        x <- x[1:candidate.num, , drop = FALSE]
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
    
    identification.table$Identification[match(names(identification.result),
                                              peak.table$name)] <-
      unlist(identification.result)
    
    identification.table$Candidate.number <-
      sapply(identification.table$Identification, function(x) {
        if (is.na(x))
          return(0)
        return(length(stringr::str_split(
          string = x, pattern = "\\{\\}"
        )[[1]]))
      })
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
              as.data.frame(matrix(c(x, rep(NA, 14)), nrow = 1),
                            stringsAsFactors = FALSE)
            
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
      dplyr::select(colnames(object@ms1.data), 
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
