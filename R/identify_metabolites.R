
#' @title Identify metabolites based on MS1 or MS/MS database
#' @description Identify metabolites based on MS1 or MS/MS database.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
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

identify_metabolites = function(
  ms1.data,
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
  threads = 3
) {

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

  if(!is(database, "databaseClass")){
    if (!all(database %in% file)) {
      stop("Database is not in this directory, please check it.\n")
    }  
  }
  
  if (is.null(ms2.data)) {
    message(crayon::yellow("You don't provide MS2 data, so only use mz and/or RT for matching."))
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


#' @title Identify peaks based on MS1 database
#' @description Identify peaks based on MS1 database.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
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
      message(crayon::yellow(
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
    
    if(!is(database, "databaseClass")){
      if (!all(database %in% file)) {
        stop("Database is not in this directory, please check it.\n")
      }  
    }
    
    #load database
    if(!is(database, "databaseClass")){
      database.name <- database
      load(file.path(path, database.name))
      database <- get(database.name) 
    }else{
      database.name = paste(database@database.info$Source, 
                            database@database.info$Version, sep = "_")
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
      message(
        crayon::yellow(
          "You set rt.match.tol > 10,000, so RT will not be used for matching."
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
    
    
    if(masstools::get_os() == "windows"){
      bpparam = BiocParallel::SnowParam(workers = threads, 
                                        progressbar = TRUE)
    }else{
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




#' @title Identify metabolites based on MS/MS database.
#' @description Identify metabolites based on MS/MS database.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
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

metIdentify = function(
  ms1.data,
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
  silence.deprecated = FALSE
) {
  
  if(!silence.deprecated){
    message(crayon::yellow(
      "`metIdentify()` is deprecated, use `identify_metabolites()`."
    ))  
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
  
  if(!is(database, "databaseClass")){
    if (!all(database %in% file)) {
      stop("Database is not in this directory, please check it.\n")
    }  
  }
  
  #load MS2 database
  if(!is(database, "databaseClass")){
    database.name <- database
    load(file.path(path, database.name))
    database <- get(database.name) 
  }else{
    database.name = paste(database@database.info$Source, 
                          database@database.info$Version, sep = "_")
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
    message(crayon::yellow("No RT information in database.\nThe weight of RT have been set as 0."))
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
      stringr::str_split(string = ms2.data.name,
                         pattern = "\\.")[[1]]
    temp.ms2.type <- temp.ms2.type[length(temp.ms2.type)]
    
    if (temp.ms2.type %in% c("mzXML", "mzML")) {
      ms2.data <-
        masstools::read_mzxml(file = file.path(path, ms2.data.name),
                   threads = threads)
    } else{
      ms2.data <- lapply(ms2.data.name, function(temp.ms2.data) {
        temp.ms2.type <- stringr::str_split(string = temp.ms2.data,
                                            pattern = "\\.")[[1]]
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
    if(length(grep("csv", ms1.data)) == 0){
      stop("Only support csv format ms1 data.\n")
    }
    
    ms1.data <-
      readr::read_csv(file = file.path(path, ms1.data),
                      col_types = readr::cols())
    
    ##check for the ms1 data
    if(ncol(ms1.data) < 3){
      stop("MS1 data should have there columns. 
           See here: \n https://tidymass.github.io/metid/articles/metabolite_annotation_using_MS1.html")
    }
    
    if(colnames(ms1.data)[1] != "name" | 
       colnames(ms1.data)[2] != "mz" | 
       colnames(ms1.data)[3] != "rt" 
    ){
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
            y <- y[seq_len(5),]
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
    
    save(match.result,
         file = file.path(intermediate_path, "match.result"),
         compress = "xz")
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




