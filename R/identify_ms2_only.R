
#' @title Identify metabolites based on MS/MS database
#' @description Identify metabolites based on MS2 data.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param ms2.data MS2 data, must be mgf, msp or mzXML format. For example, ms2.data = c("test.mgf", "test2.msp").
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


identify_ms2_only = function(ms2.data,
                             ##only msp and mgf and mz(X)ML are supported
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

  ###Check data
  if (missing(database)) {
    stop("No database is provided.\n")
  }
  
  if (missing(ms2.data)) {
    stop("No ms2.data is provided.\n")
  }
  
  ##parameter specification
  polarity <- match.arg(polarity)
  column <- match.arg(column)
  
  ##check ms1.file and ms2.file
  file <- dir(path)
  intermediate_path <- file.path(path, "intermediate_data")
  dir.create(intermediate_path, showWarnings = FALSE)
  
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
    message(
      crayon::yellow(
        "No RT information in database.\nThe weight of RT have been set as 0."
      )
    )
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
    message(crayon::red("OK."))
  }
  
  if (!missing(ms1.data)) {
    message(crayon::green("Matching peak table with MS2 spectrum..."))
    ##check ms1 data format
    
    ms1.data <-
      ms1.info
    
    ##check for the ms1 data
    if (ncol(ms1.data) < 3) {
      stop(
        "MS1 data should have there columns. See here: \n https://tidymass.github.io/metid/articles/metabolite_annotation_using_MS1.html"
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
        mz.tol = 2,
        rt.tol = 2,
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
    
    ###if one peak matches multiple peaks, select the more relibale MS2 spectrum
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