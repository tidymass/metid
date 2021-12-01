#' @title Identify metabolites using multiple databases one time
#' @description Identify metabolites using multiple databases one time.
#' \lifecycle{maturing}
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param ms1.data The name of ms1 peak table (csv format). Column 1 is "name", column 2 is
#' "mz" and column 3 is "rt" (second).
#' @param ms2.data MS2 data, must be mgf, msp or mzXML format. For example, ms2.data = c("test.mgf", "test2.msp").
#' @param parameter.list A list contains paramters for each processing.
#' The parameter must get using metIdentifyParam or mzIdentifyParam.
#' @param path Work directory.
#' @return A list containing mzIdentifyClass object.
#' @export
#' @seealso The example and demo data of this function can be found
#' https://tidymass.github.io/metid/articles/multiple_databases.html
#' @examples 
#' \dontrun{
#' ##creat a folder nameed as example
#' path <- file.path(".", "example")
#' dir.create(path = path, showWarnings = FALSE)
#' 
#' ##get MS1 peak table from metid
#' ms1_peak <- system.file("ms1_peak", package = "metid")
#' file.copy(
#'   from = file.path(ms1_peak, "ms1.peak.table.csv"),
#'   to = path,
#'   overwrite = TRUE,
#'   recursive = TRUE
#' )
#' 
#' ##get MS2 data from metid
#' ms2_data <- system.file("ms2_data", package = "metid")
#' file.copy(
#'   from = file.path(ms2_data, "QC1_MSMS_NCE25.mgf"),
#'   to = path,
#'   overwrite = TRUE,
#'   recursive = TRUE
#' )
#' 
#' ##get databases from metid
#' database <- system.file("ms2_database", package = "metid")
#' 
#' file.copy(
#'   from = file.path(
#'     database,
#'     c(
#'       "msDatabase_rplc0.0.2",
#'       "orbitrapDatabase0.0.1",
#'       "hmdbMS1Database0.0.1"
#'     )
#'   ),
#'   to = path,
#'   overwrite = TRUE,
#'   recursive = TRUE
#' )
#' param1 <-
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
#' 
#' param2 <- identify_metabolites_params(
#'   ms1.match.ppm = 15,
#'   rt.match.tol = 15,
#'   polarity = "positive",
#'   ce = "all",
#'   column = "rp",
#'   total.score.tol = 0.5,
#'   candidate.num = 3,
#'   threads = 3,
#'   database = "orbitrapDatabase0.0.1"
#' )
#' 
#' param3 <- identify_metabolites_params(
#'   ms1.match.ppm = 15,
#'   rt.match.tol = 15,
#'   polarity = "positive",
#'   ce = "all",
#'   column = "rp",
#'   total.score.tol = 0.5,
#'   candidate.num = 3,
#'   threads = 3,
#'   database = "hmdbMS1Database0.0.1"
#' )
#' result <- metIdentify4all(
#' ms1.data = "ms1.peak.table.csv",
#' ms2.data = "QC1_MSMS_NCE25.mgf",
#' parameter.list = c(param1, param2, param3),
#' path = path
#' )
#' result[[1]]
#' result[[2]]
#' result[[3]]
#' }

metIdentify4all = function(
  ms1.data,
  ms2.data,
  parameter.list,
  path = "."
) {
  
  cat(crayon::yellow(
    "`metIdentify4all()` is deprecated, use `identify_metabolite_all()`."
  ))  
  
  threads = parameter.list[[1]]$threads
  ms1.data.name <- ms1.data
  ms2.data.name <- ms2.data
  file <- dir(path)
  intermediate_path <-
    file.path(path, "intermediate_data")
  dir.create(intermediate_path, showWarnings = FALSE)
  if (all(c("ms1.info", "ms2.info") %in% file)) {
    cat(crayon::yellow("Use old data\n"))
    load(file.path(intermediate_path, "ms1.info"))
    load(file.path(intermediate_path, "ms2.info"))
  } else{
    ##read MS2 data
    cat(crayon::green("Reading MS2 data...\n"))
    temp.ms2.type <-
      stringr::str_split(string = ms2.data.name,
                         pattern = "\\.")[[1]]
    temp.ms2.type <- temp.ms2.type[length(temp.ms2.type)]
    
    if (temp.ms2.type %in% c("mzXML", "mzML")) {
      ms2.data <-
        read_mzxml(file = file.path(path, ms2.data.name),
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
          temp.ms2.data <- read_mgf(file = file.path(path, temp.ms2.data))
        }
        temp.ms2.data
      })
      
      names(ms2.data) <- ms2.data.name
      ###prepare data for metIdentification function
      cat(crayon::green("Preparing MS2 data for identification...\n"))
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
      
      if (class(ms2.data) == "matrix") {
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
          paste(x, c(1:sum(ms1.info$name == x)), sep = "_")
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
  
  database.name <-
    unlist(lapply(parameter.list, function(x) {
      x$database
    }))
  
  if (!all(database.name %in% dir(path))) {
    stop(
      "The database: ",
      paste(database.name[!database.name %in% dir(path)],  collapse = ", "),
      "\n",
      " you want to use are not in you directory: \n",
      path
    )
  }
  
  
  identification.result <-
    vector(mode = "list", length = length(database.name))
  names(identification.result) <- database.name
  for (i in 1:length(database.name)) {
    cat(crayon::yellow("-------------------------------\n"))
    cat(crayon::yellow('Database', i, ":", database.name[i], "\n"))
    cat(crayon::yellow("-------------------------------\n"))
    new.path <-
      file.path(path, paste(database.name[i], "Result", sep =  '_'))
    dir.create(new.path, showWarnings = FALSE)
    if (any(dir(new.path) == "result")) {
      load(file.path(new.path, "result"))
      identification.result[[i]] <- result
      rm(list = "result")
      next()
    }
    
    if (length(parameter.list[[i]]) < 20) {
      result <- mzIdentify(
        ms1.data = ms1.data.name,
        ms1.match.ppm = parameter.list[[i]]$ms1.match.ppm,
        polarity = parameter.list[[i]]$polarity,
        column = parameter.list[[i]]$column,
        path = path,
        candidate.num = parameter.list[[i]]$candidate.num,
        database = parameter.list[[i]]$database,
        threads = parameter.list[[i]]$threads
      )
    } else{
      result <- metIdentify(
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
        path = path,
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






#' @title Generate the parameter list for metIdentify function
#' @description Generate the parameter list for metIdentify function.
#' \lifecycle{maturing}
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
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
#' @param database MS2 database name.
#' @param threads Number of threads
#' @return A metIdentifyClass object.
#' @export
#' @seealso The example and demo data of this function can be found
#' https://tidymass.github.io/metid/articles/multiple_databases.html
#' @examples 
#'  param1 <-
#' metIdentifyParam(
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

metIdentifyParam = function(
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
  total.score.tol = 0.5,
  candidate.num = 3,
  database,
  threads = 3
){
  
  cat(crayon::yellow(
    "`metIdentifyParam()` is deprecated, use `identify_metabolites_params()`."
  ))  
  
  
  if (missing(database)) {
    stop("The database name must be provided.\n")
  }
  # if(class(database) != "metIdentifyClass"){
  #   stop("database must be metIdentifyClass object\n")
  # }
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






#' @title Generate the mzIdentify parameter list
#' @description Generate the mzIdentify parameter list.
#' \lifecycle{deprecated}
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param ms1.match.ppm Precursor match ppm tolerance.
#' @param polarity The polarity of data, "positive"or "negative".
#' @param column "hilic" (HILIC column) or "rp" (reverse phase).
#' @param candidate.num The number of candidates.
#' @param database MS1 database name.
#' @param threads Number of threads
#' @return A mzIdentifyClass object.
#' @export
#' @seealso The example and demo data of this function can be found
#' https://tidymass.github.io/metid/articles/metid.html

mzIdentifyParam = function(
  ms1.match.ppm = 25,
  polarity = c("positive", "negative"),
  column = c("rp", "hilic"),
  candidate.num = 1,
  database,
  threads = 3
){
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
