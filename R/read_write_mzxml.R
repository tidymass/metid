#' @title read_mzxml
#' @description Read mzXML data.
#' \lifecycle{experimental}
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param file The vector of names of ms2 files. MS2 file must be mzXML or mzML.
#' @param threads Thread number
#' @return Return ms2 data. This is a list.
#' @export

read_mzxml = function(file,
                     threads = 3) {
  # pbapply::pboptions(style = 1)
  cat(crayon::green("Reading MS2 data...\n"))
  # mzxml.data.list <- pbapply::pblapply(file, ListMGF)
  ms2 <-
    MSnbase::readMSData(files = file,
                        msLevel. = 2,
                        mode = "onDisk")
  cat(crayon::green("Processing...\n"))
  
  new.ms2 <- ProtGenerics::spectra(object = ms2)
  rm(list = c("ms2"))
  #
  temp.fun <- function(idx, ms2) {
    temp.ms2 <- ms2[[idx]]
    rm(list = c("ms2"))
    info <-
      data.frame(
        name = paste("mz", temp.ms2@precursorMz,
                     "rt", temp.ms2@rt, sep = ""),
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
          paste(x, c(1:sum(info$name == x)), sep = "_")
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
  
  # new.ms2 <-
  #   BiocParallel::bplapply(
  #     X = c(1:length(new.ms2)),
  #     FUN = temp.fun,
  #     BPPARAM = BiocParallel::SnowParam(workers = threads,
  #                                       progressbar = TRUE,
  #                                       type = "SOCK"),
  #     ms2 = new.ms2
  #   )
  
  new.ms2 <-
    BiocParallel::bplapply(
      X = c(1:length(new.ms2)),
      FUN = temp.fun,
      BPPARAM = BiocParallel::MulticoreParam(workers = threads,
                                             progressbar = TRUE),
      ms2 = new.ms2
    )
  
  new.ms2 <- new.ms2
}




#' @title readMZXML
#' @description Read mzXML data.
#' \lifecycle{experimental}
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param file The vector of names of ms2 files. MS2 file must be mzXML or mzML.
#' @param threads Thread number
#' @return Return ms2 data. This is a list.
#' @export

readMZXML = function(file,
                     threads = 3) {
  cat(crayon::yellow("`readMZXML()` is deprecated, use `read_mzxml()`.\n"))
  # pbapply::pboptions(style = 1)
  cat(crayon::green("Reading MS2 data...\n"))
  # mzxml.data.list <- pbapply::pblapply(file, ListMGF)
  ms2 <-
    MSnbase::readMSData(files = file,
                        msLevel. = 2,
                        mode = "onDisk")
  cat(crayon::green("Processing...\n"))
  
  new.ms2 <- ProtGenerics::spectra(object = ms2)
  rm(list = c("ms2"))
  #
  temp.fun <- function(idx, ms2) {
    temp.ms2 <- ms2[[idx]]
    rm(list = c("ms2"))
    info <-
      data.frame(
        name = paste("mz", temp.ms2@precursorMz,
                     "rt", temp.ms2@rt, sep = ""),
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
          paste(x, c(1:sum(info$name == x)), sep = "_")
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
  
  # new.ms2 <-
  #   BiocParallel::bplapply(
  #     X = c(1:length(new.ms2)),
  #     FUN = temp.fun,
  #     BPPARAM = BiocParallel::SnowParam(workers = threads,
  #                                       progressbar = TRUE,
  #                                       type = "SOCK"),
  #     ms2 = new.ms2
  #   )
  
  new.ms2 <-
    BiocParallel::bplapply(
      X = c(1:length(new.ms2)),
      FUN = temp.fun,
      BPPARAM = BiocParallel::MulticoreParam(workers = threads,
                                             progressbar = TRUE),
      ms2 = new.ms2
    )
  
  new.ms2 <- new.ms2
}