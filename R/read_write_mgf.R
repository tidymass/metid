#' @title read_mgf_gnps
#' @description Read MGF data from GNPS
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param file The vector of names of ms2 files. MS2 file must be mgf.
#' @param threads threads number.
#' @importFrom future plan
#' @importFrom furrr future_map
#' @return Return ms2 data. This is a list.
#' @export

# sxtTools::setwd_project()
# setwd("other_files/all_ms2_database/gnps/")
# x = read_mgf_gnps(file = "HMDB.mgf")

read_mgf_gnps <-
  function(file,
           threads = 3) {
    # pbapply::pboptions(style = 1)
    message(crayon::green("Reading mgf data from GNPS..."))
    
    future::plan(strategy = future::multisession, workers = threads)
    
    ms2 <- furrr::future_map(
      .x = file,
      .f = function(temp.mgf.data) {
        mgf.data <- readr::read_lines(temp.mgf.data)
        mgf.data = mgf.data[mgf.data != ""]
        mgf.data = mgf.data[mgf.data != " "]
        
        begin_idx = which(mgf.data == "BEGIN IONS")
        end_idx = which(mgf.data == "END IONS")
        
        mgf.data =
          purrr::map2(
            .x = begin_idx,
            .y = end_idx,
            .f = function(idx1, idx2) {
              temp = mgf.data[c(idx1:idx2)]
              temp = temp[temp != "BEGIN IONS"]
              temp = temp[temp != "END IONS"]
              info = temp[grep(":", temp)]
              
              info = stringr::str_split(info, ":", n = 2) %>%
                do.call(rbind, .) %>%
                as.data.frame() %>%
                dplyr::rename(info = V1, value = V2)
              
              spec = temp[-grep(":", temp)]
              spec = stringr::str_split(spec, "\t") %>%
                do.call(rbind, .) %>%
                as.data.frame() %>%
                dplyr::rename(mz = V1, intensity = V2) %>%
                dplyr::mutate(mz = as.numeric(mz),
                              intensity = as.numeric(intensity)) %>%
                as.data.frame()
              list(info = info, spec = spec)
            }
          )
      },
      .progress = TRUE
    )
    
    ms2 = Reduce(`c`, ms2)
    message(crayon::green("Done."))
    ms2
  }


#' @title read_mgf_mona
#' @description Read MGF data (from mona).
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param file The vector of names of ms2 files. MS2 file must be mgf.
#' @param threads threads number.
#' @return Return ms2 data. This is a list.
#' @export

# sxtTools::setwd_project()
# setwd("other_files/all_ms2_database/mona/2021_6_10/")
# x = read_mgf_mona(file = "MoNA-export-LC-MS-MS_Spectra.mgf")

read_mgf_mona = function(file,
                         threads = 3) {
  # pbapply::pboptions(style = 1)
  message(crayon::green("Reading mgf data from MoNA..."))
  future::plan(strategy = future::multisession, workers = threads)
  ms2 <- furrr::future_map(
    .x = file,
    .f = function(temp.mgf.data) {
      mgf.data <- readr::read_lines(temp.mgf.data)
      mgf.data = mgf.data[mgf.data != ""]
      mgf.data = mgf.data[mgf.data != " "]
      
      begin_idx = which(mgf.data == "BEGIN IONS")
      end_idx = which(mgf.data == "END IONS")
      
      mgf.data =
        purrr::map2(
          .x = begin_idx,
          .y = end_idx,
          .f = function(idx1, idx2) {
            temp = mgf.data[c(idx1:idx2)]
            temp = temp[temp != "BEGIN IONS"]
            temp = temp[temp != "END IONS"]
            info = temp[grep(":", temp)]
            
            info = stringr::str_split(info, ":", n = 2) %>%
              do.call(rbind, .) %>%
              as.data.frame() %>%
              dplyr::rename(info = V1, value = V2)
            
            spec = temp[-grep(":", temp)]
            spec = stringr::str_split(spec, " ") %>%
              do.call(rbind, .) %>%
              as.data.frame() %>%
              dplyr::rename(mz = V1, intensity = V2) %>%
              dplyr::mutate(mz = as.numeric(mz),
                            intensity = as.numeric(intensity)) %>%
              as.data.frame()
            list(info = info, spec = spec)
          }
        )
      
    },
    .progress = TRUE
  )
  
  ms2 = Reduce(`c`, ms2)
  message(crayon::green("Done."))
  ms2
}



#' @title read_mgf_experiment
#' @description Read MGF data from experiment.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param file The vector of names of ms2 files. MS2 file must be mgf.
#' @param threads threads number.
#' @importFrom magrittr %>%
#' @importFrom furrr future_map
#' @return Return ms2 data. This is a list.
#' @export

# sxtTools::setwd_project()
# setwd("example/")
# file = c("QC1_MSMS_NCE25.mgf", "QC1_MSMS_NCE25.mgf")

read_mgf_experiment = function(file,
                               threads = 3) {
  message(crayon::green("Reading mgf data..."))
  future::plan(strategy = future::multisession, workers = threads)
  ms2 <- furrr::future_map(
    .x = file,
    .f = function(temp.mgf.data) {
      mgf.data <- readr::read_lines(temp.mgf.data)
      mgf.data = mgf.data[mgf.data != ""]
      mgf.data = mgf.data[mgf.data != " "]
      
      begin_idx = which(mgf.data == "BEGIN IONS")
      end_idx = which(mgf.data == "END IONS")
      
      mgf.data =
        purrr::map2(
          .x = begin_idx,
          .y = end_idx,
          .f = function(idx1, idx2) {
            temp = mgf.data[c(idx1:idx2)]
            temp = temp[temp != "BEGIN IONS"]
            temp = temp[temp != "END IONS"]
            info = grep("=", temp, value = TRUE)
            
            info = stringr::str_split(info, "=", n = 2) %>%
              do.call(rbind, .) %>%
              as.data.frame() %>%
              dplyr::rename(info = V1, value = V2)
            
            spec = temp[-grep("=", temp, value = FALSE)]
            spec = stringr::str_split(spec, " ") %>%
              do.call(rbind, .) %>%
              as.data.frame() %>%
              dplyr::rename(mz = V1, intensity = V2) %>%
              dplyr::mutate(mz = as.numeric(mz),
                            intensity = as.numeric(intensity)) %>%
              as.data.frame()
            list(info = info, spec = spec)
          }
        )
    },
    .progress = TRUE
  )
  
  ms2 = Reduce(`c`, ms2)
  message(crayon::green("Done."))
  ms2
}



ListMGF <- function(file) {
  mgf.data <- readLines(file)
  nl.rec.new <- 1
  idx.rec <- 1
  rec.list <- list()
  for (nl in seq_along(mgf.data))
  {
    if (mgf.data[nl] == "END IONS")
    {
      rec.list[idx.rec] <- list(Compound = mgf.data[nl.rec.new:nl])
      nl.rec.new <- nl + 1
      idx.rec <- idx.rec + 1
    }
  }
  rec.list
}















#######------------------------------------------------------------------------
#' @title readMGF
#' @description Read MGF data.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param file The vector of names of ms2 files. MS2 file must be mgf.
#' @return Return ms2 data. This is a list.
#' @export

readMGF <- function(file) {
  message(crayon::yellow("`readMGF()` is deprecated, use `read_mgf()`."))
  pbapply::pboptions(style = 1)
  message(crayon::green("Reading MS2 data..."))
  # mgf.data.list <- pbapply::pblapply(file, ListMGF)
  ms2 <- pbapply::pblapply(file, function(mgf.data) {
    mgf.data <- ListMGF(mgf.data)
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
      },
      x = mgf.data,
      y = nl.spec)
    } else{
      spec <- mapply(function(x, y) {
        do.call(rbind, strsplit(x[y], split = " "))
      },
      x = mgf.data,
      y = nl.spec)
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