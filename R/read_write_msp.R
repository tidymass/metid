#---------------------------------------------------------------------------
#' @title read_msp
#' @description Read MSP data.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param file The vector of names of ms2 files. MS2 file must be msp format.
#' @param threads threads number.
#' @return Return ms2 data. This is a list.
#' @export

read_msp = function(file,
                    threads = 3) {
  
  message(crayon::green("Reading msp data..."))
  
  # future::plan(strategy = future::multisession, workers = threads)
  
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
    
    temp.idx <- data.frame(temp.idx1, temp.idx2,
                           stringsAsFactors = FALSE)
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
      temp.info <- data.frame(temp.info,
                              stringsAsFactors = FALSE)
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
        temp.spec <- data.frame(temp.spec,
                                stringsAsFactors = FALSE)
        colnames(temp.spec) <- c('mz', 'intensity')
        rownames(temp.spec) <- NULL
        temp.spec$mz <- as.numeric(as.character(temp.spec$mz))
        temp.spec$intensity <- as.numeric(temp.spec$intensity)
        temp.spec <- temp.spec[temp.spec$intensity != 0,]
      } else{
        temp.spec <- NULL
      }
      
      list('info' = temp.info,
           'spec' = temp.spec)
    })
    
    mz.idx <- grep("[Mm][Zz]", rownames(info.spec[[1]][[1]]))
    rt.idx <-
      grep("Time$|TIME$|time$|^RT|^rt|^Rt", rownames(info.spec[[1]][[1]]))
    
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


#---------------------------------------------------------------------------
#' @title read_msp_database
#' @description Read MSP data for database.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param file The vector of names of ms2 files. MS2 file must be msp format.
#' @param threads threads
#' @return Return ms2 data. This is a list.
#' @export

read_msp_database = function(file,
                             threads = 4) {
  msp.data <- readr::read_lines(file)
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
    # future::plan(strategy = future::multisession, workers = threads)
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
    
    temp.idx <- data.frame(temp.idx1, temp.idx2,
                           stringsAsFactors = FALSE)
    temp.idx <- apply(temp.idx, 1, list)
    
    temp.idx <- lapply(temp.idx, unlist)
    
    # n.spec <- which(grepl('^\\d', msp.data))
    # n.info <- seq(n.tot)[-c(n.spec, n.null)]
    
    # pbapply::pboptions(style = 1)
    info.spec =
      purrr::map(
        temp.idx,
        .f =
          function(idx) {
            temp.msp.data <- msp.data[idx[1]:idx[2]]
            
            temp.msp.data <-
              temp.msp.data[temp.msp.data != ""]
            info.idx <-
              grep("[A-Za-z]", temp.msp.data)
            temp.info <- temp.msp.data[info.idx]
            temp.info <-
              strsplit(temp.info, split = ":")
            temp.info <- do.call(rbind, temp.info)
            temp.info <- data.frame(temp.info,
                                    stringsAsFactors = FALSE)
            temp.info[, 2] <-
              stringr::str_trim(temp.info[, 2])
            colnames(temp.info) <-
              rownames(temp.info) <- NULL
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
              temp.spec <- data.frame(temp.spec,
                                      stringsAsFactors = FALSE)
              colnames(temp.spec) <-
                c('mz', 'intensity')
              rownames(temp.spec) <- NULL
              temp.spec$mz <-
                as.numeric(as.character(temp.spec$mz))
              temp.spec$intensity <-
                as.numeric(temp.spec$intensity)
              temp.spec <-
                temp.spec[temp.spec$intensity != 0,]
            } else{
              temp.spec <- NULL
            }
            
            list('info' = temp.info,
                 'spec' = temp.spec)
          }
      )
    
    mz.idx <- grep("[Mm][Zz]", rownames(info.spec[[1]][[1]]))
    rt.idx <-
      grep("Time$|TIME$|time$|^RT|^rt|^Rt", rownames(info.spec[[1]][[1]]))
    
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
      info.spec <- purrr::map(info.spec, function(x) {
        info <- x[[1]]
        mz <- as.numeric(info[mz.idx, 1])
        rt <- as.numeric(info[rt.idx, 1])
        info2 = info[-c(mz.idx, rt.idx), , drop = FALSE]
        info1 <- c(mz, rt)
        info = c(info1, info2[, 1])
        names(info) <- c("mz", "rt", rownames(info2))
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


#' ---------------------------------------------------------------------------
#' @title readMSP_MoNA
#' @description Read MSP data from MoNA.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param file The vector of names of ms2 files. MS2 file must be msp. The msp data must from MoNA.
#' @return Return ms2 data. This is a list.
#' @export

readMSP_MoNA = function(file) {
  message(crayon::yellow("`readMSP_MoNA()` is deprecated, use `read_msp_mona()`."))
  message(crayon::green("Reading MSP data..."))
  msp.data <- readr::read_lines(file)
  n.null <- which(msp.data == '')
  temp.idx1 <- c(1, n.null[-length(n.null)])
  temp.idx2 <- n.null - 1
  
  temp.idx <- data.frame(temp.idx1, temp.idx2,
                         stringsAsFactors = FALSE)
  
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
    temp.info <- data.frame(temp.info,
                            stringsAsFactors = FALSE)
    temp.info[, 2] <-
      stringr::str_trim(temp.info[, 2], side = "both")
    temp.info[, 1] <-
      stringr::str_trim(temp.info[, 1], side = "both")
    colnames(temp.info) <- rownames(temp.info) <- NULL
    #combine synons
    if (length(grep("Synon", temp.info[, 1])) > 1) {
      Synon <-
        stringr::str_c(temp.info[stringr::str_which(string = temp.info[, 1], pattern = "Synon"), 2,
                                 drop = TRUE], collapse = ";")
      temp.info[which(temp.info[, 1] == "Synon")[1], 2] <- Synon
      temp.info <- temp.info[!duplicated(temp.info[, 1]),]
      
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
      temp.spec <- data.frame(temp.spec,
                              stringsAsFactors = FALSE)
      colnames(temp.spec) <- c('mz', 'intensity')
      rownames(temp.spec) <- NULL
      temp.spec$mz <- as.numeric(as.character(temp.spec$mz))
      temp.spec$intensity <- as.numeric(temp.spec$intensity)
      temp.spec <- temp.spec[temp.spec$intensity != 0,]
    } else{
      temp.spec <- NULL
    }
    
    list('info' = temp.info,
         'spec' = temp.spec)
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



#' ---------------------------------------------------------------------------
#' @title read_msp_mona
#' @description Read MSP data from MoNA.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param file The vector of names of ms2 files. MS2 file must be msp. The msp data must from MoNA.
#' @param threads The number of threads
#' @return Return ms2 data. This is a list.
#' @export

# sxtTools::setwd_project()
# setwd("other_files/all_ms2_database/mona/2021_6_10/")
# x = read_msp_mona(file = "MoNA-export-LC-MS-MS_Spectra.msp")

read_msp_mona = function(file,
                         threads = 3) {
  message(crayon::green("Reading msp data from MoNA..."))
  future::plan(strategy = future::multisession, workers = threads)
  ms2 <- furrr::future_map(
    .x = file,
    .f = function(temp.msp.data) {
      msp.data <- readr::read_lines(temp.msp.data, progress = FALSE)
      
      if(tail(msp.data, 1) == ""){
        msp.data = msp.data[-length(msp.data)]
      }
      msp.data[msp.data == ""] = rep(c("END IONS", "BEGIN IONS"), sum(msp.data == "") /
                                       2)
      msp.data = c("BEGIN IONS", msp.data, "END IONS")
      
      begin_idx = which(msp.data == "BEGIN IONS")
      end_idx = which(msp.data == "END IONS")
      
      msp.data =
        purrr::map2(
          .x = begin_idx,
          .y = end_idx,
          .f = function(idx1, idx2) {
            temp = msp.data[c(idx1:idx2)]
            temp = temp[temp != "BEGIN IONS"]
            temp = temp[temp != "END IONS"]
            if (length(temp) == 0) {
              return(NULL)
            }
            info = temp[grep(":", temp, value = FALSE)]
            
            info = stringr::str_split(info, ":", n = 2) %>%
              do.call(rbind, .) %>%
              as.data.frame() %>%
              dplyr::rename(info = V1, value = V2)
            
            spec = temp[-grep(":", temp, value = FALSE)]
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
    }, .progress = TRUE
  )
  
  ms2 = Reduce(`c`, ms2)
  message(crayon::green("Done."))
  ms2
}


#' ---------------------------------------------------------------------------
#' @title read_msp_gnps
#' @description Read MSP data from MoNA.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param file The vector of names of ms2 files. MS2 file must be msp. The msp data must from MoNA.
#' @param threads The number of threads
#' @return Return ms2 data. This is a list.
#' @export

# sxtTools::setwd_project()
# setwd("other_files/all_ms2_database/mona/2021_6_10/")
# x = read_msp_mona(file = "MoNA-export-LC-MS-MS_Spectra.msp")

read_msp_gnps = function(file,
                         threads = 5) {
  pbapply::pboptions(style = 1)
  message(crayon::green("Reading msp data from GNPS..."))
  
  future::plan(strategy = future::multisession, workers = threads)
  ms2 <- furrr::future_map(
    .x = file,
    .f = function(temp.msp.data) {
      msp.data <- readr::read_lines(temp.msp.data)
      msp.data[msp.data == ""] = rep(c("END IONS", "BEGIN IONS"), sum(msp.data == "") /
                                       2)
      msp.data = c("BEGIN IONS", msp.data, "END IONS")
      
      begin_idx = which(msp.data == "BEGIN IONS")
      end_idx = which(msp.data == "END IONS")
      
      msp.data =
        purrr::map2(
          .x = begin_idx,
          .y = end_idx,
          .f = function(idx1, idx2) {
            # cat(paste(idx1, idx2, sep = " "))
            temp = msp.data[c(idx1:idx2)]
            temp = temp[temp != "BEGIN IONS"]
            temp = temp[temp != "END IONS"]
            if (length(temp) == 0) {
              return(NULL)
            }
            info = temp[grep(":", temp, value = FALSE)]
            
            info = stringr::str_split(info, ":", n = 2) %>%
              do.call(rbind, .) %>%
              as.data.frame() %>%
              dplyr::rename(info = V1, value = V2)
            
            spec = temp[-grep(":", temp, value = FALSE)]
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
      
    }, .progress = TRUE
  )
  
  ms2 = Reduce(`c`, ms2)
  message(crayon::green("Done."))
  ms2
}








WriteMSP = function(info, fn.pre, spec.all) {
  fn.save <- paste0(fn.pre, '_spectra.msp')
  #
  sink(fn.save)
  for (idx in seq(nrow(info))) {
    if (!is.null(spec.all[[idx]])) {
      if (nrow(spec.all[[idx]]) > 0) {
        mz <- info[idx, 'Mass']
        spec <- spec.all[[idx]]
        cat('IDX: ', idx, '\n', sep = '')
        cat('PRECURSORMZ: ', mz, '\n', sep = '')
        cat('Num Peaks: ', nrow(spec), '\n', sep = '')
        for (nr in seq(nrow(spec))) {
          cat(paste(spec[nr, ], collapse = ' '), '\n', sep = '')
        }
        cat('\n')
      }
    }
  }
  sink()
}




#---------------------------------------------------------------------------
#' @title readMSP
#' @description Read MSP data.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param file The vector of names of ms2 files. MS2 file must be msp format.
#' @param silence.deprecated Silenc the deprecated information or not.
#' @return Return ms2 data. This is a list.
#' @export

readMSP = function(file,
                   silence.deprecated = TRUE) {
  if (!silence.deprecated) {
    message(crayon::yellow("`readMSP()` is deprecated, use `read_msp()`."))
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
    
    temp.idx <- data.frame(temp.idx1, temp.idx2,
                           stringsAsFactors = FALSE)
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
      temp.info <- data.frame(temp.info,
                              stringsAsFactors = FALSE)
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
        temp.spec <- data.frame(temp.spec,
                                stringsAsFactors = FALSE)
        colnames(temp.spec) <- c('mz', 'intensity')
        rownames(temp.spec) <- NULL
        temp.spec$mz <- as.numeric(as.character(temp.spec$mz))
        temp.spec$intensity <- as.numeric(temp.spec$intensity)
        temp.spec <- temp.spec[temp.spec$intensity != 0,]
      } else{
        temp.spec <- NULL
      }
      
      list('info' = temp.info,
           'spec' = temp.spec)
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