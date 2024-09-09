#' This function reads Mass Spectrometry (MS) data in MGF format from GNPS, extracts the MS2 spectra, and returns it in a structured format.
#'
#' @param file A character vector containing file paths to the MGF files from GNPS.
#' @param threads Numeric, the number of threads to use for parallel processing. Defaults to `3`.
#'
#' @return A list where each element contains the MS2 spectra and related metadata for each entry in the MGF file. Each entry includes:
#' \item{info}{A data frame with metadata (such as m/z, retention time, etc.) for each spectrum.}
#' \item{spec}{A data frame containing the `mz` and `intensity` values of the MS2 spectrum.}
#'
#' @details
#' The function parses MGF files from GNPS, a popular platform for metabolomics analysis. It extracts both the metadata and the MS2 spectra, and organizes the data into a structured format for further processing. Parallel processing is supported to speed up reading large files.
#'
#' @examples
#' \dontrun{
#' # Read MGF data from GNPS
#' mgf_data <- read_mgf_gnps(file = c("path/to/mgf1.mgf", "path/to/mgf2.mgf"))
#' }
#'
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @export

# setwd(masstools::get_project_wd())
# setwd("other_files/all_ms2_database/gnps/")
# x = read_mgf_gnps(file = "HMDB.mgf")

read_mgf_gnps <-
  function(file, threads = 3) {
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
                dplyr::mutate(mz = as.numeric(mz), intensity = as.numeric(intensity)) %>%
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


#' Read MGF Data from MoNA
#'
#' This function reads Mass Spectrometry (MS) data in MGF format from the MoNA (MassBank of North America) database, extracts the MS2 spectra, and returns them in a structured format.
#'
#' @param file A character vector containing file paths to the MGF files from MoNA.
#' @param threads Numeric, the number of threads to use for parallel processing. Defaults to `3`.
#'
#' @return A list where each element contains the MS2 spectra and related metadata for each entry in the MGF file. Each entry includes:
#' \item{info}{A data frame with metadata (such as m/z, retention time, etc.) for each spectrum.}
#' \item{spec}{A data frame containing the `mz` and `intensity` values of the MS2 spectrum.}
#'
#' @details
#' The function parses MGF files from MoNA, a public repository of mass spectrometry data, extracts both metadata and the MS2 spectra, and organizes the data into a structured format for further processing. Parallel processing is supported to improve efficiency when handling large datasets.
#'
#' @examples
#' \dontrun{
#' # Read MGF data from MoNA
#' mgf_data <- read_mgf_mona(file = c("path/to/mgf1.mgf", "path/to/mgf2.mgf"))
#' }
#'
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @export

# setwd(masstools::get_project_wd())
# setwd("other_files/all_ms2_database/mona/2021_6_10/")
# x = read_mgf_mona(file = "MoNA-export-LC-MS-MS_Spectra.mgf")

read_mgf_mona <-
  function(file, threads = 3) {
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
                dplyr::mutate(mz = as.numeric(mz), intensity = as.numeric(intensity)) %>%
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



#' Read MGF Data from Experimental Files
#'
#' This function reads Mass Spectrometry (MS) data in MGF format from experimental files, extracts the MS2 spectra, and returns them in a structured format.
#'
#' @param file A character vector containing file paths to the MGF files.
#' @param threads Numeric, the number of threads to use for parallel processing. Defaults to `3`.
#'
#' @return A list where each element contains the MS2 spectra and related metadata for each entry in the MGF file. Each entry includes:
#' \item{info}{A data frame with metadata (such as m/z, retention time, etc.) for each spectrum.}
#' \item{spec}{A data frame containing the `mz` and `intensity` values of the MS2 spectrum.}
#'
#' @details
#' The function parses MGF files from experimental data, extracting both the metadata and the MS2 spectra, and organizes the data into a structured format for further analysis. Parallel processing is supported to improve efficiency when handling large datasets.
#'
#' @examples
#' \dontrun{
#' # Read MGF data from experimental files
#' mgf_data <- read_mgf_experiment(file = c("path/to/mgf1.mgf", "path/to/mgf2.mgf"))
#' }
#'
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @export

# setwd(masstools::get_project_wd())
# setwd("example/")
# file = c("QC1_MSMS_NCE25.mgf", "QC1_MSMS_NCE25.mgf")

read_mgf_experiment <-
  function(file, threads = 3) {
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
                dplyr::mutate(mz = as.numeric(mz), intensity = as.numeric(intensity)) %>%
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



list_mgf <- function(file) {
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





#---------------------------------------------------------------------------
#' Read MSP Data
#'
#' This function reads Mass Spectrometry (MS) data in MSP format, processes the data, and extracts MS2 spectra along with metadata (such as m/z and retention time) from the file. It supports handling MSP data from different sources and formats.
#'
#' @param file A character string specifying the file path to the MSP file.
#' @param threads Numeric, the number of threads to use for parallel processing (not yet implemented). Defaults to `3`.
#'
#' @return A list where each element contains:
#' \item{info}{A list with metadata for each spectrum, typically containing `mz` (mass-to-charge ratio) and `rt` (retention time).}
#' \item{spec}{A data frame with the `mz` and `intensity` values of the MS2 spectrum.}
#'
#' @details
#' The function parses MSP files to extract both metadata and the MS2 spectra, organizing the data into a structured format for further analysis. It handles both regular MSP formats and those generated by the MetAnalyzer software, which require some additional processing.
#'
#' @examples
#' \dontrun{
#' # Read MSP data
#' msp_data <- read_msp(file = "path/to/data.msp")
#' }
#'
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @export

read_msp <-
  function(file, threads = 3) {
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
#' Read MSP Database
#'
#' This function reads Mass Spectrometry (MS) data in MSP format from a database, processes the data, and extracts MS2 spectra along with metadata (such as m/z, retention time, and other possible descriptors) from the file. The function supports handling various MSP formats, including those generated by specific software like MetAnalyzer.
#'
#' @param file A character string specifying the file path to the MSP file.
#' @param threads Numeric, the number of threads to use for parallel processing (not yet implemented). Defaults to `4`.
#'
#' @return A list where each element contains:
#' \item{info}{A named list containing metadata such as `mz` (mass-to-charge ratio), `rt` (retention time), and other descriptors for each spectrum.}
#' \item{spec}{A data frame containing the `mz` and `intensity` values of the MS2 spectrum.}
#'
#' @details
#' The function parses MSP files, extracting both metadata and MS2 spectra, organizing the data into a structured format for further analysis. It can handle MSP data from various sources, including those with custom formats like MetAnalyzer.
#'
#' @examples
#' \dontrun{
#' # Read MSP data from a database
#' msp_data <- read_msp_database(file = "path/to/database.msp")
#' }
#'
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @export

read_msp_database <-
  function(file, threads = 4) {
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
      
      temp.idx <- data.frame(temp.idx1, temp.idx2, stringsAsFactors = FALSE)
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
              temp.info <- data.frame(temp.info, stringsAsFactors = FALSE)
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
                temp.spec <- data.frame(temp.spec, stringsAsFactors = FALSE)
                colnames(temp.spec) <-
                  c('mz', 'intensity')
                rownames(temp.spec) <- NULL
                temp.spec$mz <-
                  as.numeric(as.character(temp.spec$mz))
                temp.spec$intensity <-
                  as.numeric(temp.spec$intensity)
                temp.spec <-
                  temp.spec[temp.spec$intensity != 0, ]
              } else{
                temp.spec <- NULL
              }
              
              list('info' = temp.info, 'spec' = temp.spec)
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
#' Read MSP Data from MoNA
#'
#' This function reads and processes MSP data files from the MoNA (MassBank of North America) database. The function parses the MSP file and extracts both metadata (such as m/z and intensity) and the spectrum information for each entry. The MSP data is then returned as a list.
#'
#' @param file A character vector specifying the file path(s) to the MSP file(s).
#' @param threads Numeric, the number of threads to use for parallel processing. Defaults to `3`.
#'
#' @return A list where each element contains:
#' \item{info}{A data frame with metadata for each spectrum (typically containing identifiers like m/z, intensity, and other descriptors).}
#' \item{spec}{A data frame with the `mz` (mass-to-charge ratio) and `intensity` values of the MS2 spectrum.}
#'
#' @details
#' This function is designed to handle MSP data specifically from MoNA. The data is organized into a structured list, where each list element corresponds to a spectrum with associated metadata and peak information.
#'
#' @examples
#' \dontrun{
#' # Read MSP data from MoNA
#' msp_data <- read_msp_mona(file = "path/to/mona_data.msp")
#' }
#'
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @export

# setwd(masstools::get_project_wd())
# setwd("other_files/all_ms2_database/mona/2021_6_10/")
# x = read_msp_mona(file = "MoNA-export-LC-MS-MS_Spectra.msp")

read_msp_mona <-
  function(file, threads = 3) {
    message(crayon::green("Reading msp data from MoNA..."))
    future::plan(strategy = future::multisession, workers = threads)
    ms2 <- furrr::future_map(
      .x = file,
      .f = function(temp.msp.data) {
        msp.data <- readr::read_lines(temp.msp.data, progress = FALSE)
        
        if (tail(msp.data, 1) == "") {
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
                dplyr::mutate(mz = as.numeric(mz), intensity = as.numeric(intensity)) %>%
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


#' ---------------------------------------------------------------------------
#' Read MSP Data from GNPS
#'
#' This function reads and processes MSP data files from the GNPS (Global Natural Products Social Molecular Networking) database. The function extracts both metadata (such as m/z and intensity) and the spectrum information for each entry. The MSP data is returned as a list of parsed entries.
#'
#' @param file A character vector specifying the file path(s) to the MSP file(s).
#' @param threads Numeric, the number of threads to use for parallel processing. Defaults to `5`.
#'
#' @return A list where each element contains:
#' \item{info}{A data frame with metadata for each spectrum (typically containing identifiers like m/z, intensity, and other descriptors).}
#' \item{spec}{A data frame with the `mz` (mass-to-charge ratio) and `intensity` values of the MS2 spectrum.}
#'
#' @details
#' This function is designed to handle MSP data from the GNPS database. The data is organized into a structured list, where each list element corresponds to a spectrum with associated metadata and peak information.
#'
#' @examples
#' \dontrun{
#' # Read MSP data from GNPS
#' msp_data <- read_msp_gnps(file = "path/to/gnps_data.msp")
#' }
#'
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @export

# setwd(masstools::get_project_wd())
# setwd("other_files/all_ms2_database/mona/2021_6_10/")
# x = read_msp_mona(file = "MoNA-export-LC-MS-MS_Spectra.msp")

read_msp_gnps <-
  function(file, threads = 5) {
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
                dplyr::mutate(mz = as.numeric(mz), intensity = as.numeric(intensity)) %>%
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




#' Write MSP Data to a File
#'
#' This function writes MSP (Mass Spectral Peak) data to a `.msp` file format. It saves the spectrum data along with precursor m/z information and metadata for each entry.
#'
#' @param info A data frame containing metadata information for each spectrum, where each row represents a spectrum. The data frame must contain a column named `'Mass'` representing the precursor m/z values.
#' @param fn.pre A character string representing the prefix for the output file name. The file will be saved with the format `prefix_spectra.msp`.
#' @param spec.all A list of data frames where each element corresponds to the MS/MS spectrum of a particular compound. Each data frame must have two columns: `mz` (mass-to-charge ratio) and `intensity`.
#'
#' @return The function does not return any output but writes an MSP file to the disk.
#'
#' @details
#' This function takes metadata and spectral data and writes them in the `.msp` format, which is commonly used to store mass spectra data. The output file is generated based on the provided prefix and contains all the spectra and their associated information.
#'
#' @examples
#' \dontrun{
#' # Example of writing MSP data
#' info <- data.frame(Mass = c(300.12, 250.14), Name = c("Compound1", "Compound2"))
#' spec.all <- list(
#'   data.frame(mz = c(100, 150, 200), intensity = c(1000, 1500, 2000)),
#'   data.frame(mz = c(110, 160, 210), intensity = c(1200, 1300, 1400))
#' )
#' write_msp(info, fn.pre = "example", spec.all = spec.all)
#' }
#'
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @export

write_msp <-
  function(info, fn.pre, spec.all) {
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
