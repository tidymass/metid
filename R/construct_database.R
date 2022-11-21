##------------------------------------------------------------------------------
#' @title Construct in-house or public MS2 database for metid.
#' @description Construct MS2 spectra database according to mzXML data and compound information table (csv format).
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param path Work directory.
#' @param version The version of you database. Default is 0.0.1.
#' @param metabolite.info.name The metabolite information table name, it must be csv format.
#' The demo data can be got from the `demoData` package.
#' Please see \url{https://tidymass.github.io/metid/articles/metid.html}
#' @param source The source of your database.
#' @param link Website link of the source.
#' @param creater Creater name. For example, Xiaotao Shen.
#' @param email email address.
#' @param rt Do the metabolites have RT information or not?. If not, set it as FALSE.
#' @param mz.tol m/z tolerance for the match between metabolites and precursor m/z of MS2 spectra.
#' @param rt.tol RT tolerance for the match between metabolites and precursor m/z of MS2 spectra.
#' @param threads The number of threads
#' @importFrom magrittr %>%
#' @importFrom crayon red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @importFrom readr cols
#' @importFrom pbapply pblapply
#' @return A databaseClass object.
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/metid.html}
#' @export
#' @examples
#' \dontrun{
#' database0.0.1 <- construct_database(
#' path = ".",
#' version = "0.0.1",
#' metabolite.info.name = "metabolite.info.csv",
#' creater = "dumine",
#' email = "dumine@zju.edu.cn",
#' rt = FALSE,
#' mz.tol = 15,
#' rt.tol = 30,
#' threads = 5
#' )
#' }

construct_database <-
  function(path = ".",
           version = "0.0.1",
           metabolite.info.name = "metabolite.info.csv",
           source = "Michael Snyder Lab",
           link = "http://snyderlab.stanford.edu/",
           creater = "Xiaotao Shen",
           email = "shenxt1990@163.com",
           rt = TRUE,
           mz.tol = 15,
           rt.tol = 30,
           threads = 3) {
    ##check data first
    file <- dir(path)
    if (all(file != metabolite.info.name)) {
      message(crayon::red("No", metabolite.info.name, "in your", path, "\n"))
      return(NULL)
    }
    
    if (all(file != "POS")) {
      message(crayon::red("No POS file in your", path, "\n"))
    } else{
      file_pos <- dir(file.path(path, "POS"))
      if (length(file_pos) == 0) {
        message(crayon::red("No mzXML/mgf files in POS folder."))
      } else{
        if (sum(stringr::str_detect(file_pos, "mzXML|mgf")) == 0) {
          message(crayon::red("No mzXML/mgf files in POS folder."))
        }
      }
    }
    
    if (all(file != "NEG")) {
      message(crayon::red("No NEG file in your", path, "\n"))
    } else{
      file_neg <- dir(file.path(path, "NEG"))
      if (length(file_neg) == 0) {
        message(crayon::red("No mzXML/mgf files in NEG folder."))
      } else{
        if (sum(stringr::str_detect(file_neg, "mzXML|mgf")) == 0) {
          message(crayon::red("No mzXML/mgf files in NEG folder."))
        }
      }
    }
    
    ##read metabolite information
    message(crayon::green("Reading metabolite information..."))
    metabolite.info <-
      readTable(file = file.path(path, metabolite.info.name))
    
    message(crayon::green("Reading positive MS2 data..."))
    
    file.pos <-
      dir(file.path(path, 'POS'), full.names = TRUE)
    
    if (length(file.pos) > 0) {
      if (stringr::str_detect(file.pos, "mzXML")[1]) {
        ms2.data.pos <-
          masstools::read_mzxml(file = file.pos, threads = threads)
      }
      
      if (stringr::str_detect(file.pos, "mgf")[1]) {
        ms2.data.pos <-
          masstools::read_mgf4database(file = file.pos)
      }
      
      ms1.info.pos <- lapply(ms2.data.pos, function(x) {
        x[[1]]
      })
      
      ms1.info.pos <- do.call(rbind, ms1.info.pos) %>%
        as.data.frame()
      
      ms1.info.pos$file <- basename(ms1.info.pos$file)
      
      ms2.info.pos <- lapply(ms2.data.pos, function(x) {
        x[[2]]
      })
      
      rm(list = "ms2.data.pos")
      
      message(crayon::red("OK."))
    } else{
      ms1.info.pos <- NULL
    }
    
    message(crayon::green("Reading negative MS2 data..."))
    
    file.neg <-
      dir(file.path(path, 'NEG'), full.names = TRUE)
    
    if (length(file.neg) > 0) {
      if (stringr::str_detect(file.neg, "mzXML")[1]) {
        ms2.data.neg <-
          masstools::read_mzxml(file = file.neg, threads = threads)
      }
      
      if (stringr::str_detect(file.neg, "mgf")[1]) {
        ms2.data.neg <-
          masstools::read_mgf4database(file = file.neg)
      }
      
      ms1.info.neg <- lapply(ms2.data.neg, function(x) {
        x[[1]]
      })
      
      ms1.info.neg <- do.call(rbind, ms1.info.neg) %>%
        as.data.frame()
      
      ms1.info.neg$file <- basename(ms1.info.neg$file)
      
      ms2.info.neg <- lapply(ms2.data.neg, function(x) {
        x[[2]]
      })
      
      rm(list = "ms2.data.neg")
      message(crayon::red("OK."))
    } else{
      ms1.info.neg <- NULL
    }
    
    ###---------------------------------------------------------------------------
    message(crayon::green("Matching metabolites with MS2 spectra (positive)..."))
    
    if (!is.null(ms1.info.pos)) {
      match.result.pos <-
        masstools::mz_rt_match(
          data1 = as.data.frame(metabolite.info[, c("mz.pos", "RT")]),
          data2 = ms1.info.pos[, c(2, 3)],
          mz.tol = mz.tol,
          rt.tol = rt.tol,
          rt.error.type = "abs"
        )
      
      match.result.pos <-
        data.frame(match.result.pos,
                   "file" = ms1.info.pos$file[match.result.pos[, 2]],
                   stringsAsFactors = FALSE)
      
      unique.idx1 <- unique(match.result.pos[, 1])
      
      spectra.pos <-
        pbapply::pblapply(unique.idx1, function(idx) {
          temp.match.result.pos <-
            match.result.pos[which(match.result.pos$Index1 == idx), , drop = FALSE]
          if (nrow(temp.match.result.pos) == 0) {
            return(NULL)
          }
          temp.submitter <- metabolite.info$Submitter[idx]
          if (!is.na(temp.submitter) &
              length(grep(temp.submitter, temp.match.result.pos[, 9])) > 0) {
            temp.match.result.pos <-
              temp.match.result.pos[grep(temp.submitter, temp.match.result.pos[, 9]), ]
          }
          
          if (nrow(temp.match.result.pos) == 0) {
            return(NULL)
          }
          
          if (nrow(temp.match.result.pos) == 1) {
            temp.ms2.pos <- ms2.info.pos[temp.match.result.pos[1, 2]]
            names(temp.ms2.pos) <-
              stringr::str_extract(string = temp.match.result.pos[1, 9],
                                   pattern = "NCE[0-9]{1,3}")
            return(temp.ms2.pos)
          }
          
          unique.file.name <-
            unique(temp.match.result.pos$file)
          
          temp.ms2.pos <-
            lapply(unique.file.name, function(temp.name) {
              temp.x <-
                temp.match.result.pos[which(temp.match.result.pos$file == temp.name), , drop = FALSE]
              temp.idx <-
                which.max(unlist(lapply(ms2.info.pos[temp.x[, 2]], function(y) {
                  sum(y[, 2])
                })))
              ms2.info.pos[[temp.x[temp.idx, 2]]]
            })
          
          names(temp.ms2.pos) <-
            stringr::str_extract(string = unique.file.name,
                                 pattern = "NCE[0-9]{1,3}")
          temp.ms2.pos
        })
      
      names(spectra.pos) <-
        metabolite.info$Lab.ID[unique.idx1]
      
      spectra.pos <-
        spectra.pos[which(!unlist(lapply(spectra.pos, is.null)))]
      
      message(crayon::red("OK."))
    } else{
      spectra.pos <- NULL
    }
    
    ###---------------------------------------------------------------------------
    message(crayon::green("Matching metabolites with MS2 spectra (negative)..."))
    if (!is.null(ms1.info.neg)) {
      match.result.neg <-
        masstools::mz_rt_match(
          data1 = as.data.frame(metabolite.info[, c("mz.neg", "RT")]),
          data2 = ms1.info.neg[, c(2, 3)],
          mz.tol = mz.tol,
          rt.tol = rt.tol,
          rt.error.type = "abs"
        )
      
      match.result.neg <- data.frame(match.result.neg,
                                     "file" = ms1.info.neg$file[match.result.neg[, 2]],
                                     stringsAsFactors = FALSE)
      
      unique.idx1 <- unique(match.result.neg[, 1])
      
      spectra.neg <-
        pbapply::pblapply(unique.idx1, function(idx) {
          temp.match.result.neg <-
            match.result.neg[which(match.result.neg$Index1 == idx), , drop = FALSE]
          if (nrow(temp.match.result.neg) == 0) {
            return(NULL)
          }
          
          temp.submitter <- metabolite.info$Submitter[idx]
          if (!is.na(temp.submitter) &
              length(grep(temp.submitter, temp.match.result.neg[, 9])) > 0) {
            temp.match.result.neg <-
              temp.match.result.neg[grep(temp.submitter, temp.match.result.neg[, 9]), ]
          }
          
          if (nrow(temp.match.result.neg) == 0) {
            return(NULL)
          }
          
          if (nrow(temp.match.result.neg) == 1) {
            temp.ms2.neg <- ms2.info.neg[temp.match.result.neg[1, 2]]
            names(temp.ms2.neg) <-
              stringr::str_extract(string = temp.match.result.neg[1, 9],
                                   pattern = "NCE[0-9]{1,3}")
            return(temp.ms2.neg)
          }
          
          unique.file.name <-
            unique(temp.match.result.neg$file)
          
          temp.ms2.neg <-
            lapply(unique.file.name, function(temp.name) {
              temp.x <-
                temp.match.result.neg[which(temp.match.result.neg$file == temp.name), , drop = FALSE]
              temp.idx <-
                which.max(unlist(lapply(ms2.info.neg[temp.x[, 2]], function(y) {
                  sum(y[, 2])
                })))
              ms2.info.neg[[temp.x[temp.idx, 2]]]
            })
          
          names(temp.ms2.neg) <-
            stringr::str_extract(string = unique.file.name,
                                 pattern = "NCE[0-9]{1,3}")
          temp.ms2.neg
        })
      
      names(spectra.neg) <-
        metabolite.info$Lab.ID[unique.idx1]
      
      spectra.neg <-
        spectra.neg[which(!unlist(lapply(spectra.neg, is.null)))]
      message(crayon::red("OK."))
    } else{
      spectra.neg <- NULL
    }
    
    if (is.null(spectra.pos) & is.null(spectra.neg)) {
      Spectra <- list()
    } else{
      Spectra <- list("Spectra.positive" = spectra.pos,
                      "Spectra.negative" = spectra.neg)
    }
    
    database.info <- list(
      "Version" = version,
      "Source" = source,
      "Link" = link,
      "Creater" = creater,
      "Email" = email,
      "RT" = rt
    )
    
    spectra.info <- as.data.frame(metabolite.info)
    rm(list = "metabolite.info")
    
    msDatabase0.0.1 <- new(
      Class = "databaseClass",
      database.info = database.info,
      spectra.info = spectra.info,
      spectra.data = Spectra
    )
    
    msDatabase0.0.1@database.info$RT <-
      ifelse(all(is.na(msDatabase0.0.1@spectra.info$RT)), FALSE, TRUE)
    message(crayon::bgRed("All done!"))
    return(msDatabase0.0.1)
  }