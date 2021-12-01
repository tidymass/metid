##------------------------------------------------------------------------------
#' @title Get identification table from a metIdentifyClass object
#' @description Get identification table from a metIdentifyClass object.
#' \lifecycle{deprecated}
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param ... One or multiple metIdentifyClass objects.
#' @param candidate.num The number of candidates.
#' @param type The type of identification table.
#' @return A identification table (data.frame).
#' @export
#' @importFrom magrittr %>%
#' @importFrom crayon red yellow green bgRed
#' @importFrom stringr str_detect str_extract str_extract_all str_replace_all str_replace
#' @importFrom readr cols
#' @importFrom pbapply pblapply
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate filter everything select bind_rows
#' @importFrom plyr dlply .
#' @seealso The example and demo data of this function can be found
#' https://tidymass.github.io/metid/articles/metid.html

getIdentificationTable = function(
  ...,
  candidate.num = 3,
  type = c("old", "new")
){
  
  cat(crayon::yellow(
    "`getIdentificationTable()` is deprecated, use `get_identification_table()`."
  ))  
  
  candidate.num <- round(candidate.num)
  if (candidate.num <= 0) {
    candidate.num <- 1
  }
  
  object <- list(...)
  
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
      cat(crayon::yellow("The object is identified without MS2 spectra.\n"))
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
          x <- x[1:candidate.num, , drop = FALSE]
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
    }
    rownames(identification.table) <- NULL
    return(tibble::as_tibble(identification.table))
  } else{
    len <- lapply(object, function(x) {
      nrow(x@match.result)
    }) %>% unlist()
    
    if (any(len == 0)) {
      warning(crayon::yellow(
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
              x <- x[1:candidate.num, , drop = FALSE]
            }
            x <- data.frame(x,
                            "MS2.spectra.name" = y,
                            stringsAsFactors = FALSE)
            list(x)
          },
          x = iden.result,
          y = names(iden.result))
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
              x <- x[1:candidate.num, , drop = FALSE]
            }
            x <- data.frame(x,
                            "MS2.spectra.name" = y,
                            stringsAsFactors = FALSE)
            list(x)
          },
          x = iden.result,
          y = names(iden.result))
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
      z = 1:length(identification.result)
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
            x <- x[1:candidate.num, , drop = FALSE]
          }
          x
        }
      )
    
    # if(type == "old"){
    identification.table <-
      as.data.frame(matrix(nrow = nrow(peak.table), ncol = 3))
    colnames(identification.table) <-
      c("MS2.spectra.name",
        "Candidate.number",
        "Identification")
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
      sapply(identification.table$Identification, function(x) {
        if (is.na(x))
          return(0)
        return(length(stringr::str_split(
          string = x, pattern = "\\{\\}"
        )[[1]]))
      })
    
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
    if(type == "new"){
      identification.table <- trans_to_new_style(identification.table = identification.table)
    }
    return(tibble::as_tibble(identification.table))
  }
}






##------------------------------------------------------------------------------
#' @title Transform old style identification table to new style
#' @description Transform old style identification table to new style.
#' \lifecycle{deprecated}
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param identification.table Identification table from getIdentificationTable or getIdentificationTable2.
#' @return A identification table (data.frame).
#' @export
#' @importFrom magrittr %>%
#' @importFrom crayon red yellow green bgRed
#' @importFrom stringr str_detect str_extract str_extract_all str_replace_all str_replace
#' @importFrom readr cols
#' @importFrom pbapply pblapply
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate filter everything select bind_rows
#' @importFrom plyr dlply .
#' @seealso The example and demo data of this function can be found
#' https://tidymass.github.io/metid/articles/metid.html

trans2newStyle = function(identification.table){
  
  cat(crayon::yellow(
    "`trans2newStyle()` is deprecated, use `trans_to_new_style()`."
  ))  
  
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
            sapply(c(
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
              
            })
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


#------------------------------------------------------------------------------
#' @title Get parameters from a metIdentifyClass object
#' @description Get parameters from a metIdentifyClass object.
#' \lifecycle{deprecated}
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param object A metIdentifyClass object.
#' @return A data frame contains all the parameters of this metIdentifiyClass object.
#' @export
#' @importFrom magrittr %>%
#' @importFrom crayon red yellow green bgRed
#' @importFrom stringr str_detect str_extract str_extract_all str_replace_all str_replace
#' @importFrom readr cols
#' @importFrom pbapply pblapply
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate filter everything select bind_rows
#' @importFrom plyr dlply .

getParams = function(object){
  
  cat(crayon::yellow(
    "`getParams()` is deprecated, use `get_parameters_metid()`."
  ))  
  
  if (class(object) == "mzIdentifyClass"){
    stop(crayon::red('Please use getParams2 to get the parameters for mzIdentifyClass object.\n'))
  }
  if (class(object) != "metIdentifyClass")
    stop(crayon::red("Only for metIdentifyClass\n"))
  data.frame(
    "Parameter" = c(
      "ms1.ms2.match.mz.tol",
      "ms1.ms2.match.rt.tol",
      "ms1.match.ppm",
      "ms2.match.ppm",
      "ms2.match.tol",
      "rt.match.tol",
      "polarity",
      "ce",
      "column",
      "ms1.match.weight",
      "rt.match.weight",
      "ms2.match.weight",
      "path",
      "total.score.tol",
      "candidate.num",
      "database",
      "threads"
    ),
    "Meaning" = c(
      "MS1 features & MS spectra matching mz tolerance (ppm)",
      "MS1 features & MS spectra matching RT tolerance (s)",
      "MS1 match tolerance (ppm)",
      "MS2 fragment match tolerance (ppm)",
      "MS2 match tolerance",
      "RT match tolerance (s)",
      "Polarity",
      "Collision energy",
      "Column",
      "MS1 match weight",
      "RT match weight",
      "MS2 match weight",
      "Work directory",
      "Total score tolerance",
      "Candidate number",
      "MS2 database",
      "Thread number"
    ),
    "Value" = c(
      object@ms1.ms2.match.mz.tol,
      object@ms1.ms2.match.rt.tol,
      object@ms1.match.ppm,
      object@ms2.match.ppm,
      object@ms2.match.tol,
      object@rt.match.tol,
      object@polarity,
      object@ce,
      object@column,
      object@ms1.match.weight,
      object@rt.match.weight,
      object@ms2.match.weight,
      object@path,
      object@total.score.tol,
      object@candidate.num,
      object@database,
      object@threads
    ),
    stringsAsFactors = FALSE
  ) %>%
    tibble::as_tibble()
}


##------------------------------------------------------------------------------
#' @title Get identification information from a metIdentifyClass object
#' @description Get identification information from a metIdentifyClass object.
#' \lifecycle{deprecated}
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param object A metIdentifyClass object.
#' @param which.peak A peak name or "all". "all" means all peaks with identifications will be output.
#' @param database Database used.
#' @return A identification table (data.frame).
#' @importFrom magrittr %>%
#' @importFrom crayon red yellow green bgRed
#' @importFrom stringr str_detect str_extract str_extract_all str_replace_all str_replace
#' @importFrom readr cols
#' @importFrom pbapply pblapply
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate filter everything select bind_rows
#' @importFrom plyr dlply .
#' @export
#' @seealso The example and demo data of this function can be found
#' https://tidymass.github.io/metid/articles/metid.html


getIdenInfo = function(object,
                       which.peak,
                       database){
  
  cat(crayon::yellow(
    "`getIdenInfo()` is deprecated, use `get_iden_info()`."
  ))  
  
  if (missing(object) | missing(which.peak) | missing(database)) {
    stop("Please provide the object, which.peak and database.\n")
  }
  if (class(object) != "metIdentifyClass")
    stop("Only for metIdentifyClass\n")
  if (class(database) != "databaseClass")
    stop("Only for databaseClass\n")
  
  identification.result <- object@identification.result
  
  which.peak <- as.character(which.peak)
  
  if (!which.peak %in% object@ms1.data$name) {
    stop(which.peak, " is not in peak table, please check it.\n")
  }
  
  if (is.na(match(which.peak, object@match.result$MS1.peak.name))) {
    cat(crayon::green("The peak has no MS2 spectrum.\n"))
    return()
  }
  
  if (is.na(match(
    object@match.result$MS2.spectra.name[match(which.peak, object@match.result$MS1.peak.name)],
    names(identification.result)
  ))) {
    cat(crayon::green("The peak has no identification result.\n"))
    return()
  }
  
  temp <-
    match(object@match.result$MS2.spectra.name[match(which.peak, object@match.result$MS1.peak.name)],
          names(identification.result))
  temp <- identification.result[[temp]]
  temp <-
    data.frame(temp, database@spectra.info[match(temp$Lab.ID, database@spectra.info$Lab.ID),
                                           setdiff(colnames(database@spectra.info), colnames(temp))], stringsAsFactors = FALSE)
  temp <- tibble::as_tibble(temp)
  temp
}



#------------------------------------------------------------------------------
#' @title Get the peak names which have identifications
#' @description Get the peak names which have identifications.
#' \lifecycle{deprecated}
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param object A metIdentifyClass object.
#' @return Peak names with identifications.
#' @export
#' @importFrom magrittr %>%
#' @importFrom crayon red yellow green bgRed
#' @importFrom stringr str_detect str_extract str_extract_all str_replace_all str_replace
#' @importFrom readr cols
#' @importFrom pbapply pblapply
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate filter everything select bind_rows
#' @importFrom plyr dlply .
#' @seealso The example and demo data of this function can be found
#' https://tidymass.github.io/metid/articles/metid.html

whichHasIden = function(object){
  cat(crayon::yellow(
    "`whichHasIden()` is deprecated, use `which_has_identification()`."
  ))  
  if (class(object) != "metIdentifyClass")
    stop("Only for metIdentifyClass\n")
  
  temp <-
    object@match.result[match(names(object@identification.result),
                              object@match.result$MS2.spectra.name), c(3, 4)]
  rownames(temp) <- NULL
  temp
}



#------------------------------------------------------------------------------
#' @title Filter identifications according to m/z error, RT error, MS similarity and total score
#' @description Filter identifications according to m/z error, RT error, MS similarity and total score.
#' \lifecycle{deprecated}
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param object A metIdentifyClass object.
#' @param ms1.match.ppm MS1 match ppm.
#' @param rt.match.tol RT match tolerance.
#' @param ms2.match.tol MS2 match (MS2 similarity) tolerance.
#' @param total.score.tol Total score tolerance.
#' @return A new metIdentifyClass.
#' @importFrom magrittr %>%
#' @importFrom crayon red yellow green bgRed
#' @importFrom stringr str_detect str_extract str_extract_all str_replace_all str_replace
#' @importFrom readr cols
#' @importFrom pbapply pblapply
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate filter everything select bind_rows
#' @importFrom plyr dlply .
#' @export
#' @seealso The example and demo data of this function can be found
#' https://tidymass.github.io/metid/articles/metid.html

filterIden = function(
  object,
  ms1.match.ppm = 25,
  rt.match.tol = 30,
  ms2.match.tol = 0.5,
  total.score.tol = 0.5
){
  
  cat(crayon::yellow(
    "`filterIden()` is deprecated, use `filter_identification()`."
  ))  
  
  if (class(object) != "metIdentifyClass") {
    stop("Only for metIdentifyClass\n")
  }
  
  object@ms1.match.ppm <- ms1.match.ppm
  object@rt.match.tol <- rt.match.tol
  object@ms2.match.tol <- ms2.match.tol
  object@total.score.tol <- total.score.tol
  
  identification.result <- object@identification.result
  
  identification.result <-
    lapply(identification.result, function(x) {
      RT.error <- x$RT.error
      RT.error[is.na(RT.error)] <- rt.match.tol - 1
      SS <- x$SS
      SS[is.na(SS)] <- ms2.match.tol + 1
      SS[SS == 0] <- ms2.match.tol + 1
      x <-
        x[which(
          x$mz.error < ms1.match.ppm & RT.error < rt.match.tol &
            SS > ms2.match.tol &
            x$Total.score > total.score.tol
        ), , drop = FALSE]
    })
  
  temp.idx <-
    which(unlist(lapply(identification.result, function(x) {
      nrow(x) != 0
    })))
  
  identification.result <-
    identification.result[temp.idx]
  object@identification.result <- identification.result
  object
}




#------------------------------------------------------------------------------
#' @title Get spectra of peaks from metIdentifyClass object
#' @description Get spectra of peaks from metIdentifyClass object.
#' \lifecycle{maturing}
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param object metIdentifyClass.
#' @param peak.name Peak name.
#' @return A MS2 spectrum.
#' @importFrom magrittr %>%
#' @importFrom crayon red yellow green bgRed
#' @importFrom stringr str_detect str_extract str_extract_all str_replace_all str_replace
#' @importFrom readr cols
#' @importFrom pbapply pblapply
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate filter everything select bind_rows
#' @importFrom plyr dlply .
#' @export
#' @seealso The example and demo data of this function can be found
#' https://tidymass.github.io/metid/articles/metid.html

getMS2spectrum2Object = function(object,
                                 peak.name){
  
  cat(crayon::yellow(
    "`getMS2spectrum2Object()` is deprecated, use `get_ms2_spectrum_from_object()`."
  ))  
  
  if (class(object) != "metIdentifyClass")
    stop("Only for metIdentifyClass\n")
  if (missing(peak.name))
    stop('Please provide peak name.\n')
  object@ms2.info[[which(object@match.result$MS2.spectra.name[match(peak.name, object@match.result$MS1.peak.name)] == names(object@ms2.info))]]
}

