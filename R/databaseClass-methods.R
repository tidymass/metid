###S4 class for function metIdentification
#' An S4 class to represent MS1 or MS2 database.
#'
#' @slot database.info Database information.
#' @slot spectra.info Metabolites in database.
#' @slot spectra.data MS2 spectra data.
#' @exportClass databaseClass

setClass(
  Class = "databaseClass",
  representation(
    database.info = "list",
    spectra.info = "data.frame",
    spectra.data = "list"
  ),
  prototype = list(
    database.info = list(),
    spectra.info = data.frame(matrix(nrow = 0, ncol = 0), stringsAsFactors = FALSE),
    spectra.data = list()
  )
)

setMethod(
  f = "show",
  signature = "databaseClass",
  definition = function(object) {
    message(crayon::yellow("-----------Base information------------"))
    message("Version:", object@database.info$Version)
    message("Source:", object@database.info$Source)
    message("Link:", object@database.info$Link)
    message("Creater:",
        object@database.info$Creater,
        "(",
        object@database.info$Email,
        ")")
    message(
      ifelse(
        object@database.info$RT,
        "With RT information",
        "Without RT informtaion"
      )
    )
    message(crayon::yellow("-----------Spectral information------------"))
    message(ncol(object@spectra.info),
        " items of metabolite information:")
    
    message(paste(paste(head(
      colnames(object@spectra.info), 10
    ),
    collapse = "; "),
    ifelse(
      ncol(object@spectra.info) > 10, "(top10)", ""
    )))
    
    message(length(unique(object@spectra.info$Lab.ID)), " metabolites in total.")
    
    message(
      length(object@spectra.data$Spectra.positive),
      " metabolites with spectra in positive mode."
    )
    
    message(
      length(object@spectra.data$Spectra.negative),
      " metabolites with spectra in negative mode."
    )
    
    ce.pos <-
      unique(unlist(lapply(
        object@spectra.data$Spectra.positive, names
      )))
    
    ce.neg <-
      unique(unlist(lapply(
        object@spectra.data$Spectra.negative, names
      )))
    
    message("Collision energy in positive mode (number:):")
    message("Total number:", length(ce.pos))
    message(paste(head(ce.pos, 10), collapse = "; "), "")
    message("Collision energy in negative mode:")
    message("Total number:", length(ce.neg))
    message(paste(head(ce.neg, 10), collapse = "; "))
  }
)

#' @title Get MS2 spectra of peaks from databaseClass object
#' @description Get MS2 spectra of peaks from databaseClass object.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param lab.id The lab ID of metabolite.
#' @param database Database (databaseClass object).
#' @param polarity positive or negative.
#' @param ce Collision value.
#' @return A MS2 spectrum (data.frame).
#' @export
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/metid.html}

#' @title Get MS2 spectra of peaks from databaseClass object
#' @description Get MS2 spectra of peaks from databaseClass object.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param lab.id The lab ID of metabolite.
#' @param database Database (databaseClass object).
#' @param polarity positive or negative.
#' @param ce Collision value.
#' @return A MS2 spectrum (data.frame).
#' @export
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/metid.html}

get_ms2_spectrum = function(lab.id,
                            database,
                            polarity = c("positive", "negative"),
                            ce = "30") {
  polarity <- match.arg(polarity)
  if (!is(database, "databaseClass")) {
    stop("The database must be databaseClass object.\n")
  }
  pol <- ifelse(polarity == "positive", 1, 2)
  temp <-
    database@spectra.data[[pol]][[match(lab.id, names(database@spectra.data[[pol]]))]]
  temp[[match(ce, names(temp))]]
}


#'
#' #' setGeneric(
#' #'   name = "get_ms2_spectrum",
#' #'   def = function(object,
#' #'                  lab.id,
#' #'                  polarity = c("positive", "negative"),
#' #'                  ce = "30") standardGeneric(f = "get_ms2_spectrum"),
#' #'   signature = "object"
#' #'
#' #' )
#' #'
#'
#' #' setMethod(
#' #'   f = "get_ms2_spectrum",
#' #'   signature = "databaseClass",
#' #'   definition = function(object,
#' #'                         lab.id,
#' #'                         polarity = c("positive", "negative"),
#' #'                         ce = "30") {
#' #'     polarity <- match.arg(polarity)
#' #'     pol <- ifelse(polarity == "positive", 1, 2)
#' #'     temp <-
#' #'       object@spectra.data[[pol]][[match(lab.id, names(object@spectra.data[[pol]]))]]
#' #'     temp <- temp[[match(ce, names(temp))]]
#' #'     if(is.null(temp)){
#' #'       temp
#' #'     }else{
#' #'       tibble::as_tibble(temp)
#' #'     }
#' #'   }
#' #' )
