#' Class for Managing Spectral Database Information
#'
#' The `databaseClass` S4 class is designed to store and manage information
#' related to spectral databases. It contains slots for database metadata,
#' spectra information, and the actual spectral data.
#'
#' @slot database.info A list containing metadata or general information about the database.
#' @slot spectra.info A data frame that holds information about each spectrum,
#' such as identifiers or attributes. The data frame is initialized empty.
#' @slot spectra.data A list where each element corresponds to a spectrum's data.
#'
#' @details
#' This class is used to manage and organize spectral data within a larger analytical pipeline.
#' The `database.info` slot stores high-level metadata, while `spectra.info` and `spectra.data`
#' contain details and the actual spectral information, respectively.
#'
#' @exportClass databaseClass
#'
#' @examples
#' # Creating an instance of the databaseClass
#' db_instance <- new("databaseClass")
#'
#' # Accessing slots
#' db_instance@database.info
#' db_instance@spectra.info
#' db_instance@spectra.data
#'


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
    message(ifelse(
      object@database.info$RT,
      "With RT information",
      "Without RT informtaion"
    ))
    message(crayon::yellow("-----------Spectral information------------"))
    message(ncol(object@spectra.info),
            " items of metabolite information:")
    
    message(paste(paste(head(
      colnames(object@spectra.info), 10
    ), collapse = "; "), ifelse(
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


#' Retrieve MS2 Spectrum from a Database
#'
#' This function extracts the MS2 spectrum data for a specific metabolite from a `databaseClass` object based on the provided Lab ID, polarity, and collision energy (CE).
#'
#' @param lab.id A character string representing the Lab ID of the metabolite whose MS2 spectrum is to be retrieved.
#' @param database A `databaseClass` object that contains the MS2 spectra data.
#' @param polarity A character string specifying the polarity mode for the spectrum. It can be either `"positive"` or `"negative"`. Defaults to `"positive"`.
#' @param ce A character string specifying the collision energy (CE) to retrieve. Defaults to `"30"`.
#'
#' @return A list containing the MS2 spectrum for the specified metabolite, polarity, and collision energy.
#'
#' @details
#' The function retrieves the MS2 spectrum for a given Lab ID from the provided `databaseClass` object. The user can specify the polarity mode (positive or negative) and the collision energy (CE) value to retrieve the correct spectrum. If the `database` is not a valid `databaseClass` object, the function stops with an error.
#'
#' @examples
#' \dontrun{
#' # Assuming `db_instance` is an instance of `databaseClass`
#' spectrum <- get_ms2_spectrum(lab.id = "M123", database = db_instance, polarity = "positive", ce = "30")
#' }
#'
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/metid.html}
#' 
#' @seealso \code{\link{databaseClass}}
#'
#' @export


get_ms2_spectrum <-
  function(lab.id,
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



#' @title colnames
#' @method colnames databaseClass
#' @param x x
#' @export
#' @rdname summary-databaseClass
#' @return message

colnames.databaseClass <-
  function(x) {
    colnames(x@spectra.info)
  }


#' @title nrow
#' @method nrow databaseClass
#' @param x x
#' @export
#' @rdname summary-databaseClass
#' @return message
nrow.databaseClass <- function(x) {
  nrow(x@spectra.info)
}

#' @title ncol
#' @method ncol databaseClass
#' @param x x
#' @export
#' @rdname summary-databaseClass
#' @return message
ncol.databaseClass <- function(x) {
  ncol(x@spectra.info)
}


#' @title dim
#' @method dim databaseClass
#' @param x x
#' @export
#' @rdname summary-databaseClass
#' @return message
dim.databaseClass <- function(x) {
  dim(x@spectra.info)
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
