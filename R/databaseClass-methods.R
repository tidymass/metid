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
    cat(crayon::yellow("-----------Base information------------\n"))
    cat(crayon::green("Version:", object@database.info$Version, "\n"))
    cat(crayon::green("Source:", object@database.info$Source, "\n"))
    cat(crayon::green("Link:", object@database.info$Link, "\n"))
    cat(
      crayon::green(
        "Creater:",
        object@database.info$Creater,
        "(",
        object@database.info$Email,
        ")\n"
      )
    )
    cat(crayon::green(
      ifelse(
        object@database.info$RT,
        "With RT information\n",
        "Without RT informtaion\n"
      )
    ))
    cat(crayon::yellow("-----------Spectral information------------\n"))
    cat(crayon::green(
      "There are",
      ncol(object@spectra.info),
      "items of metabolites in database:\n"
    ))
    cat(crayon::green(paste(
      colnames(object@spectra.info), collapse = "; "
    ), "\n"))
    
    cat(crayon::green("There are", length(
      unique(object@spectra.info$Compound.name)
    ), "metabolites in total\n"))
    
    cat(
      crayon::green(
        "There are",
        length(object@spectra.data$Spectra.positive),
        "metabolites in positive mode with MS2 spectra.\n"
      )
    )
    
    cat(
      crayon::green(
        "There are",
        length(object@spectra.data$Spectra.negative),
        "metabolites in negative mode with MS2 spectra.\n"
      )
    )
    
    ce.pos <-
      unique(unlist(lapply(
        object@spectra.data$Spectra.positive, names
      )))
    
    ce.neg <-
      unique(unlist(lapply(
        object@spectra.data$Spectra.negative, names
      )))
    
    cat(crayon::green("Collision energy in positive mode (number:):\n"))
    cat(crayon::green("Total number:", length(ce.pos), "\n"))
    cat(crayon::green(paste(head(ce.pos, 10), collapse = "; "), "\n"))
    cat(crayon::green("Collision energy in negative mode:\n"))
    cat(crayon::green("Total number:", length(ce.neg), "\n"))
    cat(crayon::green(paste(head(ce.neg, 10), collapse = "; "), "\n"))
  }
)

#' @title Get MS2 spectra of peaks from databaseClass object
#' @description Get MS2 spectra of peaks from databaseClass object.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param lab.id The lab ID of metabolite.
#' @param database Database (databaseClass object).
#' @param polarity positive or negative.
#' @param ce Collision value.
#' @return A MS2 spectrum (data.frame).
#' @export
#' @seealso The example and demo data of this function can be found
#' https://tidymass.github.io/metid/articles/metid.html

#' @title Get MS2 spectra of peaks from databaseClass object
#' @description Get MS2 spectra of peaks from databaseClass object.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param lab.id The lab ID of metabolite.
#' @param database Database (databaseClass object).
#' @param polarity positive or negative.
#' @param ce Collision value.
#' @return A MS2 spectrum (data.frame).
#' @export
#' @seealso The example and demo data of this function can be found
#' https://tidymass.github.io/metid/articles/metid.html

get_ms2_spectrum = function(lab.id,
                            database,
                            polarity = c("positive", "negative"),
                            ce = "30") {
  polarity <- match.arg(polarity)
  if (class(database) != "databaseClass") {
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


