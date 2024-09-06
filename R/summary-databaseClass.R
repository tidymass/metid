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