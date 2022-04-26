#' @method filter databaseClass
#' @docType methods
#' @importFrom rlang quos !!!
#' @importFrom dplyr filter
#' @export
filter.databaseClass <-
  function(.data, ..., .preserve = FALSE) {
    dots <- rlang::quos(...)
    
    x <- .data@spectra.info
    
    x <-
      dplyr::filter(x, !!!dots, .preserve = .preserve)
    
    if(nrow(x) == 0){
      stop("No compound left.")
    }
    slot(object = .data, name = "spectra.info") = x
    
    if(length(.data@spectra.data$Spectra.positive) > 0){
      remain_idx <- 
        which(names(.data@spectra.data$Spectra.positive) %in% x$Lab.ID)
      .data@spectra.data$Spectra.positive <- 
        .data@spectra.data$Spectra.positive[remain_idx]
    }
    
    if(length(.data@spectra.data$Spectra.negative) > 0){
      remain_idx <- 
        which(names(.data@spectra.data$Spectra.negative) %in% x$Lab.ID)
      .data@spectra.data$Spectra.negative <- 
        .data@spectra.data$Spectra.negative[remain_idx]
    }
    return(.data)
  }

#' @importFrom dplyr filter
#' @export
dplyr::filter
