getExtension = function(file){
  tail(stringr::str_split(string = file, pattern = "\\.")[[1]], 1)
}

readTable = function(file, ...){
  extension <- getExtension(file = file)
  if (extension == "csv") {
    return(readr::read_csv(file = file, ...))
  }
  
  if (extension == 'xlsx') {
    return(readxl::read_xlsx(path = file, ...))
  }
  
  if (extension == "xls") {
    return(readxl::read_xls(path = file, ...))
  }
  
  if (extenstion != "csv" &
      extenstion != "xlsx" &
      extenstion != "xls") {
    cat(crayon::red("file are not csv, xlsx or xls.\n"))
  }
}





msg <- function(..., startup = FALSE) {
  if (startup) {
    if (!isTRUE(getOption("metid.quiet"))) {
      packageStartupMessage(text_col(...))
    }
  } else {
    message(text_col(...))
  }
}

text_col <- function(x) {
  # If RStudio not available, messages already printed in black
  if (!rstudioapi::isAvailable()) {
    return(x)
  }
  
  if (!rstudioapi::hasFun("getThemeInfo")) {
    return(x)
  }
  
  theme <- rstudioapi::getThemeInfo()
  
  if (isTRUE(theme$dark)) crayon::white(x) else crayon::black(x)
  
}

#' List all packages in the metid
#'
#' @param include_self Include metid in the list?
#' @export
#' @examples
#' metid_packages()
metid_packages <- function(include_self = TRUE) {
  raw <- utils::packageDescription("metid")$Imports
  imports <- strsplit(raw, ",")[[1]]
  parsed <- gsub("^\\s+|\\s+$", "", imports)
  names <- vapply(strsplit(parsed, "\\s+"), "[[", 1, FUN.VALUE = character(1))
  
  if (include_self) {
    names <- c(names, "metid")
  }
  
  names
}

invert <- function(x) {
  if (length(x) == 0) return()
  stacked <- utils::stack(x)
  tapply(as.character(stacked$ind), stacked$values, list)
}


style_grey <- function(level, ...) {
  crayon::style(
    paste0(...),
    crayon::make_style(grDevices::grey(level), grey = TRUE)
  )
}








#' @title Get MS2 spectra of peaks from databaseClass object
#' @description Get MS2 spectra of peaks from databaseClass object.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param lab.id The lab ID of metabolite.
#' @param database Database (databaseClass object).
#' @param polarity positive or negative.
#' @param ce Collision value.
#' @importFrom crayon red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @importFrom readr cols
#' @importFrom pbapply pblapply
#' @return A MS2 spectrum (data.frame).
#' @export
#' @seealso The example and demo data of this function can be found
#' https://tidymass.github.io/metid/articles/metid.html

getMS2spectrum = function(lab.id,
                          database,
                          polarity = c("positive", "negative"),
                          ce = "30") {
  cat(crayon::yellow(
    "`getMS2spectrum()` is deprecated, use `get_ms2_spectrum()`."
  ))
  polarity <- match.arg(polarity)
  if (class(database) != "databaseClass") {
    stop("The database must be databaseClass object.\n")
  }
  pol <- ifelse(polarity == "positive", 1, 2)
  temp <-
    database@spectra.data[[pol]][[match(lab.id, names(database@spectra.data[[pol]]))]]
  temp[[match(ce, names(temp))]]
}

