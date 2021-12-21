# SXTMTmatch = function(data1,
#                       data2,
#                       mz.tol,
#                       #rt.tol is relative
#                       rt.tol = 30,
#                       rt.error.type = c("relative", "abs")){
#   rt.error.type <- match.arg(rt.error.type)
#   #
#   if (nrow(data1) == 0 | nrow(data2) == 0) {
#     result <- NULL
#     return(result)
#   }
#   
#   info1 <- data1[, c(1, 2), drop = FALSE]
#   info1 <- apply(info1, 1, list)
#   
#   mz2 <- as.numeric(data2[, 1])
#   rt2 <- as.numeric(data2[, 2])
#   
#   result <- pbapply::pblapply(info1, function(x) {
#     temp.mz1 <- x[[1]][[1]]
#     temp.rt1 <- x[[1]][[2]]
#     mz.error <- abs(temp.mz1 - mz2) * 10 ^ 6 / temp.mz1
#     if (rt.error.type == "relative") {
#       rt.error <- abs(temp.rt1 - rt2) * 100 / temp.rt1
#     } else{
#       rt.error <- abs(temp.rt1 - rt2)
#     }
#     
#     j <- which(mz.error <= mz.tol & rt.error <= rt.tol)
#     if (length(j) == 0) {
#       matrix(NA, ncol = 7)
#     } else{
#       cbind(j, temp.mz1, mz2[j], mz.error[j], temp.rt1, rt2[j], rt.error[j])
#     }
#   })
#   
#   if (length(result) == 1) {
#     result <- cbind(1, result[[1]])
#   } else{
#     result <- mapply(function(x, y) {
#       list(cbind(x, y))
#     },
#     x <- 1:length(info1),
#     y = result)
#     result <- do.call(rbind, result)
#   }
#   
#   result <-
#     matrix(result[which(!apply(result, 1, function(x)
#       any(is.na(x)))), ], ncol = 8)
#   if (nrow(result) == 0)
#     return(NULL)
#   colnames(result) <-
#     c("Index1",
#       "Index2",
#       "mz1",
#       "mz2",
#       "mz error",
#       "rt1",
#       "rt2",
#       "rt error")
#   result <- result
# }
# 


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

