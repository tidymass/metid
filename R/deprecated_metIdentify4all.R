
#' @title Generate the mzIdentify parameter list
#' @description Generate the mzIdentify parameter list.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param ms1.match.ppm Precursor match ppm tolerance.
#' @param polarity The polarity of data, "positive"or "negative".
#' @param column "hilic" (HILIC column) or "rp" (reverse phase).
#' @param candidate.num The number of candidates.
#' @param database MS1 database name.
#' @param threads Number of threads
#' @return A mzIdentifyClass object.
#' @export
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/metid.html}

mzIdentifyParam = function(
  ms1.match.ppm = 25,
  polarity = c("positive", "negative"),
  column = c("rp", "hilic"),
  candidate.num = 1,
  database,
  threads = 3
){
  if (missing(database)) {
    stop("The database name must be provided.\n")
  }
  if (database != "HMDB.metabolite.data") {
    stop("database can only be HMDB.metabolite.data\n")
  }
  polarity <- match.arg(polarity)
  column <- match.arg(column)
  
  param <- list(
    ms1.match.ppm = ms1.match.ppm,
    polarity = polarity,
    column = column,
    candidate.num = candidate.num,
    database = database,
    threads = threads
  )
  list("mzIdentifyParam" = param)
}
