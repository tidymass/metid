
#' @title Identify metabolites based on MS1 or MS/MS database
#' @description Identify metabolites based on MS1 or MS/MS database.
#' \lifecycle{maturing}
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param ms1.data The name of ms1 peak table (csv format). Column 1 is "name", Column 2 is
#' "mz" and column is "rt" (second).
#' @param ms2.data MS2 data, must be mgf, msp or mzXML format. For example, ms2.data = c("test.mgf", "test2.msp").
#' @param ms1.ms2.match.mz.tol MS1 peak and MS2 spectrum matching m/z tolerance. Default is 25 pm.
#' @param ms1.ms2.match.rt.tol MS1 peak and MS2 spectrum matching RT tolerance. Default is 10 s.
#' @param ms1.match.ppm Precursor match ppm tolerance.
#' @param ms2.match.ppm Fragment ion match ppm tolerance.
#' @param mz.ppm.thr Accurate mass tolerance for m/z error calculation.
#' @param ms2.match.tol MS2 match (MS2 similarity) tolerance.
#' @param fraction.weight The weight for matched fragments.
#' @param dp.forward.weight Forward dot product weight.
#' @param dp.reverse.weight Reverse dot product weight.
#' @param rt.match.tol RT match tolerance.
#' @param polarity The polarity of data, "positive"or "negative".
#' @param ce Collision energy. Please confirm the CE values in your database. Default is "all".
#' @param column "hilic" (HILIC column) or "rp" (reverse phase).
#' @param ms1.match.weight The weight of MS1 match for total score calculation.
#' @param rt.match.weight The weight of RT match for total score calculation.
#' @param ms2.match.weight The weight of MS2 match for total score calculation.
#' @param path Work directory.
#' @param total.score.tol Total score tolerance. The total score are refering to MS-DIAL.
#' @param candidate.num The number of candidate.
#' @param database MS2 database name or MS database.
#' @param threads Number of threads
#' @return A metIdentifyClass object.
#' @importFrom crayon yellow green red bgRed
#' @importFrom magrittr %>%
#' @export
#' @seealso The example and demo data of this function can be found
#' https://tidymass.github.io/metid/articles/metid.html

identify_metabolites = function(
  ms1.data,
  ms2.data = NULL,
  ms1.ms2.match.mz.tol = 25,
  ms1.ms2.match.rt.tol = 10,
  ms1.match.ppm = 25,
  ms2.match.ppm = 30,
  mz.ppm.thr = 400,
  ms2.match.tol = 0.5,
  fraction.weight = 0.3,
  dp.forward.weight = 0.6,
  dp.reverse.weight = 0.1,
  rt.match.tol = 30,
  polarity = c("positive", "negative"),
  ce = "all",
  column = c("rp", "hilic"),
  ms1.match.weight = 0.25,
  rt.match.weight = 0.25,
  ms2.match.weight = 0.5,
  path = ".",
  total.score.tol = 0.5,
  candidate.num = 3,
  database,
  threads = 3
) {

  ###Check data
  if (missing(database)) {
    stop("No database is provided.\n")
  }

  if (missing(ms1.data)) {
    stop("Please provide MS1 data name.\n")
  }

  ##parameter specification
  polarity <- match.arg(polarity)
  column <- match.arg(column)
  ##check ms1.file and ms2.file
  file <- dir(path)

  if (!all(ms1.data %in% file)) {
    stop("MS1 data is not in the directory, please check it.\n")
  }

  if (!is.null(ms2.data)) {
    if (!all(ms2.data %in% file)) {
      stop("Some MS2 data are not in the directory, please check it.\n")
    }
  }

  if(class(database) != "databaseClass"){
    if (!all(database %in% file)) {
      stop("Database is not in this directory, please check it.\n")
    }  
  }
  
  if (is.null(ms2.data)) {
    cat(crayon::yellow("You don't provide MS2 data, so only use mz and/or RT for matching.\n"))
    mzIdentify(
      ms1.data = ms1.data,
      rt.match.tol = rt.match.tol,
      ms1.match.ppm = ms1.match.ppm,
      polarity = polarity,
      column = column,
      path = path,
      candidate.num = candidate.num,
      database = database,
      threads = threads,
      silence.deprecated = TRUE
    )
  } else{
    metIdentify(
      ms1.data = ms1.data,
      ms2.data = ms2.data,
      ##only msp and mgf and mz(X)ML are supported
      ms1.ms2.match.mz.tol = ms1.ms2.match.mz.tol,
      ms1.ms2.match.rt.tol = ms1.ms2.match.rt.tol,
      ms1.match.ppm = ms1.match.ppm,
      ms2.match.ppm = ms2.match.ppm,
      mz.ppm.thr = mz.ppm.thr,
      ms2.match.tol = ms2.match.tol,
      fraction.weight = fraction.weight,
      dp.forward.weight = dp.forward.weight,
      dp.reverse.weight = dp.reverse.weight,
      rt.match.tol = rt.match.tol,
      polarity = polarity,
      ce = ce,
      column = column,
      ms1.match.weight = ms1.match.weight,
      rt.match.weight = rt.match.weight,
      ms2.match.weight = ms2.match.weight,
      path = path,
      total.score.tol = total.score.tol,
      candidate.num = candidate.num,
      database = database,
      threads = threads,
      silence.deprecated = TRUE
    )
  }
}







