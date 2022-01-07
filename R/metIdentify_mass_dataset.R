#' @title Identify metabolites based on MS/MS database.
#' @description Identify metabolites based on MS/MS database.
#' \lifecycle{deprecated}
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param object A mass_dataset class object.
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
#' @param total.score.tol Total score tolerance. The total score are refering to MS-DIAL.
#' @param candidate.num The number of candidate.
#' @param database MS2 database name or MS2 database.
#' @param threads Number of threads
#' @return A metIdentifyClass object.
#' @importFrom crayon yellow green red bgRed
#' @export
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/metid.html}

metIdentify_mass_dataset = function(object,
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
                                    column = c("hilic", "rp"),
                                    ms1.match.weight = 0.25,
                                    rt.match.weight = 0.25,
                                    ms2.match.weight = 0.5,
                                    total.score.tol = 0.5,
                                    candidate.num = 3,
                                    database,
                                    threads = 3) {
  ###Check data
  if (missing(database)) {
    stop("No database is provided.\n")
  }
  
  ##parameter specification
  polarity <- match.arg(polarity)
  column <- match.arg(column)
  
  if (class(database) != "databaseClass") {
    stop("database should be databaseClass object.\n")
  }
  
  #load MS2 database
  database.name = paste(database@database.info$Source,
                        database@database.info$Version,
                        sep = "_")
  
  
  if (class(database) != "databaseClass") {
    stop("database must be databaseClass object\n")
  }
  
  ce.list.pos <-
    unique(unlist(lapply(
      database@spectra.data$Spectra.positive, names
    )))
  
  ce.list.neg <-
    unique(unlist(lapply(
      database@spectra.data$Spectra.negative, names
    )))
  
  ce.list <-
    ifelse(polarity == "positive", ce.list.pos, ce.list.neg)
  
  if (all(ce %in% ce.list) & ce != "all") {
    stop("All ce values you set are not in database. Please check it.\n")
    ce <- ce[ce %in% ce.list]
  }
  
  rm(list = c("ce.list.pos", "ce.list.neg", "ce.list"))
  
  ##ce values
  if (all(ce != "all")) {
    if (polarity == "positive") {
      ce.list <-
        unique(unlist(
          lapply(database@spectra.data$Spectra.positive, function(x) {
            names(x)
          })
        ))
      if (length(grep("Unknown", ce.list)) > 0) {
        ce <-
          unique(c(ce, grep(
            pattern = "Unknown", ce.list, value = TRUE
          )))
      }
    } else{
      ce.list <-
        unique(unlist(
          lapply(database@spectra.data$Spectra.negative, function(x) {
            names(x)
          })
        ))
      if (length(grep("Unknown", ce.list)) > 0) {
        ce <-
          unique(c(ce, grep(
            pattern = "Unknown", ce.list, value = TRUE
          )))
      }
    }
  }
  
  ##RT in database or not
  if (!database@database.info$RT) {
    cat(
      crayon::yellow(
        "No RT information in database.\nThe weight of RT have been set as 0.\n"
      )
    )
  }
  
  #------------------------------------------------------------------
  ##load adduct table
  if (polarity == "positive" & column == "hilic") {
    data("hilic.pos", envir = environment())
    adduct.table <- hilic.pos
  }
  
  if (polarity == "positive" & column == "rp") {
    data("rp.pos", envir = environment())
    adduct.table <- rp.pos
  }
  
  if (polarity == "negative" & column == "hilic") {
    data("hilic.neg", envir = environment())
    adduct.table <- hilic.neg
  }
  
  if (polarity == "negative" & column == "rp") {
    data("rp.neg", envir = environment())
    adduct.table <- rp.neg
  }

  if (length(object@ms2_data) == 0) {
    stop("No MS2 in you object.\n")
  }
  
  if (lapply(object@ms2_data, function(x) {
    length(x@ms2_spectra)
  }) %>%
  unlist() %>%
  sum() %>%
  `==`(0)) {
    stop("No MS2 in you object.\n")
  }
  
  #####annotaion result for each set MS2 data
  annotation_result = 
    purrr::map2(.x = names(object@ms2_data), 
                .y = object@ms2_data, 
                function(temp_ms2_data_id,
                         temp_ms2_data) {
                  cat(crayon::yellow(temp_ms2_data_id, "file:\n"))
                  cat(crayon::green(length(temp_ms2_data@ms2_spectra), "MS2 spectra.\n"))
                  
                  ms1.info = data.frame(
                    name = temp_ms2_data@ms2_spectrum_id,
                    mz = temp_ms2_data@ms2_mz,
                    rt = temp_ms2_data@ms2_rt,
                    file = temp_ms2_data@ms2_file,
                    variable_id = temp_ms2_data@variable_id
                  )
                  
                  ms2.info = temp_ms2_data@ms2_spectra
                  
                  ms2_matchresult <-
                    metIdentification(
                      ms1.info = ms1.info,
                      ms2.info = ms2.info,
                      polarity = polarity,
                      ce = ce,
                      database = database,
                      ms1.match.ppm = ms1.match.ppm,
                      ms2.match.ppm = ms2.match.ppm,
                      mz.ppm.thr = mz.ppm.thr,
                      ms2.match.tol = ms2.match.tol,
                      rt.match.tol = rt.match.tol,
                      column = column,
                      ms1.match.weight = ms1.match.weight,
                      rt.match.weight = rt.match.weight,
                      ms2.match.weight = ms2.match.weight,
                      total.score.tol = total.score.tol,
                      candidate.num = candidate.num,
                      adduct.table = adduct.table,
                      threads = threads,
                      fraction.weight = fraction.weight,
                      dp.forward.weight = dp.forward.weight,
                      dp.reverse.weight = dp.reverse.weight
                    )
                  
                  ms2_matchresult =
                    purrr::map2(
                      .x = names(ms2_matchresult),
                      .y = ms2_matchresult,
                      .f = function(temp_ms2_id, 
                                    temp_annotation_result) {
                        data.frame(ms2_files_id = temp_ms2_data_id, 
                                   ms2_spectrum_id = temp_ms2_id,
                                   temp_annotation_result
                        ) %>% 
                          dplyr::left_join(ms1.info[,c("name", "variable_id")], by = c("ms2_spectrum_id" = "name")) %>% 
                          dplyr::select(variable_id, dplyr::everything())
                        
                      }
                    ) %>% 
                    dplyr::bind_rows()
                  ms2_matchresult
                })
  
  annotation_result = 
    annotation_result %>% 
    dplyr::bind_rows() %>% 
    as.data.frame() %>% 
    dplyr::mutate(Database = database.name)
  
  cat(crayon::bgRed("All done.\n"))
  return(annotation_result)
}