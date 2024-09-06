#' @title Identify metabolites based on MS1 or MS/MS database for single peak.
#' @description Identify metabolites based on MS1 or MS/MS database.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param object A mass_dataset class obejct.
#' @param variable_id variable_id
#' @param variable_index variable_index
#' @param polarity The polarity of data, "positive"or "negative".
#' @param ms1.match.ppm Precursor match ppm tolerance.
#' @param ms2.match.ppm Fragment ion match ppm tolerance.
#' @param mz.ppm.thr Accurate mass tolerance for m/z error calculation.
#' @param database MS2 database name or MS database.
#' @param interactive_plot interactive_plot
#' @return A list of ggplot2 object.
#' @importFrom crayon yellow green red bgRed
#' @importFrom magrittr %>%
#' @importFrom masstools ms2_plot
#' @export
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/metid.html}

ms2_plot_mass_dataset <-
  function(object,
           variable_id,
           variable_index,
           polarity = c("positive", "negative"),
           ms1.match.ppm = 25,
           ms2.match.ppm = 30,
           mz.ppm.thr = 400,
           database,
           interactive_plot = FALSE) {
    polarity = match.arg(polarity)
    massdataset::check_object_class(object = object, class = "mass_dataset")
    
    if (nrow(object@annotation_table) == 0) {
      stop("No annotation in object.\n")
    } else{
      if (all(object@annotation_table$Level != 1) &
          all(object@annotation_table$Level != 2)) {
        stop("No annotations with MS2.\n")
      }
    }
    
    ###check parameters
    if (!is.numeric(ms1.match.ppm)) {
      stop("ms1.match.ppm should be numeric.\n")
    } else{
      if (ms1.match.ppm <= 0 | ms1.match.ppm >= 500) {
        stop("ms1.match.ppm should > 0 and < 500\n")
      }
    }
    
    if (!is.numeric(ms2.match.ppm)) {
      stop("ms2.match.ppm should be numeric.\n")
    } else{
      if (ms2.match.ppm <= 0 | ms2.match.ppm >= 500) {
        stop("ms2.match.ppm should > 0 and < 500\n")
      }
    }
    
    ###Check data
    if (missing(database)) {
      stop("No database is provided.\n")
    }
    
    if (!is(database, "databaseClass")) {
      stop("database should be databaseClass object.\n")
    }
    
    database.name <-
      paste(database@database.info$Source,
            database@database.info$Version,
            sep = "_")
    
    ######check variable_id and variable_index
    if (missing(variable_id) & missing(variable_index)) {
      stop("provide variable_id or variable_index.\n")
    }
    
    if (!missing(variable_id)) {
      purrr::walk(
        variable_id,
        .f = function(temp_variable_id) {
          if (!temp_variable_id %in% object@variable_info$variable_id) {
            stop(paste(temp_variable_id, "is not in variable_info.\n"))
          }
        }
      )
      variable_index = match(variable_id, object@variable_info$variable_id)
    } else{
      purrr::walk(
        variable_index,
        .f = function(temp_variable_index) {
          if (temp_variable_index <= 0 |
              temp_variable_index > nrow(object@variable_info)) {
            stop(
              "variable_index ",
              temp_variable_index,
              " should be range from 1 to ",
              nrow(object@variable_info)
            )
          }
        }
      )
    }
    
    variable_id <-
      object@variable_info$variable_id[variable_index] %>%
      unique() %>%
      `[`(1)
    
    temp_variable_id <-
      variable_id
    
    temp_annotation_table <-
      object@annotation_table %>%
      dplyr::filter(variable_id == temp_variable_id) %>%
      dplyr::filter(!is.na(SS))
    
    if (nrow(temp_annotation_table) == 0) {
      message(paste(temp_variable_id, "has no annotation with MS2."))
      return(NULL)
    }
    
    all_plot <-
      purrr::map(
        as.data.frame(t(temp_annotation_table)),
        .f = function(x) {
          temp_idx <-
            which(object@ms2_data[[x[2]]]@ms2_spectrum_id == x[3])[1]
          spectrum1 <-
            object@ms2_data[[x[2]]]@ms2_spectra[[temp_idx]]
          spectrum2 <-
            get_ms2_spectrum(
              lab.id = x[8],
              polarity = polarity,
              database = database,
              ce = x[14]
            )
          if (is.null(spectrum2)) {
            message("database may be wrong.")
            plot <-
              masstools::ms2_plot(
                spectrum1 = spectrum1,
                spectrum1_name = x[1],
                spectrum2_name = x[4],
                ppm.tol = ms1.match.ppm,
                mz.ppm.thr = ms2.match.ppm,
                interactive_plot = FALSE
              )
          } else{
            plot <-
              masstools::ms2_plot(
                spectrum1 = spectrum1,
                spectrum2 = spectrum2,
                spectrum1_name = x[1],
                spectrum2_name = x[4],
                ppm.tol = ms1.match.ppm,
                mz.ppm.thr = ms2.match.ppm,
                interactive_plot = FALSE
              )
            
          }
          temp_info = paste(colnames(temp_annotation_table),
                            x, sep = ":")
          
          plot <-
            plot +
            ggplot2::annotate(
              geom = "text",
              x = -Inf,
              y = Inf,
              label = paste(temp_info, collapse = "\n"),
              hjust = 0,
              vjust = 1
            )
          
          if (interactive_plot) {
            plot <-
              plotly::ggplotly(plot)
          }
          plot
        }
      )
    
    names(all_plot) <-
      paste(temp_annotation_table$variable_id,
            seq_len(nrow(temp_annotation_table)),
            sep = "_")
    return(all_plot)
  }
