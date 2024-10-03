#' Plot MS2 Matching Between Experimental and Library Spectra
#'
#' This function generates plots comparing MS2 spectra from experimental data with spectra from a reference database. It matches based on a set of user-specified parameters and provides visualizations with either static or interactive plots.
#'
#' @param object A `mass_dataset` object that contains MS2 data and annotations.
#' @param variable_id Optional. A vector of variable IDs to be plotted.
#' @param variable_index Optional. A vector of variable indices to be plotted. Either `variable_id` or `variable_index` must be provided.
#' @param database A `databaseClass` object containing reference MS2 spectra for matching.
#' @param polarity Character. The ionization mode, either \code{"positive"} or \code{"negative"}. Default is \code{"positive"}.
#' @param ms1.match.ppm Numeric. The mass-to-charge (m/z) tolerance in parts per million (ppm) for MS1 matching. Default is 25.
#' @param ms2.match.ppm Numeric. The mass-to-charge (m/z) tolerance in parts per million (ppm) for MS2 matching. Default is 30.
#' @param mz.ppm.thr Numeric. The m/z threshold for fragment matching. Default is 400.
#' @param remove_fragment_intensity_cutoff Numeric. The intensity cutoff for removing low-intensity fragments. Default is 0.
#' @param interactive_plot Logical. If \code{TRUE}, an interactive plot using plotly is returned. Default is \code{FALSE}.
#'
#' @details
#' The function compares MS2 spectra from the experimental data (provided by the `mass_dataset` object) to those in the reference `database`. It generates matching plots for each of the specified variables and annotations. Users can specify either the `variable_id` or `variable_index` for plotting. The function supports both static and interactive visualizations.
#'
#' @return A list of plots (either static or interactive) comparing experimental and reference MS2 spectra.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' all_plots <- plot_ms2_matching(
#'   object = my_mass_dataset,
#'   variable_id = c("pRPLC_1112", "pRPLC_1860"),
#'   database = my_database
#' )
#' }
#'
#' @export

plot_ms2_matching <-
  function(object,
           variable_id,
           variable_index,
           database,
           polarity,
           ms1.match.ppm = 25,
           ms2.match.ppm = 30,
           mz.ppm.thr = 400,
           remove_fragment_intensity_cutoff = 0,
           interactive_plot = FALSE) {
    
    if(missing(polarity)){
      stop("polarity should be provided.\n")
    }
    
    if(polarity %in% c("positive", "negative")){
      polarity <- polarity
    } else{
      stop("polarity should be either 'positive' or 'negative'.\n")
    }
    
    if (!is(object, "mass_dataset")) {
      stop("object should be mass_dataset object.\n")
    }
    
    ###Check data
    if (missing(database)) {
      stop("No database is provided.\n")
    }
    
    if (!is(database, "databaseClass")) {
      stop("database should be databaseClass object.\n")
    }
    
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
    
    database.name <-
      extract_database_name(database)
    
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
          temp_info = paste(colnames(temp_annotation_table), x, sep = ":")
          
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
      paste(temp_annotation_table$variable_id, seq_len(nrow(temp_annotation_table)), sep = "_")
    return(all_plot)
  }




#' Plot MS2 Spectra for a Single Peak in a mass_dataset Object
#'
#' This function generates MS2 spectra comparison plots for a single peak in a `mass_dataset` object by comparing experimental MS2 data with reference MS2 data from a spectral database.
#'
#' @param object A `mass_dataset` object containing the peak data.
#' @param variable_id The ID of the peak to plot. Either `variable_id` or `variable_index` must be provided.
#' @param variable_index The index of the peak to plot. Either `variable_id` or `variable_index` must be provided.
#' @param polarity Character, ionization mode, either `"positive"` or `"negative"`. Defaults to `"positive"`.
#' @param ms1.match.ppm Numeric, mass accuracy threshold for MS1 matching in parts per million (ppm). Defaults to `25`.
#' @param ms2.match.ppm Numeric, mass accuracy threshold for MS2 matching in ppm. Defaults to `30`.
#' @param mz.ppm.thr Numeric, m/z threshold in ppm for matching MS1 and MS2. Defaults to `400`.
#' @param database A `databaseClass` object containing the reference spectral database for MS2 data.
#' @param interactive_plot Logical, if `TRUE`, generates an interactive plot using `plotly`. Defaults to `FALSE`.
#'
#' @return A list of MS2 spectra comparison plots for the specified peak, with one plot per matched annotation. If `interactive_plot = TRUE`, the plots are returned as interactive `plotly` plots.
#'
#' @details
#' This function retrieves the MS2 spectra for a specified peak and compares them to the reference MS2 spectra from a provided database. It generates a plot for each matched annotation, showing the experimental spectrum and the reference spectrum side by side.
#'
#' @examples
#' \dontrun{
#' # Plot MS2 spectra for a peak
#' ms2_plots <- ms2_plot_mass_dataset(
#'   object = mass_object,
#'   variable_id = "P001",
#'   database = reference_database
#' )
#' }
#'
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @export

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
    polarity <- match.arg(polarity)
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
          temp_info = paste(colnames(temp_annotation_table), x, sep = ":")
          
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
      paste(temp_annotation_table$variable_id, seq_len(nrow(temp_annotation_table)), sep = "_")
    return(all_plot)
  }
