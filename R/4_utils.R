get_extension <- function(file) {
  tail(stringr::str_split(string = file, pattern = "\\.")[[1]], 1)
}

readTable <- function(file, ...) {
  extension <- get_extension(file = file)
  if (extension == "csv") {
    return(readr::read_csv(file = file, show_col_types = FALSE, ...))
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
    message(crayon::red("file is not csv, xlsx or xls."))
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
  
  if (isTRUE(theme$dark))
    crayon::white(x)
  else
    crayon::black(x)
  
}

#' List all packages in the metid
#'
#' @param include_self Include metid in the list?
#' @export
#' @return metid_packages
#' @examples
#' metid_packages()
metid_packages <- function(include_self = TRUE) {
  raw <- utils::packageDescription("metid")$Imports
  imports <- strsplit(raw, ",")[[1]]
  parsed <- gsub("^\\s+|\\s+$", "", imports)
  names <-
    vapply(strsplit(parsed, "\\s+"), "[[", 1, FUN.VALUE = character(1))
  
  if (include_self) {
    names <- c(names, "metid")
  }
  
  names
}

invert <- function(x) {
  if (length(x) == 0)
    return()
  stacked <- utils::stack(x)
  tapply(as.character(stacked$ind), stacked$values, list)
}


style_grey <- function(level, ...) {
  crayon::style(paste0(...),
                crayon::make_style(grDevices::grey(level), grey = TRUE))
}





#' Summary of Annotation Table with Plot
#'
#' This function generates a summary of an annotation table from a `mass_dataset` object or a 
#' user-provided object. The function filters annotation data based on a specified `Level` and
#' creates a polar bar plot to visualize the distribution of annotations across levels and databases.
#'
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @param object Either a `mass_dataset` object or a data frame containing annotation information.
#' @param level A numeric vector specifying the annotation levels to include in the summary. Defaults to levels 1, 2, 3, and 4.
#'
#' @return A ggplot2 object representing a polar bar plot of the annotation summary, 
#' showing the number and percentage of annotations by level and database.
#'
#' @details
#' If the input is a `mass_dataset` object, the function extracts the annotation information 
#' using the `massdataset::extract_variable_info` function. For non-`mass_dataset` objects, 
#' it expects the input to be in a specific format and selects the most relevant annotations 
#' based on the `Level`, `SS`, and `Total.score` columns.
#'
#' The function then filters out annotations with missing `Compound.name` and 
#' those that are not in the specified `level`. It calculates the number and percentage 
#' of annotations by level and database, which are visualized in a polar bar plot.
#'
#' @examples
#' \dontrun{
#' # For a mass_dataset object:
#' summary_annotation_table(mass_object, level = c(1, 2, 3))
#'
#' # For a custom annotation table:
#' custom_annotation <- data.frame(variable_id = ..., Level = ..., SS = ..., Total.score = ..., Compound.name = ...)
#' summary_annotation_table(custom_annotation, level = c(1, 2))
#' }
#'
#' @importFrom plyr . dlply
#' @importFrom purrr map
#' @importFrom dplyr filter select mutate group_by distinct ungroup arrange
#' @importFrom dplyr pull
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual theme_void guides
#' @importFrom ggplot2 guide_legend coord_polar
#' @importFrom massdataset extract_variable_info
#'
#' @seealso \code{\link{massdataset::extract_variable_info}}
#'
#' @export

summary_annotation_table <-
  function(object,
           level = c(1, 2, 3, 4)) {
    if (is(object, class2 = "mass_dataset")) {
      annotation_table <-
        massdataset::extract_variable_info(object)
    } else{
      annotation_table <-
        object %>%
        plyr::dlply(.variables = .(variable_id)) %>%
        purrr::map(function(y) {
          if (nrow(y) == 1) {
            return(y)
          } else{
            y %>%
              dplyr::filter(Level == min(Level)) %>%
              dplyr::filter(SS == max(SS, na.rm = TRUE)) %>%
              dplyr::filter(Total.score == max(Total.score, na.rm = TRUE)) %>%
              head(1)
          }
        }) %>%
        dplyr::bind_rows()
    }
    
    annotation_table <-
      annotation_table %>%
      dplyr::filter(!is.na(Compound.name)) %>%
      dplyr::filter(Level %in% level)
    
    temp_data <-
      rbind(
        annotation_table %>%
          dplyr::select(fill = Level) %>%
          dplyr::mutate(class = "level") %>%
          dplyr::group_by(fill) %>%
          dplyr::mutate(n = dplyr::n()) %>%
          dplyr::distinct() %>%
          dplyr::mutate(fill = as.character(fill)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(rate = round(n * 100 / sum(n), 2)),
        annotation_table %>%
          dplyr::select(fill = Database) %>%
          dplyr::mutate(class = "database") %>%
          dplyr::group_by(fill) %>%
          dplyr::mutate(n = dplyr::n()) %>%
          dplyr::distinct() %>%
          dplyr::ungroup() %>%
          dplyr::mutate(rate = round(n * 100 / sum(n), 2))
      ) %>%
      dplyr::mutate(class = factor(class, levels = c("level", "database")))
    
    class_color <-
      annotation_table %>%
      dplyr::select(Level, Database) %>%
      dplyr::distinct() %>%
      dplyr::arrange(Level)
    
    temp_data$fill <-
      factor(temp_data$fill,
             levels = unique(c(
               class_color$Level, class_color$Database
             )))
    
    temp_color <-
      c(
        "#3B4992FF",
        "#EE0000FF",
        "#008B45FF",
        "#631879FF",
        "#008280FF",
        "#BB0021FF",
        "#5F559BFF",
        "#A20056FF",
        "#808180FF",
        "#1B1919FF"
      )
    
    level_color <-
      temp_color[1:length(unique(class_color$Level))]
    
    names(level_color) <-
      unique(class_color$Level)
    
    database_color <-
      purrr::map(unique(class_color$Level), function(x) {
        temp =
          class_color %>%
          dplyr::filter(Level == x) %>%
          dplyr::pull(Database)
        alpha = seq(1, 0.1, length.out = length(temp))
        database_color = purrr::map(alpha, function(y) {
          ggplot2::alpha(level_color[x], y)
        }) %>%
          unlist()
        names(database_color) = temp
        database_color
      }) %>%
      unlist()
    
    temp_data$new_name <-
      paste(temp_data$fill,
            " (",
            temp_data$n,
            ",",
            temp_data$rate,
            "%)",
            sep = "")
    
    names(level_color) <-
      temp_data$new_name[match(names(level_color), temp_data$fill)]
    
    names(database_color) <-
      temp_data$new_name[match(names(database_color), temp_data$fill)]
    
    new_level <-
      temp_data$new_name[match(levels(temp_data$fill), temp_data$fill)]
    
    temp_data$fill <-
      temp_data$new_name
    
    temp_data$fill <-
      factor(temp_data$fill, levels = new_level)
    
    ggplot(temp_data,
           aes(x = class, y = n, fill = fill)) +
      geom_col(color = "black",
               position = 'stack') +
      # scale_x_discrete(limits = c(" ", "level","database")) +
      scale_fill_manual(values = c(level_color, database_color)) +
      theme_void() +
      guides(fill = guide_legend(title = "", ncol = 1)) +
      coord_polar("y")
  }
