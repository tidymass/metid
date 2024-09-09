

#' Validate Parameters for Metabolite Annotation
#'
#' This function validates the input parameters required for the metabolite annotation process. It checks the `mass_dataset` object, the annotation database, the adduct table, and various matching and weighting parameters to ensure they are correctly specified.
#'
#' @param object A `mass_dataset` object containing the MS1, RT, and/or MS2 data for annotation.
#' @param database A `databaseClass` object containing the reference database for metabolite annotation.
#' @param based_on Character vector. Specifies the criteria to base the annotation on. Can include `"ms1"`, `"rt"`, and/or `"ms2"`. Default is `c("ms1", "rt", "ms2")`.
#' @param polarity Character. The ionization mode, either `"positive"` or `"negative"`. Default is `"positive"`.
#' @param column Character. The chromatographic column type, either `"hilic"` or `"rp"` (reversed-phase). Default is `"hilic"`.
#' @param adduct.table A data frame containing the adduct table to be used in the annotation process.
#' @param ce Character or vector. Collision energy values to be used in MS2 matching. Default is `"all"`.
#' @param ms1.match.ppm Numeric. The mass tolerance in parts per million (ppm) for MS1 peak matching. Default is 25.
#' @param ms2.match.ppm Numeric. The mass tolerance in ppm for MS2 peak matching. Default is 30.
#' @param mz.ppm.thr Numeric. m/z threshold for ppm calculation. Default is 400.
#' @param ms2.match.tol Numeric. The score tolerance for MS2 matches. Default is 0.5.
#' @param fraction.weight Numeric. Weight for the fraction of matched fragments in the total MS2 score calculation. Default is 0.3.
#' @param dp.forward.weight Numeric. Weight for the forward dot product score in MS2 matching. Default is 0.6.
#' @param dp.reverse.weight Numeric. Weight for the reverse dot product score in MS2 matching. Default is 0.1.
#' @param rt.match.tol Numeric. Retention time matching tolerance in seconds. Default is 30.
#' @param ms1.match.weight Numeric. Weight for MS1 matching in the total score. Must be between 0 and 1. Default is 0.25.
#' @param rt.match.weight Numeric. Weight for RT matching in the total score. Must be between 0 and 1. Default is 0.25.
#' @param ms2.match.weight Numeric. Weight for MS2 matching in the total score. Must be between 0 and 1. Default is 0.5.
#' @param total.score.tol Numeric. The threshold for the total score. Must be between 0 and 1. Default is 0.5.
#' @param candidate.num Numeric. The number of top candidate annotations to retain. Must be greater than 0. Default is 3.
#' @param remove_fragment_intensity_cutoff Numeric. The intensity cutoff for removing low-intensity MS2 fragments. Default is 0.
#' @param threads Numeric. The number of threads to use for parallel processing. Must be greater than 0. Default is 3.
#'
#' @return The function does not return a value but throws an error if any validation checks fail.
#'
#' @details
#' The function checks the following:
#' * Ensures the `object` is a valid `mass_dataset` object and that the `database` is a valid `databaseClass` object.
#' * Verifies the presence of required MS1, RT, and MS2 data based on the `based_on` argument.
#' * Validates the adduct table and checks if collision energy (CE) values match with those in the database.
#' * Ensures that the matching parameters (`ms1.match.ppm`, `ms2.match.ppm`, `mz.ppm.thr`, `ms2.match.tol`) are numeric and within valid ranges.
#' * Validates the weights (`ms1.match.weight`, `rt.match.weight`, `ms2.match.weight`) and ensures they sum to 1.
#' * Ensures that the number of candidates to retain (`candidate.num`) and the number of threads are greater than 0.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' check_parameters4annotate_metabolites(
#'   object = my_dataset,
#'   database = my_database,
#'   based_on = c("ms1", "ms2"),
#'   polarity = "positive",
#'   column = "rp",
#'   ms1.match.ppm = 20,
#'   ms2.match.ppm = 25,
#'   ms1.match.weight = 0.3,
#'   rt.match.weight = 0.3,
#'   ms2.match.weight = 0.4,
#'   candidate.num = 5,
#'   threads = 4
#' )
#' }
#'
#'
#' @export


check_parameters4annotate_metabolites <-
  function(object,
           database,
           based_on = c("ms1", "rt", "ms2"),
           polarity = c("positive", "negative"),
           column = c("hilic", "rp"),
           adduct.table = adduct.table,
           ce = "all",
           ms1.match.ppm = 25,
           ms2.match.ppm = 30,
           mz.ppm.thr = 400,
           ms2.match.tol = 0.5,
           fraction.weight = 0.3,
           dp.forward.weight = 0.6,
           dp.reverse.weight = 0.1,
           rt.match.tol = 30,
           ms1.match.weight = 0.25,
           rt.match.weight = 0.25,
           ms2.match.weight = 0.5,
           total.score.tol = 0.5,
           candidate.num = 3,
           remove_fragment_intensity_cutoff = 0,
           threads = 3) {
    based_on <-
      match.arg(based_on, several.ok = TRUE)
    
    ###object
    check_mass_dataset(object = object, based_on = based_on)
    
    ###database
    check_database(database = database, based_on = based_on)
    
    ###check adduct.table
    check_adduct_table(adduct.table)
    
    polarity <- match.arg(polarity)
    column <- match.arg(column)
    
    ##CE values
    ## Check if CE values (collision energies) match with the database
    
    ###ms1.match.ppm
    if (!is.numeric(ms1.match.ppm)) {
      stop("ms1.match.ppm should be numeric.\n")
    } else{
      if (ms1.match.ppm <= 0 | ms1.match.ppm >= 500) {
        stop("ms1.match.ppm should > 0 and < 500\n")
      }
    }
    
    # ms2.match.ppm
    if (!is.numeric(ms2.match.ppm)) {
      stop("ms2.match.ppm should be numeric.\n")
    } else{
      if (ms2.match.ppm <= 0 | ms2.match.ppm >= 500) {
        stop("ms2.match.ppm should > 0 and < 500\n")
      }
    }
    
    
    # ms1.match.weight
    if (!is.numeric(ms1.match.weight)) {
      stop("ms1.match.weight should be numeric.\n")
    } else{
      if (ms1.match.weight < 0 | ms1.match.weight > 1) {
        stop("ms1.match.weight should > 0 and < 1\n")
      }
    }
    
    # rt.match.weight
    if (!is.numeric(rt.match.weight)) {
      stop("rt.match.weight should be numeric.\n")
    } else{
      if (rt.match.weight < 0 | rt.match.weight > 1) {
        stop("rt.match.weight should > 0 and < 1\n")
      }
    }
    
    # ms2.match.weight
    if (!is.numeric(ms2.match.weight)) {
      stop("ms2.match.weight should be numeric.\n")
    } else{
      if (ms2.match.weight < 0 | ms2.match.weight > 1) {
        stop("ms2.match.weight should > 0 and < 1\n")
      }
    }
    
    if (ms1.match.weight + rt.match.weight + ms2.match.weight != 1) {
      stop("ms1.match.weight + rt.match.weight + ms2.match.weight should be 1.\n")
    }
    
    # total.score.tol
    if (!is.numeric(total.score.tol)) {
      stop("total.score.tol should be numeric.\n")
    } else{
      if (total.score.tol < 0 | total.score.tol > 1) {
        stop("total.score.tol should > 0 and < 1\n")
      }
    }
    
    ###candidate.num
    if (!is.numeric(candidate.num)) {
      stop("candidate.num should be numeric.\n")
    } else{
      if (candidate.num <= 0) {
        stop("candidate.num should > 0.\n")
      }
    }
    
    ###threads
    if (!is.numeric(threads)) {
      stop("threads should be numeric.\n")
    } else{
      if (threads <= 0) {
        stop("threads should > 0.\n")
      }
    }
    
  }


#' Validate the mass_dataset Object for Metabolite Annotation
#'
#' This function checks if the provided `mass_dataset` object is valid for metabolite annotation based on the specified criteria (`ms1`, `rt`, and/or `ms2`). It ensures that the object contains the necessary data for the chosen annotation method.
#'
#' @param object A `mass_dataset` object containing MS1, RT, and/or MS2 data.
#' @param based_on Character vector. Specifies the criteria to base the validation on. Can include `"ms1"`, `"rt"`, and `"ms2"`. Default is `c("ms1", "rt", "ms2")`.
#'
#' @return The function does not return a value but throws an error if any validation checks fail.
#'
#' @details
#' The function checks the following:
#' * Ensures that `object` is a valid `mass_dataset` object.
#' * If `"ms2"` is included in `based_on`, it checks that the `mass_dataset` object contains MS2 data (`object@ms2_data`). If no MS2 data is present, the function throws an error.
#'
#' If any of these conditions are not met, the function will stop with an appropriate error message.
#'
#' @examples
#' \dontrun{
#' # Load a sample mass dataset
#' my_data <- load_mass_dataset("path/to/dataset")
#'
#' # Validate the dataset for MS1 and MS2 information
#' check_mass_dataset(object = my_data, based_on = c("ms1", "ms2"))
#' }
#'
#'
#' @export


check_mass_dataset <-
  function(object, based_on = c("ms1", "rt", "ms2")) {
    based_on <- match.arg(based_on, several.ok = TRUE)
    if (!is(object, "mass_dataset")) {
      stop("object should be mass_dataset object.\n")
    }
    
    if ("ms2" %in% based_on) {
      if (length(object@ms2_data) == 0) {
        stop("The object should have MS2 information.\n")
      }
    }
  }


#' Validate a Database Object for Metabolite Annotation
#'
#' This function validates a `databaseClass` object to ensure it contains the necessary information for metabolite annotation based on specified criteria. It checks for the presence of retention time (RT) information and MS2 spectra, verifies collision energy (CE) values, and ensures that the database structure aligns with the selected polarity and annotation basis.
#'
#' @param database A `databaseClass` object containing MS2 spectra and related information.
#' @param polarity Character. The ionization mode, either `"positive"` or `"negative"`. Default is `"positive"`.
#' @param ce Character or vector. The collision energy (CE) values to validate against the database. Set to `"all"` to include all CE values present in the database. Default is `"all"`.
#' @param based_on Character vector. Specifies the criteria for validation, which can include `"ms1"`, `"rt"`, and `"ms2"`. Multiple criteria can be provided. Default is `c("ms1", "rt", "ms2")`.
#'
#' @return The function does not return a value but throws an error if any validation checks fail.
#'
#' @details
#' The function performs the following validations on the `database` object:
#'
#' * **Polarity Check**: Ensures that the `polarity` argument is either `"positive"` or `"negative"`.
#' * **Database Class Check**: Confirms that the `database` object is of class `databaseClass`.
#' * **Retention Time (RT) Validation**:
#'   * If `"rt"` is included in the `based_on` argument, the function checks that the `database@spectra.info` data frame contains an `RT` column.
#'   * Ensures that the `RT` column does not consist entirely of `NA` values.
#' * **MS2 Information Validation**:
#'   * If `"ms2"` is included in the `based_on` argument, the function verifies that both the `database` object and the associated data contain MS2 information.
#'   * Checks that the `database@spectra.data` contains MS2 spectra corresponding to the specified `polarity`.
#'   * Validates that the provided `ce` values are present in the database's collision energy list unless `"all"` is specified.
#'
#' If any of these conditions are not met, the function will terminate execution and provide an informative error message indicating the nature of the issue.
#'
#' @examples
#' \dontrun{
#' # Assuming you have a databaseClass object named my_database
#' # Validate the database with default settings
#' check_database(database = my_database)
#'
#' # Validate the database for negative polarity and specific CE values
#' check_database(
#'   database = my_database,
#'   polarity = "negative",
#'   ce = c("10", "20"),
#'   based_on = c("ms1", "ms2")
#' )
#'
#' # Validate the database using only MS1 and RT criteria
#' check_database(
#'   database = my_database,
#'   based_on = c("ms1", "rt")
#' )
#' }
#'
#'
#' @export


check_database <-
  function(database,
           polarity = c("positive", "negative"),
           ce = "all",
           based_on = c("ms1", "rt", "ms2")) {
    polarity <- match.arg(polarity)
    based_on <- match.arg(based_on, several.ok = TRUE)
    
    if (!is(database, "databaseClass")) {
      stop("database should be a databaseClass object.\n")
    }
    
    ###check RT information
    if ("rt" %in% based_on) {
      if (!"RT" %in% colnames(database@spectra.info)) {
        stop("The database should have RT information.\n")
      } else{
        if (all(is.na(database@spectra.info$RT))) {
          stop("The database RT information should not all be NA.\n")
        }
      }
    }
    
    ##check MS2 information
    ####if based_on has ms2, the object and database should have ms2 information
    if ("ms2" %in% based_on) {
      ##check if object have ms2 information
      if (length(object@ms2_data) == 0) {
        stop("The object should have MS2 information.\n")
      }
      
      ##Check if database have MS2 spectra or not
      spectra_pos_number <-
        database@spectra.data[['Spectra.positive']] %>%
        length()
      
      spectra_neg_number <-
        database@spectra.data[['Spectra.negative']] %>%
        length()
      
      spectra_number <-
        ifelse(polarity == "positive",
               spectra_pos_number,
               spectra_neg_number)
      
      if (spectra_number == 0) {
        stop("The database should have MS2 spectra information.\n")
      }
      
      ##CE values
      ## Check if CE values (collision energies) match with the database
      if (spectra_number > 0) {
        ce.list.pos <- unique(unlist(lapply(
          database@spectra.data$Spectra.positive, names
        )))
        ce.list.neg <- unique(unlist(lapply(
          database@spectra.data$Spectra.negative, names
        )))
        
        ce.list <- ifelse(polarity == "positive", ce.list.pos, ce.list.neg)
        
        ## Stop execution if provided CE values are not in the database
        if (all(ce %in% ce.list) & !("all" %in% ce)) {
          message(
            "The CE values you can chose are: ",
            paste(ce.list, collapse = ", "),
            "\n",
            "or just set it as ",
            "'all'"
          )
          stop("All ce values you set are not in the database. Please check it.\n")
        }
      }
      
    }
    
  }


#' Validate the Adduct Table
#'
#' This function checks whether the input adduct table is valid by ensuring it is a data frame and contains the required columns (`adduct` and `mz`). It also validates the data types of the columns and checks for any missing values.
#'
#' @param adduct.table A data frame containing information about adducts. The table must have two columns: `adduct` (character) and `mz` (numeric).
#'
#' @return The function does not return a value but throws an error if the adduct table is invalid.
#'
#' @details
#' The function performs several checks on the `adduct.table`:
#' * Ensures that the input is a data frame.
#' * Confirms the presence of the required columns: `adduct` (character) and `mz` (numeric).
#' * Verifies that the `adduct` column is of character type and the `mz` column is numeric.
#' * Checks for missing values (`NA`) in either column.
#'
#' If any of these conditions are not met, the function throws an error with a message explaining the issue.
#'
#' @examples
#' \dontrun{
#' # Example of a valid adduct table
#' adduct_table <- data.frame(
#'   adduct = c("(M+H)+", "(M+Na)+"),
#'   mz = c(1.0073, 22.9893)
#' )
#'
#' # Validate the adduct table
#' check_adduct_table(adduct_table)
#'
#' # Example of an invalid adduct table (missing 'mz' column)
#' invalid_table <- data.frame(
#'   adduct = c("(M+H)+", "(M+Na)+")
#' )
#'
#' check_adduct_table(invalid_table) # This will throw an error
#' }
#'
#' @export


check_adduct_table <-
  function(adduct.table) {
    # Check if adduct.table is a data frame
    if (!is.data.frame(adduct.table)) {
      stop("The input adduct.table is not a data frame.")
    }
    
    # Check if the required columns are present
    required_columns <- c("adduct", "mz")
    if (!all(required_columns %in% colnames(adduct.table))) {
      stop("The adduct.table must have two columns: 'adduct' and 'mz'.")
    }
    
    # Check if the 'adduct' column is of character type
    if (!is.character(adduct.table$adduct)) {
      stop("The 'adduct' column must be of type character.")
    }
    
    # Check if the 'mz' column is of numeric type
    if (!is.numeric(adduct.table$mz)) {
      stop("The 'mz' column must be of type numeric.")
    }
    
    # Check if there are any missing values in either column
    if (any(is.na(adduct.table$adduct))) {
      stop("The 'adduct' column contains missing values.")
    }
    
    if (any(is.na(adduct.table$mz))) {
      stop("The 'mz' column contains missing values.")
    }
    
  }


#' Validate MS1 and MS2 Information
#'
#' This function checks the integrity and validity of MS1 and MS2 data used for metabolite annotation. It ensures that the required columns and formats are present and correctly typed, and that the MS1 and MS2 information is consistent with the specified matching criteria.
#'
#' @param ms1.info A data frame containing MS1 peak information, including `mz`, `rt`, and `variable_id`. If `based_on` includes `"ms2"`, it should also include `ms2_spectrum_id` and `ms2_files_id`.
#' @param ms2.info A list of MS2 spectra, where each element is a matrix with two columns (`mz` and `intensity`). The names of `ms2.info` should match `ms2_spectrum_id` in `ms1.info`.
#' @param based_on Character vector. Specifies the criteria to be used for validation, which can include `"ms1"`, `"rt"`, and `"ms2"`. Default is `c("ms1", "rt", "ms2")`.
#'
#' @return The function does not return a value but throws an error if any validation checks fail.
#'
#' @details
#' The function performs the following checks:
#' * Ensures that `ms1.info` is a data frame and contains the required columns, depending on the `based_on` argument. For `"ms1"` and `"rt"`, it checks for the presence of `mz`, `rt`, and `variable_id`. For `"ms2"`, it checks for `ms2_spectrum_id` and `ms2_files_id`.
#' * Validates the data types of each column in `ms1.info`, ensuring that `mz` and `rt` are numeric and `variable_id`, `ms2_spectrum_id`, and `ms2_files_id` are character strings.
#' * Verifies that `ms2.info` is a list where each element is a matrix with exactly two columns (`mz` and `intensity`), and that both columns contain numeric values.
#' * Checks for missing values (`NA`) in the required columns.
#' * Ensures that the names of `ms2.info` match the `ms2_spectrum_id` column in `ms1.info`, if `ms2` data is being used.
#'
#' @examples
#' \dontrun{
#' # Example MS1 info data frame
#' ms1_info <- data.frame(
#'   variable_id = c("id1", "id2"),
#'   mz = c(150.08, 180.12),
#'   rt = c(12.5, 14.7),
#'   ms2_spectrum_id = c("spectrum1", "spectrum2"),
#'   ms2_files_id = c("file1", "file2")
#' )
#'
#' # Example MS2 info
#' ms2_info <- list(
#'   spectrum1 = matrix(c(75, 1000, 80, 2000), 
#'   ncol = 2, byrow = TRUE, 
#'   dimnames = list(NULL, c("mz", "intensity"))),
#'   spectrum2 = matrix(c(85, 3000, 90, 1500), 
#'   ncol = 2, byrow = TRUE, 
#'   dimnames = list(NULL, c("mz", "intensity")))
#' )
#'
#' # Validate MS1 and MS2 info
#' check_ms1_ms2_info(ms1_info, ms2_info, based_on = c("ms1", "rt", "ms2"))
#' }
#'
#'
#' @export

check_ms1_ms2_info <-
  function(ms1.info = NULL,
           ms2.info = NULL,
           based_on = c("ms1", "rt", "ms2")) {
    # Ensure based_on is a valid vector of arguments
    based_on <- match.arg(based_on, several.ok = TRUE)
    
    # If "ms1" or "rt" is in based_on, ms1.info must not be NULL and must be checked
    if ("ms1" %in% based_on || "rt" %in% based_on) {
      if (is.null(ms1.info)) {
        stop("ms1.info is required when 'ms1' or 'rt' is specified in based_on.")
      }
      
      # Check if ms1.info is a data frame
      if (!is.data.frame(ms1.info)) {
        stop("ms1.info should be a data frame.")
      }
      
      # Ensure required columns are present
      required_columns <- c("mz", "rt", "variable_id")
      if ("ms2" %in% based_on) {
        required_columns <- c(required_columns, "ms2_spectrum_id", "ms2_files_id")
      }
      
      if (!all(required_columns %in% colnames(ms1.info))) {
        stop(paste(
          "ms1.info must contain the following columns:",
          paste(required_columns, collapse = ", ")
        ))
      }
      
      # Check for missing values (NA) in required columns
      if (any(is.na(ms1.info[required_columns]))) {
        stop("ms1.info contains missing values (NA) in required columns.")
      }
      
      # Check column types
      if (!is.numeric(ms1.info$mz)) {
        stop("The 'mz' column in ms1.info must be numeric.")
      }
      if (!is.numeric(ms1.info$rt)) {
        stop("The 'rt' column in ms1.info must be numeric.")
      }
      if (!is.character(ms1.info$variable_id)) {
        stop("The 'variable_id' column in ms1.info must be a character.")
      }
      
      if ("ms2" %in% based_on) {
        if (!is.character(ms1.info$ms2_spectrum_id)) {
          stop("The 'ms2_spectrum_id' column in ms1.info must be a character.")
        }
        if (!is.character(ms1.info$ms2_files_id)) {
          stop("The 'ms2_files_id' column in ms1.info must be a character.")
        }
      }
    }
    
    # If "ms2" is in based_on, ms2.info must not be NULL and must be checked
    if ("ms2" %in% based_on) {
      if (is.null(ms2.info)) {
        stop("ms2.info is required when 'ms2' is specified in based_on.")
      }
      
      # Check that ms2.info is a list
      if (!is.list(ms2.info)) {
        stop("ms2.info should be a list.")
      }
      
      # Check that all elements of ms2.info are matrices with two columns (mz, intensity)
      for (name in names(ms2.info)) {
        spectrum <- ms2.info[[name]]
        
        if (!is.matrix(spectrum)) {
          stop(
            paste(
              "Each element of ms2.info must be a matrix. Element '",
              name,
              "' is not a matrix."
            )
          )
        }
        
        if (ncol(spectrum) != 2) {
          stop(
            paste(
              "Each matrix in ms2.info must have exactly two columns ('mz' and 'intensity'). Element '",
              name,
              "' does not have two columns."
            )
          )
        }
        
        if (!all(colnames(spectrum) == c("mz", "intensity"))) {
          stop(
            paste(
              "The matrix in element '",
              name,
              "' must have column names 'mz' and 'intensity'."
            )
          )
        }
        
        if (!is.numeric(spectrum[, "mz"]) ||
            !is.numeric(spectrum[, "intensity"])) {
          stop(
            paste(
              "Both 'mz' and 'intensity' columns in element '",
              name,
              "' must be numeric."
            )
          )
        }
        
        if (any(is.na(spectrum))) {
          stop(paste(
            "The matrix in element '",
            name,
            "' contains missing values (NA)."
          ))
        }
      }
      
      # Ensure all names of ms2.info are in ms1.info$ms2_spectrum_id
      if (!is.null(ms1.info)) {
        if (!all(names(ms2.info) %in% ms1.info$ms2_spectrum_id)) {
          stop(
            "All names of ms2.info must be present in the 'ms2_spectrum_id' column of ms1.info."
          )
        }
      }
      
    }
  }

#' Validate Parameters for Calculating Total Score
#'
#' This function validates the input parameters required for calculating the total score based on MS1, RT, and MS2 matching. It ensures the weights sum to 1 and that the matching scores are vectors of the same length.
#'
#' @param mz.match.score Numeric vector. The matching scores for MS1 peak matching (m/z).
#' @param RT.match.score Numeric vector. The matching scores for retention time (RT) matching.
#' @param SS Numeric vector. The matching scores for MS2 spectral matching (similarity score).
#' @param ms1.match.weight Numeric. Weight assigned to MS1 matching score in the total score calculation. Default is 0.25.
#' @param rt.match.weight Numeric. Weight assigned to RT matching score in the total score calculation. Default is 0.25.
#' @param ms2.match.weight Numeric. Weight assigned to MS2 matching score in the total score calculation. Default is 0.5.
#'
#' @return The function does not return a value but throws an error if any validation checks fail.
#'
#' @details
#' The function checks the following:
#' * Ensures that `mz.match.score`, `RT.match.score`, and `SS` are vectors of the same length.
#' * Checks that `ms1.match.weight`, `rt.match.weight`, and `ms2.match.weight` are numeric and that their sum equals 1.
#' * Ensures that all required arguments (`mz.match.score`, `RT.match.score`, `SS`, and the weights) are provided.
#'
#' If any of these conditions are not met, the function throws an error with a message explaining the issue.
#'
#' @examples
#' \dontrun{
#' # Example of valid inputs
#' mz_match <- c(0.9, 0.8, 0.7)
#' rt_match <- c(0.8, 0.7, 0.9)
#' ss_match <- c(0.7, 0.9, 0.8)
#'
#' # Validate the parameters
#' check_parameters4calculate_total_score(
#'   mz.match.score = mz_match,
#'   RT.match.score = rt_match,
#'   SS = ss_match,
#'   ms1.match.weight = 0.3,
#'   rt.match.weight = 0.3,
#'   ms2.match.weight = 0.4
#' )
#' }
#'
#'
#' @export

check_parameters4calculate_total_score <-
  function(mz.match.score,
           RT.match.score,
           SS,
           ms1.match.weight = 0.25,
           rt.match.weight = 0.25,
           ms2.match.weight = 0.5) {
    # Check if any argument is missing
    if (missing(mz.match.score) ||
        missing(RT.match.score) || missing(SS) ||
        missing(ms1.match.weight) ||
        missing(rt.match.weight) || missing(ms2.match.weight)) {
      stop("All arguments must be provided.")
    }
    
    # # Check if mz.match.score, RT.match.score, and SS are numeric vectors
    # if (!is.numeric(mz.match.score) ||
    #     !is.numeric(RT.match.score) || !is.numeric(SS)) {
    #   stop("mz.match.score, RT.match.score, and SS must be numeric vectors.")
    # }
    
    # Check if mz.match.score, RT.match.score, and SS are of the same length
    if (length(mz.match.score) != length(RT.match.score) ||
        length(RT.match.score) != length(SS)) {
      stop("mz.match.score, RT.match.score, and SS must be vectors of the same length.")
    }
    
    # Check if ms1.match.weight, rt.match.weight, and ms2.match.weight are numeric and sum to 1
    if (!is.numeric(ms1.match.weight) ||
        !is.numeric(rt.match.weight) ||
        !is.numeric(ms2.match.weight)) {
      stop("ms1.match.weight, rt.match.weight, and ms2.match.weight must be numeric.")
    }
    
    if (ms1.match.weight + rt.match.weight + ms2.match.weight != 1) {
      stop("ms1.match.weight, rt.match.weight, and ms2.match.weight must sum to 1.")
    }
    
  }