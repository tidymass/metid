##------------------------------------------------------------------------------
#' @title get_identification_table_all
#' @description Extract the identifications from multiple results of `identify_metabolite_all()`.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param ... One or multiple results from `identify_metabolite_all()`.
#' @param candidate.num candidate.num
#' @param level_condition Condition for level assign.
#' @return A identification table (data.frame).
#' @export
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom magrittr %>%
#' @importFrom stringr str_split
#' @importFrom purrr map map2
#' @importFrom crayon green red yellow
#' @importFrom dplyr select mutate everything filter left_join
#' @seealso The example and demo data of this function can be found
#' \url{https://tidymass.github.io/metid/articles/metid.html}

get_identification_table_all =
  function(...,
           candidate.num = 1,
           level_condition = c(
             "mz&rt&ms2" = 1,
             "mz&rt" = 2,
             "mz&ms2" = 2,
             "mz" = 3
           )) {
    options(warn = -1)
    result = list(...)
    ##rename result
    result =
      purrr::map2(
        .x = result,
        .y = seq_along(result),
        .f = function(x, y) {
          names(x) = paste(y, names(x), sep = "_")
          x
        }
      )
    
    result = unlist(result)
    
    database_level =
      result %>%
      purrr::map(function(y) {
        database = y@database
        mz_match = ifelse(is.na(y@identification.result[[1]]$mz.error[1]), "", "mz")
        rt_match = ifelse(is.na(y@identification.result[[1]]$RT.error[1]), "", "rt")
        ms2_match = ifelse(all(y@identification.result[[1]]$SS == 0), "", "ms2")
        final = paste(mz_match, rt_match, ms2_match, sep = "&")
        final = stringr::str_split(string = final, pattern = "&")[[1]]
        final = paste(final[final != ""], collapse = "&")
        level = level_condition[match(final, names(level_condition))]
        c(database, level)
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()

    colnames(database_level)[seq_len(2)] = c("database", "level")
    database_level =
      database_level %>% 
      tibble::rownames_to_column(var = "name")
    
    database_level$level[is.na(database_level$level)] = 5
    
    ####annotation table level 1
    if (any(database_level$level == 1)) {
      message(crayon::green("Level 1 table..."))
      annotation_table1 <-
        get_identification_table(result[which(database_level$level == 1)],
                                 type = "new",
                                 candidate.num = candidate.num)
      annotation_table1 <-
        data.frame(annotation_table1,
                   Level = 1,
                   stringsAsFactors = FALSE)
      message(crayon::red("OK"))
    } else{
      message(crayon::yellow("No level 1."))
      annotation_table1 = NULL
    }
    
    ####annotation table level 2
    if (any(database_level$level == 2)) {
      message(crayon::green("Level 2 table..."))
      annotation_table2 <-
        get_identification_table(result[which(database_level$level == 2)],
                                 type = "new",
                                 candidate.num = candidate.num)
      annotation_table2 <-
        data.frame(annotation_table2,
                   Level = 2,
                   stringsAsFactors = FALSE)
      message(crayon::red("OK"))
    } else{
      message(crayon::yellow("No level 2."))
      annotation_table2 = NULL
    }
    
    ####annotation table level 3
    if (any(database_level$level == 3)) {
      message(crayon::green("Level 3 table..."))
      annotation_table3 <-
        get_identification_table(result[which(database_level$level == 3)],
                                 type = "new",
                                 candidate.num = candidate.num)
      annotation_table3 <-
        data.frame(annotation_table3,
                   Level = 3,
                   stringsAsFactors = FALSE)
      message(crayon::red("OK"))
    } else{
      message(crayon::yellow("No level 3."))
      annotation_table3 = NULL
    }
    
    ###remove the redundant result
    if(!is.null(annotation_table1)){
      annotation_table1 =
        annotation_table1 %>%
        dplyr::filter(!is.na(Compound.name))
    }
    
    if(!is.null(annotation_table2)){
      annotation_table2 =
        annotation_table2 %>%
        dplyr::filter(!is.na(Compound.name))
      if(!is.null(annotation_table1)){
        annotation_table2 =
          annotation_table2 %>%
          dplyr::filter(!name %in% annotation_table1$name)
      }
    }
    
    if(!is.null(annotation_table3)){
      annotation_table3 =
        annotation_table3 %>%
        dplyr::filter(!is.na(Compound.name))
      
      if(!is.null(annotation_table1)){
        annotation_table3 =
          annotation_table3 %>%
          dplyr::filter(!name %in% annotation_table1$name)
      }
      
      if(!is.null(annotation_table2)){
        annotation_table3 =
          annotation_table3 %>%
          dplyr::filter(!name %in% annotation_table2$name)
      }
    }
    
    annotation_table =
      rbind(annotation_table1,
            annotation_table2,
            annotation_table3)

    ms1_peak = result[[1]]@ms1.data

    diff_name = 
    setdiff(ms1_peak$name, annotation_table$name)
    
    if(length(diff_name) > 0){
      new_matrix = 
        ms1_peak[match(diff_name, ms1_peak$name),]
      new_matrix = 
      new_matrix %>% 
        dplyr::mutate(MS2.spectra.name = NA,
                      Candidate.number = NA, 
                      Compound.name = NA,   
                      CAS.ID = NA,
                      HMDB.ID = NA,
                      KEGG.ID = NA,        
                      Lab.ID = NA,
                      Adduct = NA,
                      mz.error = NA,        
                      mz.match.score = NA,
                      RT.error = NA,
                      RT.match.score = NA,  
                      CE = NA,
                      SS = NA,
                      Total.score = NA,     
                      Database = NA,
                      Level = NA
                      )
      annotation_table = 
        rbind(annotation_table, 
              new_matrix) %>% 
        as.data.frame()
    }
    
    ms1_peak_name = colnames(ms1_peak)
    ms1_peak_name = ms1_peak_name[ms1_peak_name != "name"]
  
    annotation_table = 
    annotation_table %>% 
      dplyr::select(-ms1_peak_name) %>% 
      dplyr::left_join(ms1_peak, by = "name") %>% 
      dplyr::select(name, ms1_peak_name, dplyr::everything())

    if(candidate.num == 1){
      annotation_table = 
      annotation_table[match(result[[1]]@ms1.data$name, annotation_table$name),]
    }
    
    message(crayon::red("All done."))
    return(tibble::as_tibble(annotation_table))
  }