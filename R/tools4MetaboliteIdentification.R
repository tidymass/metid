#---------------------------------------------------------------------------
metIdentification <-
  function(ms1.info,
           ms2.info,
           polarity = c("positive", "negative"),
           ce = '30',
           database,
           ms1.match.ppm = 25,
           ms2.match.ppm = 30,
           mz.ppm.thr = 400,
           ms2.match.tol = 0.5,
           rt.match.tol = 30,
           column = "rp",
           ms1.match.weight = 0.25,
           rt.match.weight = 0.25,
           ms2.match.weight = 0.5,
           total.score.tol = 0.5,
           candidate.num = 3,
           adduct.table,
           threads = 3,
           fraction.weight = 0.3,
           dp.forward.weight = 0.6,
           dp.reverse.weight = 0.1,
           remove_fragment_intensity_cutoff = 0) {
    polarity <- match.arg(polarity)
    ms1.info$mz <- as.numeric(ms1.info$mz)
    ms1.info$rt <- as.numeric(ms1.info$rt)
    ##filter the database for using
    ##polarity
    if (polarity == "positive") {
      spectra.data <- database@spectra.data$Spectra.positive
    } else{
      spectra.data <- database@spectra.data$Spectra.negative
    }
    
    ##get the MS2 spectra within the CE values
    if (any(ce == "all")) {
      message(crayon::yellow("Use all CE values."))
      ce <- unique(unlist(lapply(spectra.data, function(x) {
        names(x)
      })))
    } else{
      spectra.data <- lapply(spectra.data, function(x) {
        x <- x[which(names(x) %in% ce)]
        if (length(x) == 0)
          return(NULL)
        return(x)
      })
    }
    
    ##remove some metabolites which have no spectra
    spectra.data <-
      spectra.data[which(!unlist(lapply(spectra.data, is.null)))]
    if (length(spectra.data) == 0) {
      stop("No spectra with CE: ",
           paste(ce, collapse = ", "),
           " in you database.\n")
    }
    
    spectra.info <- database@spectra.info
    
    spectra.info <-
      spectra.info[which(spectra.info$Lab.ID %in% names(spectra.data)),]
    
    rm(list = c("database"))
    message("\n")
    message(crayon::green('Identifing metabolites with MS/MS database...'))
    
    if (masstools::get_os() == "windows") {
      bpparam = BiocParallel::SnowParam(workers = threads,
                                        progressbar = TRUE)
    } else{
      bpparam = BiocParallel::MulticoreParam(workers = threads,
                                             progressbar = TRUE)
    }

    #####bug fixing
    # for(i in  seq_len(nrow(ms1.info))){
    #   cat(i, " ")
    #   identify_peak(idx = i,
    #                ms1.info = ms1.info,
    #                ms2.info = ms2.info,
    #                spectra.info = spectra.info,
    #                spectra.data = spectra.data,
    #                ppm.ms1match = ms1.match.ppm,
    #                ppm.ms2match = ms2.match.ppm,
    #                mz.ppm.thr = mz.ppm.thr,
    #                ms2.match.tol = ms2.match.tol,
    #                rt.match.tol = rt.match.tol,
    #                ms1.match.weight = ms1.match.weight,
    #                rt.match.weight = rt.match.weight,
    #                ms2.match.weight = ms2.match.weight,
    #                total.score.tol = total.score.tol,
    #                adduct.table = adduct.table,
    #                candidate.num = candidate.num,
    #                fraction.weight = fraction.weight,
    #                dp.forward.weight = dp.forward.weight,
    #                dp.reverse.weight = dp.reverse.weight,
    #                remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff
    #                )
    # }
        
    identification.result <-
      suppressMessages(
        BiocParallel::bplapply(
          seq_len(nrow(ms1.info)),
          FUN = identify_peak,
          BPPARAM = bpparam,
          ms1.info = ms1.info,
          ms2.info = ms2.info,
          spectra.info = spectra.info,
          spectra.data = spectra.data,
          ppm.ms1match = ms1.match.ppm,
          ppm.ms2match = ms2.match.ppm,
          mz.ppm.thr = mz.ppm.thr,
          ms2.match.tol = ms2.match.tol,
          rt.match.tol = rt.match.tol,
          ms1.match.weight = ms1.match.weight,
          rt.match.weight = rt.match.weight,
          ms2.match.weight = ms2.match.weight,
          total.score.tol = total.score.tol,
          adduct.table = adduct.table,
          candidate.num = candidate.num,
          fraction.weight = fraction.weight,
          dp.forward.weight = dp.forward.weight,
          dp.reverse.weight = dp.reverse.weight,
          remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff
        )
      )
    
    names(identification.result) <- ms1.info$name
    identification.result <-
      identification.result[which(!unlist(lapply(identification.result, function(x)
        all(is.na(x)))))]
    if (length(identification.result) == 0) {
      return(list(NULL))
    }
    return(identification.result)
  }


#---------------------------------------------------------------------------

#' @title Identify metabolites based on MS1 or MS/MS database
#' @description Identify metabolites based on MS1 or MS/MS database.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param idx idx
#' @param ms1.info ms1.info
#' @param ms2.info ms2.info
#' @param spectra.info spectra.info
#' @param spectra.data spectra.data
#' @param ppm.ms1match ppm.ms1match
#' @param ppm.ms2match ppm.ms2match
#' @param mz.ppm.thr Accurate mass tolerance for m/z error calculation.
#' @param ms2.match.tol MS2 match (MS2 similarity) tolerance.
#' @param rt.match.tol The weight for matched fragments.
#' @param ms1.match.weight Forward dot product weight.
#' @param rt.match.weight Reverse dot product weight.
#' @param ms2.match.weight RT match tolerance.
#' @param total.score.tol The polarity of data, "positive"or "negative".
#' @param adduct.table Collision energy. Please confirm the CE values in your database. Default is "all".
#' @param candidate.num "hilic" (HILIC column) or "rp" (reverse phase).
#' @param fraction.weight The weight of MS1 match for total score calculation.
#' @param dp.forward.weight The weight of RT match for total score calculation.
#' @param dp.reverse.weight The weight of MS2 match for total score calculation.
#' @param remove_fragment_intensity_cutoff remove_fragment_intensity_cutoff
#' @param ... other parameters
#' @return A metIdentifyClass object.


identify_peak <-
  function(idx,
           ms1.info,
           ms2.info,
           spectra.info,
           spectra.data,
           ppm.ms1match = 25,
           ppm.ms2match = 30,
           mz.ppm.thr = 400,
           ms2.match.tol = 0.5,
           rt.match.tol = 30,
           ms1.match.weight = 0.25,
           rt.match.weight = 0.25,
           ms2.match.weight = 0.5,
           total.score.tol = 0.5,
           adduct.table,
           candidate.num = 3,
           fraction.weight = 0.3,
           dp.forward.weight = 0.6,
           dp.reverse.weight = 0.1,
           remove_fragment_intensity_cutoff = 0,
           ...) {
    pk.precursor <- ms1.info[idx, , drop = FALSE]
    rm(list = c("ms1.info"))
    pk.mz <- pk.precursor$mz
    pk.rt <- pk.precursor$rt
    pk.spec <- ms2.info[[idx]]
    rm(list = c("ms2.info"))
    
    if (length(pk.spec) == 0) {
      return(NA)
    }
    
    spectra.mz <-
      apply(adduct.table, 1, function(x) {
        temp.n <-
          stringr::str_extract(string = as.character(x[1]), pattern = "[0-9]{1}M")
        temp.n <-
          as.numeric(stringr::str_replace(
            string = temp.n,
            pattern = "M",
            replacement = ""
          ))
        temp.n[is.na(temp.n)] <- 1
        as.numeric(x[2]) + temp.n * as.numeric(spectra.info$mz)
      })
    
    colnames(spectra.mz) <- adduct.table[, 1]
    rownames(spectra.mz) <- spectra.info$Lab.ID
    
    ###mz match
    # match.idx <-
    #   apply(spectra.mz, 1, function(x) {
    #     temp.mz.error <-
    #       abs(x - pk.mz) * 10 ^ 6 / ifelse(pk.mz < 400, 400, pk.mz)
    #     temp.mz.match.score <-
    #       exp(-0.5 * (temp.mz.error / (ppm.ms1match)) ^ 2)
    #     data.frame(
    #       "addcut" = names(temp.mz.error)[which(temp.mz.error < ppm.ms1match)],
    #       "mz.error" = temp.mz.error[which(temp.mz.error < ppm.ms1match)],
    #       "mz.match.score" = temp.mz.match.score[which(temp.mz.error < ppm.ms1match)],
    #       stringsAsFactors = FALSE
    #     )
    #   })
    
    temp.mz.error <-
      abs(spectra.mz - pk.mz) * 10 ^ 6 / ifelse(pk.mz < 400, 400, pk.mz)
    
    temp.mz.error[which(is.na(temp.mz.error), arr.ind = TRUE)] <-
      ppm.ms1match + 1
    
    temp.mz.match.score <-
      exp(-0.5 * (temp.mz.error / (ppm.ms1match)) ^ 2)
    
    if (sum(temp.mz.error < ppm.ms1match) == 0) {
      return(NA)
    }
    
    match.idx <-
      which(temp.mz.error < ppm.ms1match, arr.ind = TRUE) %>%
      as.data.frame()
    
    match.idx <-
      seq_len(nrow(match.idx)) %>%
      purrr::map(function(i) {
        data.frame(
          Lab.ID = rownames(spectra.mz)[match.idx$row[i]],
          Adduct = colnames(spectra.mz)[match.idx$col[i]],
          mz.error = as.numeric(temp.mz.error[match.idx$row[i], match.idx$col[i]]),
          mz.match.score = as.numeric(temp.mz.match.score[match.idx$row[i], match.idx$col[i]])
        )
      }) %>%
      dplyr::bind_rows()
    
    rownames(match.idx) <- NULL
    
    rm(list = c("spectra.mz", "adduct.table"))
    
    ###RT match
    if (rt.match.tol > 10000) {
      temp.rt.error <-
        temp.rt.match.score <-
        rep(NA, nrow(match.idx))
    } else{
      temp.rt.error <-
        abs(pk.rt - spectra.info$RT[match(match.idx$Lab.ID, spectra.info$Lab.ID)])
      
      temp.rt.match.score <-
        exp(-0.5 * (temp.rt.error / (rt.match.tol)) ^ 2)
    }
    
    RT.error <-
      data.frame(RT.error = temp.rt.error,
                 RT.match.score = temp.rt.match.score)
    
    match.idx <-
      data.frame(match.idx, RT.error, stringsAsFactors = FALSE)
    
    rm(list = c("RT.error"))
    
    if (any(!is.na(match.idx$RT.error))) {
      RT.error <- match.idx$RT.error
      RT.error[is.na(RT.error)] <- rt.match.tol - 1
      match.idx <-
        match.idx[RT.error < rt.match.tol, , drop = FALSE]
    }
    
    if (nrow(match.idx) == 0) {
      return(NA)
    }
    
    ###MS2 spectra match
    ms2.score <-
      seq_len(nrow(match.idx)) %>% 
      purrr::map(function(i){
        lib.spec <- spectra.data[[match.idx$Lab.ID[i]]]
        dp <- lapply(lib.spec, function(y) {
          y <- as.data.frame(y)
          y$mz <- as.numeric(y$mz)
          y$intensity <- as.numeric(y$intensity)
          masstools::get_spectra_match_score(
            exp.spectrum = as.data.frame(pk.spec),
            lib.spectrum = y,
            ppm.tol = ppm.ms2match,
            mz.ppm.thr = mz.ppm.thr,
            fraction.weight = fraction.weight,
            dp.forward.weight = dp.forward.weight,
            dp.reverse.weight = dp.reverse.weight,
            remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff
          )
        })
        dp <- dp[which.max(unlist(dp))]
        dp <- unlist(dp)
        data.frame("CE" = names(dp),
                   "SS" = dp,
                   stringsAsFactors = FALSE)
      }) %>% 
      dplyr::bind_rows()
    
    rownames(ms2.score) <- NULL
    
    match.idx <-
      data.frame(match.idx, ms2.score, stringsAsFactors = FALSE)
    match.idx <-
      match.idx[which(match.idx$SS > ms2.match.tol), , drop = FALSE]
    rm(list = c("ms2.score"))
    if (nrow(match.idx) == 0) {
      return(NA)
    }
    
    ###total score
    total.score <- 
      seq_len(nrow(match.idx)) %>% 
      purrr::map(function(i){
        if(is.na(match.idx$RT.error[i])){
          as.numeric(match.idx$mz.match.score[i]) * (ms1.match.weight + rt.match.weight /
                                                       2) +
            as.numeric(match.idx$SS[i]) * (ms2.match.weight + rt.match.weight /
                                             2)
        }else{
          as.numeric(match.idx$mz.match.score[i]) * ms1.match.weight +
            as.numeric(match.idx$RT.match.score[i]) * rt.match.weight +
            as.numeric(match.idx$SS[i]) * ms2.match.weight
        }
      }) %>% 
      unlist() %>% 
      as.numeric()
    
    match.idx <- data.frame(match.idx,
                            "Total.score" = total.score,
                            stringsAsFactors = FALSE)
    rm(list = c("total.score"))
    match.idx <-
      match.idx[match.idx$Total.score > total.score.tol, , drop = FALSE]
    
    if (nrow(match.idx) == 0) {
      return(NA)
    }
    
    match.idx <-
      match.idx %>% 
      dplyr::arrange(dplyr::desc(Total.score))
    
    if (nrow(match.idx) > candidate.num) {
      match.idx <- match.idx[seq_len(candidate.num), ]
    }
    ##add other information
    match.idx <-
      data.frame(spectra.info[match(match.idx$Lab.ID, spectra.info$Lab.ID),
                              c("Compound.name", "CAS.ID", "HMDB.ID", "KEGG.ID")], 
                 match.idx,
                 stringsAsFactors = FALSE)
    
    return(match.idx)
  }




# identify_peak <-
#   function(idx,
#            ms1.info,
#            ms2.info,
#            spectra.info,
#            spectra.data,
#            ppm.ms1match = 25,
#            ppm.ms2match = 30,
#            mz.ppm.thr = 400,
#            ms2.match.tol = 0.5,
#            rt.match.tol = 30,
#            ms1.match.weight = 0.25,
#            rt.match.weight = 0.25,
#            ms2.match.weight = 0.5,
#            total.score.tol = 0.5,
#            adduct.table,
#            candidate.num = 3,
#            fraction.weight = 0.3,
#            dp.forward.weight = 0.6,
#            dp.reverse.weight = 0.1,
#            remove_fragment_intensity_cutoff = 0,
#            ...) {
#     pk.precursor <- ms1.info[idx, ,drop = FALSE]
#     rm(list = c("ms1.info"))
#     pk.mz <- pk.precursor$mz
#     pk.rt <- pk.precursor$rt
#     pk.spec <- ms2.info[[idx]]
#     rm(list = c("ms2.info"))
#     
#     if (length(pk.spec) == 0) {
#       return(NA)
#     }
#     
#     spectra.mz <-
#       apply(adduct.table, 1, function(x) {
#         temp.n <-
#           stringr::str_extract(string = as.character(x[1]), pattern = "[0-9]{1}M")
#         temp.n <-
#           as.numeric(stringr::str_replace(
#             string = temp.n,
#             pattern = "M",
#             replacement = ""
#           ))
#         temp.n[is.na(temp.n)] <- 1
#         as.numeric(x[2]) + temp.n * as.numeric(spectra.info$mz)
#       })
#     
#     colnames(spectra.mz) <- adduct.table[, 1]
#     rownames(spectra.mz) <- spectra.info$Lab.ID
#     
#     ###mz match
#     match.idx <-
#       apply(spectra.mz, 1, function(x) {
#         temp.mz.error <-
#           abs(x - pk.mz) * 10 ^ 6 / ifelse(pk.mz < 400, 400, pk.mz)
#         temp.mz.match.score <-
#           exp(-0.5 * (temp.mz.error / (ppm.ms1match)) ^ 2)
#         data.frame(
#           "addcut" = names(temp.mz.error)[which(temp.mz.error < ppm.ms1match)],
#           "mz.error" = temp.mz.error[which(temp.mz.error < ppm.ms1match)],
#           "mz.match.score" = temp.mz.match.score[which(temp.mz.error < ppm.ms1match)],
#           stringsAsFactors = FALSE
#         )
#       })
#     
#     # temp.mz.error <-
#     #   abs(spectra.mz - pk.mz) * 10 ^ 6 / ifelse(pk.mz < 400, 400, pk.mz)
#     # 
#     # temp.mz.error[which(is.na(temp.mz.error), arr.ind = TRUE)] <- ppm.ms1match + 1
#     # 
#     # temp.mz.match.score <-
#     #   exp(-0.5 * (temp.mz.error / (ppm.ms1match)) ^ 2)
#     # 
#     # if(sum(temp.mz.error < ppm.ms1match) == 0){
#     #   return(NA)
#     # }
#     
#     # match.idx <-
#     #   spectra.mz %>%
#     #   t() %>%
#     #   as.data.frame() %>%
#     #   purrr::map(function(x) {
#     #     temp.mz.error <-
#     #       abs(x - pk.mz) * 10 ^ 6 / ifelse(pk.mz < 400, 400, pk.mz)
#     #     temp.mz.match.score <-
#     #       exp(-0.5 * (temp.mz.error / (ppm.ms1match)) ^ 2)
#     #     
#     #     data.frame(
#     #       "addcut" = colnames(spectra.mz),
#     #       "mz.error" = temp.mz.error,
#     #       "mz.match.score" = temp.mz.match.score,
#     #       stringsAsFactors = FALSE
#     #     ) %>%
#     #       dplyr::filter(temp.mz.error < ppm.ms1match)
#     #   })
#     
#     rm(list = c("spectra.mz", "adduct.table"))
#     
#     ##remove some none matched
#     match.idx <-
#       match.idx[which(unlist(lapply(match.idx, function(x) {
#         nrow(x)
#       })) != 0)]
#     
#     if (length(match.idx) == 0) {
#       return(NA)
#     }
#     
#     match.idx <- 
#       mapply(function(x, y) {
#       list(data.frame("Lab.ID" = y, x, stringsAsFactors = FALSE))
#     },
#     x = match.idx,
#     y = names(match.idx))
#     
#     match.idx <- do.call(rbind, match.idx)
#     # match.idx <- data.frame(rownames(match.idx), match.idx, stringsAsFactors = FALSE)
#     colnames(match.idx) <-
#       c("Lab.ID", "Adduct", "mz.error", "mz.match.score")
#     rownames(match.idx) <- NULL
#     
#     ###RT match
#     RT.error <- t(apply(match.idx, 1, function(x) {
#       temp.rt.error <-
#         abs(pk.rt - spectra.info$RT[match(x[1], spectra.info$Lab.ID)])
#       temp.rt.match.score <-
#         exp(-0.5 * (temp.rt.error / (rt.match.tol)) ^ 2)
#       
#       ####if user set rt.match.tol as FALSE, set RT.error as NA
#       if (rt.match.tol > 10000) {
#         temp.rt.error = NA
#         temp.rt.match.score = NA
#       }
#       c(temp.rt.error, temp.rt.match.score)
#     }))
#     colnames(RT.error) <- c("RT.error", "RT.match.score")
#     match.idx <-
#       data.frame(match.idx, RT.error, stringsAsFactors = FALSE)
#     rm(list = c("RT.error"))
#     
#     if (any(!is.na(match.idx$RT.error))) {
#       RT.error <- match.idx$RT.error
#       RT.error[is.na(RT.error)] <- rt.match.tol - 1
#       match.idx <-
#         match.idx[RT.error < rt.match.tol, , drop = FALSE]
#     }
#     
#     if (nrow(match.idx) == 0) {
#       return(NA)
#     }
#     
#     ###MS2 spectra match
#     ms2.score <-
#       apply(match.idx, 1, function(x) {
#         x <- as.character(x)
#         lib.spec <- spectra.data[[x[1]]]
#         dp <- lapply(lib.spec, function(y) {
#           masstools::get_spectra_match_score(
#             exp.spectrum = as.data.frame(pk.spec),
#             lib.spectrum = y,
#             ppm.tol = ppm.ms2match,
#             mz.ppm.thr = mz.ppm.thr,
#             fraction.weight = fraction.weight,
#             dp.forward.weight = dp.forward.weight,
#             dp.reverse.weight = dp.reverse.weight,
#             remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff
#           )
#         })
#         dp <- dp[which.max(unlist(dp))]
#         dp <- unlist(dp)
#         data.frame("CE" = names(dp),
#                    "SS" = dp,
#                    stringsAsFactors = FALSE)
#       })
#     
#     ms2.score <- do.call(rbind, ms2.score)
#     rownames(ms2.score) <- NULL
#     
#     match.idx <-
#       data.frame(match.idx, ms2.score, stringsAsFactors = FALSE)
#     match.idx <-
#       match.idx[which(match.idx$SS > ms2.match.tol), , drop = FALSE]
#     rm(list = c("ms2.score"))
#     if (nrow(match.idx) == 0) {
#       return(NA)
#     }
#     
#     ###total score
#     total.score <- apply(match.idx, 1, function(x) {
#       if (is.na(x["RT.match.score"])) {
#         as.numeric(x["mz.match.score"]) * (ms1.match.weight + rt.match.weight /
#                                              2) +
#           as.numeric(x['SS']) * (ms2.match.weight + rt.match.weight /
#                                    2)
#       } else{
#         as.numeric(x['mz.match.score']) * ms1.match.weight +
#           as.numeric(x['RT.match.score']) * rt.match.weight +
#           as.numeric(x["SS"]) * ms2.match.weight
#       }
#     })
#     
#     match.idx <- data.frame(match.idx,
#                             "Total.score" = total.score,
#                             stringsAsFactors = FALSE)
#     rm(list = c("total.score"))
#     match.idx <-
#       match.idx[match.idx$Total.score > total.score.tol, , drop = FALSE]
#     
#     if (nrow(match.idx) == 0) {
#       return(NA)
#     }
#     
#     match.idx <-
#       match.idx[order(match.idx$Total.score, decreasing = TRUE),]
#     if (nrow(match.idx) > candidate.num) {
#       match.idx <- match.idx[seq_len(candidate.num),]
#     }
#     ##add other information
#     match.idx <-
#       data.frame(spectra.info[match(match.idx$Lab.ID, spectra.info$Lab.ID),
#                               c("Compound.name", "CAS.ID", "HMDB.ID", "KEGG.ID")], match.idx,
#                  stringsAsFactors = FALSE)
#     match.idx <-
#       match.idx[order(match.idx$Total.score, decreasing = TRUE), , drop = FALSE]
#     
#     return(match.idx)
#   }



# plotMS2match(matched.info = temp.matched.info, exp.spectrum = exp.spectrum,
#              lib.spectrum = lib.spectrum, database = database)
#
# matched.info <- temp.matched.info
# range.mz = range.mz
# exp.spectrum = exp.spectrum
# lib.spectrum = lib.spectrum
# col.lib = col.lib
# col.exp = col.exp
# ce = ce
# polarity = polarity
# database = database
plotMS2match = function(matched.info,
                        range.mz,
                        ppm.tol = 30,
                        mz.ppm.thr = 400,
                        exp.spectrum,
                        lib.spectrum,
                        polarity = c("positive", "negative"),
                        xlab = "Mass to charge ratio (m/z)",
                        ylab = "Relative intensity",
                        col.lib = "red",
                        col.exp = "black",
                        ce = "30",
                        title.size = 15,
                        lab.size = 15,
                        axis.text.size = 15,
                        legend.title.size = 15,
                        legend.text.size = 15,
                        database) {
  polarity <- match.arg(polarity)
  exp.spectrum[, 1] <- as.numeric(exp.spectrum[, 1])
  exp.spectrum[, 2] <- as.numeric(exp.spectrum[, 2])
  
  lib.spectrum[, 1] <- as.numeric(lib.spectrum[, 1])
  lib.spectrum[, 2] <- as.numeric(lib.spectrum[, 2])
  
  exp.spectrum[, 2] <-
    exp.spectrum[, 2] / max(exp.spectrum[, 2])
  lib.spectrum[, 2] <-
    lib.spectrum[, 2] / max(lib.spectrum[, 2])
  
  exp.spectrum <- as.data.frame(exp.spectrum)
  lib.spectrum <- as.data.frame(lib.spectrum)
  if (missing(range.mz)) {
    range.mz <- c(min(exp.spectrum[, 1], lib.spectrum[, 1]),
                  max(exp.spectrum[, 1], lib.spectrum[, 1]))
    
  }
  
  matched.spec <-
    masstools::ms2_match(
      exp.spectrum = exp.spectrum,
      lib.spectrum = lib.spectrum,
      ppm.tol = ppm.tol,
      mz.ppm.thr = mz.ppm.thr
    )
  matched.idx <-
    which(matched.spec[, "Lib.intensity"] > 0 &
            matched.spec[, "Exp.intensity"] > 0)
  
  plot <- ggplot2::ggplot(data = matched.spec) +
    ggplot2::geom_segment(
      mapping = ggplot2::aes(
        x = Exp.mz,
        y = Exp.intensity - Exp.intensity,
        xend = Exp.mz,
        yend = Exp.intensity
      ),
      colour = col.exp
    ) +
    ggplot2::geom_point(
      data = matched.spec[matched.idx, , drop = FALSE],
      mapping = ggplot2::aes(x = Exp.mz, y = Exp.intensity),
      colour = col.exp
    ) +
    ggplot2::xlim(range.mz[1], range.mz[2]) +
    ggplot2::ylim(-1, 1) +
    ggplot2::labs(x = xlab,
                  y = ylab,
                  title = as.character(matched.info["Compound.name"])) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      # axis.line = ggplot2::element_line(arrow = ggplot2::arrow()),
      plot.title = ggplot2::element_text(
        color = "black",
        size = title.size,
        face = "plain",
        hjust = 0.5
      ),
      axis.title = ggplot2::element_text(
        color = "black",
        size = lab.size,
        face = "plain"
      ),
      axis.text = ggplot2::element_text(
        color = "black",
        size = axis.text.size,
        face = "plain"
      ),
      legend.title = ggplot2::element_text(
        color = "black",
        size = legend.title.size,
        face = "plain"
      ),
      legend.text = ggplot2::element_text(
        color = "black",
        size = legend.text.size,
        face = "plain"
      )
    )
  
  temp.info <-
    unlist(database@spectra.info[match(matched.info["Lab.ID"], database@spectra.info$Lab.ID), , drop = TRUE])
  temp.info <- temp.info[!is.na(temp.info)]
  temp.info <-
    temp.info[unlist(lapply(temp.info, stringr::str_count)) < 50]
  temp.info <-
    paste(names(temp.info), temp.info, sep = ": ")
  
  plot <- plot +
    ggplot2::annotate(
      geom = "text",
      x = -Inf,
      y = Inf,
      label = paste(temp.info, collapse = "\n"),
      hjust = 0,
      vjust = 1
    )
  
  temp.info2 <-
    matched.info[c("mz.error", "RT.error", "SS", "Total.score", "Adduct", "CE")]
  temp.info2 <- temp.info2[!is.na(temp.info2)]
  
  temp.info2 <-
    paste(names(temp.info2), temp.info2, sep = ": ")
  plot <- plot +
    ggplot2::annotate(
      geom = "text",
      x = -Inf,
      y = -Inf,
      label = paste(temp.info2, collapse = "\n"),
      hjust = 0,
      vjust = 0
    )
  
  plot <- plot +
    ggplot2::annotate(
      geom = "text",
      x = Inf,
      y = Inf,
      label = "Experiment MS2 spectrum",
      color = col.exp,
      hjust = 1,
      vjust = 1
    ) +
    ggplot2::annotate(
      geom = "text",
      x = Inf,
      y = -Inf,
      label = "Database MS2 spectrum",
      color = col.lib,
      hjust = 1,
      vjust = -1
    )
  
  plot <- plot +
    ggplot2::geom_segment(
      data = matched.spec,
      mapping = ggplot2::aes(
        x = Lib.mz,
        y = Lib.intensity - Lib.intensity,
        xend = Lib.mz,
        yend = -Lib.intensity
      ),
      colour = col.lib
    ) +
    ggplot2::geom_point(
      data = matched.spec[matched.idx, , drop = FALSE],
      mapping = ggplot2::aes(x = Lib.mz, y = -Lib.intensity),
      colour = col.lib
    )
  plot
}
