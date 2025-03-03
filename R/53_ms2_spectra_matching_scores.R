#' Calculate MS2 Spectra Matching Score
#'
#' This function computes the matching score between experimental and library MS2 spectra using the dot product method. The function allows for fragment matching based on mass-to-charge ratio (m/z) tolerance and intensity comparison, and calculates the score based on forward and reverse dot product and the fraction of matched fragments.
#'
#' @param experimental.spectrum A data frame representing the experimental spectrum, with columns `mz` for mass-to-charge ratios and `intensity` for corresponding intensities.
#' @param library.spectrum A data frame representing the library spectrum, with columns `mz` for mass-to-charge ratios and `intensity` for corresponding intensities.
#' @param ms2.match.ppm Numeric. The mass tolerance in parts per million (ppm) for matching MS2 fragments. Default is 30.
#' @param mz.ppm.thr Numeric. Threshold for m/z to use in ppm calculation. Default is 400.
#' @param method Character. The method to calculate the matching score. Currently, only "dotproduct" is supported. Default is "dotproduct".
#' @param fraction.weight Numeric. The weight of the fraction of matched fragments in the overall score. Default is 0.2.
#' @param dp.forward.weight Numeric. The weight for the forward dot product score in the overall score. Default is 0.7.
#' @param dp.reverse.weight Numeric. The weight for the reverse dot product score in the overall score. Default is 0.1.
#' @param remove_fragment_intensity_cutoff Numeric. A cutoff value to remove low-intensity fragments from both experimental and library spectra. Default is 0.
#'
#' @return A numeric value representing the MS2 matching score, which is a weighted combination of the forward dot product score, reverse dot product score, and the fraction of matched fragments.
#'
#' @details
#' The function normalizes the intensities of both the experimental and library spectra and filters out fragments below the specified intensity cutoff. It then calculates the matching score based on the dot product method. The forward dot product score is calculated using all matched fragments, while the reverse dot product score is computed using only the fragments that have non-zero intensity in the library spectrum.
#'
#' @examples
#' experimental.spectrum <- data.frame(mz = 1:10, intensity = 1:10)
#' library.spectrum <- data.frame(mz = 1:10, intensity = 1:10)
#' get_spectra_match_score(experimental.spectrum, library.spectrum)
#'
#' @author
#' Xiaotao Shen \email{xiaotao.shen@outlook.com}
#'
#' @export

calculate_ms2_matching_score <-
  function(experimental.spectrum,
           library.spectrum,
           ms2.match.ppm = 30,
           mz.ppm.thr = 400,
           method = c("dotproduct"),
           fraction.weight = 0.2,
           dp.forward.weight = 0.7,
           dp.reverse.weight = 0.1,
           remove_fragment_intensity_cutoff = 0) {
    method <- match.arg(method)
    experimental.spectrum <- as.data.frame(experimental.spectrum)
    library.spectrum <- as.data.frame(library.spectrum)
    
    experimental.spectrum$intensity <-
      experimental.spectrum$intensity / max(experimental.spectrum$intensity)
    library.spectrum$intensity <-
      library.spectrum$intensity / max(library.spectrum$intensity)
    
    experimental.spectrum <-
      experimental.spectrum %>%
      dplyr::filter(intensity > remove_fragment_intensity_cutoff)
    
    library.spectrum <-
      library.spectrum %>%
      dplyr::filter(intensity > remove_fragment_intensity_cutoff)
    
    match.matrix <-
      match_ms2_fragments(
        experimental.spectrum = experimental.spectrum,
        library.spectrum = library.spectrum,
        ms2.match.ppm = ms2.match.ppm,
        mz.ppm.thr = mz.ppm.thr
      )
    
    if (method == "dotproduct") {
      fraction <-
        sum(!is.na(match.matrix$Lib.index) &
              !is.na(match.matrix$Exp.index)) / nrow(match.matrix)
      
      dp_forward <- calculate_dotproduct(exp.int = match.matrix$Exp.intensity,
                                         lib.int = match.matrix$Lib.intensity)
      dp_reverse <-
        calculate_dotproduct(exp.int =
                               match.matrix$Exp.intensity[which(match.matrix$Lib.intensity > 0)],
                             lib.int =
                               match.matrix$Lib.intensity[which(match.matrix$Lib.intensity > 0)])
      dp_forward[is.na(dp_forward)] <- 0
      dp_reverse[is.na(dp_reverse)] <- 0
      ms2_match_score <-
        dp_forward * dp.forward.weight + dp_reverse * dp.reverse.weight +
        fraction * fraction.weight
    }
    
    return(ms2_match_score)
  }


#' Calculate Dot Product Score for MS2 Spectra
#'
#' This function calculates the dot product score between experimental and library MS2 spectra intensities. It weights the intensities and computes the similarity score using a modified dot product formula.
#'
#' @param exp.int A numeric vector of intensities from the experimental spectrum.
#' @param lib.int A numeric vector of intensities from the library spectrum.
#'
#' @return A numeric value representing the dot product similarity score between the experimental and library spectra.
#'
#' @details
#' The function applies a custom weighting scheme to the experimental and library intensities, where the weight for each intensity value is computed based on its proportion relative to the total sum of intensities in its respective spectrum. The weighted intensities are then used to calculate the dot product score, which is normalized to give a similarity score between the two spectra.
#'
#' @examples
#' calculate_dotproduct(exp.int = 1:10, lib.int = 1:10)
#'
#' @author
#' Xiaotao Shen \email{xiaotao.shen@outlook.com}
#'
#' @export

calculate_dotproduct <-
  function(exp.int, lib.int) {
    exp.weight <- lapply(exp.int, function(x) {
      1 / (1 + x / (sum(exp.int) - 0.5))
    }) %>%
      unlist()
    
    lib.weight <- lapply(lib.int, function(x) {
      1 / (1 + x / (sum(lib.int) - 0.5))
    }) %>%
      unlist()
    
    x <- exp.weight * exp.int
    y <- lib.weight * lib.int
    return(sum(x * y) ^ 2 / (sum(x ^ 2) * sum(y ^ 2)))
  }


#' @title Match MS2 Fragments
#' @description Matches MS2 fragment ions between an experimental spectrum and a library spectrum.
#' @author Xiaotao Shen \email{xiaotao.shen@outlook.com}
#' 
#' @param experimental.spectrum A data frame with `mz` and `intensity` columns representing the experimental MS2 spectrum.
#' @param library.spectrum A data frame with `mz` and `intensity` columns representing the reference MS2 spectrum.
#' @param ms2.match.ppm Numeric, mass tolerance in parts per million (ppm) for fragment matching. Default is 30 ppm.
#' @param mz.ppm.thr Numeric, minimum m/z threshold for ppm-based error calculation. Default is 400.
#' @param direction Character, either `"reverse"` (library to experimental matching) or `"forward"` (experimental to library matching). Default is `"reverse"`.
#' @param remove.noise Logical, whether to remove low-intensity fragment ions before matching. Default is `TRUE`.
#'
#' @details This function aligns the fragment ions from an experimental MS2 spectrum with those from a library spectrum based on m/z similarity within a given ppm tolerance.
#' 
#' - If `direction = "reverse"`, the function looks for matches of each **library fragment ion** in the **experimental spectrum**.
#' - If `direction = "forward"`, the function looks for matches of each **experimental fragment ion** in the **library spectrum**.
#' 
#' If multiple matches are found for a given fragment, the most intense one is selected.
#'
#' @return A data frame containing the matched fragment ions with the following columns:
#' \item{Lib.index}{Index of the fragment in the library spectrum.}
#' \item{Exp.index}{Index of the fragment in the experimental spectrum.}
#' \item{Lib.mz}{m/z value of the library fragment ion.}
#' \item{Lib.intensity}{Intensity of the library fragment ion.}
#' \item{Exp.mz}{m/z value of the experimental fragment ion.}
#' \item{Exp.intensity}{Intensity of the experimental fragment ion.}
#'
#' Unmatched peaks in either spectrum are included with missing values (`NA`) in the corresponding columns.
#' 
#' @examples
#' # Example data for experimental and library MS2 spectra
#' experimental.spectrum <- data.frame(mz = c(100.1, 150.2, 200.3), intensity = c(300, 500, 200))
#' library.spectrum <- data.frame(mz = c(100.09, 150.25, 210.4), intensity = c(250, 600, 150))
#' 
#' # Perform MS2 fragment matching
#' matched_fragments <- match_ms2_fragments(
#'   experimental.spectrum = experimental.spectrum,
#'   library.spectrum = library.spectrum,
#'   ms2.match.ppm = 30,
#'   direction = "reverse"
#' )
#'
#' print(matched_fragments)
#'
#' @export

match_ms2_fragments <-
  function(experimental.spectrum,
           library.spectrum,
           ms2.match.ppm = 30,
           mz.ppm.thr = 400,
           direction = c("reverse", "forward"),
           remove.noise = TRUE) {
    direction <-
      match.arg(direction)
    ## remove noisy fragments
    if (remove.noise) {
      experimental.spectrum <-
        remove_noise(
          spectrum = experimental.spectrum,
          ms2.match.ppm = ms2.match.ppm,
          mz.ppm.thr = mz.ppm.thr
        )
      library.spectrum <-
        remove_noise(
          spectrum = library.spectrum,
          ms2.match.ppm = ms2.match.ppm,
          mz.ppm.thr = mz.ppm.thr
        )
    }
    
    ## for each fragment in library.spectrum,
    ## its matched fragments index in experimental.spectrum
    if (direction == "reverse") {
      match.idx <-
        lapply(library.spectrum$mz, function(x) {
          diff.mz <- abs(x - experimental.spectrum$mz)
          x[x < mz.ppm.thr] <- mz.ppm.thr
          mz.error <- diff.mz * 10 ^ 6 / x
          temp.idx <- which(mz.error < ms2.match.ppm)
          if (length(temp.idx) == 0) {
            return(NA)
          }
          if (length(temp.idx) > 1) {
            return(temp.idx[which.max(experimental.spectrum$intensity[temp.idx])])
          }
          return(temp.idx)
        })
      
      match.idx <- do.call(rbind, match.idx)
      match.idx <- cbind(seq_len(nrow(match.idx)), match.idx)
      colnames(match.idx) <- c("Lib", "Exp")
      
      non.idx2 <-
        setdiff(c(seq_len(nrow(
          experimental.spectrum
        ))), match.idx[, 2][!is.na(match.idx[, 2])])
      
      if (length(non.idx2) != 0) {
        match.idx2 <- data.frame(NA, non.idx2, stringsAsFactors = FALSE)
        colnames(match.idx2) <- c("Lib", "Exp")
      } else {
        match.idx2 <- NULL
      }
      
      match.matrix <-
        as.data.frame(rbind(match.idx, match.idx2), stringsAsFactors = FALSE)
      
      match.matrix <-
        data.frame(match.matrix, library.spectrum[match.matrix$Lib, c(1, 2)], experimental.spectrum[match.matrix$Exp, c(1, 2)])
      colnames(match.matrix) <-
        c("Lib.index",
          "Exp.index",
          "Lib.mz",
          "Lib.intensity",
          "Exp.mz",
          "Exp.intensity")
    } else{
      match.idx <-
        lapply(experimental.spectrum$mz, function(x) {
          diff.mz <- abs(x - library.spectrum$mz)
          x[x < mz.ppm.thr] <- mz.ppm.thr
          mz.error <- diff.mz * 10 ^ 6 / x
          temp.idx <- which(mz.error < ms2.match.ppm)
          if (length(temp.idx) == 0) {
            return(NA)
          }
          if (length(temp.idx) > 1) {
            return(temp.idx[which.max(library.spectrum$intensity[temp.idx])])
          }
          return(temp.idx)
        })
      
      match.idx <- do.call(rbind, match.idx)
      match.idx <- cbind(seq_len(nrow(match.idx)), match.idx)
      colnames(match.idx) <- c("Exp", "Lib")
      
      non.idx2 <-
        setdiff(c(seq_len(nrow(
          library.spectrum
        ))), match.idx[, 2][!is.na(match.idx[, 2])])
      
      if (length(non.idx2) != 0) {
        match.idx2 <- data.frame(NA, non.idx2, stringsAsFactors = FALSE)
        colnames(match.idx2) <- c("Exp", "Lib")
      } else {
        match.idx2 <- NULL
      }
      
      match.matrix <-
        as.data.frame(rbind(match.idx, match.idx2), stringsAsFactors = FALSE)
      
      match.matrix <-
        data.frame(match.matrix, library.spectrum[match.matrix$Lib, c(1, 2)], experimental.spectrum[match.matrix$Exp, c(1, 2)])
      colnames(match.matrix) <-
        c("Exp.index",
          "Lib.index",
          "Lib.mz",
          "Lib.intensity",
          "Exp.mz",
          "Exp.intensity")
    }
    
    match.matrix$Lib.intensity[is.na(match.matrix$Lib.intensity)] <- 0
    match.matrix$Exp.intensity[is.na(match.matrix$Exp.intensity)] <- 0
    rownames(match.matrix) <- NULL
    match.matrix
  }


#' Remove Noise from MS2 Spectrum
#'
#' This function removes noisy fragments from an MS2 spectrum based on mass-to-charge ratio (m/z) tolerance and intensity comparison, ensuring that only the most relevant peaks are retained.
#'
#' @param spectrum A data frame representing the MS2 spectrum, with columns `mz` for mass-to-charge ratios and `intensity` for corresponding intensities.
#' @param ms2.match.ppm Numeric. The mass tolerance in parts per million (ppm) for determining whether peaks are considered duplicates/noise. Default is 30.
#' @param mz.ppm.thr Numeric. A threshold for m/z values to use in ppm calculation. If the m/z difference is below this value, it is set to this threshold. Default is 400.
#'
#' @return A data frame with the noisy fragments removed. The data frame retains the columns `mz` and `intensity` for the cleaned spectrum.
#'
#' @details
#' The function sorts the spectrum by m/z values and calculates the differences between adjacent fragments. If the difference between any two fragments' m/z values is less than the specified ppm threshold, the fragment with the lower intensity is removed. This helps clean the spectrum by filtering out low-intensity fragments that are close to higher-intensity ones.
#'
#' @examples
#' \dontrun{
#' # Example spectrum
#' spectrum <- data.frame(mz = c(100, 100.001, 150, 200, 250), intensity = c(10, 5, 50, 100, 30))
#'
#' # Remove noise from the spectrum
#' clean_spectrum <- remove_noise(spectrum = spec, ms2.match.ppm = 30, mz.ppm.thr = 400)
#'
#' print(clean_spec)
#' }
#'
#' @author
#' Xiaotao Shen \email{xiaotao.shen@outlook.com}
#'
#' @export

remove_noise <-
  function(spectrum,
           ms2.match.ppm = 30,
           mz.ppm.thr = 400) {
    if (nrow(spectrum) == 1) {
      return(spectrum)
    }
    spectrum <- spectrum[order(spectrum[, 1]), ]
    mz <- spectrum[, 1]
    mz <- mz[-1]
    diff.mz <- diff(spectrum[, 1])
    mz[which(mz < mz.ppm.thr)] <- mz.ppm.thr
    mz.error <- diff.mz * 10 ^ 6 / mz
    temp.idx <- which(mz.error < ms2.match.ppm)
    if (length(temp.idx) > 0) {
      remove.idx <- lapply(temp.idx, function(idx) {
        c(idx, idx + 1)[which.min(spectrum[c(idx, idx + 1), 2])]
      })
      
      remove.idx <- unique(unlist(remove.idx))
      spectrum <- spectrum[-remove.idx, , drop = FALSE]
    } else {
      return(spectrum)
    }
  }
