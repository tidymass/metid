#' @title Show the base information of metid pacakge
#' @description Show the base information of metid pacakge.
#' \lifecycle{maturing}
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @return A ASCII log of metid
#' @importFrom magrittr %>%
#' @importFrom crayon red yellow green bgRed
#' @importFrom stringr str_detect str_extract str_extract_all
#' @importFrom stringr str_replace_all str_replace str_trim str_c str_count
#' @importFrom readr cols read_csv
#' @importFrom pbapply pblapply pboptions
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom dplyr mutate filter everything select bind_rows left_join pull
#' @importFrom plyr dlply .
#' @importFrom BiocParallel bplapply MulticoreParam SnowParam
#' @importFrom readxl read_xlsx read_xls
#' @importFrom purrr map map2
#' @importFrom ggplot2 aes ggplot geom_point geom_line geom_smooth theme annotate
#' @importFrom ggplot2 geom_abline theme_bw ggsave geom_segment xlim ylim labs
#' @importFrom ggplot2 element_line element_text
#' @importFrom MSnbase readMSData
#' @importFrom ProtGenerics spectra
#' @importFrom tinytools ms2Match get_os mz_rt_match
#' @importFrom stats lm loess predict
#' @importFrom plotly ggplotly
#' @importFrom future plan multisession 
#' @importFrom furrr future_map2 future_map
#' @import lifecycle
#' @import RColorBrewer
#' @import utils
#' @import ggplot2
#' @import methods
#' @import graphics
#' @import grDevices
#' @import utils
#' @importClassesFrom massdataset mass_dataset
#' @export
#' @examples
#' metid_logo()

metid_logo <- function() {
  cat(crayon::green("Thank you for using metid!\n"))
  cat(crayon::green("Version 1.1.0 (20210702)\n"))
  cat(
    crayon::green(
      "More information can be found at https://tidymass.github.io/metid/\n"
    )
  )
  cat(crayon::green(
    c(
      "                _    _____  ___ ",
      " _ __ ___   ___| |_  \\_   \\/   \\",
      "| '_ ` _ \\ / _ \\ __|  / /\\/ /\\ /",
      "| | | | | |  __/ |_/\\/ /_/ /_// ",
      "|_| |_| |_|\\___|\\__\\____/___,'  ",
      "                                "
    )
    
  ), sep = "\n")
}






# library(cowsay)
# # https://onlineasciitools.com/convert-text-to-ascii-art
# # writeLines(capture.output(say("Hello"), type = "message"), con = "ascii_art.txt")
# art <- readLines("logo.txt")
# dput(art)
# metid_logo <-
#   c("                _    _____  ___ ", " _ __ ___   ___| |_  \\_   \\/   \\",
#     "| '_ ` _ \\ / _ \\ __|  / /\\/ /\\ /", "| | | | | |  __/ |_/\\/ /_/ /_// ",
#     "|_| |_| |_|\\___|\\__\\____/___,'  ", "                                "
#   )
# cat(metid_logo, sep = "\n")
