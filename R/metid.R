#' @title Show the base information of metid pacakge
#' @description Show the base information of metid pacakge.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @return A ASCII log of metid
#' @export
#' @examples 
#' metid()

metid <- function() {
  message(crayon::yellow(
    "`metid()` is deprecated, please use `metid_logo()`."
  ))  
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
