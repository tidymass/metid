.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    crayon::green(
      "metid,
More information can be found at https://tidymass.github.io/metid/
Authors: Xiaotao Shen (shenxt1990@163.com)
Maintainer: Xiaotao Shen."
    )
# cat(crayon::green(
#   c(
#     "                _    _____  ___ ",
#     " _ __ ___   ___| |_  \\_   \\/   \\",
#     "| '_ ` _ \\ / _ \\ __|  / /\\/ /\\ /",
#     "| | | | | |  __/ |_/\\/ /_/ /_// ",
#     "|_| |_| |_|\\___|\\__\\____/___,'  ",
#     "                                "
#   )
# ), sep = "\n")
  )
}
