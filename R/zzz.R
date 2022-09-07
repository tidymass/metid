.onAttach <- function(libname, pkgname) {
  # needed <- core[!is_attached(core)]
  # if (length(needed) == 0)
  #   return()
  #
  crayon::num_colors(TRUE)
  metid_attach()
  #
  # if (!"package:conflicted" %in% search()) {
  #   x <- metid_conflicts()
  #   msg(metid_conflict_message(x), startup = TRUE)
  # }
  msg(paste0("Version ", metid_version, " (", update_date, ')'))
}

is_attached <- function(x) {
  paste0("package:", x) %in% search()
  
}
