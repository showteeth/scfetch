.onLoad <- function(libname, pkgname) {
  # set global variables in order to avoid CHECK notes
  utils::globalVariables("PanglaoDBComposition")
  utils::globalVariables("PanglaoDBMeta")
  invisible()
}
