
.onAttach <- function(libname, pkgname) {
  if (interactive()) {
    packageStartupMessage('ARTP2 ', packageVersion('ARTP2'), '  For help type ?rARTP, ?sARTP, or ?warm.start. ')
    packageStartupMessage('The most updated version can be downloaded from https://github.com/zhangh12/ARTP2')
  }
}
