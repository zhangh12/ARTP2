
.onAttach <- function(libname, pkgname) {
  if (interactive()) {
    packageStartupMessage('ARTP2 ', packageVersion('ARTP2'), '  For help type ?rARTP, ?sARTP, ?warm.start, or ?meta. ')
    packageStartupMessage('The most frequently updated version can be downloaded from https://github.com/zhangh12/ARTP2')
  }
}
