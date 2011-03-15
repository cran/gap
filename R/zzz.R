.onLoad <- function(lib, pkg) print("R/gap is loaded")
.noGenerics <- TRUE
.onUnload <- function(libpath) library.dynam.unload("gap", libpath)
