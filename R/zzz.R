
#' @import TMB
#' @useDynLib motus

.onUnload <- function(lib) {
    library.dynam.unload("motus", lib)
}
