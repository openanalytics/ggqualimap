.onAttach <- function(libname, pkgname) {
    # Runs when attached to search() path such as by library() or require()
    if (interactive()) {
        packageStartupMessage('v', as.character(packageVersion("ggqualimap")), 
        ', type vignette("ggqualimap-vignette", package="ggqualimap") ", 
        "to get started.')
    }
}
