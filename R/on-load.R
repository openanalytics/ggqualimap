.onLoad <- function(libname, pkgname) {
    ggqualimap_verbose = "FALSE"
    # global options
    opts = c(
              "ggqualimap_verbose" = ggqualimap_verbose
            )
    for (i in setdiff(names(opts), names(options())) ) {
        text = paste('options(', i, '=', opts[i], ')', sep="")
        eval(parse(text=text))
    }
    invisible()
}
