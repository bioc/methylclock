.onLoad <- function(libname, pkgname) {
    # Populate the registry with the built-in clocks on load.
    register_builtin_clocks()
}
