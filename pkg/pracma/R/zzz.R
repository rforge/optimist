##
##  z z z . R
##


# Define 'ans' as in Matlab
# But: "Package namespaces are locked when loaded!"
# makeActiveBinding("ans", function() .Last.value, .GlobalEnv)

.pracmaEnv <- new.env()
assign("elapsedTime", 0, envir = .pracmaEnv)

.onLoad <- function(libname, pkgname) {
    # require(some_packages)

    # Load dynamic libraries
    # library.dynam(pkg, pkg, lib)

    environment(.pracmaEnv) <- asNamespace("pracma")

    packageStartupMessage(
    paste("PRACMA 1.2.0 (Sept-2012 HwB): Practical Numerical Math Functions.",    
          "'pracma::matlab()' will redefine function names to Matlab idioms.",
          sep = "\n")
    )
}