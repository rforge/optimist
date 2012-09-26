##
##  m a t l a b . R  Matlab Idioms
##


matlab <- function() {
    assign("diag", pracma::Diag, envir = .GlobalEnv)
    assign("find", pracma::finds, envir = .GlobalEnv)
    assign("fix", pracma::ml.fix, envir = .GlobalEnv)
    assign("imag", pracma::Imag, envir = .GlobalEnv)
    assign("lcm", pracma::ml.lcm, envir = .GlobalEnv)
    assign("mode", pracma::Mode, envir = .GlobalEnv)
    assign("norm", pracma::Norm, envir = .GlobalEnv)
    assign("poly", pracma::Poly, envir = .GlobalEnv)
    assign("real", pracma::Real, envir = .GlobalEnv)
    assign("reshape", pracma::Reshape, envir = .GlobalEnv)
    assign("strtrim", pracma::strTrim, envir = .GlobalEnv)
    assign("toeplitz", pracma::Toeplitz, envir = .GlobalEnv)
    # assign("", pracma::, envir = .GlobalEnv)
    invisible(NULL)
}
