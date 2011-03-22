##
##  s t r f i n d . R  Test suite
##


strfind <- pracma::strfind
#findstr <- pracma::findstr

identical(strfind("", "aba"), NULL)
identical(strfind("ab", "aba"), NULL)
identical(strfind("aba", "aba"), 1)
identical(strfind("ababa", "aba"), c(1, 3))
identical(strfind("ababa", "aba", overlap=FALSE), 1)
