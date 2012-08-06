###
### wilkinson.R  +++ Test suite +++
###

Wilkinson <- pracma::Wilkinson

identical(Wilkinson(0), NULL)
identical(Wilkinson(1), matrix(0, nrow=1, ncol=1))
identical(Wilkinson(3), matrix(c(1,1,0, 1,0,1, 0,1,1), 3, 3))
