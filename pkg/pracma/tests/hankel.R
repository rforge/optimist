##
##  hankel.R  Test
##

Hankel <- pracma::Hankel

identical(Hankel(2), matrix(2, nrow=1, ncol=1))
identical(Hankel(1:3), matrix(c(1,2,3,2,3,0,3,0,0), 3, 3))
identical(Hankel(1:3, 3:1), matrix(c(1,2,3,2,3,2,3,2,1), 3, 3))
identical(Hankel(1:3, 2:1), matrix(c(1,2,3,2,3,1), 3, 2))
identical(Hankel(1:2, 3:1), matrix(c(1,2,2,2,2,1), 2, 3))
