##
##  i n t e g r a l 2 . R
##


integral2 <- function(fun, xmin, xmax, ymin, ymax, sector = FALSE,
                        reltol = 1e-6, abstol = reltol, maxlist = 5000,
                        vectorized = TRUE, singular = FALSE, ...) {
    # check input parameters
    nlist <- floor(maxlist/10)

    # check function and vectorization
    fun <- match.fun(fun)
    FUN <- function(x, y) fun(x, y, ...)

    phiBvar <- function(x) ymin * ones(size(x)[1], size(x)[2])
    phiTvar <- function(x) ymax * ones(size(x)[1], size(x)[2])

    # check borders and redefine
    thetaL <- xmin; thetaR <- xmax
    phiB <- 0; phiT <- 1
    area <- (thetaR - thetaL) * (phiT - phiB)

    # initial quadrature
    Qs   <- .tensor(thetaL, thetaR, phiB, phiT, FUN, phiBvar, phiTvar)
    Qsub <- Qs$qsub; esub <- Qs$esub
    Q    <- sum(Qsub)

    # some more parameters and main list
    eps <- .Machine$double.eps
    tol <- 100 * eps * abs(Q)
    err_ok <- 0
    adjust <- 1
    mainList <- zeros(nlist, 7)
    nList <- 0

    # save info on rectangles in main list
    s2l <- .save2list(mainList, nList, Qsub, esub, thetaL, thetaR, phiB, phiT,
                        tol, area, adjust, err_ok)
    mainList <- s2l$mlist
    nList    <- s2l$nlist
    errbnd   <- s2l$errbnd

    if (nList == 0 || errbnd <= tol)
        return(list(Q = Q, error = errbnd))

    while (TRUE) {
        ne <- .nextEntry(mainList, nList)
        mainList <- ne$mlist
        nList <- ne$nlist

        temp <- ne$entry
        q <- temp[1]; e <- temp[2]
        thetaL <- temp[3]; thetaR <- temp[4]
        phiB <- temp[5]; phiT <- temp[6]

        # Approximate integral over four subrectangles
        Qs <- .tensor(thetaL, thetaR, phiB, phiT, FUN, phiBvar, phiTvar)
        Qsub <- Qs$qsub; esub <- Qs$esub

        newq <- sum(Qsub)
        adjust <- min(1, abs(q - newq)/e)
        Q <- Q + (newq - q)
        tol <- max(abstol, reltol * abs(Q)) / 8  
        tol <- max(tol, 100 * eps * abs(Q))

        s2l <- .save2list(mainList, nlist, Qsub, esub, thetaL, thetaR, phiB, phiT,
                            tol, area, adjust, err_ok)
        mainList <- s2l$mlist
        nList    <- s2l$nlist
        errbnd   <- s2l$errbnd

        if (nList == 0 || errbnd <= tol) {
            break
        } else if (nList > maxlist) {
            if (errbnd > max(abstol, max(100*eps, reltol) * abs(Q))) {
                stop("Maximum number of subintervals: w/o convergence.")
            } else {
                warning("Maximum number of subintervals: maybe low accuracy.")
            }
            break
        }
    }
    # accuracy or max. number of subdivisions rached
    return(list(Q = Q, error = errbnd))
}


.tensor <- function(thetaL, thetaR, phiB, phiT, FUN, phiBvar, phiTvar)
{
    # Gauss-Kronrod (3,7) pair with degrees of precision 5 and 11
    nodes <- c( -0.9604912687080202, -0.7745966692414834, -0.4342437493468026,
                0, 0.4342437493468026, 0.7745966692414834, 0.9604912687080202) 
    nnodes <- length(nodes)    
    onevec <- ones(2*nnodes,1)
    narray <- 0.25 * cbind(nodes, nodes)
    wt3    <- c(0, 5/9, 0, 8/9, 0, 5/9, 0)
    wt7    <- c(0.1046562260264672, 0.2684880898683334, 0.4013974147759622,
                0.4509165386584744, 0.4013974147759622, 0.2684880898683334,
                0.1046562260264672)

    Qsub <- zeros(4,1)
    esub <- zeros(4,1)
    
    dtheta <- thetaR - thetaL
    etheta <- thetaL + dtheta * c(0.25,0.75)
    theta  <- c(dtheta*narray + rep(etheta, each = nnodes))
    x <- theta
    X <- onevec %*% x

    dphi <- phiT - phiB;
    ephi <- phiB + dphi * c(0.25,0.75)
    phi  <- c(dphi*narray + rep(ephi, each = nnodes))
    phi  <- as.matrix(c(phi))

    top    <- phiTvar(x)
    bottom <- phiBvar(x)
    dydt <- top - bottom
    t <- phi  
    Y <- onevec %*% bottom + t %*% dydt

    Z <- FUN(X, Y)
    temp <- onevec %*% dydt
    Z <- Z * temp

    # Tensor product: Gauss 3 and 7 points formulae
    esub[1] <- wt3 %*% t(wt3 %*% Z[1:nnodes,1:nnodes])
    esub[2] <- wt3 %*% t(wt3 %*% Z[1:nnodes,(nnodes+1):(2*nnodes)])  
    esub[3] <- wt3 %*% t(wt3 %*% Z[(nnodes+1):(2*nnodes),1:nnodes])
    esub[4] <- wt3 %*% t(wt3 %*% Z[(nnodes+1):(2*nnodes),(nnodes+1):(2*nnodes)])
    esub <- (esub/4)*(dtheta/2)*(dphi/2) 

    Qsub[1] <- wt7 %*% t(wt7 %*% Z[1:nnodes,1:nnodes])
    Qsub[2] <- wt7 %*% t(wt7 %*% Z[1:nnodes,(nnodes+1):(2*nnodes)])
    Qsub[3] <- wt7 %*% t(wt7 %*% Z[(nnodes+1):(2*nnodes),1:nnodes])
    Qsub[4] <- wt7 %*% t(wt7 %*% Z[(nnodes+1):(2*nnodes),(nnodes+1):(2*nnodes)])   

    Qsub <- (Qsub/4)*(dtheta/2)*(dphi/2);  
    esub <- abs(esub - Qsub);

    return(list(qsub = Qsub, esub = esub))
}


.save2list <- function(mainList, nList,
                      Qsub, esub, thetaL, thetaR, phiB, phiT,
                      tol, area, adjust, err_ok) {

    eps <- .Machine$double.eps
    dtheta <- thetaR - thetaL
    dphi <- phiT - phiB
    localtol <- tol * (dtheta/2) * (dphi/2) / area
    localtol <- max(localtol, 100*eps*abs(sum(Qsub)))

    adjerr <- adjust * esub
    if (nList+4 > size(mainList,1))
        mainList <- rbind(mainList, zeros(100, 7))
    
    if (adjerr[1] > localtol) {
        nList <- nList + 1
        mainList[nList, ] <- c(Qsub[1], esub[1], thetaL, thetaL + dtheta/2,
                                             phiB ,phiB + dphi/2, adjerr[1])
    } else {
      err_ok <- err_ok + adjerr[1]
    }

    if (adjerr[2] > localtol) {
      nList <- nList + 1
      mainList[nList, ] <- c(Qsub[2], esub[2], thetaL + dtheta/2, thetaR,
                                           phiB, phiB + dphi/2, adjerr[2])
    } else {
      err_ok <- err_ok + adjerr[2]
    }

    if (adjerr[3] > localtol) {
      nList <- nList + 1
      mainList[nList, ] <- c(Qsub[3], esub[3], thetaL, thetaL + dtheta/2,
                                           phiB + dphi/2, phiT, adjerr[3])
    } else {
      err_ok <- err_ok + adjerr[3]
    }

    if (adjerr[4] > localtol) {
      nList <- nList + 1
      mainList[nList, ] <- c(Qsub[4], esub[4], thetaL + dtheta/2, thetaR,
                                           phiB + dphi/2, phiT, adjerr[4])
    } else {
      err_ok <- err_ok + adjerr[4]
    }   
    
    errbnd <- err_ok + sum(mainList[, 7])
    return(list(mlist = mainList, nlist = nList, errbnd = errbnd))
}


.nextEntry <- function(mainList, nList) {
    indx <- which.max(abs(mainList[1:nList, 7]))
    temp <- mainList[ indx, ]
    mainList <- mainList[-indx, ]
    nList <- nList - 1
    return(list(mlist = mainList, nlist = nList, entry = temp))
}

#-- --------------------------------------------------------------------------

.integral3 <- function(fun, xmin, xmax, ymin, ymax, zmin, zmax) {
    fz <- function(z) {
        fxy <- function(x, y) fun(x, y, z)
        integral2(fxy, xmin, xmax, ymin, ymax)$Q    # singular = TRUE
    }
    romberg(fz, zmin, zmax)$value
}
