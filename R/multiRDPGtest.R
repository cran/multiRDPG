#' \code{nullestimation} calculates the estimation under the null hypothesis
#'
#' @param A List of symmetric A matrices
#' @param d Dimension of the latent space
#'
#' @return Returns a list of the following\cr
#'  \tabular{ll}{
#'           \code{U}        \tab The common latent space vectors. U in R^{n x d}\cr
#'           \code{Lambda}   \tab List of Lambdas. Each is a positive diagonal matrix of size d x d.
#' }
#'
#' @examples
#' #simulate data
#' U <- matrix(0, nrow=20, ncol=3)
#' U[,1] <- 1/sqrt(20)
#' U[,2] <- rep(c(1,-1), 10)/sqrt(20)
#' U[,3] <- rep(c(1,1,-1,-1), 5)/sqrt(20)
#'
#' L<-list(diag(c(11,6,2)),diag(c(15,4,1)))
#' A <- list()
#' for(i in 1:2){
#'   P <- U%*%L[[i]]%*%t(U)
#'   A[[i]] <-apply(P,c(1,2),function(x){rbinom(1,1,x)})
#'   A[[i]][lower.tri(A[[i]])]<-t(A[[i]])[lower.tri(A[[i]])]
#' }
#'
#' #fit model
#' nullestimation(A,3)
#'
#'
#' @author Agnes Martine Nielsen (agni@@dtu.dk)
#'
#' @seealso \code{\link{multiRDPG}}
#'
#' @export
nullestimation <-function(A,d){
  # Positive semi definite part of A
  Aplus <- lapply(A, aplusfun)

  #Eigendecomposition of average
  meanA <- Reduce("+", Aplus) / length(Aplus)
  eigA <- eigen(meanA)

  U <- eigA$vectors[,1:d]
  Lambda <- diag(eigA$values[1:d])
  out <- list(U = U, Lambda = Lambda)
  return(out)
}

#'\code{calculateLRstat} calculates the likelihood ratio statistic
#'

#' @param A List of symmetric A matrices
#' @param null Model fitted under the null
#' @param alt Model fitted under the alternative
#'
#' @return Value of the likelihood ratio statistic
#'
#'
#' @author Agnes Martine Nielsen (agni@@dtu.dk)
#'
#' @keywords internal
calculateLRstat <-function(A,null,alt){
  Aplus <- lapply(A, aplusfun)

  nullpart <- sum(sapply(Aplus, function(X){norm(X - null$U %*% null$Lambda %*% t(null$U),'F')^2}))
  altpart <-  sum(mapply(function(X, altLambda){norm(X - alt$U %*% altLambda %*% t(alt$U),'F')^2}, Aplus, alt$Lambda))
  Tval <- nullpart - altpart
  return(Tval)
}

#' \code{swapA} swaps the elements of A at random
#'
#' @param A List of symmetric A matrices
#'
#' @return Returns list of matrices with elements swapped
#'
#'
#' @author Agnes Martine Nielsen (agni@@dtu.dk)
#'
#' @importFrom stats runif
#'
#' @keywords internal
swapA <-function(A){
  nn<-dim(A[[1]])
  Ab <-list()
  index <- matrix(runif(prod(nn)) > 0.5, nn[1], nn[2])
  index[lower.tri(index)] <- t(index)[lower.tri(index)]

  Ab[[1]] <- matrix(NA, nn[1], nn[2])
  Ab[[1]][index] <- A[[1]][index]
  Ab[[1]][!index] <- A[[2]][!index]

  Ab[[2]] <- matrix(NA, nn[1], nn[2])
  Ab[[2]][!index] <- A[[1]][!index]
  Ab[[2]][index] <- A[[2]][index]

  return(Ab)
}

#' Performs test based on Multiple Random Dot Product Graph
#'
#' \code{multiRDPG_test} calculates the likelihood ratio test for whether a set of graphs
#' comes from the same disribution.
#'
#' @param A List of symmetric A matrices
#' @param d Dimension of the latent space
#' @param maxiter Maximum number of iterations in the fit of multiRDPG. Default is 100.
#' @param tol Tolerance for the step in the objective function in multiRDPG. Default is 1e-6.
#' @param B Number of permutation iterations. Default is 1000.
#'
#' @return Returns a list of the following elements:
#' \tabular{ll}{
#'            \code{pvalue}     \tab Estimated p-values \cr
#'            \code{Tval}       \tab Value of the test statistic \cr
#'            \code{Tstar}      \tab Vector of the test statistic for each permutation iteration \cr
#'            \code{nullmodel}  \tab Model fit under the null \cr
#'            \code{altmodel}   \tab Modelfit under the alternative
#'}
#'
#' @examples
#' #simulate data
#' U <- matrix(0, nrow=20, ncol=3)
#' U[,1] <- 1/sqrt(20)
#' U[,2] <- rep(c(1,-1), 10)/sqrt(20)
#' U[,3] <- rep(c(1,1,-1,-1), 5)/sqrt(20)
#'
#' L<-list(diag(c(11,6,2)),diag(c(15,4,1)))
#' A <- list()
#' for(i in 1:2){
#'   P <- U%*%L[[i]]%*%t(U)
#'   A[[i]] <-apply(P,c(1,2),function(x){rbinom(1,1,x)})
#'   A[[i]][lower.tri(A[[i]])]<-t(A[[i]])[lower.tri(A[[i]])]
#' }
#'
#' #perform test
#' multiRDPG_test(A,3,B=100)
#'
#'
#' @author Agnes Martine Nielsen (agni@@dtu.dk)
#'
#' @seealso \code{\link{multiRDPG}}
#'
#'@export
multiRDPG_test<-function(A,d,maxiter=100, tol=1e-6,B = 1000){

  # Checks of A
  if(!all(sapply(A, function(x){identical(dim(x), dim(A[[1]]))}))) stop("All matrices must be of the same dimensions")
  if(!all(sapply(A, isSymmetric))) stop("All matrices must be symmetric")
  if(!all(sapply(A, function(x){min(x%in%c(0,1))}))) stop("All matrix elements must be 0 or 1")

  alt <- multiRDPG(A, d, maxiter=maxiter, tol=tol)
  null <- nullestimation(A, d)

  # calculate LR stat
  Tval <-calculateLRstat(A, null, alt)

  # bootstraps
  Tstar <- rep(NA, B)
  for(b in 1:B){

    #Create bootstrap sample
    Ab <- swapA(A)

    #calculate LR stat
    altb <- multiRDPG(Ab, d, maxiter=maxiter, tol=tol)
    nullb <- nullestimation(Ab, d)

    Tstar[b] <- calculateLRstat(Ab, nullb, altb)

  }

  # calculate p value
  pvalue <- sum(Tstar >= Tval) / B

  out<-list(pvalue = pvalue,Tval = Tval, Tstar = Tstar, nullmodel = null, altmodel = alt,
            data.name = deparse(substitute(A)),d = d)

  class(out) <- "multiRDPGtest"
  return(out)
}
