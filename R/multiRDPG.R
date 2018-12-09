#' Calculates positive semi-definite part of a matrix
#'
#' \code{aplusfun} calculates the positive semidefinite part of a matrix
#'
#' @param Amat A symmetric matrix
#'
#' @return Returns positive semidefinite part of Amat
#'
#' @author Agnes Martine Nielsen (agni@@dtu.dk)
#'
#' @keywords internal
aplusfun<-function(Amat){
  eigA <- eigen(Amat)
  Aplus <- eigA$vectors %*% diag(pmax(eigA$values, 0)) %*% t(eigA$vectors)
  return(Aplus)
}

#' \code{updateU} calculates the U update
#'
#' @param AA  All Aplus matrices column binded
#' @param Lambda List of Lambdas
#' @param Uold U found in past iteration
#'
#' @return Returns updated U
#'
#' @author Agnes Martine Nielsen (agni@@dtu.dk)
#'
#' @keywords internal
updateU<-function(AA,Lambda,Uold){
  LUold <- do.call(cbind, lapply(Lambda, function(x){x %*% t(Uold)}))
  svdsol <- svd(AA %*% t(LUold))
  U<- svdsol$u %*% t(svdsol$v)
  return(U)
}

#' \code{updateL} updates all Lambda
#'
#' @param Aplus List of positive definite parts of A
#' @param U Current U matrix
#' @param d Dimension of latent space
#'
#' @return Returns list of updated Lambdas
#'
#'
#' @author Agnes Martine Nielsen (agni@@dtu.dk)
#'
#' @keywords internal
updateL<-function(Aplus,U,d){
  Z <- lapply(Aplus, function(x){t(U) %*% x %*% U})
  Lambda <- lapply(Z, function(x){diag(pmax(0, diag(x)))}) #max because of non-negative constraint

  return(Lambda)
}

#' \code{calculateobjfun} calculates the objective function
#'
#' @param Aplus List of positive semi definite part of A matrices
#' @param Lambda List of Lambdas
#' @param U U matrix that contains the othogonal vectors
#'
#' @return Returns the value of the objective function, sum(||Aplus-U Lambda U^T||_F^2)
#'
#' @author Agnes Martine Nielsen (agni@@dtu.dk)
#'
#' @keywords internal
calculateobjfun <-function(Aplus,Lambda,U){
  obj <- sum(mapply(function(X,Y){norm(X - U %*% Y %*% t(U),"F")^2}, Aplus, Lambda))
  return(obj)
}

#' Fitting Multiple Random Dot Product Graphs
#'
#'\code{multiRDPG} is used to fit Multiple Random Dot Product Graphs from a set of adjacency matrices.
#'
#'
#' @param A List of adjacency matrices representing graphs. Each matrix must be symmetric. All matrices of the same size n x n.
#' @param d Dimension of latent space. d<= n.
#' @param maxiter Maximal number of iterations. Default is 100.
#' @param tol Tolerance for update of the objective function. Default is 1e-6.
#'
#' @return
#'
#' Returns a list of the following:\cr
#' \tabular{ll}{
#'             \code{U}         \tab Matrix of the joint vectors. n x d. \cr
#'             \code{Lambda}    \tab List of diagonal matrices. One for each graph. d x d.\cr
#'             \code{Converged} \tab Represent of the algorithm converged. 1 if converged, 0 if not.\cr
#'             \code{iter}      \tab Number of iterations \cr
#'             \code{maxiter}   \tab Maximal number of iterations.Default is 100.\cr
#'             \code{objfun}    \tab Value of the objective function. sum_k ||A^k - U Lambda U^T||_F^2
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
#' #fit model
#' multiRDPG(A,3)
#'
#'
#' @author Agnes Martine Nielsen (agni@@dtu.dk)
#'
#' @seealso \code{\link{multiRDPG_test}}
#'
#' @export
multiRDPG<- function(A,d,maxiter=100, tol=1e-6){

  # Checks of A
  if(!all(sapply(A, function(x){identical(dim(x), dim(A[[1]]))}))) stop("All matrices must be of the same dimensions")
  if(!all(sapply(A, isSymmetric))) stop("All matrices must be symmetric")
  if(!all(sapply(A, function(x){min(x%in%c(0,1))}))) stop("All matrix elements must be 0 or 1")

  n <- dim(A[[1]])[[1]]
  if(d>n) stop("d must be less than or equal to n")

  k <- length(A)

  nn <- nullestimation(A, d)
  Lambda <- list()
  for(ll in 1:k) Lambda[[ll]] <- nn$Lambda
  U <- nn$U

  # Initialize for checking convergence
  converged <- 0
  iter <- 0
  objfun <- rep(NA, maxiter)
  objfun[1] <- Inf

  # 2) Calculate positive semi definte part of A
  Aplus <- lapply(A,aplusfun)

  # Column bind A for U update (to only do this once)
  AA <- do.call(cbind, Aplus)

  while(!converged && iter < maxiter){
    # 3) Solve for U (procrustes) => min ||[A1, A2] - Uold [LV1 LV2]||_F
    Uold <- U
    U<- updateU(AA, Lambda, Uold)

    # 4) Solve for Lambda1 and Lambda2 (quadratic)
    Lambda <- updateL(Aplus, U, d)

    # 5) Go to 3) until convergence
    objfun[iter+2] <-calculateobjfun(Aplus, Lambda, U)
    if(abs(objfun[iter+1] - objfun[iter+2]) < tol) converged = 1

    iter <- iter + 1
  }

  # Check objective function. If non-decreasing then something is wrong.
  objfun = objfun[!is.na(objfun)]
  if(max(diff(objfun)) > 1e-10) stop("Objective is non-decreasing")

  # 6) Output U, Lambda
  out<-list(U = U,Lambda = Lambda, converged = converged, iter = iter,
            maxiter = maxiter, objfun = objfun, call = match.call(),tol = tol)

  class(out) <- "multiRDPGfit"
  return(out)
}
