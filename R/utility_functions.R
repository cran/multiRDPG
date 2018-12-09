#' Plots object from \code{multiRDPG}
#'
#' @param x multiRDPGfit object from function \code{multiRDPG}
#' @param ... further arguments passed to or from other methods
#'
#'
#' @author Agnes Martine Nielsen (agni@@dtu.dk)
#'
#' @seealso \code{\link{multiRDPG}}
#'
#' @export
#' @importFrom graphics plot
plot.multiRDPGfit <- function(x,...){
 plot(x$objfun[-1],xlab = "Iterations", ylab = "Objective Function",type = "b",main = x$call,...)
}

#' Print object from \code{multiRDPG}
#'
#' @param x multiRDPGfit object from function \code{multiRDPG}
#' @param ... further arguments passed to or from other methods
#'
#'
#' @author Agnes Martine Nielsen (agni@@dtu.dk)
#'
#' @seealso \code{\link{multiRDPG}}
#'
#' @export
#' @importFrom utils head
print.multiRDPGfit <- function(x,...){
  cat("Call:\n")
  print(x$call,...)

  cat("\nNumber of graphs: ", length(x$Lambda),"\n")
  cat("Converged: ")
  if(x$converged == 1) cat("Yes\n") else cat("No\n")
  cat("Number of iterations: ",x$iter,"\n")

  cat("\nEstimated U:\n")
  headU<-head(x$U)
  print(headU,...)
  if(dim(x$U)[1] - dim(headU)[1] > 0) cat("[",dim(x$U)[1] - dim(headU)[1]," rows omitted ]\n")

  cat("\nEstimated Lambda: \n")
  print(x$Lambda[1:2],...)
  if(length(x$Lambda)-2 > 0) cat("[",length(x$Lambda)-2," matrices omitted ]\n")
}


#' Plots object from \code{multiRDPG_test}
#'
#' Plots histogram of permutation test statistics and indicates test statistic value with red line.
#'
#' Red line indicates the value of the test statistics with a red line.
#'
#' @param x multiRDPGtest object from function \code{multiRDPG_test}
#' @param ... further arguments passed to or from other methods
#'
#'
#' @author Agnes Martine Nielsen (agni@@dtu.dk)
#'
#' @seealso \code{\link{multiRDPG_test}}
#'
#' @export
#' @importFrom graphics hist lines
plot.multiRDPGtest <- function(x,...){
  hist(x$Tstar, xlab = "Tstar",xlim = c(min(x$Tstar,x$Tval)*0.9,max(x$Tstar,x$Tval)*1.1),main = x$call,...)
  lines(c(x$Tval,x$Tval),c(0,length(x$Tstar)),col = "red")
}


#' Print object from \code{multiRDPG_test}
#'
#' @param x multiRDPGtest object from function \code{multiRDPG_test}
#' @param ... further arguments passed to or from other methods
#'
#'
#' @author Agnes Martine Nielsen (agni@@dtu.dk)
#'
#' @seealso \code{\link{multiRDPG_test}}
#'
#' @export
print.multiRDPGtest <- function(x,...){
  cat("\n")
  cat("\t MultiRDPG Graph Hypotesis Test\n")
  cat("\n")
  cat("data:  ",x$data.name, "\n",...)
  cat("t = ",x$Tval,", p-value = ",x$pvalue,", d = ",x$d ,"\n",...)
  cat("alternative hypothesis: Lambdas are not equal\n")
}
