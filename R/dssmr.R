#' Double Structured Sparse Multi-task Regression
#'
#' This is the function that used for the Double Structured Sparse Multi-task Regression (DSSMR) model.
#' The core algorithm is based on the FISTA optimization method. The model deal with the multi-task regression model
#' with \code{py} tasks, and for each task, the data has \code{px} features and sample size \code{n}.
#'
#'
#' @param x The list of the design matrix with length \code{py}, i.e. the number of tasks. Each element is
#'  the design matrix with dimention \code{n} by \code{px}, where \code{n} is the sample size, and \code{px} is the number of features
#'  in each task.
#' @param y The matrix of the response variables, with dimension \code{n} by \code{py}, where \code{n} is the sample size and
#'  \code{py} is the number of task, and each column corresponse to each task.
#' @param Lpx The penalty matrix corresponds to the independent variables for the second order term with dimension \code{px} by \code{px}. Default is identity matrix.
#' @param Lpy The penalty matrix corresponds to the response (dependent) variables for the second order term with dimension \code{py} by \code{py}. Default is identity matrix.
#' @param lambda The parameter for the l1 term.
#' @param mu The parameter for the (two) l2 terms.
#' @param alpha The parameter for combining the \code{Lpx} (1-\code{alpha}) and \code{Lpy} (\code{alpha}).
#' @param stepsize Parameter for the FISTA algorithm, either "fixed" or "backtracking".
#' @param control The parameter list for the computational algorithm. L: Lipschitz constant of the gradient;
#'  use.gram: whether the Gram matrix t(X) %*% X should be computed stored. maxiter: an upper bound on the number of iterations.
#'  tol: If successive iterates differ less than tol (w.r.t. ell_2-norm), the algorithm stops. init: starting value. Default is the zero vector.
#'  sigma: parameter for back-tracking.
#' @return A list, \code{beta} is the estimated coefficient matrix and \code{fit} is the list of the fitting results:
#'  beta: The estimate as a long vector. iter: number of iterations run. obj: objective values of the iterates. objbetahat: the objective value of the estimate betahat. kkt: KKT optimality. tol: difference of the last two successive iterates (w.r.t. ell_2-norm). L: Lipschitz constant of the gradient or an estimate thereof (if stepsize was set to 'backtracking').
#' @examples
#' n <- 600
#' px <- 10
#' py <- 10
#' alpha <- 0.5
#' Lpx <- crossprod(sl(px))
#' Lpy <- crossprod(sl(py))
#' Lx <- kronecker(Lpx,diag(rep(1,py)));
#' Ly <- kronecker(Lpy,diag(rep(1,px)))
#' P <- matrix(0,px*py,px*py)
#' for(i in 1:px){ for(j in 1:py){ P[(i-1)*py+j,(j-1)*px+i] <- 1 } } # Permutation matrix
#'
#' set.seed(1)
#' X <- matrix(runif(n*px),n,px)
#' beta0 <- 1:max(px,py)
#' beta <- toeplitz(beta0*(beta0>(max(px,py)/2)))[1:px,1:py]/max(px,py)
#' Y <- X%*%beta+matrix(rnorm(n*py,0,0.01),n,py)
#' x <- list()
#' for(i in 1:py){
#'   x[[i]] <- X
#' }
#' Lambda <- alpha*Ly+(1-alpha)*t(P)%*%Lx%*%P
#'
#' fit1 <- dssmr(x=x,y=c(Y),Lpx=Lpx,Lpy=Lpy,lambda=10,mu=0.1,alpha=0.5,stepsize="backtracking",
#'               control=list(L = NULL, use.gram = TRUE, maxiter = 5000, tol = 1e-5, init = NULL, sigma = 0.9))
#' norm(fit1$beta-beta,"F")
#' @export

dssmr=function(x,y,Lpx=NULL,Lpy=NULL,lambda=10,mu=0.1,alpha=0.5,stepsize="fixed",
               control=list(L = NULL, use.gram = TRUE, maxiter = 5000, tol = 1e-3, init = NULL, sigma = 0.9)){
  py=length(x)
  px=dim(x[[1]])[2]
  if(is.null(Lpx)){
    Lpx=diag(rep(1,px))
  }
  if(is.null(Lpy)){
    Lpy=diag(rep(1,py))
  }
  Lx=kronecker(Lpx,diag(rep(1,py))); Ly=kronecker(Lpy,diag(rep(1,px)))
  P=matrix(0,px*py,px*py)
  for(i in 1:px){ for(j in 1:py){ P[(i-1)*py+j,(j-1)*px+i]=1 } } # Permutation matrix
  Lambda=alpha*Ly+(1-alpha)*t(P)%*%Lx%*%P

  datax=as.matrix(Matrix::bdiag(x))
  datay=c(y)

  out = senet_fista(datax, datay, Lambda, lossfun=l2loss, gradfun=l2grad, lambda1=lambda, lambda2=mu, stepsize=stepsize, control = control)
  Bhat=matrix(out$beta,byrow = FALSE,nrow = px, ncol = py)

  return(list(fit=out,beta=Bhat))
}



