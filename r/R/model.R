# Convention: H[t] = C C' + A e[t-1] e[t-1]' A' + B H[t-1] B'.
.layout <- function(k, type) {
  ci <- which(lower.tri(matrix(0,k,k),diag=TRUE),arr.ind=TRUE)
  ai <- if(type == "diagonal") cbind(seq_len(k),seq_len(k)) else which(matrix(TRUE,k,k),arr.ind=TRUE)
  list(ci=ci,ai=ai,nc=nrow(ci),na=nrow(ai),
       names=c(sprintf("C[%d,%d]",ci[,1],ci[,2]),sprintf("A[%d,%d]",ai[,1],ai[,2]),sprintf("B[%d,%d]",ai[,1],ai[,2])))
}
.unpack <- function(p,k,l) {
  C <- A <- B <- matrix(0,k,k)
  C[l$ci] <- p[seq_len(l$nc)]
  diag(C) <- exp(diag(C))
  A[l$ai] <- p[l$nc+seq_len(l$na)]
  B[l$ai] <- p[l$nc+l$na+seq_len(l$na)]
  list(C=C,A=A,B=B)
}
.radius <- function(A,B) max(Mod(eigen(kronecker(A,A)+kronecker(B,B),only.values=TRUE)$values))
.filter <- function(p,e,S,l,keep=FALSE) {
  k <- ncol(e); n <- nrow(e); m <- .unpack(p,k,l)
  rho <- .radius(m$A,m$B)
  if(!is.finite(rho) || rho >= .9999) return(NULL)
  H <- S; ll <- numeric(n)
  if(keep) { hs <- array(NA_real_,c(k,k,n)); z <- matrix(NA_real_,n,k) }
  for(t in seq_len(n)) {
    if(t>1) H <- tcrossprod(m$C)+tcrossprod(m$A %*% e[t-1,])+m$B %*% H %*% t(m$B)
    H <- (H+t(H))/2
    U <- tryCatch(chol(H),error=function(e) NULL)
    if(is.null(U) || any(!is.finite(U))) return(NULL)
    zz <- backsolve(U,e[t,],transpose=TRUE)
    ll[t] <- -.5*(k*log(2*pi)+2*sum(log(diag(U)))+sum(zz^2))
    if(keep) { hs[,,t] <- H; z[t,] <- zz }
  }
  if(keep) list(loglik=ll,H=hs,z=z,matrices=m,radius=rho) else ll
}
#' Fit a stationary BEKK(1,1) model
bekk_fit <- function(data, type=c("full","diagonal"), demean=TRUE,
                     se=c("robust","hessian","none"), starts=3L,
                     control=list(maxit=1500,reltol=1e-8)) {
  type <- match.arg(type); se <- match.arg(se)
  x <- as.matrix(data)
  if(!is.numeric(x) || ncol(x)<2L || nrow(x)<30L || any(!is.finite(x)))
    stop("data must be a finite numeric matrix with >= 2 columns and >= 30 rows")
  if(length(starts)!=1L || !is.finite(starts) || starts<1 || starts!=as.integer(starts)) stop("starts must be a positive integer")
  k <- ncol(x); n <- nrow(x); nm <- colnames(x)
  if(is.null(nm)) nm <- paste0("Series ",seq_len(k))
  mu <- if(demean) colMeans(x) else rep(0,k)
  e <- sweep(x,2,mu); S <- crossprod(e)/n
  if(min(eigen(S,symmetric=TRUE,only.values=TRUE)$values)<=0) stop("data covariance must be positive definite")
  l <- .layout(k,type); np <- length(l$names)
  if(n<=np+1) stop("too few observations for this model")
  objective <- function(p) { v <- tryCatch(.filter(p,e,S,l),error=function(e) NULL); if(is.null(v) || any(!is.finite(v))) 1e20 else -sum(v) }
  fits <- vector("list",starts)
  for(j in seq_len(starts)) {
    a <- .12+.08*((j-1)%%4); b <- .85-.05*((j-1)%%4)
    C <- t(chol(S*(1-a^2-b^2))); diag(C) <- log(diag(C))
    A <- diag(a,k); B <- diag(b,k)
    p <- c(C[l$ci],A[l$ai],B[l$ai])
    fits[[j]] <- stats::optim(p,objective,method="BFGS",control=control)
  }
  values <- vapply(fits,function(f) f$value,numeric(1)); best <- fits[[which.min(values)]]
  f <- .filter(best$par,e,S,l,TRUE)
  if(is.null(f)) stop("optimization failed to find a stationary model")
  if(best$convergence!=0) warning("optimizer did not converge; inspect convergence and do not rely on inference")
  V <- matrix(NA_real_,np,np); inference <- "not requested"
  if(se!="none") {
    h <- stats::optimHess(best$par,objective)
    ev <- tryCatch(eigen(h,symmetric=TRUE,only.values=TRUE)$values,error=function(e) NA_real_)
    if(all(is.finite(ev)) && min(ev)>max(ev)*1e-10 && best$convergence==0) {
      bread <- solve(h); V <- bread
      if(se=="robust") {
        scores <- matrix(0,n,np)
        for(i in seq_len(np)) {
          step <- 1e-5*(abs(best$par[i])+1); pp <- pm <- best$par
          pp[i] <- pp[i]+step; pm[i] <- pm[i]-step
          lp <- .filter(pp,e,S,l); lm <- .filter(pm,e,S,l)
          if(is.null(lp)||is.null(lm)) { scores[,] <- NA_real_; break }
          scores[,i] <- (lp-lm)/(2*step)
        }
        V <- bread %*% crossprod(scores) %*% bread
      }
      # Delta method from log C diagonal to natural parameter scale.
      d <- rep(1,np); idx <- which(l$ci[,1]==l$ci[,2]); d[idx] <- diag(f$matrices$C)
      V <- V*outer(d,d); inference <- if(all(is.finite(V))) se else "unavailable: boundary scores"
    } else inference <- "unavailable: nonconvergence or non-positive/ill-conditioned Hessian"
    if(!all(is.finite(V))) warning(inference)
  }
  estimates <- c(f$matrices$C[l$ci],f$matrices$A[l$ai],f$matrices$B[l$ai])
  names(estimates) <- l$names; dimnames(V) <- list(l$names,l$names); colnames(f$z) <- nm
  structure(list(coefficients=estimates,vcov=V,matrices=f$matrices,H=f$H,
    standardized_residuals=f$z,residuals=e,mean=mu,series=nm,type=type,
    loglik=sum(f$loglik),nobs=n,npar=np+if(demean) k else 0L,
    convergence=best$convergence,message=best$message,radius=f$radius,
    inference=inference,objective=best$value,call=match.call()),class="bekk_fit")
}
coef.bekk_fit <- function(object,...) object$coefficients
vcov.bekk_fit <- function(object,...) object$vcov
logLik.bekk_fit <- function(object,...) structure(object$loglik,df=object$npar,nobs=object$nobs,class="logLik")
residuals.bekk_fit <- function(object,standardized=TRUE,...) if(standardized) object$standardized_residuals else object$residuals
print.bekk_fit <- function(x,...) {
  cat("BEKK(1,1):",x$type,"| n =",x$nobs,"| logLik =",round(x$loglik,3),"\n")
  cat("Convergence:",x$convergence,"| covariance persistence:",round(x$radius,5),"| SE:",x$inference,"\n")
  print(bekk_table(x)); invisible(x)
}
#' Simulate a stationary BEKK(1,1) process
bekk_simulate <- function(n,C,A,B,burn=500L,seed=NULL) {
  if(!is.null(seed)) {
    existed <- exists(".Random.seed",.GlobalEnv,inherits=FALSE)
    if(existed) old <- get(".Random.seed",.GlobalEnv)
    on.exit(if(existed) assign(".Random.seed",old,.GlobalEnv) else if(exists(".Random.seed",.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv))
    set.seed(seed)
  }
  C <- as.matrix(C); A <- as.matrix(A); B <- as.matrix(B); k <- nrow(C)
  if(k<2 || any(vapply(list(C,A,B),function(m) !is.numeric(m)||!identical(dim(m),c(k,k))||any(!is.finite(m)),logical(1)))) stop("C, A, B must be finite square matrices of identical size, >= 2")
  if(length(n)!=1 || !is.finite(n) || n<1 || n!=as.integer(n) || length(burn)!=1 || !is.finite(burn) || burn<0 || burn!=as.integer(burn)) stop("n and burn must be integer counts")
  if(.radius(A,B)>=1) stop("model is not covariance stationary")
  Q <- tcrossprod(C)
  H <- matrix(solve(diag(k*k)-kronecker(A,A)-kronecker(B,B),as.vector(Q)),k,k)
  if(inherits(try(chol(H),silent=TRUE),"try-error")) stop("unconditional covariance must be positive definite")
  out <- matrix(0,n+burn,k); e <- rep(0,k)
  for(i in seq_len(n+burn)) {
    if(i>1) H <- Q+tcrossprod(A%*%e)+B%*%H%*%t(B)
    e <- drop(t(chol((H+t(H))/2))%*%stats::rnorm(k)); out[i,] <- e
  }
  out[burn+seq_len(n),,drop=FALSE]
}
