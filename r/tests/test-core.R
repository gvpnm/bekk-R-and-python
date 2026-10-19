library(bekkpub)
ns <- asNamespace("bekkpub")
layout <- get(".layout",ns); filter <- get(".filter",ns); radius <- get(".radius",ns)
k <- 2; C <- diag(c(.4,.3)); A <- diag(.2,2); B <- diag(.8,2)
stopifnot(abs(radius(A,B)-.68)<1e-12)
set.seed(8); e <- matrix(rnorm(200),100,2); S <- crossprod(e)/100
l <- layout(2,"diagonal"); cc <- C; diag(cc) <- log(diag(C))
p <- c(cc[l$ci],diag(A),diag(B)); f <- filter(p,e,S,l,TRUE)
expected <- tcrossprod(C)+tcrossprod(A%*%e[1,])+B%*%S%*%t(B)
stopifnot(max(abs(f$H[,,2]-expected))<1e-12)
for(t in c(1,2,50)) {
 H <- f$H[,,t]
 direct <- -.5*(2*log(2*pi)+as.numeric(determinant(H,logarithm=TRUE)$modulus)+sum(e[t,]*(solve(H)%*%e[t,])))
 stopifnot(abs(direct-f$loglik[t])<1e-10,
           max(abs(drop(t(chol(H))%*%f$z[t,])-e[t,]))<1e-10)
}
x <- bekk_simulate(400,C,A,B,seed=99)
stopifnot(identical(x,bekk_simulate(400,C,A,B,seed=99)))
set.seed(19); before <- .Random.seed
invisible(bekk_simulate(30,C,A,B,seed=2)); stopifnot(identical(before,.Random.seed))
fit <- bekk_fit(x,type="diagonal",starts=2,se="robust")
stopifnot(fit$convergence==0,fit$radius<1,all(is.finite(coef(fit))),all(is.finite(vcov(fit))),
          all(diag(vcov(fit))>0),length(coef(fit))==7)
tab <- bekk_table(fit); d <- bekk_diagnostics(fit,5)
stopifnot(nrow(tab)==7,nrow(d)==10,all(d$p_value>=0 & d$p_value<=1),
          all(d$p_adjusted>=d$p_value),all(tab$conf_low<tab$conf_high))
b <- stats::Box.test(residuals(fit)[,1],lag=5,type="Ljung-Box")
stopifnot(abs(d$statistic[1]-b$statistic)<1e-10)
html <- bekk_report(fit,format="html",lags=5)
tex <- bekk_report(fit,format="latex",lags=5)
stopifnot(grepl("<table>",html,fixed=TRUE),grepl("\\toprule",tex,fixed=TRUE),grepl("\\\\",tex,fixed=TRUE))
escape <- get(".escape",ns)
stopifnot(escape("A&B_<", "html")=="A&amp;B_&lt;",escape("A&B_", "latex")=="A\\&B\\_")
stopifnot(inherits(try(bekk_fit(matrix(NA_real_,40,2)),silent=TRUE),"try-error"),
          inherits(try(bekk_diagnostics(fit,100),silent=TRUE),"try-error"))
# Exercise full-matrix parameter indexing and fit.
full <- bekk_fit(x,type="full",starts=1,se="none",control=list(maxit=1200,reltol=1e-7))
stopifnot(length(coef(full))==11,full$radius<1,all(is.finite(full$H)))
