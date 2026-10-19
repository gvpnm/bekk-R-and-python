#' Residual diagnostic tests
bekk_diagnostics <- function(object,lags=10L) {
  if(!inherits(object,"bekk_fit")) stop("object must be a bekk_fit")
  z <- residuals(object); n <- nrow(z); k <- ncol(z)
  if(length(lags)!=1 || !is.finite(lags) || lags<1 || lags!=as.integer(lags) || lags>=n/5) stop("lags must be a positive integer smaller than n/5")
  rows <- list()
  add <- function(test,series,stat,df,p) {
    rows[[length(rows)+1L]] <<- data.frame(test=test,series=series,statistic=unname(stat),df=df,p_value=unname(p))
  }
  for(i in seq_len(k)) {
    v <- z[,i]
    for(squared in c(FALSE,TRUE)) {
      b <- stats::Box.test(if(squared) v^2 else v,lag=lags,type="Ljung-Box")
      add(if(squared) "Ljung-Box squared residuals" else "Ljung-Box residuals",object$series[i],b$statistic,lags,b$p.value)
    }
    y <- embed(v^2,lags+1); reg <- stats::lm(y[,1]~y[,-1,drop=FALSE])
    lmstat <- nrow(y)*summary(reg)$r.squared
    add("ARCH LM",object$series[i],lmstat,lags,stats::pchisq(lmstat,lags,lower.tail=FALSE))
    u <- v-mean(v); m2 <- mean(u^2); skew <- mean(u^3)/m2^1.5; kurt <- mean(u^4)/m2^2
    jb <- n*(skew^2/6+(kurt-3)^2/24)
    add("Jarque-Bera",object$series[i],jb,2,stats::pchisq(jb,2,lower.tail=FALSE))
  }
  # Hosking (1980) multivariate portmanteau with finite-sample weights.
  port <- function(v) {
    v <- scale(v,center=TRUE,scale=FALSE); nn <- nrow(v); kk <- ncol(v)
    s <- crossprod(v)/nn
    if(qr(s)$rank<kk) return(NA_real_)
    inv <- solve(s); q <- 0
    for(h in seq_len(lags)) {
      ch <- crossprod(v[(h+1):nn,,drop=FALSE],v[seq_len(nn-h),,drop=FALSE])/nn
      q <- q+sum(diag(ch%*%inv%*%t(ch)%*%inv))/(nn-h)
    }
    nn^2*q
  }
  for(squared in c(FALSE,TRUE)) {
    q <- port(if(squared) z^2 else z); df <- k*k*lags
    add(if(squared) "Multivariate portmanteau squared" else "Multivariate portmanteau", "Joint",q,df,stats::pchisq(q,df,lower.tail=FALSE))
  }
  result <- do.call(rbind,rows)
  result$p_adjusted <- stats::p.adjust(result$p_value,method="holm")
  attr(result,"notes") <- paste("Asymptotic diagnostic p-values use nominal degrees of freedom without correction for estimated BEKK parameters; interpret as screening tests.",
    "Holm adjustment covers all tests in this table. Whitening uses lower Cholesky factors and is series-order dependent.",
    "Normality rejection concerns the Gaussian distribution assumption, not necessarily covariance adequacy.")
  result
}
