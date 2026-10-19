#' Coefficient estimates and asymptotic inference
bekk_table <- function(object,level=.95) {
  if(!inherits(object,"bekk_fit")) stop("object must be a bekk_fit")
  if(length(level)!=1 || !is.finite(level) || level<=0 || level>=1) stop("level must lie between zero and one")
  est <- coef(object); se <- sqrt(pmax(diag(vcov(object)),0)); z <- est/se
  p <- 2*stats::pnorm(-abs(z)); critical <- stats::qnorm((1+level)/2)
  data.frame(term=names(est),estimate=unname(est),std_error=unname(se),z=unname(z),
    p_value=unname(p),conf_low=unname(est-critical*se),conf_high=unname(est+critical*se),
    significance=ifelse(is.na(p),"",ifelse(p<.001,"***",ifelse(p<.01,"**",ifelse(p<.05,"*","")))),row.names=NULL)
}
.escape <- function(x,format) {
  if(format=="html") {
    x <- gsub("&","&amp;",x,fixed=TRUE); x <- gsub("<","&lt;",x,fixed=TRUE); return(gsub(">","&gt;",x,fixed=TRUE))
  }
  chars <- strsplit(x,"",fixed=TRUE)
  vapply(chars,function(ch) paste0(vapply(ch,function(c) switch(c,
    "\\"="\\textbackslash{}","&"="\\&","%"="\\%","$"="\\$","#"="\\#",
    "_"="\\_","{"="\\{","}"="\\}","~"="\\textasciitilde{}","^"="\\textasciicircum{}",c),character(1)),collapse=""),character(1))
}
.render <- function(d,caption,format,digits) {
  d[] <- lapply(d,function(v) if(is.numeric(v)) ifelse(is.na(v),"NA",formatC(v,digits=digits,format="f")) else as.character(v))
  headers <- .escape(names(d),format); cap <- .escape(caption,format)
  if(format=="html") {
    rows <- apply(d,1,function(v) paste0("<tr>",paste0("<td>",.escape(v,format),"</td>",collapse=""),"</tr>"))
    paste0("<table><caption>",cap,"</caption><thead><tr>",paste0("<th scope=\"col\">",headers,"</th>",collapse=""),"</tr></thead><tbody>",paste(rows,collapse="\n"),"</tbody></table>")
  } else {
    rows <- apply(d,1,function(v) paste(.escape(v,format),collapse=" & "))
    paste0("\\begin{table}[htbp]\n\\centering\n\\caption{",cap,"}\n\\begin{tabular}{l",paste(rep("r",ncol(d)-1),collapse=""),"}\n\\toprule\n",paste(headers,collapse=" & ")," \\\\\n\\midrule\n",paste0(rows," \\\\",collapse="\n"),"\n\\bottomrule\n\\end{tabular}\n\\end{table}")
  }
}
#' Export a complete publication report
bekk_report <- function(object,file=NULL,format=c("html","latex"),lags=10L,digits=4L) {
  format <- match.arg(format)
  if(length(digits)!=1 || !is.finite(digits) || digits<1 || digits>10 || digits!=as.integer(digits)) stop("digits must be an integer from 1 to 10")
  ct <- bekk_table(object); dg <- bekk_diagnostics(object,lags)
  fit <- data.frame(metric=c("Observations","Parameters (including estimated means)","Log likelihood","AIC","BIC","Covariance persistence","Optimizer convergence code"),
    value=c(object$nobs,object$npar,object$loglik,stats::AIC(object),stats::BIC(object),object$radius,object$convergence))
  notes <- paste0("BEKK(1,1), ",object$type,". Standard errors: ",object$inference,
    ". C is lower triangular with positive diagonal. H[t] = C C' + A e[t-1] e[t-1]' A' + B H[t-1] B'. ",
    "Stationarity requires spectral radius of A tensor A + B tensor B < 1. ",
    "Initial covariance is the sample residual second moment. Means are estimated separately; reported covariance-parameter SEs condition on those means. ",
    "Intervals are 95% asymptotic Wald intervals. Stars: * p<0.05, ** p<0.01, *** p<0.001. ",
    "Displayed p-values rounded to zero are not exactly zero. ",attr(dg,"notes"))
  parts <- c(.render(ct,"BEKK parameter estimates",format,digits),.render(fit,"Model fit and stationarity",format,digits),.render(dg,"Residual diagnostics",format,digits))
  if(format=="html") output <- paste0('<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>BEKK estimation report</title><style>body{font-family:Georgia,serif;max-width:1200px;margin:40px auto;padding:20px}table{border-collapse:collapse;width:100%;margin:30px 0;font-variant-numeric:tabular-nums}caption{font-weight:bold;text-align:left;margin-bottom:10px}th,td{padding:8px 12px;text-align:right}th:first-child,td:first-child{text-align:left}thead{border-top:2px solid;border-bottom:1px solid}tbody{border-bottom:2px solid}p{line-height:1.5}</style></head><body><h1>BEKK estimation report</h1>',paste(parts,collapse="\n"),"<p>",.escape(notes,format),"</p></body></html>")
  else output <- paste0("% Requires \\usepackage{booktabs}; use landscape or resizebox for wide tables.\n",paste(parts,collapse="\n\n"),"\n\\par\\small ",.escape(notes,format),"\n")
  if(!is.null(file)) writeLines(output,file,useBytes=TRUE)
  invisible(output)
}
