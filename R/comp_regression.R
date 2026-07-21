#' @title Compare regression coefficients
#' @author Kara J.M. Taylor
#' @description Use manual T-test calculation to test differences in regression coefficients
#' @param m1 model 1
#' @param m2 model 2
#' @param coef Either 'slope' or 'intercept'
#' @export


compare_regressions <- function(m1, m2, coef){
  
  if (coef == "slope"){
    coef1 <- m1[[4]][2]
    coef2 <- m2[[4]][2]
  }
  if (coef == "intercept"){
    coef1 <- m1[[4]][1]
    coef2 <- m2[[4]][1]
  }
  
  rr1 <- m1[[9]]
  rr2 <- m2[[9]]
  
  df1 <- m1[[10]][3]
  df2 <- m2[[10]][3]
  
  v1 <- (coef1^2*(1-rr1)) / ((df1*rr1))
  v2 <- (coef2^2*(1-rr2)) / ((df2*rr2))
  
  ttest <- abs((coef1/sqrt(rr1)) - (coef2/sqrt(rr2))) / sqrt(v1+v2)
  
  tp <- 2*pt(ttest, df=df1+df2, lower=F)
  
  cat(paste0("t= ", round(ttest,3), ", p= ",round(tp,4)))
}