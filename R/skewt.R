#' Get differential expression p-values after fitting with a skew-t distribution. 
#'
#' This function fits the intensity values of a set of genes with a 
#' skew-t distribution and calculates the p-values of upregulation.
#'
#' @param dat numeric vector which holds the intensity values of all genes
#' @export
#' @return p-values
skewt.pval <- function(dat) {
    fit <- st.mle(y=dat)
    #xfit <- seq(min(dat), max(dat), length=100)
    #yfit <- dst(xfit, dp=fit$dp)
    pval <- 1 - pst(dat, dp=fit$dp)
    pval
}
