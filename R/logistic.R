fit.logistic <- function(tpoint, ts) {
    require(minpack.lm)
    xfit <- seq(0,max(tpoint),length=100)
    ## model based on a list of parameters
    logistic <- function(par, xx) par$a + par$b / (1 + exp(-par$c * (xx - par$d)))
    ## residual function
    residFun <- function(p, observed, xx) observed - logistic(p,xx)
    ## starting values for parameters
    parStart <- list(a=0, b=1, c=1, d=1)
    ## perform fit
    nls.out <- nls.lm(par=parStart, fn = residFun, observed = ts, xx = tpoint)
                                     
    ## plot model evaluated at final parameter estimates
    logistic(as.list(coef(nls.out)), xfit)
}
