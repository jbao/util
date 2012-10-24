#' Calculate the p-value of the likelihood-ratio test between the full and reduced
#' model. 
#'
#' This function calculates the likelihood-ratio test p-value of comparing the
#' full and reduced model fitting to the time series data.
#'
#' @param core the data vector, consists of time series under both conditions
#' @param ctrl.tpoint time points of the control measurement
#' @param treated.tpoint time points of the treatment measurement
#' @export
#' @return p-value
lrtest <- function(core, ctrl.tpoint, treated.tpoint) {
    design.ctrl <- cbind(rep(1,length(ctrl.tpoint)), ctrl.tpoint, 
        ctrl.tpoint^2, ctrl.tpoint^3, matrix(0,nrow=length(ctrl.tpoint),ncol=4))
    design.treated <- cbind(matrix(0,nrow=length(treated.tpoint),ncol=4), 
        rep(1,length(treated.tpoint)), treated.tpoint, treated.tpoint^2, 
        treated.tpoint^3)
    design.full <- rbind(design.ctrl, design.treated)
    design.red <- cbind(rep(1,length(ctrl.tpoint)+length(treated.tpoint)), 
        c(ctrl.tpoint,treated.tpoint), c(ctrl.tpoint^2,treated.tpoint^2), 
        c(ctrl.tpoint^3,treated.tpoint^3))
    colnames(design.full) <- c(paste("a", 0:3, sep=""), paste("b", 0:3, sep=""))
    colnames(design.red) <- paste("d", 0:3, sep="")
    mod.full <- lm(core ~ a0 + a1 + a2 + a3 + b0 + b1 + b2 + b3 - 1,
        dat=data.frame(design.full))
    mod.red <- lm(core ~ d0 + d1 + d2 + d3 - 1, dat=data.frame(design.red))
    test <- anova(mod.red, mod.full)
    test$'Pr(>F)'[2]
}

#' Calculate the p-value of the likelihood-ratio test between the full and reduced
#' model. 
#'
#' This function calculates the likelihood-ratio test p-value of comparing the
#' full and reduced model fitting to the time series data.
#'
#' @param core the data vector, consists of time series under both conditions
#' @param ctrl.tpoint time points of the control measurement
#' @param treated.tpoint time points of the treatment measurement
#' @export
#' @return p-value
lrtest.replicate <- function(df) {
    idx.ctrl <- which(as.numeric(df$cond)==1)
    design.ctrl <- cbind(rep(1,length(idx.ctrl)), df$time[idx.ctrl], 
        df$time[idx.ctrl]^2, df$time[idx.ctrl]^3, 
        matrix(0,nrow=length(idx.ctrl),ncol=4))
    idx.treated <- which(as.numeric(df$cond)==2)
    design.treated <- cbind(matrix(0,nrow=length(idx.treated),ncol=4), 
        rep(1,length(idx.treated)), df$time[idx.treated], df$time[idx.treated]^2, 
        df$time[idx.treated]^3)
    design.full <- rbind(design.ctrl, design.treated)
    design.red <- cbind(rep(1,nrow(df)), df$time, df$time^2, df$time^3)
    colnames(design.full) <- c(paste("a", 0:3, sep=""), paste("b", 0:3, sep=""))
    colnames(design.red) <- paste("d", 0:3, sep="")
    mod.full <- lm(df$val ~ a0 + a1 + a2 + a3 + b0 + b1 + b2 + b3 - 1,
        dat=data.frame(design.full))
    mod.red <- lm(df$val ~ d0 + d1 + d2 + d3 - 1, dat=data.frame(design.red))
    test <- anova(mod.red, mod.full)
    test$'Pr(>F)'[2]
}

#' Fit the time series to a reduced model assuming no differential expression
#' between conditions.
#'
#' This function fits the given time series to a reduced cubic polynomial model
#' by assuming there is no differential expression between conditions
#'
#' @param core the data vector, consists of time series under both conditions
#' @param ctrl.tpoint time points of the control measurement
#' @param treated.tpoint time points of the treatment measurement
#' @export
#' @return the fitted curve
fit.red <- function(core, ctrl.tpoint, treated.tpoint) {
    design.red <- cbind(rep(1,length(ctrl.tpoint)+length(treated.tpoint)), 
        c(ctrl.tpoint,treated.tpoint), c(ctrl.tpoint^2,treated.tpoint^2), 
        c(ctrl.tpoint^3,treated.tpoint^3))
    colnames(design.red) <- paste("d", 0:3, sep="")
    mod.red <- lm(core ~ d0 + d1 + d2 + d3 - 1, dat=data.frame(design.red))
    predict(mod.red)
}

#' Fit the time series to a full model assuming differential expression
#' between conditions.
#'
#' This function fits the given time series to a full cubic polynomial model
#' by assuming there is differential expression between conditions
#'
#' @param core the data vector, consists of time series under both conditions
#' @param ctrl.tpoint time points of the control measurement
#' @param treated.tpoint time points of the treatment measurement
#' @export
#' @return the fitted curve
fit.full <- function(core, ctrl.tpoint, treated.tpoint) {
    design.ctrl <- cbind(rep(1,length(ctrl.tpoint)), ctrl.tpoint, 
        ctrl.tpoint^2, ctrl.tpoint^3, matrix(0,nrow=length(ctrl.tpoint),ncol=4))
    design.treated <- cbind(matrix(0,nrow=length(treated.tpoint),ncol=4), 
        rep(1,length(treated.tpoint)), treated.tpoint, treated.tpoint^2, 
        treated.tpoint^3)
    design.full <- rbind(design.ctrl, design.treated)
    colnames(design.full) <- c(paste("a", 0:3, sep=""), paste("b", 0:3, sep=""))
    mod.full <- lm(core ~ a0 + a1 + a2 + a3 + b0 + b1 + b2 + b3 - 1,
        dat=data.frame(design.full))
    predict(mod.full)
}

#' Calculate the p-value of the likelihood-ratio test between the full and reduced
#' model. 
#'
#' This function calculates the likelihood-ratio test p-value of comparing the
#' full and reduced model fitting to the time series data.
#'
#' @param core the data vector, consists of time series under both conditions
#' @param ctrl.tpoint time points of the control measurement
#' @param treated.tpoint time points of the treatment measurement
#' @export
#' @return p-value
fit.full.replicate <- function(df) {
    idx.ctrl <- which(as.numeric(df$cond)==1)
    design.ctrl <- cbind(rep(1,length(idx.ctrl)), df$time[idx.ctrl], 
        df$time[idx.ctrl]^2, df$time[idx.ctrl]^3, 
        matrix(0,nrow=length(idx.ctrl),ncol=4))
    idx.treated <- which(as.numeric(df$cond)==2)
    design.treated <- cbind(matrix(0,nrow=length(idx.treated),ncol=4), 
        rep(1,length(idx.treated)), df$time[idx.treated], df$time[idx.treated]^2, 
        df$time[idx.treated]^3)
    design.full <- rbind(design.ctrl, design.treated)
    colnames(design.full) <- c(paste("a", 0:3, sep=""), paste("b", 0:3, sep=""))
    mod.full <- lm(df$val ~ a0 + a1 + a2 + a3 + b0 + b1 + b2 + b3 - 1,
        dat=data.frame(design.full))
    err <- stats::predict(mod.full, se = TRUE)
    res <- data.frame(time=df$time, cond=df$cond, replicate=df$replicate, 
        fit=err$fit)
    res$ucl <- err$fit + 1.96 * err$se.fit
    res$lcl <- err$fit - 1.96 * err$se.fit
    res
}

