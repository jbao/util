# $Id: full_reduced.R 313 2012-05-10 08:26:08Z jbao $
# requires predefinition of ctrl.tpoint and treated.tpoint
lrtest <- function(core) {
    #core <- c(as.numeric(ctrl), as.numeric(treated))
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

fit.red <- function(core) {
    design.red <- cbind(rep(1,length(ctrl.tpoint)+length(treated.tpoint)), 
        c(ctrl.tpoint,treated.tpoint), c(ctrl.tpoint^2,treated.tpoint^2), 
        c(ctrl.tpoint^3,treated.tpoint^3))
    colnames(design.red) <- paste("d", 0:3, sep="")
    mod.red <- lm(core ~ d0 + d1 + d2 + d3 - 1, dat=data.frame(design.red))
    predict(mod.red)
}

fit.full <- function(core) {
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

