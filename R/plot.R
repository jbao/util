#' Plot expression time series for the given genes. 
#'
#' This function plots the expression profiles of the given genes.
#'
#' @param ctrl list of time series data under the control condition, contains
#' the field dat, cond, tpoint
#' @param treatd list of time series data under the treatment condition, contains
#' the field dat, cond, tpoint
#' @param filename filename
#' @export
plot.full.reduced <- function(ctrl, treated, filename) {
    pdf(filename,8,11)
    for (i in 1:ceiling(nrow(ctrl$dat)/12)) {
        symbol <- ctrl$dat$Gene
        slice <- seq((i-1)*12+1,i*12)
        slice <- slice[!(slice > nrow(ctrl$dat))]
        df.ctrl <- melt(ctrl$dat[slice,])
        df.ctrl <- transform(df.ctrl, cond=ctrl$cond)
        df.ctrl <- transform(df.ctrl, time=rep(ctrl$tpoint,each=length(slice)))
        df.treated <- melt(treated$dat[slice,])
        df.treated <- transform(df.treated, cond=treated$cond)
        df.treated <- transform(df.treated, time=rep(treated$tpoint,
            each=length(slice)))
        df <- rbind(df.ctrl, df.treated)
        df <- within(df, Gene <- factor(Gene, levels = symbol[slice]))
        with(df, levels(Gene))
        df <- within(df, cond <- factor(cond, levels = c(treated$cond,ctrl$cond)))
        with(df, levels(cond))
        
        # full/reduced model fit
        idx.ctrl <- which(names(ctrl$dat)=='Gene')
        idx.treated <- which(names(treated$dat)=='Gene')
        core <- cbind(ctrl$dat[slice,-idx.ctrl],treated$dat[slice,-idx.treated])
        red <- apply(core,1,function(x) fit.red(x,ctrl$tpoint,treated$tpoint))
        fit.reduced <- as.data.frame(red[-(1:length(ctrl$tpoint)),])
        names(fit.reduced) <- symbol[slice]
        full <- apply(core,1,function(x) fit.full(x,ctrl$tpoint,treated$tpoint))
        fit.ctrl <- as.data.frame(full[1:length(ctrl$tpoint),])
        names(fit.ctrl) <- symbol[slice]
        fit.treated <- as.data.frame(full[-(1:length(ctrl$tpoint)),])
        names(fit.treated) <- symbol[slice]

        # get the p-value
        df.param <- apply(core,1,function(x) lrtest(x,ctrl$tpoint,treated$tpoint))
        df.param <- melt(df.param)
        df.param <- transform(df.param, Gene=symbol[slice], 
            label=sprintf("p==%.1e",df.param$value))

        # rearrange the data frame to plot the fitting curve
        fit.reduced <- transform(fit.reduced, time=treated$tpoint)
        fit.reduced <- melt(fit.reduced, id.vars='time')
        fit.reduced <- transform(fit.reduced, Gene=rep(symbol[slice],
            each=length(treated$tpoint)))
        fit.reduced <- within(fit.reduced, Gene <- factor(Gene, 
            levels = symbol[slice]))
        with(fit.reduced, levels(Gene))
        
        fit.ctrl <- transform(fit.ctrl, time=ctrl$tpoint)
        fit.ctrl <- melt(fit.ctrl, id.vars='time')
        fit.ctrl <- transform(fit.ctrl, Gene=rep(symbol[slice],
            each=length(ctrl$tpoint)))
        fit.ctrl <- within(fit.ctrl, Gene <- factor(Gene, levels = symbol[slice]))
        with(fit.ctrl, levels(Gene))
        
        fit.treated <- transform(fit.treated, time=treated$tpoint)
        fit.treated <- melt(fit.treated, id.vars='time')
        fit.treated <- transform(fit.treated, Gene=rep(symbol[slice],
            each=length(treated$tpoint)))
        fit.treated <- within(fit.treated, Gene <- factor(Gene, 
            levels = symbol[slice]))
        with(fit.treated, levels(Gene))
        
        print(ggplot(df, aes(time,value)) + geom_point(aes(colour=cond)) + 
            geom_line(aes(colour=cond)) + 
            geom_line(data=fit.reduced,aes(time,value),alpha=0.2,
                linetype='dotdash') +
            geom_line(data=fit.ctrl,aes(time,value),alpha=0.2,linetype='dashed') +
            geom_line(data=fit.treated,aes(time,value),alpha=0.2,
                linetype='dashed') +
            #geom_text(data=df.param, aes(x=3,y=max(df$value),label=label),size=3,
            #    parse=T) +
            facet_wrap(~Gene, ncol=3) + 
            scale_colour_brewer('',palette='Set1') +
            #ylim(floor(min(df$value)),ceiling(max(df$value))) +
            xlab('time (h)') + ylab(expression(paste(log[2]~'fold expression'))) +
            opts(legend.position=c(0.9,0.9),legend.key=theme_blank())
        )
    }
    dev.off()
    df
}

#' Plot expression time series for the given genes. 
#'
#' This function plots the expression profiles of the given genes.
#'
#' @param ctrl list of time series data under the control condition, contains
#' the field dat, cond, tpoint
#' @param treatd list of time series data under the treatment condition, contains
#' the field dat, cond, tpoint
#' @param filename filename
#' @export
plot.full.reduced.example <- function(ctrl, treated, title.str) {
    symbol <- ctrl$dat$Gene
    df.ctrl <- melt(ctrl$dat)
    df.ctrl <- transform(df.ctrl, cond=ctrl$cond)
    df.ctrl <- transform(df.ctrl, time=ctrl$tpoint)
    df.treated <- melt(treated$dat)
    df.treated <- transform(df.treated, cond=treated$cond)
    df.treated <- transform(df.treated, time=treated$tpoint)
    df <- rbind(df.ctrl, df.treated)
    df <- within(df, cond <- factor(cond, levels = c(treated$cond,ctrl$cond)))
    with(df, levels(cond))
    
    # full/reduced model fit
    idx.ctrl <- which(names(ctrl$dat)=='Gene')
    idx.treated <- which(names(treated$dat)=='Gene')
    core <- cbind(ctrl$dat[,-idx.ctrl],treated$dat[,-idx.treated])
    red <- apply(core,1,function(x) fit.red(x,ctrl$tpoint,treated$tpoint))
    fit.reduced <- as.data.frame(red[-(1:length(ctrl$tpoint)),])
    names(fit.reduced) <- symbol
    full <- apply(core,1,function(x) fit.full(x,ctrl$tpoint,treated$tpoint))
    fit.ctrl <- as.data.frame(full[1:length(ctrl$tpoint),])
    names(fit.ctrl) <- symbol
    fit.treated <- as.data.frame(full[-(1:length(ctrl$tpoint)),])
    names(fit.treated) <- symbol

    # get the p-value
    df.param <- apply(core,1,function(x) lrtest(x,ctrl$tpoint,treated$tpoint))
    df.param <- melt(df.param)
    df.param <- transform(df.param, Gene=symbol, 
        label=sprintf("p==%.1e",df.param$value))

    # rearrange the data frame to plot the fitting curve
    fit.reduced <- transform(fit.reduced, time=treated$tpoint)
    fit.reduced <- melt(fit.reduced, id.vars='time')
    fit.reduced <- transform(fit.reduced, Gene=rep(symbol,
        each=length(treated$tpoint)))
    
    fit.ctrl <- transform(fit.ctrl, time=ctrl$tpoint)
    fit.ctrl <- melt(fit.ctrl, id.vars='time')
    fit.ctrl <- transform(fit.ctrl, Gene=rep(symbol,
        each=length(ctrl$tpoint)))
    
    fit.treated <- transform(fit.treated, time=treated$tpoint)
    fit.treated <- melt(fit.treated, id.vars='time')
    fit.treated <- transform(fit.treated, Gene=rep(symbol,
        each=length(treated$tpoint)))
    
    print(ggplot(df, aes(time,value)) + geom_point(aes(colour=cond)) + 
        geom_line(aes(colour=cond)) + 
        geom_line(data=fit.reduced,aes(time,value),alpha=0.2,
            linetype='dotdash') +
        geom_line(data=fit.ctrl,aes(time,value),alpha=0.2,linetype='dashed') +
        geom_line(data=fit.treated,aes(time,value),alpha=0.2,
            linetype='dashed') +
        geom_text(data=df.param, aes(x=3,y=1,label=label),size=3,
            parse=T) +
        scale_colour_brewer('',palette='Set1') +
        #ylim(floor(min(df$value)),ceiling(max(df$value))) +
        labs(x='time (h)',y=expression(paste(log[2]~'fold expression')),
            title=paste(title.str,'\n')) +
        opts(legend.position=c(0.9,0.9),legend.key=theme_blank())
    )
    df
}

#' Plot expression time series for the given genes. 
#'
#' This function plots the expression profiles of the given genes.
#'
#' @param ctrl list of time series data under the control condition, contains
#' the field dat, cond, tpoint
#' @param treatd list of time series data under the treatment condition, contains
#' the field dat, cond, tpoint
#' @param filename filename
#' @export
fit.impulse <- function(tpoint, ts) {
    require(minpack.lm)
    xfit <- seq(0,max(tpoint),length=100)
    ## model based on a list of parameters
    #logistic <- function(par, xx) par$a + par$b / (1 + exp(-par$c * (xx - par$d)))
    logistic <- function(par, xx) 1/par$h1 * (par$h0+(par$h1-par$h0)/(1+exp(par$beta*(par$t1-xx)))) * (par$h2+(par$h1-par$h2)/(1+exp(par$beta*(xx-par$t2))))
    ## residual function
    residFun <- function(p, observed, xx) observed - logistic(p,xx)
    ## starting values for parameters
    #parStart <- list(a=1, b=1, c=1, d=1)
    #if (abs(max(ts[1:8])) > abs(min(ts[1:8])))
    #    parStart <- list(h1=2,h0=0,beta=-1,t1=1,h2=1,t2=2)
    #else
        parStart <- list(h1=-0.5,h0=0,beta=1,t1=1,h2=-1,t2=2)
    ## perform fit
    nls.out <- nls.lm(par=parStart, fn = residFun, observed = ts, xx = tpoint)
                                             
    ## plot model evaluated at final parameter estimates
    logistic(as.list(coef(nls.out)), xfit)
}

#' Plot expression time series for the given genes. 
#'
#' This function plots the expression profiles of the given genes.
#'
#' @param ctrl list of time series data under the control condition, contains
#' the field dat, cond, tpoint
#' @param treatd list of time series data under the treatment condition, contains
#' the field dat, cond, tpoint
#' @param filename filename
#' @export
fit.impulse.param <- function(tpoint, ts) {
    ## model based on a list of parameters
    #logistic <- function(par, xx) par$a + par$b / (1 + exp(-par$c * (xx - par$d)))
    logistic <- function(par, xx) 1/par$h1 * (par$h0+(par$h1-par$h0)/(1+exp(par$beta1*(par$t1-xx)))) * (par$h2+(par$h1-par$h2)/(1+exp(par$beta2*(xx-par$t2))))
    ## residual function
    residFun <- function(p, observed, xx) observed - logistic(p,xx)
    ## starting values for parameters
    #parStart <- list(a=1, b=1, c=1, d=1)
    parStart <- list(h1=1,h0=0,beta1=1,beta2=1,t1=1,h2=1,t2=1)
    ## perform fit
    nls.out <- nls.lm(par=parStart, fn = residFun, observed = ts, xx = time.points)
                                             
    ## plot model evaluated at final parameter estimates
    param <- coef(nls.out)
}

#' Plot expression time series for the given genes. 
#'
#' This function plots the expression profiles of the given genes.
#'
#' @param ctrl list of time series data under the control condition, contains
#' the field dat, cond, tpoint
#' @param treatd list of time series data under the treatment condition, contains
#' the field dat, cond, tpoint
#' @param filename filename
#' @export
plot.impulse <- function(treated, filename) {
    pdf(filename,8,11)
    for (i in 1:ceiling(nrow(treated$dat)/12)) {
    #for (i in 1:10) {
        symbol <- treated$dat$Gene
        slice <- seq((i-1)*12+1,i*12)
        slice <- slice[!(slice > nrow(treated$dat))]
        #df.ctrl <- melt(ctrl$dat[slice,])
        #df.ctrl <- transform(df.ctrl, cond=ctrl$cond)
        #df.ctrl <- transform(df.ctrl, time=rep(ctrl$tpoint,each=length(slice)))
        df.treated <- melt(treated$dat[slice,])
        df.treated <- transform(df.treated, cond=treated$cond)
        df.treated <- transform(df.treated, time=rep(treated$tpoint,
            each=length(slice)))
        #df <- rbind(df.ctrl, df.treated)
        df <- df.treated
        df <- within(df, Gene <- factor(Gene, levels = symbol[slice]))
        with(df, levels(Gene))
        #df <- within(df, cond <- factor(cond, levels = c(treated$cond,ctrl$cond)))
        #with(df, levels(cond))
        
        # full/reduced model fit
        #idx.ctrl <- which(names(ctrl$dat)=='Gene')
        idx.treated <- which(names(treated$dat)=='Gene')
        #core <- cbind(ctrl$dat[slice,-idx.ctrl],treated$dat[slice,-idx.treated])
        #red <- apply(core,1,function(x) fit.red(x,ctrl$tpoint,treated$tpoint))
        #fit.reduced <- as.data.frame(red[-(1:length(ctrl$tpoint)),])
        #names(fit.reduced) <- symbol[slice]
        #full <- apply(core,1,function(x) fit.full(x,ctrl$tpoint,treated$tpoint))
        #fit.ctrl <- as.data.frame(full[1:length(ctrl$tpoint),])
        #names(fit.ctrl) <- symbol[slice]
        #fit.treated <- as.data.frame(full[-(1:length(ctrl$tpoint)),])
        #names(fit.treated) <- symbol[slice]
        fit.treated <- apply(treated$dat[slice,-idx.treated],1,function(x) fit.impulse(treated$tpoint,x))
        names(fit.treated) <- symbol[slice]

        # get the p-value
        #df.param <- apply(core,1,function(x) lrtest(x,ctrl$tpoint,treated$tpoint))
        #df.param <- melt(df.param)
        #df.param <- transform(df.param, Gene=symbol[slice], 
        #    label=sprintf("p==%.1e",df.param$value))

        # rearrange the data frame to plot the fitting curve
        #fit.reduced <- transform(fit.reduced, time=treated$tpoint)
        #fit.reduced <- melt(fit.reduced, id.vars='time')
        #fit.reduced <- transform(fit.reduced, Gene=rep(symbol[slice],
        #    each=length(treated$tpoint)))
        #fit.reduced <- within(fit.reduced, Gene <- factor(Gene, 
        #    levels = symbol[slice]))
        #with(fit.reduced, levels(Gene))
        #
        #fit.ctrl <- transform(fit.ctrl, time=ctrl$tpoint)
        #fit.ctrl <- melt(fit.ctrl, id.vars='time')
        #fit.ctrl <- transform(fit.ctrl, Gene=rep(symbol[slice],
        #    each=length(ctrl$tpoint)))
        #fit.ctrl <- within(fit.ctrl, Gene <- factor(Gene, levels = symbol[slice]))
        #with(fit.ctrl, levels(Gene))
        
        nfit <- 100
        xfit <- seq(0,max(treated$tpoint),length=nfit)
        fit.treated <- transform(fit.treated, time=xfit)
        fit.treated <- melt(fit.treated, id.vars='time')
        fit.treated <- transform(fit.treated, Gene=rep(symbol[slice], each=nfit))
        fit.treated <- within(fit.treated, Gene <- factor(Gene, 
            levels = symbol[slice]))
        with(fit.treated, levels(Gene))
        
        print(ggplot(df, aes(time,value)) + geom_point(aes(colour=cond)) + 
            geom_line(aes(colour=cond)) + 
            #geom_line(data=fit.reduced,aes(time,value),alpha=0.2,
            #    linetype='dotdash') +
            #geom_line(data=fit.ctrl,aes(time,value),alpha=0.2,linetype='dashed') +
            geom_line(data=fit.treated,aes(time,value),alpha=0.2,
                linetype='dashed') +
            #geom_text(data=df.param, aes(x=3,y=max(df$value),label=label),size=3,
            #    parse=T) +
            facet_wrap(~Gene, ncol=3) + 
            scale_colour_brewer('',palette='Set1') +
            #ylim(floor(min(df$value)),ceiling(max(df$value))) +
            xlab('time (h)') + ylab(expression(paste(log[2]~'fold expression'))) +
            opts(legend.position=c(0.9,0.9),legend.key=theme_blank())
        )
    }
    dev.off()
    df
}

#' Plot expression time series for the given genes. 
#'
#' This function plots the expression profiles of the given genes.
#'
#' @param ctrl list of time series data under the control condition, contains
#' the field dat, cond, tpoint
#' @param treatd list of time series data under the treatment condition, contains
#' the field dat, cond, tpoint
#' @param filename filename
#' @export
plot.logistic <- function(treated, filename) {
    require(reshape)
    pdf(filename,8,11)
    for (i in 1:ceiling(nrow(treated$dat)/12)) {
    #for (i in 1:10) {
        symbol <- treated$dat$Gene
        slice <- seq((i-1)*12+1,i*12)
        slice <- slice[!(slice > nrow(treated$dat))]
        df.treated <- melt(treated$dat[slice,])
        df.treated <- transform(df.treated, cond=treated$cond)
        df.treated <- transform(df.treated, time=rep(treated$tpoint,
            each=length(slice)))
        df <- df.treated
        df <- within(df, Gene <- factor(Gene, levels = symbol[slice]))
        with(df, levels(Gene))
        
        # logistic fit
        idx.treated <- which(names(treated$dat)=='Gene')
        fit.treated <- apply(treated$dat[slice,-idx.treated],1,function(x) fit.logistic(treated$tpoint,x))
        colnames(fit.treated) <- symbol[slice]

        # get upregulation time
        nfit <- 100
        xfit <- seq(0,max(treated$tpoint),length=nfit)
        df.param <- apply(fit.treated,2,function(x)xfit[which.max(abs(diff(x)/diff(xfit)))])
        df.param <- melt(df.param)
        df.param <- transform(df.param, Gene=symbol[slice], 
            label=sprintf("t==%.1e",df.param$value))

        fit.treated <- transform(fit.treated, time=xfit)
        fit.treated <- melt(fit.treated, id.vars='time')
        fit.treated <- transform(fit.treated, Gene=rep(symbol[slice], each=nfit))
        fit.treated <- within(fit.treated, Gene <- factor(Gene, 
            levels = symbol[slice]))
        with(fit.treated, levels(Gene))
        
        print(ggplot(df, aes(time,value)) + geom_point(aes(colour=cond)) + 
            geom_line(aes(colour=cond)) + 
            #geom_line(data=fit.reduced,aes(time,value),alpha=0.2,
            #    linetype='dotdash') +
            #geom_line(data=fit.ctrl,aes(time,value),alpha=0.2,linetype='dashed') +
            geom_line(data=fit.treated,aes(time,value),alpha=0.2,
                linetype='dashed') +
            geom_text(data=df.param, aes(x=3,y=5,label=label),size=3,parse=T) +
            facet_wrap(~Gene, ncol=3) + 
            scale_colour_brewer('',palette='Set1') +
            #ylim(floor(min(df$value)),ceiling(max(df$value))) +
            xlab('time (h)') + ylab(expression(paste(log[2]~'fold expression'))) +
            opts(legend.position=c(0.9,0.9),legend.key=theme_blank())
        )
    }
    dev.off()
    df
}

#' Plot expression time series for the given genes. 
#'
#' This function plots the expression profiles of the given genes.
#'
#' @param ctrl list of time series data under the control condition, contains
#' the field dat, cond, tpoint
#' @param treatd list of time series data under the treatment condition, contains
#' the field dat, cond, tpoint
#' @param filename filename
#' @export
plot.logistic.example <- function(treated, title.str) {
    require(reshape)
    symbol <- treated$dat$Gene
    df.treated <- melt(treated$dat)
    df.treated <- transform(df.treated, cond=treated$cond)
    df.treated <- transform(df.treated, time=treated$tpoint)
    df <- df.treated
    
    # logistic fit
    idx.treated <- which(names(treated$dat)=='Gene')
    fit.treated <- apply(treated$dat[,-idx.treated],1,function(x) fit.logistic(treated$tpoint,x))
    colnames(fit.treated) <- symbol

    # get upregulation time
    nfit <- 100
    xfit <- seq(0,max(treated$tpoint),length=nfit)
    df.param <- apply(fit.treated,2,function(x)xfit[which.max(abs(diff(x)/diff(xfit)))])
    df.param <- melt(df.param)
    df.param <- transform(df.param, Gene=symbol, 
        label=sprintf("t==%.1e",df.param$value))

    fit.treated <- transform(fit.treated, time=xfit)
    fit.treated <- melt(fit.treated, id.vars='time')
    fit.treated <- transform(fit.treated, Gene=rep(symbol, each=nfit))
    
    print(ggplot(df, aes(time,value)) + geom_point(aes(colour=cond)) + 
        geom_line(aes(colour=cond)) + 
        geom_line(data=fit.treated,aes(time,value),linetype='dashed') +
        geom_text(data=df.param, aes(x=1,y=3,label=label),size=5,parse=T) +
        scale_colour_brewer('',palette='Set1') +
        #ylim(floor(min(df$value)),ceiling(max(df$value))) +
        labs(x='time (h)',y=expression(paste(log[2]~'fold expression')),
            title=paste(title.str,'\n')) +
        opts(legend.position=c(0.9,0.9),legend.key=theme_blank())
    )
    df
}

