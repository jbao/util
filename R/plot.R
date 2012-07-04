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
plot.timeseries <- function(ctrl, treated, filename) {
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
