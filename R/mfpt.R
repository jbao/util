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
get.mfpt.sp <- function(mfpt.mat, idmap, net, source.entrez, target.entrez) {
    #require(doMC)
    #registerDoMC()
    mfpt.sp <- foreach(i=1:length(source.entrez), .combine='rbind') %:%
        foreach(j=1:length(target.entrez), .combine='rbind') %dopar% {
            ii <- idmap$vertex[which(source.entrez[i] == idmap$entrez)] + 1
            jj <- idmap$vertex[which(target.entrez[j] == idmap$entrez)] + 1
            source.idx <- match(source.entrez[i], V(net)$gene)
            target.idx <- match(target.entrez[j], V(net)$gene)
            sp <- shortest.paths(net, source.idx, target.idx, weights=NA) # BFS
            data.frame(mfpt=mfpt.mat[ii,jj], sp=as.numeric(sp),
                #source=get(source.entrez[i], org.Hs.egSYMBOL),
                source=idmap$official[which(source.entrez[j]==idmap$entrez)],
                target=idmap$official[which(target.entrez[j]==idmap$entrez)])
    }
}
