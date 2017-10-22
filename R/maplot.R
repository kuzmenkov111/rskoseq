#' maplot
#' @description maplot
#' @usage maplot(dat, fdr, gp1, gp2)
#' @param dat dataframe: result of DEG analysis using 'tmm_cnt'
#' @param fdr q-value default 0.05
#' @param gp1,gp2 group name
#' @importFrom grDevices adjustcolor
#' @importFrom ggplot2 ggplot aes geom_point labs scale_color_manual theme_bw facet_wrap theme element_rect
#' @export
maplot <-function(dat, fdr = 0.05, gp1, gp2){
  # dat <- deg_list[[1]]; fdr = 0.01; gtitle=""; gp1="ctrl"; gp2 = "c2";
  # dat <- deg_list

  # argument check: dat ----
  vcols <- grDevices::adjustcolor(col = c("magenta", "gray30"), alpha.f = 0.5)
  A <- NULL; M <- NULL; a.value <- NULL; m.value <- NULL

  # maplot function ----
  maplt <- function(dat, fdr, gp1, gp2){
    # graph title
    comp <- paste0("gp1: ", gp1, "   gp2: ", gp2)
    # number of deg or non-deg ----
    nsig <- sum(dat$q.value < fdr); nnsig <- sum(dat$q.value >= fdr)
    lab_sig <- paste0("DEG:", nsig); lab_nsig <- paste0("non-DEG:", nnsig)
    # colour for significant genes
    sig <- factor(ifelse(dat$q.value < fdr, lab_sig, lab_nsig),
                  levels = c(lab_sig, lab_nsig))
    # ggplot
    magg <- ggplot2::ggplot(dat, aes(x = a.value, y = m.value, colour=sig )) +
      ggplot2::geom_point(size=0.3) +
      ggplot2::labs(x="A=(log2(gp2)+log2(gp1))/2", y= "M=log2(gp2)-log2(gp1)",
                    colour="", title = comp) +
      ggplot2::scale_color_manual(values = vcols)+
      ggplot2::theme_bw(base_size = 15) +
      ggplot2::theme(legend.position = c(.85, .85),
                     legend.text = element_text(size = 15),
                     legend.key.size = grid::unit(3.0, "lines"))

    return(magg)
  }

  if (class(dat)=="list"){
    # gg_list <-lapply(dat, function(x)maplt(dat = x, fdr = fdr, gp1 = gp1, gp2 = gp2))
    # do.call(gridExtra::grid.arrange, c(gg_list, list(ncol = 3)))

    # maplot with facet ----
    tmp <- do.call(rbind, dat)
    names(tmp)[names(tmp) =="a.value"] <- "A"; names(tmp)[names(tmp)=="m.value"] <- "M"
    sig <-  factor(ifelse(tmp$q.value < fdr, 1, 0), levels = c(1, 0))

    if (is.null(names(tmp))){
      comp <- rep(seq(length(dat)), c(sapply(dat, nrow)))
    }else{
      comp <- factor(rep(names(dat), c(sapply(dat, nrow))), levels = names(dat))
    }

    tmp <- cbind(tmp, comp, sig)
    magg <- ggplot2::ggplot(tmp, ggplot2::aes(x = A, y = M, colour = sig)) +
      ggplot2::geom_point(size=0.3) +
      ggplot2::labs(x="A=(log2(a)+log2(b))/2", y= "M=log2(b)-log2(a)",
                    colour="", title = "") +
      ggplot2::scale_color_manual(values = vcols) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = "FDR < 0.01") +
      ggplot2::theme(strip.background = element_rect(colour="white", fill="white"),
                     legend.position="none") +
      ggplot2::facet_wrap(~comp, ncol=3)


  } else if (class(dat)=="data.frame"){
    res_magg <- maplt(dat = dat, fdr = fdr, gp1 = gp1, gp2 = gp2)
    res_magg

    # base plot ----
    # plot(dat$a.value, dat$m.value,
    #      xlab = "A=(log2(gp2)+log2(gp1))/2",
    #      ylab = "M=log2(gp2)-log2(gp1)",
    #      pch=16, cex=0.3, main=title
    #      )
    # mtext(paste("gp1:", gp1, "   gp2:", gp2, sep = ""), cex = 0.8)
    # #sigdat <- dat %>% filter(estimatedDEG ==1)
    # sigdat <- dat %>% dplyr::filter(q.value < fdr)
    #
    # points(sigdat$a.value, sigdat$m.value, pch=16, cex=0.3, col="magenta")
    # abline(h = 0, v = 0, lty=3)
    # ndeg <- nrow(sigdat); nnondeg <- nrow(dat)-ndeg
    # legend("topright", bg="transparent",
    #        col=c("black", "magenta"), pch=16,box.lty = 0,
    #        cex= 0.8,
    #        c(paste0("non_DEG: ", nnondeg),
    #          paste("DEG: ", ndeg, "\n",  "(FDR<", fdr, ")", sep=""))
    # )
  }

  return(magg)
}
