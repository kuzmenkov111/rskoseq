#' Differential Expression analysis using all combination
#' @description split count data with index and Differential expression analysis of all comination groups.
#' @usage degall(dat, idx, normalize)
#' @param dat dataframe: RNA-seq count table.  row: samples, column: genes
#' @param idx factor: group of data. E.g. factor(c(1,1,2,2)); factor(c("A","A","B","B"))
#' @param normalize integer: If dat is already normalized count data, 'norm' is NULL. The default value is NULL
#'     If dat is not normalized, select a pipe line number from 1:"DEGES/TbT", 2:"DEGES/edgeR", 3:"DEGES/DESeq".
#' @return ggplot object which containing result of 'TCC::estimateDE', without normalized count data.
#' @examples ## normalized fpkm
#' # nfpkm <- rskodat::nfpkm[c(1:3, 10:12, 19:21, 28:30),]
#' # index <- factor(rep(1:4, each=3))
#' # res <- rskoseq::degall(dat=nfpkm[-1:-4], idx=index, normalize=NULL)
#' ## get result of estimateDE
#' # head(res$data)
#' ## ggplot
#' # ggplus::facet_multiple(res, facets="comp", ncol = 3, nrow = 5)
#'
#' #
#' comp1 <- nfpkm[nfpkm$days =="d3" & nfpkm$runs %in% c("S1", "S2"),]
#' comp2 <- nfpkm[nfpkm$days =="d6" & nfpkm$runs %in% c("S1", "S2"),]
#' comp3 <- nfpkm[nfpkm$days =="d12" & nfpkm$runs %in% c("S1", "S2"),]
#' index <- factor(rep(1:2, each=3))
#' res1 <- rskoseq::degall(dat=comp1[-1:-4], idx=index, normalize=NULL)
#' res2 <- rskoseq::degall(dat=comp2[-1:-4], idx=index, normalize=NULL)
#' res3 <- rskoseq::degall(dat=comp3[-1:-4], idx=index, normalize=NULL)
#'
#' @import TCC
#' @importFrom dplyr %>% select mutate
#' @importFrom ggplot2 ggplot aes scale_color_manual theme_minimal labs facet_wrap theme geom_text
#' @importFrom graphics legend
#' @importFrom grDevices adjustcolor
#' @importFrom utils combn
#' @importFrom plyr .
#' @importFrom stats setNames
#' @importFrom tidyr gather
#' @export
degall <- function(dat, idx, normalize = NULL){
  # argument check: dat ----
  if (!is.data.frame(dat) & !is.matrix(dat)){
    stop("dat is a dataframe or matrix object")
  } else {
    if( nrow(dat) != length(idx)){
      stop("dat contains rows as samples, and columns as genes.")
    }
  }

  # split dataframe using index ----
  gplist <- split(seq(nrow(dat)), idx)
  dats <- setNames(lapply(seq_along(gplist), function(i){
    tmp <- data.frame(t(dat[gplist[[i]],])) %>%
    setNames(., paste(levels(idx)[i], 1:length(gplist[[i]]), sep="_"))
    }), levels(idx))


  # Round-robin combination ----
  idx_comb <- combn(levels(idx), 2)
  dat_comb <- lapply(seq(ncol(idx_comb)), function(i) {
    data.frame(dats[[idx_comb[1,i]]], dats[[idx_comb[2,i]]], check.names = F)
  })

  # combination names ----
  comb_name <- paste0("G1(",idx_comb[1,], ") vs G2(", idx_comb[2,], ")")
  names(dat_comb) <- comb_name

  # all combinatins group ----
  gmatch <- function(p, v){
    l <- list()
    for(i in 1:length(p)){ l[[i]] <- which(v %in% p[i]) }
    return(l)
  }

  gps <- vector("list", length=ncol(idx_comb))
  for ( i in 1:ncol(idx_comb)){
    gpi <- gmatch(idx_comb[,i], idx)
    gps[[i]] <- unlist(lapply(seq_along(gpi), function(j){
      rep(j, length(gpi[[j]]))
    }))
  }
  # mapply aruguments check ----
  gp1s <- idx_comb[1,]; gp2s <- idx_comb[2,]

  # normalize and DE function ----
  nde <- function(d, gp, fdr, gp1, gp2, method =2){
    # create tcc object ----
    tcc <- methods::new("TCC", d, gp)
    # calcNormFactors & estimateDE ----
    if (method == 2){ # 2:"DEGES/edgeR"
      tcc <- TCC::calcNormFactors(tcc, norm.method = "tmm",
                                  test.method = "edger",
                                  iteration = 3,
                                  FDR = 0.1,
                                  floorPDEG = 0.05)
      tcc <- TCC::estimateDE(tcc, test.method = "edger", FDR = fdr)

      # get result of estimateDE ----
      res_de <- TCC::getResult(tcc, sort = TRUE)
      return(res_de)
    }
  }
  # DE function ----
  de <- function(d, gp, fdr, gp1, gp2){
    # create tcc object ----
    tcc <- methods::new("TCC", d, gp)
    # estimateDE ----
    tcc <- TCC::estimateDE(tcc, test.method = "edger", FDR = fdr)
    # get result of estimateDE ----
    res_de <- TCC::getResult(tcc, sort = TRUE)
    return(res_de)
  }

  # normalization and DE ----
  if ( !is.null(normalize) ){
    meth.norm <- as.integer(normalize)
    if (!any(c(1,2,3) %in% meth.norm )){ #nde
      stop(' Select a pipeline number froa as follows, 1:"DEGES/TbT", 2:"DEGES/edgeR", 3:"DEGES/DESeq".')
    } else {
      cat (paste(c("iDEGES/TbT", "iDEGES/edgeR", "iDEGES/DESeq")[meth.norm], "was selected.\n" ))
    }
    res <- mapply(FUN = nde, d = dat_comb, gp = gps, gp1 = gp1s, gp2 = gp2s, MoreArgs = list(fdr=0.05, method =meth.norm), SIMPLIFY = FALSE)

  } else if (is.null(normalize)){ # de
    res <- mapply(FUN = de, d = dat_comb, gp = gps, gp1 = gp1s, gp2 = gp2s, MoreArgs = list(fdr=0.05), SIMPLIFY = FALSE)
  }

  gene_id = NULL; m.value=NULL; a.value=NULL; q.value=NULL; estimatedDEG=NULL; fct=NULL

  de_list <- lapply(seq_along(res), function(i){
    res[[i]] %>%
      dplyr::select(c(gene_id, m.value, a.value, q.value, rank, estimatedDEG)) %>%
      dplyr::mutate(fct = factor(ifelse(estimatedDEG==1 & m.value < 0, "down",
                                        ifelse(estimatedDEG==1, "up", "non-DEG")),
                                 levels = c("non-DEG", "up", "down")),
                    comp = factor(rep(comb_name[i], nrow(.))))
  })

  # rbind all deg table ----
  tmp <- do.call(rbind, de_list)

  # dataframe for geom_text using facet ----
  ## text position ----
  x = NULL; y=NULL; key=NULL; value=NULL
  pos_up <- round(max(tmp$a.value) -(max(tmp$a.value)-min(tmp$a.value))/2)-2
  pos_dw <- round(max(tmp$a.value) -(max(tmp$a.value)-min(tmp$a.value))/2)+2
  pos_y <- ceiling(max(tmp$m.value))

  tmp.add <- tapply(tmp$fct, tmp$comp, table)
  add_de <- data.frame(comp=names(tmp.add), up= sapply(tmp.add, "[", 2), down= sapply(tmp.add, "[", 3), row.names = NULL) %>%
    tidyr::gather(key="key", value="value", -1) %>%
    dplyr::mutate(x=ifelse(key =="up", pos_up, pos_dw), y = rep(pos_y, nrow(.)))


  ## ggplot -----
  magg <- ggplot2::ggplot(tmp, ggplot2::aes(x = a.value, y=m.value, colour=fct )) +
    ggplot2::geom_point(size=0.3) +
    ggplot2::scale_color_manual(values =grDevices::adjustcolor(c(8L,2L,4L), alpha.f = 0.5)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "A=(log2(G2)+log2(G1))/2", y = "M=log2(G2)-log2(G1)", colour="") +
    ggplot2::theme(legend.position="top") +
    ggplot2::geom_text(data=add_de, ggplot2::aes(x=x,y=y, label=value, colour=key), size=3, show.legend  = F) +
    ggplot2::facet_wrap(~comp, ncol=3)
  print(magg)
  return(magg)
}

