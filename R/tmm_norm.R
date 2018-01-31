#' Normarization of RNA-seq count data using tmm
#' @description This functions returns the normarized count data of multiple groups
#' @usage tmm_norm(dat, column, gp, method)
#' @param dat A data frame, matrix. The 'dat' contains columns which names genes and rows contains samples.
#' @param column Count data columns without id column. The default is 1:ncol(dat)
#' @param gp replicate group
#' @param method  a pipe line select from 1:"iDEGES/TbT", 2:"iDEGES/edgeR", 3:"iDEGES/DESeq".
#' @return ggobject which containing resulto of 'TCC::getNormalizedData' ,
#' @examples
#' # # sample data of rna-seq
#' # fpkm <- rskodat::fpkm
#' # gp <- sort(as.integer(factor(sapply(strsplit(names(fpkm), "_"),
#' #           function(x) paste(head(x,2), collapse = "")))))
#' # nfpkm <- rskoseq::tmm_norm(dat = fpkm, column = 1:ncol(fpkm), gp = gp, method = 2)
#' # fpkm[1:4,1:4]; nfpkm[1:4,1:4]
#'
#' @import TCC
#' @importFrom methods new
#' @export
tmm_norm <- function(dat, column = 1:ncol(dat), gp, method =2){

  # argument check: dat ----
  if(class(dat) == "data.frame"){
    if (!all(sapply(dat[column], mode) == "numeric")){
      stop("'dat[column]' just contains count data.")
    } else {
      d <- dat[column]
    }
  } else {
    stop("dat must to be data.frame")
  }

  # argument check: nrow(dat) and length(gp) ----
  if(ncol(dat) != length(gp)){
    stop("'dat' contains columns which names samples, and rows contains genes")
  }
  # create TCC object ----
  tcc <- methods::new("TCC", d, gp)

  # calcNormFactors ----
  if (method == 1){ # 1:"iDEGES/TbT"
    samplesize <- 10000
    tcc <- TCC::calcNormFactors(tcc, norm.method = "tmm",
                                test.method = "bayseq",
                                iteration = 1,
                                samplesize = samplesize)

  } else if (method == 2){ # 2:"iDEGES/edgeR"
    tcc <- TCC::calcNormFactors(tcc, norm.method = "tmm",
                                test.method ="edger",
                                iteration = 3)

  } else if (method == 3){  # 3:"iDEGES/DESeq"
    tcc <- TCC::calcNormFactors(tcc, norm.method = "deseq",
                                test.method = "deseq",
                                iteration = 1)
  } else {
    stop('Select from a number of method  1:"iDEGES/TbT", 2:"iDEGES/edgeR", 3:"iDEGES/DESeq" for multi group normarization. ')
  }

  # getNormalizedData ----
  dat[column] <- TCC::getNormalizedData(tcc)
  return(dat)
}
