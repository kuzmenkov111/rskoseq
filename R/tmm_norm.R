#' Normarization of RNA-seq count data using tmm
#' @description This functions returns the normarized count data of multi group
#' @usage tmm_norm(dat, column, gp, method)
#' @param dat A data frame, matrix. The 'dat' contains columns which names genes and rows contains samples.
#' @param column Count data columns without id column. The default is 1:ncol(dat)
#' @param gp replicate group
#' @param method  a pipe line select from 1:"DEGES/TbT", 2:"DEGES/edgeR", 3:"DEGES/DESeq".
#' @return dataframe of normalized count data
#' @examples
#' # sample data of rna-seq
#' data(fpkm2)
#' ## normalized count
#' # group <- fpkm2$reps
#' # nfpkm <- tmm_norm(dat = fpkm2, column = 5:ncol(fpkm2), gp = group, method = 2)
#' @import TCC
#' @importFrom methods new
#' @export
tmm_norm <- function(dat, column = 1:ncol(dat), gp, method =2){
  # dat = fpkm2; column = 5:ncol(fpkm2); gp = rep(1:4, each=3); method =2
  # argument check: nrow dat  and length gp ----
  if(nrow(dat) != length(gp)){
    stop("'dat' contains columns which names genes, and rows contains samples")
  }

  # argument check: dat ----
  if(class(dat) == "data.frame" & !all(sapply(dat[column], mode) == "numeric")){
    stop("'dat[column]' just contains count data.")
  } else if (class(dat) == "matrix" & mode(dat) !="numeric"){
    stop("'dat[column]' just contains count data.")
  }

  # transpose dat for TCC ----
  if(class(dat)=="data.frame"){
    d <- t(dat[column])
  }else if (class(dat)=="matrix"){
    d <- t(data.frame(dat)[column])
  }

  # create TCC object ----
  tcc <- methods::new("TCC", d, gp)

  # calcNormFactors ----
  if (method == 1){ # 1:"DEGES/TbT"
    samplesize <- 10000
    tcc <- TCC::calcNormFactors(tcc, norm.method = "tmm",
                                test.method = "bayseq",
                                iteration = 1,
                                samplesize = samplesize)

  } else if (method == 2){ # 2:"DEGES/edgeR"
    tcc <- TCC::calcNormFactors(tcc, norm.method = "tmm",
                                test.method ="edger",
                                iteration = 3)

  } else if (method == 3){  # 3:"DEGES/DESeq"
    tcc <- TCC::calcNormFactors(tcc, norm.method = "deseq",
                                test.method = "deseq",
                                iteration = 1)
  } else {
    stop('Select from a number of method  1:"DEGES/TbT", 2:"DEGES/edgeR", 3:"DEGES/DESeq" for multi group normarization. ')
  }

  # getNormalizedData ----
  dat[column] <- t(TCC::getNormalizedData(tcc))
  return(dat)

}
