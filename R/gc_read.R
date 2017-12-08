#' Quality assesment with ShortRead package and return ggplot object.
#' @description This functions returns ggplot objects and create data.frame of GC content per read
#' @usage gc_read(fqdir, suffix, prefix, facet_col)
#' @param fqdir A vector of file path of fastq files, or dir path, containing fastq files
#' @param suffix A pattern of fastq file suffix. The default is ".fastq.gz"
#' @param prefix A vector of samples name. The default values are names of fastq files containing in fqdir, which substitute 'suffix' character.
#' @param facet_col facet of sequence content and quality score per sample plot.
#' @return  The list of ggplot objects returns and write.table for quality assesment with 'qa'.
#' @examples
#' ## arguments
#' # p <- system.file("extdata/E-MTAB-1147", package = "ShortRead")
#' # sffx <- ".fastq.gz"
#' # fqs <- list.files(p, sffx)
#' # prfx <- sapply(strsplit(fqs, "/"), function(x)sub(sffx, "", tail(x, n=1)))
#'
#' ## execution
#' # res_qa <- gc_read(fqdir=p, suffix=sffx, prefix=prfx, facet_col=2)
#' # do.call(gridExtra::grid.arrange, c(res_qa, list(ncol=2)))
#'
#' @importFrom ShortRead FastqSampler yield sread
#' @importFrom dplyr %>%
#' @importFrom ggplot2 ggplot aes geom_histogram theme_bw labs facet_wrap
#' @importFrom Biostrings letterFrequency
#' @importFrom grDevices dev.off pdf
#' @importFrom tidyr gather
#' @export

gc_read <- function(fqdir,
                    suffix=".fastq.gz",
                    prefix =sub(suffix, "", list.files(fqdir, suffix)),
                    facet_col){
  # file exists or not(full path) ----
  if (!all(file.exists(fqdir))){
    stop("I cannot find these all files.")
  }

  # fastq files ----
  fqs <- list.files(fqdir, suffix, full.names = T)
  if (identical(fqs, character(0))){
    stop(paste0("There is not fastq files in ", fqdir))

  }

  # function: gc sampler ----
  gc_sampler <- function(file_path, n=NULL){
    if (is.null(n)) n <- 1e+6
    sampler_fq <- ShortRead::FastqSampler(file_path, n)
    fqs <- ShortRead::yield(sampler_fq)
    close(sampler_fq)
    # gc content per read
    v_gc <- Biostrings::letterFrequency(ShortRead::sread(fqs), "GC", as.prob=TRUE)*100
    return(v_gc)
  }
  # function: cbind.fillna ----
  cbind.fillna <- function(l, fill){
    mlen <- max(sapply(l, length))
    d <- data.frame(sapply(l,  function(x){ c(x, rep(NA, mlen-length(x)))}), check.names = F)
    return(d)
  }

  # data.frame of GC content per read
  GC=NULL; count=NULL; nOccurrences=NULL; nReads=NULL; `nReads(%)`=NULL

  gc_dat <- lapply(fqs, function(x)gc_sampler(file_path = x)) %>%
    cbind.fillna(.) %>%
    setNames(., prefix) %>%
    tidyr::gather(., key="sample", value="GC")

  # ggplot ----
  gggc <- ggplot2::ggplot(gc_dat, ggplot2::aes(x=GC)) +
    ggplot2::geom_histogram(bins = 20, alpha=0.5) +
    ggplot2::theme_bw(base_size = 15)+
    ggplot2::labs(x="GC(%)")+
    ggplot2::facet_wrap(~sample, ncol=facet_col)

  # return ----
  return(gggc)
}



