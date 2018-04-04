#' Quality assesment with ShortRead package and return ggplot object.
#' @description This functions returns ggplot objects and create report of qa. The qa report is created at above of the fastq directory.
#' @usage gg_qa(fqdir, suffix, prefix, facet_col, outdir)
#' @param fqdir A vector of file path of fastq files, or dir path, containing fastq files
#' @param suffix A pattern of fastq file suffix. The default is ".fastq.gz"
#' @param prefix A vector of samples name. The default values are names of fastq files containing in fqdir, which substitute 'suffix' character.
#' @param facet_col facet of sequence content and quality score per sample plot.
#' @param outdir output directory. The default values is same at fqdir.
#' @return  The list of ggplot objects returns and write.table for quality assesment with 'qa'.
#' @examples
#' ## arguments
#' # p <- system.file("extdata/E-MTAB-1147", package = "ShortRead")
#' # sffx <- ".fastq.gz"
#' # prfx <- sub(sffx,"",list.files(p, sffx))
#'
#' ## execution
#' # res_qa <- rskoseq::gg_qa(fqdir=p, suffix=sffx, prefix=prfx, facet_col=2)
#' # do.call(gridExtra::grid.arrange, c(res_qa, list(ncol=2)))
#'
#' @importFrom ShortRead qa report
#' @importFrom dplyr %>% arrange group_by mutate ungroup filter
#' @importFrom ggplot2 ggplot aes geom_bar theme_bw theme element_text
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics plot
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather
#' @importFrom stats median
#' @importFrom utils head
#' @export
gg_qa <- function(fqdir,
                  suffix = ".fastq.gz",
                  prefix = sub(suffix, "", list.files(fqdir, suffix)),
                  facet_col,
                  outdir = paste0(dirname(fqdir),"/qa")){
  # fqdir <- "~/pub/sampledata/rnaseq/project1/fastq";
  # suffix <- ".fastq.gz";
  # prefix =sub(suffix, "", list.files(fqdir, suffix))
  # facet_col=2

  # file exists or not(full path) ----
  if (!all(file.exists(fqdir))){
    stop("I cannot find these all files.")
  }

  # quarity assesment
  ## outputdir ----
  if (!file.exists(outdir)){
    dir.create(outdir)
  } else if (file.exists(outdir) & !identical(list.files(outdir), character(0))){
    stop(paste("There is some file at ", outdir))
  }

  ## execute qa and reporting ----
  qadat <- ShortRead::qa(dirPath = fqdir)
  ShortRead::report(qadat, dest = paste0(outdir, "/report"))


  # initialize of tbl object columns ----
  read <- NULL; Base <- NULL; Count <- NULL; Cycle <- NULL; Qscore <- NULL;
  `Score Sequence Content(%)` <- NULL;  lane <- NULL; value <- NULL;
  Score <- NULL; `Sequence Content(%)` <- NULL;

  # lane name ----
  if(length(rownames(qadat[["readCounts"]])) != length(prefix) | identical(prefix, character(0))){
    stop("The number of prefix and samples must to be the same.")
  } else {
    smp <- prefix
  }

  # data frame: read count ----
  readat <- qadat[["readCounts"]] %>%
    tibble::rownames_to_column("sample") %>%
    dplyr::mutate(sample = sub(suffix,"", sample))

  ## ggplot: read count ----
  read_gg <-
    ggplot2::ggplot(data = readat, ggplot2::aes(x = sample, y = read), fill = sample) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::labs(x = "") +
    ggplot2::geom_text(ggplot2::aes(label=read), color="white", angle=90, size=3.5) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))



  # sequence content ----
  ## tbl object: baseCall ----
  scrate <- qadat[["perCycle"]]$baseCall %>%
    dplyr::mutate(lane=sub(suffix, "", lane)) %>%
    dplyr::arrange(Base) %>%
    dplyr::mutate(Base=factor(Base, levels=c("A","T","G","C","N"))) %>%
    dplyr::mutate(lane=sub(suffix, "", lane)) %>%
    dplyr::group_by(Cycle, lane) %>%
    dplyr::mutate(`Sequence Content(%)` = (Count/sum(Count))*100) %>%
    dplyr::ungroup()

  ## ggplot: sequence content per sample ----
  sc_gg <-
    ggplot2::ggplot(scrate, ggplot2::aes(x=Cycle, y=`Sequence Content(%)`, group=Base, colour=Base))+
    ggplot2::geom_line() +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::facet_wrap(~lane, ncol=facet_col)

  ## ggplot: sequence content per nuc ----
  sc_gg2 <-
    ggplot2::ggplot(scrate, ggplot2::aes(x=Cycle, y = `Sequence Content(%)`, colour=Base, group=lane))+
    ggplot2::geom_line(alpha=0.5) +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::facet_wrap(~Base, ncol=facet_col)

  ## ggplot: fraquency of 'N' per sample ----
  nrgg_smp <- scrate %>%
    dplyr::filter(Base=="N") %>%
    ggplot2::ggplot(ggplot2::aes(x=Cycle, y=`Sequence Content(%)`)) +
    ggplot2::geom_point() +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::labs(y="N Content(%)") +
    ggplot2::facet_wrap(~lane, ncol=facet_col)



  # Quality Scores ----
  ## tbl object: quality socre per sample ----
  qsdat <- qadat[["perCycle"]]$quality %>%
    mutate(lane = sub(suffix, "", lane)) %>%
    dplyr::group_by(Cycle, lane) %>%
    dplyr::summarise(median=median(rep(Score, Count)), mean=mean(rep(Score, Count))) %>%
    dplyr::ungroup() %>%
    tidyr::gather(key = Qscore, value = value, median, mean)

  ## ggplot: quality socre per sample ----
  qs_gg <- ggplot2::ggplot(qsdat, ggplot2::aes(x=Cycle, y=value, group=Qscore, colour=Qscore))+
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::facet_wrap(~lane, ncol=facet_col)

  ## ggplot: quality score all samples -----
  nsmp <- length(unique(qsdat$lane))
  gttle <- paste0("All ", nsmp, " samples")
  qs_gg2 <- qsdat %>%
    dplyr::filter(Qscore=="mean") %>%
    ggplot2::ggplot(ggplot2::aes(x = Cycle, y = value))+
    ggplot2::geom_point(size = 0.7, alpha = 0.5) +
    ggplot2::theme_light(base_size = 15) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(title = gttle, y = "Mean Quality Score")


  # frequentSequences ----
  ## initialize of tbl object columns ----
  count=NULL; nReads=NULL; nOccurrences=NULL; `nReads(%)`=NULL

  ## tbl object: top10 frequentSequences ----
  freqSeq_top10 <- qadat[["frequentSequences"]] %>%
    dplyr::mutate(lane = sub(suffix, "", lane)) %>%
    dplyr::group_by(lane) %>%
    dplyr::top_n(count, n=10)

  ## fasta object: top10 frequentSequences ----
  freqSeq_fna <- setNames(Biostrings::DNAStringSet(freqSeq_top10$sequence),
                          paste(freqSeq_top10$lane, freqSeq_top10$count, sep="|"))
  Biostrings::writeXStringSet(x = freqSeq_fna, filepath = paste0(outdir, "/freqSeq.fna"))

  ## tbl object: top10 sequenceDistribution ----
  freqSeq_dat <- qadat[["sequenceDistribution"]] %>%
    dplyr::group_by(lane) %>%
    dplyr::mutate(`nReads(%)`=nReads/sum(nReads)*100) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(lane = sub(suffix, "", lane))

  ## freqSeq ggplot -----
  gg_frqseq <- ggplot2::ggplot(freqSeq_dat, ggplot2::aes(x=nOccurrences, y=`nReads(%)`))+
    ggplot2::scale_color_manual(values = "red") +
    ggplot2::geom_bar(stat="identity", fill="red") +
    ggplot2::facet_wrap(~lane, ncol=facet_col, scale="free") +
    ggplot2::theme_minimal(base_size = 15)

  # return list of ggplot obj and write.table
  # list of data frame and ggplot objects   ----
  df_list <- list(readat, scrate, qsdat, freqSeq_top10)
  gg_list <- list (read = read_gg, seq_cnt_smp = sc_gg, seq_cnt_nuc = sc_gg2,
                   n_gg = nrgg_smp, qsc_smp = qs_gg, qsc_all = qs_gg2,
                   freqSeq = gg_frqseq)


  # write.table ----
  fls <- paste0(outdir, c("/readat.txt", "/scrate.txt", "/qsdat.txt", "/freqSeq.txt"))
  invisible(lapply(seq_along(fls), function(i){
    write.table(df_list[[i]], fls[[i]], sep="\t", quote = F, col.names = T, row.names = F)
  }))

  # save pdf ----
  ## read ----
  grDevices::pdf(paste0(outdir, "/", "read.pdf"), width = 10, height = 5)
  print(read_gg)
  grDevices::dev.off()

  ## sc per sample ----
  if (nrow(readat) > 48){
    grDevices::pdf(paste0(outdir, "/", "sc.pdf"), width = 10, height = 5)
    sc <- ggplus::facet_multiple(plot=sc_gg, facets="lane", ncol = 8, nrow = 6)
    grDevices::dev.off()

  } else {
    grDevices::pdf(paste0(outdir, "/", "sc_s.pdf"), width = 10, height = 5)
    print(sc_gg)
    grDevices::dev.off()
  }

  ## sc per base ----
  grDevices::pdf(paste0(outdir, "/", "sc_b.pdf"), width = 10, height = 5)
  print(sc_gg2)
  grDevices::dev.off()

  ## qs per sample ----
  if (nrow(readat) > 48){
    grDevices::pdf(paste0(outdir, "/", "qs_s.pdf"), width = 10, height = 5)
    sc <- ggplus::facet_multiple(plot=qs_gg, facets="lane", ncol = 8, nrow = 6)
    grDevices::dev.off()

  } else {
    grDevices::pdf(paste0(outdir, "/", "qs_s.pdf"), width = 10, height = 5)
    print(qs_gg)
    grDevices::dev.off()
  }

  ## qs all ----
  grDevices::pdf(paste0(outdir, "/", "qs_al.pdf"), width = 10, height = 5)
  print(qs_gg2)
  grDevices::dev.off()

  ## freqSeq ----
  if (nrow(readat) > 48){
    grDevices::pdf(paste0(outdir, "/", "freqSeq.pdf"), width = 10, height = 5)
    sc <- ggplus::facet_multiple(plot=gg_frqseq, facets="lane", ncol = 8, nrow = 6)
    grDevices::dev.off()

  } else {
    grDevices::pdf(paste0(outdir, "/", "freqSeq.pdf"), width = 10, height = 5)
    print(gg_frqseq)
    grDevices::dev.off()
  }

  # return ----
  return(gg_list)

}
