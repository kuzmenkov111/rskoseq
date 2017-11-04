#' Quality assesment with ShortRead package and return ggplot object.
#' @description This functions returns ggplot objects and create report of qa. The qa report is created at above of the fastq directory.
#' @usage gg_qa(fqdir, suffix, prefix, facet_col)
#' @param fqdir A vector of file path of fastq files, or dir path, containing fastq files
#' @param suffix A pattern of fastq file suffix. The default is ".fastq.gz"
#' @param prefix A vector of samples name. The default values are names of fastq files containing in fqdir, which substitute 'suffix' character.
#' @param facet_col facet of sequence content and quality score per sample plot.
#' @return  list of data frame for quality assesment with 'qa'
#'     and list ofggplot objects
#' @examples
#' ## arguments
#' # p <- system.file("extdata/E-MTAB-1147", package = "ShortRead")
#' # sffx <- ".fastq.gz"
#' # fqs <- list.files(p, sffx)
#' # prfx <- sapply(strsplit(fqs, "/"), function(x)sub(sffx, "", tail(x, n=1)))
#'
#' ## execution
#' # res_qa <- gg_qa(fqdir=p, suffix=sffx, prefix=prfx, facet_col=2)
#' # do.call(gridExtra::grid.arrange, c(res_qa, list(ncol=2)))
#'
#' @importFrom ShortRead qa report
#' @importFrom dplyr %>% arrange group_by mutate ungroup filter
#' @importFrom ggplot2 ggplot aes geom_bar theme_bw theme element_text
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics plot
#' @importFrom rsko mstrings
#' @importFrom tidyr gather
#' @importFrom stats median
#' @importFrom utils head
#' @export
gg_qa <- function(fqdir,
                  suffix=".fastq.gz",
                  prefix =sub(suffix, "", list.files(fqdir, suffix)),
                  facet_col){
  # fqdir <- "~/pub/dat/sampledata/rnaseq/project1/fastq";
  # suffix <- ".fastq.gz";
  # facet_col=2

  # file exists or not(full path) ----
  if (!all(file.exists(fqdir))){
    stop("I cannot find these all files.")
  }

  # quarity assesment
  ## outputdir ----
  qadir <-
    paste0(
      paste(lapply(strsplit(fqdir, "/"), function(x)head(x, length(x)-1))[[1]],
                  collapse = "/"),
      "/qa")
  if (!file.exists(qadir)){
    dir.create(qadir)
  } else if (file.exists(qadir) & !identical(list.files(qadir), character(0))){
    stop(paste("There is some file at ", qadir))
  }

  ## execute qa and reporting ----
  qadat <- ShortRead::qa(dirPath = fqdir)
  ShortRead::report(qadat, dest = paste0(qadir, "/report"))

  ## read number ----
  read <- NULL; Base <- NULL; Count <- NULL; Cycle <- NULL; Qscore <- NULL;
  `Score Sequence Content(%)` <- NULL;  lane <- NULL; value <- NULL;
  Score <- NULL; `Sequence Content(%)` <- NULL;

  ## lane name ----
  if(length(rownames(qadat[["readCounts"]])) != length(prefix) | identical(prefix, character(0))){
    stop("The number of prefix and samples must to be the same.")
  } else {
    smp <- prefix
  }


  ## read data frame ----
  readat <- data.frame(sample = smp,
                       read = qadat[["readCounts"]]$read)

  ## ggplot object
  read_gg <-
    ggplot2::ggplot(data = readat, ggplot2::aes(x = sample, y = read), fill = sample) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::theme_bw()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))


  # sequence content
  scdat <- qadat[["perCycle"]]$baseCall
  scdat$lane <- as.character(scdat$lane)

  ## lane name ----
  scdat$lane <- rsko::mstrings(patterns = unique(as.character(scdat$lane)),
                               target = as.character(scdat$lane),
                               action = "mmsub",
                               replacements = prefix )

  ## sequence content ----
  scrate <- scdat %>%
    dplyr::arrange(Base) %>%
    dplyr::mutate(Base=factor(Base, levels=c("A","T","G","C","N"))) %>%
    dplyr::mutate(lane=sub(suffix, "", lane)) %>%
    dplyr::group_by(Cycle, lane) %>%
    dplyr::mutate(`Sequence Content(%)` = (Count/sum(Count))*100) %>%
    dplyr::ungroup()

  ## sequence content per sample ----
  sc_gg <-
    ggplot2::ggplot(scrate, ggplot2::aes(x=Cycle, y=`Sequence Content(%)`, group=Base, colour=Base))+
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~lane, ncol=facet_col)

  ## sequence content per nuc ----
  sc_gg2 <-
    ggplot2::ggplot(scrate, ggplot2::aes(x=Cycle, y = `Sequence Content(%)`, colour=Base, group=lane))+
    ggplot2::geom_line(alpha=0.5) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~Base, ncol=facet_col)

  ## fraquency of 'N' per sample ----
  nrgg_smp <- scrate %>%
    dplyr::filter(Base=="N") %>%
    ggplot2::ggplot(ggplot2::aes(x=Cycle, y=`Sequence Content(%)`)) +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::labs(y="N Content(%)") +
    ggplot2::facet_wrap(~lane, ncol=2)

  # quality socre per sample
  qsdat <- qadat[["perCycle"]]$quality
  qsdat$lane <- as.character(qsdat$lane)

  ## replace the lane name by prefix ----
  qsdat$lane <- rsko::mstrings(patterns = unique(as.character(qsdat$lane)),
                               target = as.character(qsdat$lane),
                               action = "mmsub",
                               replacements = prefix )

  qsds <- qsdat %>%
    dplyr::group_by(Cycle, lane) %>%
    dplyr::summarise(median=median(rep(Score, Count)), mean=mean(rep(Score, Count))) %>%
    dplyr::ungroup() %>%
    tidyr::gather(key = Qscore, value = value, median, mean)

  qs_gg <- ggplot2::ggplot(qsds, ggplot2::aes(x=Cycle, y=value, group=Qscore, colour=Qscore))+
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~lane, ncol=facet_col)

  # quality score all sample -----
  nsmp <- length(unique(qsds$lane))
  gttle <- paste0("All ", nsmp, " samples")
  qs_gg2 <- qsds %>%
    dplyr::filter(Qscore=="mean") %>%
    ggplot2::ggplot(ggplot2::aes(x = Cycle, y = value))+
    ggplot2::geom_point(size = 0.7, alpha = 0.5) +
    ggplot2::theme_light() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(title = gttle, y = "Mean Quality Score")


  # return ----
  ## return list of data frames
  df_list <- list(Read = readat, QScore = qsdat, SContent=scdat)
  gg_list <- list (read = read_gg, seq_cnt_smp = sc_gg, seq_cnt_nuc = sc_gg2,
                   n_gg = nrgg_smp, qsc_smp = qs_gg, qsc_all = qs_gg2)

  ## save pdf ----
  grDevices::pdf(paste0(qadir, "/", "qa_plot.pdf"))
  invisible(lapply(gg_list, graphics::plot))
  grDevices::dev.off()

  return(gg_list)


}
