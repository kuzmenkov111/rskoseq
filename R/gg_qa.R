#' Quality assesment with ShortRead package and return ggplot object.
#' @description This functions returns ggplot objects
#' @usage gg_qa(filepath, suffix, prefix, facet_col)
#' @param filepath A vector of file path of fastq files, or dir path, containing fastq files
#' @param suffix A pattern of fastq file suffix. eg. ".fastq.gz"
#' @param prefix A vector of samples name(optional).
#' @param facet_col facet
#' @return  list of data frame for quality assesment with 'qa'
#'     and list ofggplot objects
#' @examples
#' ## arguments
#' # p <- system.file("extdata/E-MTAB-1147", package = "ShortRead")
#' # sffx <- ".fastq.gz"
#' # fls <- list.files(p, sffx, full.names = TRUE)
#' # prfx <- sapply(strsplit(fls, "/"), function(x)sub(sffx, "", tail(x, n=1)))
#' ## execution
#' # res_qa <- gg_qa(filepath=fls, suffix=sffx, prefix=prfx, facet_col=2)
#' # do.call(gridExtra::grid.arrange, c(res_qa, list(ncol=2)))
#'
#' @importFrom ShortRead qa
#' @importFrom dplyr %>% arrange group_by mutate ungroup filter
#' @importFrom ggplot2 ggplot aes geom_bar theme_bw theme element_text
#' @importFrom tidyr gather
#' @importFrom stats median
#' @export
gg_qa <- function(filepath, suffix, prefix, facet_col){
  # file exists or not(full path)
  if (!all(file.exists(filepath))){
    stop("I cannot find these all files.")
  }

  # quarity assesment ----
  qadat <- ShortRead::qa(dirPath = filepath)

  # read number ----
  read <- NULL; Base <- NULL; Count <- NULL; Cycle <- NULL; Qscore <- NULL;
  `Score Sequence Content(%)` <- NULL;  lane <- NULL; value <- NULL;
  Score <- NULL; `Sequence Content(%)` <- NULL;

  ## lane name
  if(!missing(prefix)){
    smp <- sub(suffix, "", rownames(qadat[["readCounts"]]))
  }else{
    smp <- prefix
  }

  ## read data frame
  readat <- data.frame(sample = smp,
                       read = qadat[["readCounts"]]$read)

  ## ggplot object
  read_gg <-
    ggplot2::ggplot(data = readat, ggplot2::aes(x = sample, y = read), fill = sample) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::theme_bw()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))


  # sequence content ----
  scdat <- qadat[["perCycle"]]$baseCall
  scdat$lane <- as.character(scdat$lane)

  ## lane name
  ggrep <- function(pattern, target){
    lapply(pattern, function(x){
      grep(x, target)
    })
  }
  if(!missing(prefix)){
    smp <- sub(suffix, "", as.character(scdat$lane))
    scdat$lane <- sub(suffix,"",as.character(scdat$lane))

  }else{
    pos <- ggrep(prefix, as.character(scdat$lane))
    invisible(lapply(seq_along(pos),
                     function(i){scdat$lane[pos[[i]]] <<- prefix[i]}))
  }

  ## sequence content
  scrate <- scdat %>%
    dplyr::arrange(Base) %>%
    dplyr::mutate(Base=factor(Base, levels=c("A","T","G","C","N"))) %>%
    dplyr::mutate(lane=sub(suffix, "", lane)) %>%
    dplyr::group_by(Cycle, lane) %>%
    dplyr::mutate(`Sequence Content(%)` = (Count/sum(Count))*100) %>%
    dplyr::ungroup()

  ## sequence content per sample
  sc_gg <-
    ggplot2::ggplot(scrate, ggplot2::aes(x=Cycle, y=`Sequence Content(%)`, group=Base, colour=Base))+
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~lane, ncol=facet_col)

  ## sequence content per nuc
  sc_gg2 <-
    ggplot2::ggplot(scrate, ggplot2::aes(x=Cycle, y = `Sequence Content(%)`, colour=Base, group=lane))+
    ggplot2::geom_line(alpha=0.5, size=1) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~Base, ncol=facet_col)

  ## fraquency of 'N' per sample
  nrgg_smp <- scrate %>%
    dplyr::filter(Base=="N") %>%
    ggplot2::ggplot(ggplot2::aes(x=Cycle, y=`Sequence Content(%)`)) +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::labs(y="N Content(%)") +
    ggplot2::facet_wrap(~lane, ncol=2)

  # quality socre per sample -----
  qsdat <- qadat[["perCycle"]]$quality
  qsdat$lane <- as.character(qsdat$lane)

  ## lane name
  if(!missing(prefix)){
    smp <- sub(suffix, "", as.character(qsdat$lane))
    qsdat$lane <- sub(suffix,"",as.character(qsdat$lane))

  }else{
    pos <- ggrep(prefix, as.character(qsdat$lane))
    invisible(lapply(seq_along(pos),
                     function(i){qsdat$lane[pos[[i]]] <<- prefix[i]}))
  }

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

  return(gg_list)

}
