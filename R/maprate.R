#' Parse hisat2 log file
#' @description Parse hisat2 log file and extract mapping rate of all samples.
#' @usage maprate(fp, lab)
#' @param fp file path of bowtie log
#' @param lab file name as sample name
#' @return data frame of mapping rate and ggplot object
#' @examples
#' h2.paired.log <- system.file("extdata", "hisat2_log.txt", package = "rskoseq")
#' label <- c("A","B","C","D")
#' res <- maprate(fp = h2.paired.log, lab = label)
#' @export

maprate <- function(fp, lab){

  # argument check: hista2 log file
  if(!file.exists(fp)){
    stop("This file does not exist.")
  }

  # argument check: labels length
  lines <- readLines(fp)
  nsample <- length(grep("overall alignment rate", lines))
  if(!length(lab)==nsample){
    stop("number of samples and label's length are different")
  }

  # read log file ----
  lines <- gsub("^[ ]+|[ ]+$", "",  lines)
  rn <- sub(" reads; of these:", "", lines[grep("reads; of these", lines)])

  prn_lines <- grep(" were paired; of these:", lines, value = T)
  prn <- sapply(strsplit(prn_lines, " "), "[", 1)
  prn_rate <- sub(".*\\((.*)\\%\\).*", "\\1", prn_lines)


  if (!identical(prn_lines, character(0))){
    # unmap ----
    unm_lines <- grep(") aligned concordantly 0 times", lines, value = T)
    unm <- sapply(strsplit(unm_lines, " "), "[", 1)
    unm_rate <- sub("\\((.*)\\%\\)", "\\1", sapply(strsplit(unm_lines, " "), "[", 2))

    # unique ----
    uqm_lines <- grep(" aligned concordantly exactly 1 time", lines, value = T)
    uqm <- sapply(strsplit(uqm_lines, " "), "[", 1)
    uqm_rate <- sub("\\((.*)\\%\\)", "\\1", sapply(strsplit(uqm_lines, " "), "[", 2))

    # multiple ----
    mm_lines <- grep(" aligned concordantly >1 times", lines, value = T)
    mm <- sapply(strsplit(mm_lines, " "), "[", 1)
    mm_rate <- sub("\\((.*)\\%\\)", "\\1", sapply(strsplit(mm_lines, " "), "[", 2))

    # data.frame of paired-reads mapping rate ----
    mrate <- data.frame(id=lab,
                        nreads=rn, npreads = prn, nunmreads = unm,
                        nuqreads = uqm, nmpreads = mm,
                        rpreads = prn_rate,
                        unmapped = as.numeric(unm_rate),
                        multiple = as.numeric(mm_rate),
                        unique = as.numeric(uqm_rate),
                        stringsAsFactors = F
    )
  } else {
    # unmap ----
    unm_lines <- grep(") aligned 0 times", lines, value = T)
    unm <- sapply(strsplit(unm_lines, " "), "[", 1)
    unm_rate <- sub("\\((.*)\\%\\)", "\\1", sapply(strsplit(unm_lines, " "), "[", 2))

    # unique ----
    uqm_lines <- grep(") aligned exactly 1 time", lines, value = T)
    uqm <- sapply(strsplit(uqm_lines, " "), "[", 1)
    uqm_rate <- sub("\\((.*)\\%\\)", "\\1", sapply(strsplit(uqm_lines, " "), "[", 2))

    # multiple ----
    mm_lines <- grep(") aligned >1 times", lines, value = T)
    mm <- sapply(strsplit(mm_lines, " "), "[", 1)
    mm_rate <- sub("\\((.*)\\%\\)", "\\1", sapply(strsplit(mm_lines, " "), "[", 2))

    # data.frame of single-reads mapping rate ----
    mrate <- data.frame(id=lab,
                        nreads=rn, nunmreads = unm,
                        nuqreads = uqm, nmpreads = mm,
                        unmapped = as.numeric(unm_rate),
                        multiple = as.numeric(mm_rate),
                        unique = as.numeric(uqm_rate),
                        stringsAsFactors = F
    )
  }


  id <- NULL; unmapped <- NULL; multiple <- NULL;  key <- NULL; value <- NULL; pos <- NULL
  ggmrate <- mrate %>%
    dplyr::select(id, unmapped, multiple, unique) %>%
    tidyr::gather(key="key", value="value", -1) %>%
    dplyr::mutate(key = factor(key, levels=c("unmapped", "multiple", "unique"))) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(pos = rev(cumsum(rev(value)) - rev(value) + rev(value)/2)) %>%
    ggplot2::ggplot(ggplot2::aes (x = id, y = value, fill = key)) +
    ggplot2::theme_minimal() +
    ggplot2::geom_bar(stat = "identity") +
    #ggplot2::geom_text(ggplot2::aes(label=value), col = "white",size = 3, hjust = 0.5, vjust = 3, position = "stack") +
    ggplot2::geom_text(ggplot2::aes(label=value, y=pos), col = "white") +
    ggplot2::theme(axis.text.x =ggplot2::element_text(angle=90, hjust=1)) +
    ggplot2::labs(x="", y="mapping rate(%)", fill="" )

  return(list(mrate, ggmrate))
}

