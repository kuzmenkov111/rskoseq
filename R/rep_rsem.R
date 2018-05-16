#' Continuous execution of sequence alignment using RSEM
#' @description Description
#' @usage rep_rsem (alndir, fqdir, paired, idx_name, suffix_fq)
#' @param alndir alignment directory path, created by rskoseq::project_rnsq.
#' @param fqdir character: the fully path of fastq files. The default is 'paste0(dirname(alndir), "/fastq")'.
#'     If this directory containing still analyzed fastq and additional fastq files,
#' @param paired logical: paired or single read. If paired end, the names of fastq file must be "_R1" and "_R2".
#' @param idx_name rsem index file path
#' @param suffix_fq suffix of fastq files. The default value is ".fastq.gz"
#' @examples ##
#' ### rep rsem
#' # single ----
#' # idx <- "~/db/index/rsem_idx/cge25207.add"
#' # alnd <- "~/pub/sampledata/rnaseq/project1/test.rsm"
#' # rskoseq::rep_rsem(alndir=alnd, idx_name=idx)
#' # paired ----
#' # alnd <- "~/pub/sampledata/rnaseq/project1/test2.rsm"
#' # fqd <- "~/pub/sampledata/rnaseq/project1/paired_fastq"
#' # rskoseq::rep_rsem(alndir = alnd, fqdir = fqd, paired =T)
#' @importFrom dplyr select arrange
#' @importFrom parallel detectCores
#' @importFrom readr read_delim
#' @export
#'
rep_rsem <- function(alndir,
                     fqdir = paste0(dirname(alndir), "/fastq"),
                     paired = FALSE,
                     idx_name,
                     suffix_fq=".fastq.gz"){

  # argument check: collect PATH of fastq files in fastq directory ########################
  ## fastq files or directory exist or not ----
  path_fq <- list.files(fqdir, suffix_fq, full.names = TRUE)
  if (identical(path_fq, character(0))){
    stop(paste("There is not fastq files in", fqdir, ", or the suffix of fastq is different from", suffix_fq, "."))
  }

  ## collect fastq files pqth ----
  if (paired ==T){
    r1fqs <- grep("R1", path_fq, value = T)
    r2fqs <- grep("R2", path_fq, value = T)

  } else {
    r1fqs <- path_fq
  }

  ## prefix of fastq files ----
  if (all(grepl("R1", r1fqs))){
    prefix <- sub("_R1", "",
                  sub(suffix_fq, "",
                      sapply(strsplit(r1fqs, "\\/"), function(x)tail(x, 1))))
  } else {
    prefix <- sub(suffix_fq, "", sapply(strsplit(r1fqs, "/"), function(x)tail(x, 1)))
  }

  ## if alignment directory is not exists  ----
  if (!file.exists(alndir)){dir.create(path = alndir, recursive = T)}
  wd <- getwd()
  setwd(alndir)

  ## log files output under the alignment directory ----
  path_maplog <- paste0(alndir, "/", prefix, "_map_log.txt") # default

  ## index file exists or not ----
  if (!file.exists(paste0(idx_name, ".1.bt2"))){
    stop("Can't find index file.")
  }

  ## rsem-calculate-expression path ----
  rsem <- suppressWarnings(system("which rsem-calculate-expression", intern = T))
  if (!length(rsem)){
    stop("Ther is no rsem-calculate-expression, or PATH environmental variable")
  }

  ## check version ----
  # system(paste0(system("which bowtie2", intern=T), " --version"), intern = T)[[1]]
  # system(paste0(system("which samtools", intern=T), " --version"), intern = T)[[1]]

  ## project name and alignment directory name----
  alnd <- basename(alndir)
  prjn <- basename(dirname(alndir))

  ## command log file ----
  datestrings <- gsub(":", ".", gsub(" ", "_", date()))
  path_comlog <- list.files(alndir, "log.txt", full.names = T)
  if (identical(path_comlog, character(0))){
    path_comlog <- paste0(alndir, "/", prjn, "_", alnd, "_", datestrings, "_log.txt")
    file.create(path_comlog)
  }

  ## detect cores ----
  cores <- parallel::detectCores()

  # rsem execution ############################
  ## open connection of command-log file ----
  con <- file(path_comlog, "a")
  writeLines("# rsem", con)

  ## execute command ----
  if (paired == F){
    for(i in seq_along(r1fqs)){
      com <- paste(rsem, "--bowtie2 --sort-bam-by-coordinate -p",
                   cores,
                   r1fqs[i], idx_name, prefix[i], ">>", path_maplog[i], "2>&1",
                   sep=" ")
      system(com, wait = T)
      cat(paste0(com, " \n"))
      writeLines(com, con)
    }
  } else {
    for(i in seq_along(r1fqs)){
      com <- paste(rsem, "--bowtie2 --sort-bam-by-coordinate -p",
                   cores,
                   "--paired-end",
                   r1fqs[i], r2fqs[i],
                   idx_name,
                   prefix[i], ">>", path_maplog[i], "2>&1")
      system(com, wait = T)
      cat(paste0(com, " \n"))
      writeLines(com, con)
    }
  }

  ## close connection of command-log file ----
  close(con)

  # file manipulation ############################
  ## merge results files  ----
  transcript_id <- NULL; gene_id <- NULL;
  gen <- list.files(".", ".genes.results")
  iso <- list.files(".", ".isoforms.results")
  ifpkms <- lapply(seq_along(iso), function(i) {
    readr::read_delim(iso[i], delim="\t", col_names=T) %>%
      dplyr::select("transcript_id", "gene_id", "FPKM") %>%
      setNames(., c("transcript_id", "gene_id", sub(".isoforms.results", "", iso[i]))) %>%
      dplyr::arrange(transcript_id)
  })
  ifpkm <- data.frame(transcript_id=ifpkms[[1]]$transcript_id,
                      gene_id = ifpkms[[1]]$gene_id,
                      do.call(cbind, lapply(ifpkms, function(x)x[,3])),
                      check.names=F, stringsAsFactors=F)
  ## ecount ----
  iecnts <- lapply(seq_along(iso), function(i) {
    readr::read_delim(iso[i], delim="\t", col_names=T) %>%
      dplyr::select("transcript_id", "gene_id","expected_count") %>%
      setNames(., c("transcript_id", "gene_id", sub(".isoforms.results", "", iso[i]))) %>%
      dplyr::arrange(transcript_id)
  })
  iecnt <- data.frame(transcript_id=iecnts[[1]]$transcript_id,
                      gene_id = iecnts[[1]]$gene_id,
                      do.call(cbind, lapply(iecnts, function(x)x[,3])),
                      check.names=F, stringsAsFactors=F)

  write.table(iecnt, paste0("./",prjn, ".", alnd, ".iecnt.txt"),
              sep="\t", quote = F, col.names = T, row.names = F)
  write.table(ifpkm, paste0("./",prjn, ".", alnd, ".ifpkm.txt"),
              sep="\t", quote = F, col.names = T, row.names = F)

  ## gfpkm ----
  `transcript_id(s)` <- NULL
  gfpkms <- lapply(seq_along(iso), function(i) {
    readr::read_delim(gen[i], delim="\t", col_names=T) %>%
      dplyr::select("gene_id", `transcript_id(s)`, "FPKM") %>%
      setNames(., c("gene_id", "transcript_ids", sub(".genes.results", "", gen[i]))) %>%
      dplyr::arrange(gene_id)
  })
  gecnts <- lapply(seq_along(iso), function(i) {
    readr::read_delim(gen[i], delim="\t", col_names=T) %>%
      dplyr::select("gene_id", `transcript_id(s)`, "expected_count") %>%
      setNames(., c("gene_id", "transcript_ids", sub(".genes.results", "", gen[i]))) %>%
      dplyr::arrange(gene_id)
  })

  gfpkm <- data.frame(gene_id = gfpkms[[1]]$gene_id,
                      transcript_ids = gfpkms[[1]]$transcript_ids,
                      do.call(cbind, lapply(gfpkms, function(x)x[,3])),
                      check.names = F, stringsAsFactors = F)

  gecnt <- data.frame(gene_id = gecnts[[1]]$gene_id,
                      transcript_ids = gfpkms[[1]]$transcript_ids,
                      do.call(cbind, lapply(gecnts, function(x)x[,3])),
                      check.names = F, stringsAsFactors = F)

  write.table(gecnt, paste0("./",prjn, ".", alnd, ".gecnt.txt"),
              sep="\t", quote = F, col.names = T, row.names = F)
  write.table(gfpkm, paste0("./",prjn, ".", alnd, ".gfpkm.txt"),
              sep="\t", quote = F, col.names = T, row.names = F)


  # result files remove to ./resdir ----
  fc1 <- sapply(c(gen, iso), function(x)file.rename(x, paste0(alndir, "/resdir/", x)))

  # all stat directories removed to ./stats ----
  fc2 <- sapply(list.files(alndir, ".stat"), function(x)file.rename(x, paste0(alndir, "/stats/", x)))

  # sorted bam files removed to ./sortedbam ----
  fc3 <- sapply(list.files(alndir, ".sorted.bam"), function(x)file.rename(x, paste0(alndir, "/sortbam/", x)))

  if (all(fc3)){
    fr3 <- file.remove(list.files(alndir, "\\.bam"))
    if (all(fr3)){
      print("all sorted bam files removed")
    } else {
      print("Could not all sorted bam files removed")
    }
  }

  # mapping rate ----
  mread <- lapply(seq_along(path_maplog), function(i){
    lines <- readLines(path_maplog[i]) %>% gsub("^[ ]+|[ ]+$", "",  .)
    if (paired == F){
      as.numeric(c(
        sapply(strsplit(grep("aligned 0 times", lines, value = T), " "), "[", 1),
        sapply(strsplit(grep("aligned exactly 1 time", lines, value = T), " "), "[", 1),
        sapply(strsplit(grep("aligned >1 times", lines, value = T), " "), "[", 1)
      ))

    } else if (paired == T){
      as.numeric(c(
        sapply(strsplit(grep("aligned concordantly 0 times", lines, value = T), " "), "[", 1),
        sapply(strsplit(grep("aligned concordantly exactly 1 time", lines, value = T), " "), "[", 1),
        sapply(strsplit(grep("aligned concordantly >1 times", lines, value = T), " "), "[", 1)
      ))

    }

  })
  mread <- setNames(mread, prefix)


  m <- NULL; rate <- NULL; read <- NULL; smp <- NULL
  mdat <- data.frame(m =factor(c("unmapped", "mapped", "multiple"),
                               levels = c("unmapped", "multiple","mapped")),
                               mread, check.names = F) %>%
    tidyr::gather(., "smp", "read", -1) %>%
    group_by(smp) %>%
    mutate(rate = read/sum(read)*100)

  write.table(mdat, "./mapping_rate.txt", sep="\t", quote = F, col.names = T, row.names = F)

  ggmrate <- ggplot2::ggplot(mdat, ggplot2::aes(smp, rate, colour=m)) +
    ggplot2::geom_bar(ggplot2::aes(fill=m, colour = NULL), stat = "identity") +
    ggplot2::scale_fill_manual(values = c("#7c3daf", "#ff9ee7", "#94a80d")) +
    ggplot2::labs(x = NULL, fill = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  ggplot2::ggsave(filename = "mrate.pdf", plot = ggmrate)


  # restore working directory ----
  setwd(wd)

}
