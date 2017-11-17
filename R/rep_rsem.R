#' Continuous execution of sequence alignment using RSEM
#' @description Description
#' @usage rep_rsem (path_aln, idx_name, suffix_fq)
#' @param path_aln alignment directory path, created by rskoseq::project_rnsq.
#' @param idx_name rsem index file path
#' @param suffix_fq suffix of fastq files. The default value is ".fastq.gz"
#' @examples ##
#' ### rep rsem
#' # idx <- "~/db/index/rsem_idx/cge25207.add"
#' # alnd <- "~/pub/sampledata/rnaseq/project1/rsem1/"
#' # prj <- "~/pub/sampledata/rnaseq/project1"
#' # rep_rsem(path_prj = prj, alnd = alnd, idx_name=idx)
#' @importFrom parallel detectCores
#' @importFrom readr read_delim
#' @export
#'
# ## rsem prepare ref ----
# listf <- "~/db/index/rsem_idx/CHO-K1-v1.refseq_2014.transcripts.list"
# fas <- "~/db/index/rsem_idx/CHO-K1-v1.refseq_2014.transcripts.fa"
# idxname <- "CHO-K1-v1.refseq_2014.transcripts"
# rep_rsem(prj = getwd(), alnd="rsem171104")
rep_rsem <- function(path_aln, idx_name, suffix_fq=".fastq.gz"){
  # default setting arguments
  ## project name and alignment directory name----
  alnd <- basename(path_aln)
  path_prj <- dirname(path_aln)
  prjn <- basename(path_prj)


  ## move to result path ----
  setwd(path_aln)

  ## fastq files path and prefix and suffix fqs ----
  path_fq <- list.files(paste0(path_prj, "/fastq"), suffix_fq, full.names = T) # default
  prefix_fq <- sub(suffix_fq, "", sapply(strsplit(path_fq, "/"), function(x)tail(x, 1))) # default
  path_maplog <- paste0(path_aln, "/",prefix_fq, "_map_log.txt") # default

  # index file exists or not ----
  if (!file.exists(paste0(idx_name, ".1.bt2"))){
    stop("Can't find index file.")
  }

  ## rsem-calculate-expression path ----
  rsem <- suppressWarnings(system("which rsem-calculate-expression", intern = T))
  if (!length(rsem)){
    stop("Ther is no rsem-calculate-expression, or PATH environmental variable")
  }

  # ## check version
  # system(paste0(system("which bowtie2", intern=T), " --version"), intern = T)[[1]]
  # system(paste0(system("which samtools", intern=T), " --version"), intern = T)[[1]]

  ## command log file ----
  path_comlog <- list.files(path_aln, "log.txt", full.names = T)
  if (identical(path_comlog, character(0))){
    path_comlog <- paste0(path_aln, "/", prjn, "_", alnd, "_log.txt")
    file.create(path_comlog)
  }

  ## detect cores ----
  cores <- parallel::detectCores()

  # rsem execution
  ## open connection of command-log file ----
  con <- file(path_comlog, "a")
  writeLines("# rsem", con)

  for(i in seq_along(path_fq)){
    com <- paste(rsem, "--bowtie2 --sort-bam-by-coordinate -p",
                 cores,
                 path_fq[i], idx_name, prefix_fq[i], ">>", path_maplog[i], "2>&1",
                 sep=" ")
    system(com, wait = T)
    cat(paste0(com, " \n"))
    writeLines(com, con)
  }
  close(con)

  # file manipulation
  ## result files removed to ./resdir ----
  ## fpkm
  gen <- list.files(".", ".genes.results")
  iso <- list.files(".", ".isoforms.results")
  ifpkms <- lapply(seq_along(iso), function(i) {
    readr::read_delim(iso[i], delim="\t", col_names=T) %>%
      select("transcript_id", "FPKM") %>%
      setNames(., c("transcript_id", sub("_R1.isoforms.results", "", iso[i]))) %>%
      arrange(transcript_id)
  })
  ifpkm <- data.frame(transcript_id=ifpkms[[1]]$transcript_id,
                      do.call(cbind, lapply(ifpkms, function(x)x[,2])),
                      check.names=F, stringsAsFactors=F)
  ## ecount
  iecnts <- lapply(seq_along(iso), function(i) {
    readr::read_delim(iso[i], delim="\t", col_names=T) %>%
      select("transcript_id", "expected_count") %>%
      setNames(., c("transcript_id", sub("_R1.isoforms.results", "", iso[i])))
  })
  iecnt <- data.frame(transcript_id=iecnts[[1]]$transcript_id,
                      do.call(cbind, lapply(iecnts, function(x)x[,2])),
                      check.names=F, stringsAsFactors=F)
  write.table(iecnt, "./iecnt.txt", sep="\t", quote = F, col.names = T, row.names = F)
  write.table(ifpkm, "./ifpkm.txt", sep="\t", quote = F, col.names = T, row.names = F)
  ## gfpkm
  gfpkms <- lapply(seq_along(iso), function(i) {
    readr::read_delim(gen[i], delim="\t", col_names=T) %>%
      select("gene_id", "FPKM") %>%
      setNames(., c("gene_id", sub("_R1.genes.results", "", gen[i])))
  })
  gecnts <- lapply(seq_along(iso), function(i) {
    readr::read_delim(gen[i], delim="\t", col_names=T) %>%
      select("gene_id", "expected_count") %>%
      setNames(., c("gene_id", sub("_R1.genes.results", "", gen[i])))
  })


  fc1 <- sapply(c(gen, iso), function(x)file.copy(x, paste0(path_aln, "/resdir")))
  if (all(fc1)){
    fr1 <- file.remove(c(iso, gen))
    if (all(fr1)){
      print("all result files removed")
    } else {
      print("Could not all result files removed")
    }
  }

  ## all stat directories removed to ./stats ----
  fc2 <- sapply(list.files(path_aln, ".stat"), function(x)file.copy(x, paste0(path_aln, "/stats"), recursive = T))
  if (all(fc2)){
    unlink(list.files(path_aln, ".stat"), recursive=T)
    print("all stat directory removed")
  }

  ## sorted bam files removed to ./sortedbam ----
  fc3 <- sapply(list.files(path_aln, ".sorted.bam"), function(x)file.copy(x, paste0(path_aln, "/sortbam")))
  if (all(fc3)){
    fr3 <- file.remove(list.files(path_aln, "\\.bam"))
    if (all(fr3)){
      print("all sorted bam files removed")
    } else {
      print("Could not all sorted bam files removed")
    }
  }

}
