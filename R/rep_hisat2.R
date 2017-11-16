#' Replicate execution of read alignment using hisat2
#' @description The NGS read alignment using hisat2 for multiple samples.
#'     The input and output directory must be created by 'rskoseq::project_rnsq'.
#' @usage rep_hisat2(alndir, idx, fqdir, paired, suffix_fq, ...)
#' @param alndir character: the name of alignment directory. results output
#' @param idx character: the fully path of hisat2 index name.
#' @param project logical:default is TRUE
#' @param fqdir character: the fully path of fastq files. The default is 'paste0(dirname(alndir), "/fastq")'.
#'     If this directory containing still analyzed fastq and additional fastq files,
#' @param paired logical: paired or single read. If paired end, the names of fastq file must be "_R1" and "_R2".
#' @param suffix_fq character: suffix of fastq files. The default is ".fastq.gz"
#' @param ... additional hisat2 options. E.g.  "--no-spliced-alignment"
#' @importFrom utils tail
#' @examples # # arguments ----
#' # alndir <- "~/pub/sampledata/rnaseq/project1/alignment1"
#' # idx <- "~/db/index/hisat2_idx"
#' # # excution ----
#' # rep_hisat2(prjd = prj, idx = idx, paired=TRUE, ...="--no-spliced-alignment")
#'
#' # # if all result create under the alignment directory
#' # rep_hisat2(prjd = prj, idx = idx, paired=FALSE, alnd = "alignment2")
#' @export
# memo delete --
# # project
# alndir <- "~/pub/sampledata/rnaseq/project1/alignment1" # result output directory
# idx <- "~/db/index/hisat2_idx/CriGri_1.0"
# fqdir <- paste0(dirname(alndir), "/fastq")
# project <- TRUE # default
# paired <- FALSE # default
# suffix_fq <- ".fastq.gz" # default
#
# # no project
# alndir <- "~/pub/sampledata/rnaseq/project1/aln_171112" # result output directory create
# idx <- "~/db/index/hisat2_idx/CriGri_1.0" # hisat2 index
# fqdir <- "~/pub/sampledata/rnaseq/project1/fastq" # fastq containing dir
# project <- FALSE # default
# paired <- FALSE # default
# suffix_fq <- ".fastq.gz" # default

# if all result output under the 'alndir', project argument must to be the FALSE
# if paired is TRUE, suffix of read must to be _R1.fastq.gz,  _R2.fastq.gz

rep_hisat2 <- function(alndir, idx, project = TRUE,
                       fqdir=paste0(dirname(alndir), "/fastq"),
                       paired=FALSE, suffix_fq=".fastq.gz", ...){
  # system command check: hisat2 and samtools program PATH ----
  hs2c <- suppressWarnings(system("which hisat2", intern = T))
  samc <- suppressWarnings(system("which samtools", intern = T))
  if (hs2c == 1){
    stop("There is not hisat2 program, or the PATH does not found.")
  }
  if (samc == 1){
    stop("There is not samtools program, or the PATH does not found.")
  }

  # argument check: hisat2 index exists or not  ----
  if (!file.exists(paste0(idx, ".1.ht2"))){
    stop(paste0("Thres is not hisat2 index named as ", idx, " ."))
  }
  # argument check: collect PATH of fastq files in fastq directory  ----
  ## fastq files or directory exist or not ----
  path_fq <- list.files(fqdir, suffix_fq, full.names = TRUE)
  if (identical(path_fq, character(0))){
    stop(paste("There is not fastq files in", fqdir, ", or the suffix of fastq is different from", suffix_fq, "."))
  }

  ## collect fastq files pqth ----
  if (paired ==T){
    r1fqs <- grep(paste0("_R1",suffix_fq), path_fq, value = T)
    r2fqs <- grep(paste0("_R2",suffix_fq), path_fq, value = T)
  } else {
    r1fqs <- path_fq
  }

  ## prefix of fastq files ----
  if (all(grepl(paste0("_R1", suffix_fq), r1fqs))){
    prefix <- sub(paste0("_R1", suffix_fq), "",
                  sapply(strsplit(r1fqs, "\\/"), function(x)tail(x, 1)))
  } else {
    prefix <- sub(suffix_fq, "", sapply(strsplit(r1fqs, "/"), function(x)tail(x, 1)))
  }

  # log files output under the alignment directory ----
  ## if alignment directory is not exists  ----
  if (!file.exists(alndir)){dir.create(path = alndir, recursive = T)}
  ## commando log file create ----
  aln <- basename(alndir)
  path_prj <- dirname(alndir)
  prjn <- basename(path_prj)
  path_comlog <- paste0(alndir, "/", prjn, "_", aln,"_", "log.txt")
  file.create(path_comlog)

  ## open connection of command log file ----
  con <- file(path_comlog, "a")
  writeLines(date(), con)

  ## hisat2 log file(created through hisat2 program) ----
  com_h2log <- paste0(" 2>> ", alndir, "/hisat2_log_", gsub(" ", "_", date()), ".txt ")

  ## detect cores ----
  cores <- parallel::detectCores()


  # hisat2 additional options ----
  if (!missing(...)){
    add_op <- paste0(" ", ..., " ")
  }else{
    add_op <- ""
  }


  # project or not ----
  ## if project truth, all result directory are still exists .
  ## Otherwise project is FALSE, all result files and directory are under the alignment directory
  if (project == TRUE){
    # argument check: alignment directory  exists or not ----
    if (!file.exists(alndir)){
      stop(paste0(" There is not alignment directory '", alndir, "'. \n"))
    }

    ## bam, h2_log.txt, and met files must be in the same directory. ----
    for(i in seq_along(r1fqs)){
      ## fail align ----
      failaln <- paste0(alndir, "/res_hisat2/failalign/", prefix[i], ".failalign.fq.gz")
      ## hisat2 meta file ----
      metfile <- paste0(alndir, "/res_hisat2/",
                        prefix[i], ".hisat.met.txt ")
      ## samtools sam -> bam -> sort ----
      bamdir <- paste0(alndir, "/res_hisat2/", prefix[i], ".sort.bam")
      bampfx <- paste0(alndir, "/res_hisat2/", prefix[i], "_sort")
      com_smtools <-
        paste0( samc, " view -bS -@ ", cores, " - | ",
                samc, " sort -O bam -o ", bamdir," -T ", bampfx, " -@ ",cores )

      # execute command ----
      if (paired == TRUE){ ## paired end
        com <- paste0(hs2c, " -p ", cores, " -x ", idx,
                      " --un-conc-gz ", failaln,
                      " --met-file ", metfile,
                      " --dta",
                      add_op,
                      " -1 ", r1fqs[i],
                      " -2 ", r2fqs[i],
                      com_h2log, " | ",
                      com_smtools)
        cat(paste0(com, " \n"))
        return_com <- system(com, intern = T)

        ## if only hisat2 result has statement 1, I don't know whether error occurance. ----
        if (!identical(return_com, character(0))){
          stop("hisat2 or returned error.  ")
        }

        ## write command to command_log.txt ----
        writeLines(com, con)

      } else if (paired == FALSE){ # single
        com <- paste0(hs2c, " -p ", cores, " -x ", idx,
                      " --un-conc-gz ", failaln,
                      " --met-file ", metfile,
                      add_op,
                      "-U ", r1fqs[i],
                      com_h2log, " | ",
                      com_smtools)
        cat(paste0(com, " \n"))
        return_com <- system(com, intern = T)
        ## if only hisat2 result has statement 1, I don't know whether error occurance. ----
        if (!identical(return_com, character(0))){
          stop("hisat2 or returned error.  ")
        }

        ## write command to command_log.txt ----
        writeLines(com, con)
      }
    }
  } else if (project==FALSE){
    ## create failalign dir ----
    dir.create(path = paste0(alndir, "/failalign"), recursive = T)

    ## bam, h2_log.txt, and met files must be in the same directory. ----
    for(i in seq_along(r1fqs)){
      ## fail align ----
      failaln <- paste0(alndir, "/failalign/", prefix[i], ".failalign.fq.gz")
      ## hisat2 meta file ----
      metfile <- paste0(alndir, "/", prefix[i], ".hisat.met.txt ")
      ## samtools sam -> bam -> sort ----
      bamdir <- paste0(alndir, "/", prefix[i], ".sort.bam")
      bampfx <- paste0(alndir, "/", prefix[i], "_sort")
      com_smtools <-
        paste0( samc, " view -bS -@ ", cores, " - | ",
                samc, " sort -O bam -o ", bamdir," -T ", bampfx, " -@ ",cores )

      # execute command ----
      if (paired == TRUE){ ## paired end
        com <- paste0(hs2c, " -p ", cores, " -x ", idx,
                      " --un-conc-gz ", failaln,
                      " --met-file ", metfile,
                      " --dta",
                      add_op,
                      " -1 ", r1fqs[i],
                      " -2 ", r2fqs[i],
                      com_h2log, " | ",
                      com_smtools)
        cat(paste0(com, " \n"))
        return_com <- system(com, intern = T)

        ## if only hisat2 result has statement 1, I don't know whether error occurance. ----
        if (!identical(return_com, character(0))){
          stop("hisat2 or returned error.  ")
        }

        ## write command to command_log.txt ----
        writeLines(com, con)

      } else if (paired == FALSE){ # single
        com <- paste0(hs2c, " -p ", cores, " -x ", idx,
                      " --un-conc-gz ", failaln,
                      " --met-file ", metfile,
                      add_op,
                      "-U ", r1fqs[i],
                      com_h2log, " | ",
                      com_smtools)
        cat(paste0(com, " \n"))
        return_com <- system(com, intern = T)
        ## if only hisat2 result has statement 1, I don't know whether error occurance. ----
        if (!identical(return_com, character(0))){
          stop("hisat2 or returned error.  ")
        }

        ## write command to command_log.txt ----
        writeLines(com, con)
      }
    }

  } else {
    stop("fail")
  }

}



