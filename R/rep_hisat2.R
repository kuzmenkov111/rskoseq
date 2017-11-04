#' Replicate execution of read alignment using hisat2
#' @description The NGS read alignment using hisat2 for multiple samples.
#'     The input and output directory must be created by 'rskoseq::project_rnsq'.
#' @usage rep_hisat2(prjd, idx, paired, alnd, fqd, suffix_fq, ...)
#' @param prjd character: project directory. create by 'rskoseq::dir_template'
#' @param idx character: the fully path of hisat2 index name
#' @param paired logical: paired or sinle read
#' @param alnd character: the name of alignment directory. The default is "alignment1"
#' @param fqd character: the fully path of fastq files. The default is 'paste0(prjd, "/fastq")'
#' @param suffix_fq character: suffix of fastq files. The default is ".fastq.gz"
#' @param ... additional hisat2 options. E.g.  "--no-spliced-alignment"
#' @importFrom utils tail
#' @examples # # arguments ----
#' # prj <- "~/pub/sampledata/rnaseq/project1"
#' # idx <- "~/db/index/hisat2_idx"
#'
#' # # excution ----
#' # rep_hisat2(prjd = prj, idx = idx, paired=TRUE, ...="--no-spliced-alignment")
#'
#' # # if excute different alignment setting/
#' # rep_hisat2(prjd = prj, idx = idx, paired=FALSE, alnd = "alignment2")
#'
#' # # now directory component
#' # project1/
#' # ├── alignment1
#' # │   ├── project1_alignment1_log.txt
#' # │   ├── res_hisat2
#' # │   │   ├── failalign
#' # │   │   │   ├── test1.failalign.fq.1.gz
#' # │   │   │   ├── test1.failalign.fq.2.gz
#' # │   │   │   ├── test2.failalign.fq.1.gz
#' # │   │   │   └── test2.failalign.fq.2.gz
#' # │   │   ├── h2_log.txt
#' # │   │   ├── test1.hisat.met.txt
#' # │   │   ├── test1.sort.bam
#' # │   │   ├── test2.hisat.met.txt
#' # │   │   └── test2.sort.bam
#' # │   └── res_stringtie
#' # │       ├── ballgown
#' # │       ├── gff
#' # │       ├── mgff
#' # │       └── tab
#' # ├── alignment2
#' # │   ├── res_hisat2
#' # │   │   └── failalign
#' # │   └── res_stringtie
#' # │       ├── ballgown
#' # │       ├── gff
#' # │       ├── mgff
#' # │       └── tab
#' # ├── fastq
#' # │   ├── test1_R1.fastq.gz
#' # │   ├── test1_R2.fastq.gz
#' # │   ├── test2_R1.fastq.gz
#' # │   └── test2_R2.fastq.gz
#' # └── qa
#' #
#' # 18 directories, 14 files
#' @export
## memo delete
## prjd="~/pub/dat/sampledata/rnaseq/project1";
## idx <- "~/db/index/hisat2_idx/SD0218_11_a_contig01"; paired=T;alnd = "alignment2"
## fqd = paste0(prjd, "/fastq"); suffix_fq=".fastq.gz"
rep_hisat2 <- function(prjd, idx, paired,
                       alnd = "alignment1",
                       fqd = paste0(prjd, "/fastq"),
                       suffix_fq=".fastq.gz", ...){

  # hisat2 and samtools program PATH ----
  hs2c <- suppressWarnings(system("which hisat2", intern = T))
  samc <- suppressWarnings(system("which samtools", intern = T))
  if (hs2c == 1){
    stop("There is not hisat2 program, or the PATH does not found.")
  }
  if (samc == 1){
    stop("There is not samtools program, or the PATH does not found.")
  }

  # argument check: project directory ----
  if(!file.exists(prjd)){
    stop(paste0("\'", prjd, "\'", " does not found."))
  }
  # argument check: hisat2 index exists or not  ----
  idx1 <- paste0(idx, ".1.ht2")
  if (!file.exists(idx1)){
    stop(paste0("Thres is not hisat2 index named as ", idx, " ."))
  }

  # argument check: alignment directory  exists or not ----
  alnp <- paste0(prjd, "/", alnd)
  if (!file.exists(alnp)){
    stop(paste0(" There is not alignment directory '", alnd, "', which is at directly under the '", prjd, "'. \n"))
  }

  # argument check: alignment data still exists or not ----
  bamf <- paste0(prjd, "/", alnd, "/res_hisat2")
  if ( length(list.files(bamf, "bam")) ){
    stop(paste0(" There is already bam files in the '", bamf, "' ."))
  }

  # argument check: path of fastq files and collect fastq files fqd and sample name
  # argument check: fastq directory  ----
  if (file.exists(fqd)){
    fqpath <- list.files(fqd, suffix_fq, full.names = TRUE)
    if (identical(fqpath, character(0))) {
      stop(paste0("There is no fastq files, or the suffix of these fastq files is different from ",  suffix_fq))
    }
  } else {
    stop(paste0("There is not '",fqd, "'  directory."))
  }

  # collect PATH of fastq files in fastq directory  ----
  r1fqs <- grep("R1", list.files(fqd, suffix_fq, full.names = T), value = T)
  r2fqs <- grep("R2", list.files(fqd, suffix_fq, full.names = T), value = T)
  prefix <- sub(paste0("_R1", suffix_fq), "",
                sapply(strsplit(r1fqs, "\\/"), function(x)tail(x, 1)))

  # create command and execution ----
  ## hisat2 log file ----
  logout <- paste0(" 2>> ", prjd, "/", alnd, "/res_hisat2/h2_log.txt")

  ## detect cores ----
  cores <- parallel::detectCores()

  ## command log maked by rskoseq::project_rnsq or newly creaion. ----
  logfile  <-  list.files(paste0(prjd, "/", alnd), "log.txt", full.names = T)
  if (identical(logfile, character(0))){
    prjn <- sapply(strsplit(prjd, "/"), function(x)tail(x,1))
    logfile <- paste0(prjd, "/",alnd, "/", prjn, "_", alnd, "_", "log.txt")
    file.create(logfile)
  }

  ## open connection of command-log file ----
  con <- file(logfile, "a")
  writeLines("# hisat2", con)

  ## bam, h2_log.txt, and met files must be in the same directory. ----
  for(i in seq_along(r1fqs)){
    ## fail align ----
    failaln <- paste0(prjd, "/", alnd, "/res_hisat2/failalign/",
                      prefix[i], ".failalign.fq.gz")
    ## hisat2 meta file ----
    metfile <- paste0(prjd, "/", alnd, "/res_hisat2/",
                      prefix[i], ".hisat.met.txt ")
    ## hisat2 additional options ----
    if (!missing(...)){
      add_op <- paste0(" ", ..., " ")
    }else{
      add_op <- ""
    }

    ## samtools sam -> bam -> sort ----
    bamdir <- paste0(prjd, "/", alnd, "/res_hisat2/", prefix[i], ".sort.bam")
    bampfx <- paste0(prjd, "/", alnd, "/res_hisat2/", prefix[i], "_sort")
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
                    logout, " | ",
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
                    logout, " | ",
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

  # close command_log connection ----
  close(con)
}

