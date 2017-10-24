#' Replicate execution of read mapping using hisat2
#' @description Description
#' @usage rep_hisat2(prj, fqdir, suffix_fq, idx, paired, ...)
#' @param prj character: project directory. create by 'rskoseq::dir_template'
#' @param fqdir character: path of fastq files. default: paste0(prj, "/fastq")
#' @param suffix_fq character: suffix of fastq files. The default is ".fastq.gz"
#' @param idx character: hisat2 index name
#' @param paired logical: paired or sinle read
#' @param ... additional hisat2 options. E.g.  "--no-spliced-alignment"
#' @examples #
#' # # arguments ----
#' # prj <- "~/pub/dat/sampledata/rnaseq/project1"
#' # # fqdir <- paste0(prj, "/fastq")
#' # # suffix_fq <- ".fastq.gz"
#' # idx <- "~/db/index/hisat2_idx/SD0218_11_a_contig01"
#' # paired <- FALSE
#' # # excution ----
#' # rep_hisat2(prj=prj, idx=idx, paired=TRUE, ...="--no-spliced-alignment")
#' @importFrom utils tail
#' @export
rep_hisat2 <- function(prj, fqdir=paste0(prj, "/fastq"), suffix_fq=".fastq.gz", idx, paired, ...){
  # argument check ----
  ## argument check: hisat2 program PATH ----
  if (!any(grep("hisat2", unlist(strsplit(Sys.getenv("PATH"), ":"))))){
    stop("There is not hisat2 program, or the PATH does not found.")
  }
  ## argument check: project directory ----
  if(!file.exists(prj)){
    stop(paste0("\'", prj, "\'", " does not fount."))
  }
  ## argument check: hisat2 index  ----
  idx1 <- paste0(idx, ".1.ht2")
  if (!file.exists(idx1)){
    stop("create hisat2 index.")
  }

  # argument check: path of fastq files and collect fastq files fqdir and sample name----
  if (file.exists(fqdir)){
    fqpath <- list.files(fqdir, suffix_fq, full.names = TRUE)
  }
  if (identical(fqpath, character(0))){
    stop(paste0("There is no fastq files or the suffix of these fastq files is different from ",  suffix_fq))
  }


  # fastq files path ----
  r1fqs <- grep("R1", list.files(fqdir, suffix_fq, full.names = T), value = T)
  r2fqs <- grep("R2", list.files(fqdir, suffix_fq, full.names = T), value = T)
  prefix <- sub(paste0("_R1", suffix_fq), "",
                sapply(strsplit(r1fqs, "\\/"), function(x)tail(x, 1)))
  # hisat2 log file ----
  logout <- paste0(" 2>> ", prj, "/res_hisat2/log.txt")

  # execute hisat2 ----
  ## command log ----
  com_log <-  paste0(prj, "/command_log.txt")
  con <- file(com_log, "a")

  ## bam, log, and met files must be in the same directory. ----
  for(i in seq_along(r1fqs)){
    # Synthesis of command parts
    ## fail align ----
    failaln <- paste0(prj, "/res_hisat2/failalign/", prefix[i], ".failalign.fq.gz")
    ## hisat2 meta file ----
    metfile <- paste0(prj, "/res_hisat2/", prefix[i], ".hisat.met.txt ")
    ## hisat2 additional options ----
    if (!missing(...)){
      add_op <- paste0(" ", ..., " ")
    }else{
      add_op <- ""
    }

    ## sam -> bam -> sort ----
    bamdir <- paste0(prj, "/res_hisat2/", prefix[i], ".sort.bam")
    bampfx <- paste0(prefix[i], "_sort")
    com_smtools <-
      paste0("samtools view -bS -@ 4 - | samtools sort -O bam -o ", bamdir,
             " -T ", bampfx, " -@ 4" )

    # execute command ----
    if (paired == TRUE){ ## paired end
      com <- paste0("hisat2 -p 4 -x ", idx,
                    " --un-conc-gz ", failaln,
                    " --met-file ", metfile,
                    " --dta ",
                    add_op,
                    " -1 ", r1fqs[i],
                    " -2 ", r2fqs[i],
                    logout, " | ",
                    com_smtools)
      cat(paste0(com, " \n"))
      system(com)

      ## write command to command_log.txt ----
      writeLines(com, con)

    } else if (paired == FALSE){ # single
      com <- paste0("hisat2 -p 4 -x ", idx,
                    " --un-conc-gz ", failaln,
                    " --met-file ", metfile,
                    add_op,
                    "-U ", r1fqs[i],
                    logout, " | ",
                    com_smtools)
      cat(paste0(com, " \n"))
      system(com, wait = TRUE)

      ## write command to command_log.txt ----
      writeLines(com, con)
    }

  }

  # close command_log connection ----
  close(con)

  # move files to the respective directory. ----
  # mvbam <- paste0("mv ", paste0(prj, "/res_hisat2/*.bam ", prj, "/res_hisat2/sortbam") )
  # system(mvbam)

}

