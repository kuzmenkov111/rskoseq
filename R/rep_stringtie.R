#' Consecutive execution of Transcript assembly and quantification for RNA-Seq using stringtie
#' @description Consecutive processiong of stringtie for multiple samples, then FPKM and cov data table are created.
#'    Execution following from 'rskoseq::project_rnsq', 'rskoseq::rep_hisat2'.
#' @usage rep_stringtie(bamdir,suffix_bam, guide_gff, res_dir, ...)
#' @param bamdir character: the name of alignment directory, sorted-bam files searched from under this directory.
#'     If designate your own sorted-bam containing directory, give the path of the directory
#'     and 'res_dir' must be directory name under the 'bamdir'.
#' @param suffix_bam character: The default is ".sort.bam".
#' @param guide_gff The file path of guide gff.
#' @param res_dir output directory path, the default is 'paste0(dirname(bamdir), "/res_stringtie")'
#' @param ... additional options of stringtie. E.g. "-e"
#' @examples #
#' # bamdir <- "~/pub/sampledata/rnaseq/project1/h2.171117/res_hisat2"
#' # guide "~/db/index/hisat2_idx/CriGri_1.0.gff3"
#' # rep_stringtie(bamdir=bamdir, guide_gff=guide, ... = "-p 8")
#' @importFrom utils tail write.table read.table
#' @export
rep_stringtie <- function(bamdir,
                          suffix_bam=".sort.bam", guide_gff,
                          res_dir=paste0(dirname(bamdir), "/res_stringtie"), ...){
  # argument check: stringtie program PATH ----
  if (!any(grep("stringtie", unlist(strsplit(Sys.getenv("PATH"), ":"))))){
    stop("There is not stringtie program, or the PATH does not found.")
  }

  # argument check: project directory and project name ----
  if(!file.exists(bamdir)){
    stop(paste0("\'", bamdir, "\'", " does not found."))
  }

  # argument check: guide_gff ----
  if (!file.exists(guide_gff)){
    stop(paste0("There is not ", guide_gff))
  }

  # argument check: path of sorted bam files and get all samples name----
  if (file.exists(bamdir)){
    bamfls <- list.files(bamdir, suffix_bam, full.names = T)
    if (identical(bamfls, character(0))){
      stop(paste0("There is not '.srot.bam' files in ", bamdir,
                  ", or the suffix of these bam files is different from '",  suffix_bam, "'."))
    }
  } else {
    stop(paste0("There is not ", bamdir))
  }

  # collect sample names from bam files, and create path of result gff files. ----
  smps <- sub(".sort.bam", "",
              sapply(strsplit(bamfls, "/"), function(x)tail(x, 1)))
  res_gff <- paste0(res_dir, "/gff/", smps, ".gff")

  # stringtie additional options ----
  if (!missing(...)){
    add_op <- paste0(" ", ..., " ")
  }else{
    add_op <- ""
  }

  # stringtie execution
  ## commando log file create ----
  datestrings <- gsub(":", ".", gsub(" ", "_", date()))
  path_comlog <- paste0(dirname(bamdir),"/stringtie_", datestrings,"_log.txt")
  file.create(path_comlog)

  ## open path_comlog of ----
  con <- file(path_comlog, "a")
  writeLines(date(), con)

  ## detect cores ----
  cores <- parallel::detectCores()

  ## execute ----
  for(i in seq_along(bamfls)){
    com <- paste("stringtie -p", cores, bamfls[i], "-G", guide_gff, "-o", res_gff[i], add_op, sep = " ")
    cat(paste0(com, " \n"))
    system(com, wait = TRUE)

    ## write command log ----
    writeLines(com, con)
  }

  # stringtie --merge execution ----
  ## create gff list files (it must be full path) ----
  resgff <- list.files(paste0(res_dir, "/gff"), ".gff", full.names = T)
  resgff_dat <- data.frame(resgff)
  out_f <- paste0(res_dir, "/resgff.txt")
  write.table(resgff_dat, out_f, sep="\t", quote=F, row.names = F, col.names = F)

  ## stringtie --merge ----
  com2 <- paste0("stringtie --merge -G ", guide_gff, " -o ", res_dir, "/merged.gff ", out_f)
  system(com2, wait = T)

  ## over write command log ----
  writeLines(paste0("# stringtie --merge \n", com2), con)

  # stringtie execution using merged gff
  ## define output files ----
  mgff <- paste0(res_dir, "/merged.gff")
  res_mgff <- paste0(res_dir, "/mgff/", smps, "_m.gff")
  res_ballgown <- paste0(res_dir, "/ballgown/", smps)
  res_tab <- paste0(res_dir, "/tab/", smps, ".tab")
  ## replicate execution ----
  writeLines("# stringtie -eb", con)
  for(i in seq_along(bamfls)){
    ## create ballgown directory
    if(!file.exists(res_ballgown[i])){
      dir.create(res_ballgown[i])
    }
    ## command of 'stringtie -eb'
    com3 <-
      paste0("stringtie -p ", cores,
             " ", bamfls[i],
             " -e ",
             " -A ", res_tab[i],
             " -b ", res_ballgown[i],
             " -G ", mgff,
             " -o ", res_mgff[i])
    cat(paste0(com3, " \n"))
    system(com3)

    ## write command log ----
    writeLines(com, con)
  }
  close(con)

  # create fpkm and cov table from t_data.ctab files ----
  ballgown_smp <- paste0(res_dir, "/ballgown/", smps)
  t_dats <- lapply(ballgown_smp, function(x){
    smp <- sapply(strsplit(x, "/"), function(x)tail(x, 1))
    read.table(paste0(x,"/t_data.ctab"), sep="\t", header = T, stringsAsFactors = F)
  })
  cov_list <- lapply(t_dats, function(x)x[c("t_name", "cov")])
  fpkm_list <- lapply(t_dats, function(x)x[c("t_name", "FPKM")])
  invisible(lapply(seq_along(cov_list), function(i) names(cov_list[[i]])[[2]] <<- smps[i]))
  invisible(lapply(seq_along(fpkm_list), function(i) names(fpkm_list[[i]])[[2]] <<- smps[i]))

  f <- function(x, y)dplyr::full_join(x, y, by="t_name")
  fpkm <- Reduce(f, fpkm_list)
  cov <- Reduce(f, cov_list)

  # output file ----
  prjn <- basename(dirname((dirname(bamdir))))
  alnd <- basename(dirname(bamdir))
  fpkmout <- paste0(res_dir, "/", prjn, "_", alnd, "_FPKM.txt")
  covout <- paste0(res_dir, "/", prjn, "_", alnd, "_cov.txt")

  write.table(fpkm, fpkmout, quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(cov, covout, quote = F, sep = "\t", row.names = F, col.names = T)

}


