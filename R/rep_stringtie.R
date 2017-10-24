#' Consecutive execution of Transcript assembly and quantification for RNA-Seq using stringtie
#' @description Consecutive processiong of stringtie for multi samples, then FPKM and cov data table are created.
#' @usage rep_stringtie(prj, bam_dir, suffix_bam, guide_gff, res_dir, ...)
#' @param prj character: project directory path, created by 'rskoseq::dir_template'
#' @param bam_dir character: The directory path of soreted bam files. The default is paste0(prj, "/res_hisat2"),
#'     created by 'rskoseq::dir_template'.
#' @param suffix_bam character: The default is ".sort.bam".
#' @param guide_gff The file path of guide gff.
#' @param res_dir output directory path, the default is 'paste0(prj, "/res_stringtie")'
#' @param ... additional options of stringtie. E.g. "-e"
#' @examples #
#' # project_name <- "~/pub/dat/sampledata/rnaseq/project1"
#' # guide <- "~/db/index/hisat2_idx/SD0218_11_a_contig01.gff"
#' # rep_stringtie(prj=project_name, guide_gff=guide)
#' @importFrom utils tail write.table read.table
#' @export
rep_stringtie <- function(prj, bam_dir=paste0(prj, "/res_hisat2"), suffix_bam=".sort.bam", guide_gff, res_dir=paste0(prj, "/res_stringtie"), ...){
  # argument check: stringtie program PATH ----
  if (!any(grep("stringtie", unlist(strsplit(Sys.getenv("PATH"), ":"))))){
    stop("There is not stringtie program, or the PATH does not found.")
  }

  # argument check: project directory ----
  if(!file.exists(prj)){
    stop(paste0("\'", prj, "\'", " does not found."))
  }

  # argument check: guide_gff ----
  if (!file.exists(guide_gff)){
    stop(paste0("There is not ", guide_gff))
  }

  # argument check: path of sorted bam files and get all samples name----
  if (file.exists(bam_dir)){
    bamfls <- list.files(bam_dir, suffix_bam, full.names = T)
  }
  if (identical(bamfls, character(0))){
    stop(paste0("There is not sorted bam files, or the suffix of these bam files is different from ",  suffix_bam))
  }

  # collect sample names and result gff files ----
  smps <- sub(".sort.bam", "",
              sapply(strsplit(bamfls, "/"), function(x)tail(x, 1)))
  res_gff <- paste0(res_dir, "/gff/", smps, ".gff")

  # stringtie additional options ----
  if (!missing(...)){
    add_op <- paste0(" ", ..., " ")
  }else{
    add_op <- ""
  }

  # stringtie execution ----
  ## command log ----
  com_log <-  paste0(prj, "/command_log.txt")
  con <- file(com_log, "a")
  writeLines("# stringtie", con)
  ## execute ----
  for(i in seq_along(bamfls)){
    com <- paste0("stringtie -p 4 ", bamfls[i], " -G ", guide_gff, " -o ", res_gff[i], add_op)
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
  system(com2)

  ## over write command log ----
  writeLines(paste0("# stringtie --merge \n", com2), con)

  # stringtie execution using merged gff ----
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
      paste0("stringtie -p 4 ", bamfls[i],
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
  write.table(fpkm, paste0(prj, "/res_stringtie/FPKM.txt"),
              quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(cov, paste0(prj, "/res_stringtie/cov.txt"),
              quote = F, sep = "\t", row.names = F, col.names = T)


}


