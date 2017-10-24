#' Replicate execution of Transcript assembly and quantification for RNA-Seq using stringtie
#' @description Description
#' @usage rep_stringtie(prj, bam_dir, suffix_bam, guide_gff, res_dir)
#' @param prj character: project directory path, created by 'rskoseq::dir_template'
#' @param bam_dir character: The directory path of soreted bam files. The default is paste0(prj, "/res_hisat2"),
#'     created by 'rskoseq::dir_template'.
#' @param suffix_bam character: The default is ".sort.bam".
#' @param guide_gff The file path of guide gff.
#' @param res_dir output directory path, the default is 'paste0(prj, "/res_stringtie")'
#' @param ... additional options of stringtie.
#' @examples #
#' # project_name <- "~/pub/dat/sampledata/rnaseq/project1"
#' # guide <- "~/db/index/hisat2_idx/SD0218_11_a_contig01.gff"
#' # rep_stringtie(prj=project_name, guide_gff=guide)
#' @importFrom utils tail write.table
#' @export
rep_stringtie <- function(prj, bam_dir=paste0(prj, "/res_hisat2"), suffix_bam=".sort.bam", guide_gff, res_dir=paste0(prj, "/res_stringtie")){
  # argument check: stringtie program PATH ----
  if (!any(grep("stringtie", unlist(strsplit(Sys.getenv("PATH"), ":"))))){
    stop("There is not stringtie program, or the PATH does not found.")
  }

  # argument check: project directory ----
  if(!file.exists(prj)){
    stop(paste0("\'", prj, "\'", " does not fount."))
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

  # stringtie execution ----
  ## command log ----
  com_log <-  paste0(prj, "/command_log.txt")
  con <- file(com_log, "a")
  ## execute ----
  for(i in seq_along(bamfls)){
    com <- paste0("stringtie -p 4 ", bamfls[i], " -G ", guide_gff, " -o ", res_gff[i])
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
  writeLines(com, con)

  # stringtie execution using merged gff ----
  ## define output files ----
  mgff <- paste0(res_dir, "/merged.gff")
  res_mgff <- paste0(res_dir, "/mgff/", smps, "_m.gff")
  res_ballgown <- paste0(res_dir, "/ballgown/", smps)
  ## replicate execution ----
  for(i in seq_along(bamfls)){
    dir.create(res_ballgown[i])
    com3 <-
      paste0("stringtie -p 4 ", bamfls[i],
             " -b ", res_ballgown[i],
             " -G ",  mgff,
             " -o ", res_mgff[i])
    cat(paste0(com3, " \n"))
    system(com3)

    ## write command log ----
    writeLines(com, con)

  }

  close(con)

}


