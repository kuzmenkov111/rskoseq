#' Create directory for RNA-seq pipeline
#' @description Create Project directories for RNA-seq pipeline using hisat2-stringtie.
#'    Then fastq files move to './project/fastq'　directory, followning execution of 'rskoseq::rep_hisat2', for read mapping using hisat2.
#' @usage dir_template(prj)
#' @param prj character: path of a project name
#' @examples #
#' # project_name <- "~/pub/dat/sampledata/rnaseq/project1"
#' # dir_template(prj = project_name)
#' # system(paste0("tree ", project_name))
#' # project1
#' # ├── command_log.txt
#' # ├── fastq
#' # ├── qa
#' # ├── res_hisat2
#' # │   └── failalign
#' # └── res_stringtie
#' #     ├── ballgown
#' #     ├── gff
#' #     ├── mgff
#' #     └── tab
#' # 9 directories, 1 file
#'
#' @export
dir_template <- function(prj){
  if (file.exists(prj)){
    stop(paste0("The directory ", prj, " still exists."))
  } else {
    com <- paste0(
      "mkdir ", prj, " \\\n",
      prj, "/fastq ", " \\\n",
      prj, "/qa ", " \\\n",
      prj, "/res_hisat2", " \\\n",
      prj, "/res_hisat2/failalign ", " \\\n",
      prj, "/res_stringtie ", " \\\n",
      prj, "/res_stringtie/ballgown ", " \\\n",
      prj, "/res_stringtie/gff ", " \\\n",
      prj, "/res_stringtie/mgff", " \\\n",
      prj, "/res_stringtie/tab", " \\\n"
    )
    system(com, wait = TRUE)
    cat(com)

    # command log ----
    com_log <-  paste0(prj, "/command_log.txt")
    fc <- file(com_log, "w")
    writeLines(com, fc)
    close(fc)
  }

}
