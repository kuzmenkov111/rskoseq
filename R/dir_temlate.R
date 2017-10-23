#' Create directory for RNA-seq pipeline
#' @description Create Project directories for RNA-seq pipeline using hisat2-stringtie.
#'    Then fastq files move to './project/fastq'　directory, followning execution of 'rskoseq::rep_hisat2', for read mapping using hisat2.
#' @usage dir_template(prj)
#' @param prj character: path of a project name
#' @examples #
#' dir_template(prj = "~/pub/dat/sampledata/rnseq/project1")
#' # project1
#' # ├── fastq
#' # ├── qa
#' # ├── res_hisat2
#' # │.. ├── failalign
#' # │.. ├── mat
#' # │.. └── sortbam
#' # └── res_stringtie
#' #    ├── ballgown
#' #    └── gff
#' #
#' # 9 directories, 0 files
#' @export
dir_template <- function(prj){
  com <- paste0(
    "mkdir ", prj, " \\\n",
    prj, "/fastq ", " \\\n",
    prj, "/qa ", " \\\n",
    prj, "/res_hisat2", " \\\n",
    prj, "/res_hisat2/sortbam ", " \\\n",
    prj, "/res_hisat2/failalign ", " \\\n",
    prj, "/res_hisat2/mat ", " \\\n",
    prj, "/res_stringtie ", " \\\n",
    prj, "/res_stringtie/ballgown ", " \\\n",
    prj, "/res_stringtie/gff ", " \\\n"
  )
  cat(com)
  system(com)
}
