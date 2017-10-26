#' Create directory for RNA-seq pipeline
#' @description Create Project directories for RNA-seq pipeline using hisat2-stringtie.
#'    Then fastq files move to './project/fastq'　directory, followning execution of 'rskoseq::rep_hisat2', for read mapping using hisat2.
#' @usage project_rnsq(prjd, alnd)
#' @param prjd character: path of a project directory name
#' @param alnd character: path of a alignment directory name, the default is "alignment1"
#' @examples
#' # # create new project
#' # prj <- "~/pub/dat/sampledata/rnaseq/project1"
#' # project_rnsq(prjd = prj)
#' # system(paste("tree", prj))
#' # # now there is a project directory, create another alignment directory.
#' # project_rnsq(prjd = prj, alnd = "alignment2")
#' # ./project1
#' # ├── alignment
#' # │   ├── command_log.txt
#' # │   ├── res_hisat2
#' # │   │   └── failalign
#' # │   └── res_stringtie
#' # │       ├── ballgown
#' # │       ├── gff
#' # │       ├── mgff
#' # │       └── tab
#' # ├── alignment2
#' # │   ├── res_hisat2
#' # │   │   └── failalign
#' # │   └── res_stringtie
#' # │       ├── ballgown
#' # │       ├── gff
#' # │       ├── mgff
#' # │       └── tab
#' # ├── fastq
#' # └── qa
#'
#' @export
project_rnsq <- function(prjd, alnd="alignment1"){
  # system command ----
  if (file.exists(prjd) & !any(grepl(alnd, list.files(prjd)))){
    # create another alignment directory
    # res_hisat2 dir ----
    h2dir <- paste0(prjd, "/", alnd, "/res_hisat2/failalign")
    dir.create(path = h2dir, recursive = T)

    # res_stringtie dir ----
    bgdir <- paste0(prjd, "/", alnd, "/res_stringtie/ballgown")
    tabdir <- paste0(prjd, "/", alnd, "/res_stringtie/tab")
    gffdir <- paste0(prjd, "/", alnd, "/res_stringtie/gff")
    mgffdir <- paste0(prjd, "/", alnd, "/res_stringtie/mgff")

    dir.create(path = bgdir, recursive = T)
    dir.create(path = tabdir, recursive = T)
    dir.create(path = gffdir, recursive = T)
    dir.create(path = mgffdir, recursive = T)

  } else if (file.exists(prjd) & any(grepl(alnd, list.files(prjd)))) {
    stop(paste0("The directory ", prjd, " still exists, and give different name of alignment directory from '", alnd, "'."))
  } else if (!file.exists(prjd)){
    # project directory ----
    dir.create(prjd)
    dir.create(paste(prjd, "fastq", sep="/"))
    dir.create(paste(prjd, "qa", sep="/"))


    # res_hisat2 dir ----
    h2dir <- paste0(prjd, "/", alnd, "/res_hisat2/failalign")
    dir.create(path = h2dir, recursive = T)

    # res_stringtie dir ----
    bgdir <- paste0(prjd, "/", alnd, "/res_stringtie/ballgown")
    tabdir <- paste0(prjd, "/", alnd, "/res_stringtie/tab")
    gffdir <- paste0(prjd, "/", alnd, "/res_stringtie/gff")
    mgffdir <- paste0(prjd, "/", alnd, "/res_stringtie/mgff")

    dir.create(path = bgdir, recursive = T)
    dir.create(path = tabdir, recursive = T)
    dir.create(path = gffdir, recursive = T)
    dir.create(path = mgffdir, recursive = T)

    # command log file ----
    file.create(paste(prjd, alnd, "command_log.txt", sep="/"))

  }
}






