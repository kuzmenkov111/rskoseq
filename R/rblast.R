#' Execution of local blast search usign system command and output result as data frame
#' @description This functions returns a data frame as BLAST output, if multiple databese gave as an arugment, they are all in one,
#'     which fomat is outfmt 6 and additional column.
#' @usage rblast(in_f, out_f, program, db, ...)
#' @param in_f character: input file path of multifasta format.
#' @param out_f character: output file path or stdout as "-".
#' @param program character: select a blast search program "blastn", "blastp", "blastx", "tblastn", "tblastx"
#' @param db blast database name or file path of fasta: if you gives fasta file path, 'makeblastdb' was executed.
#'  it is not file path. if multiple database using, corresponding out put files path gave as 'out_f'.
#' @param ... additional parameter as character strings E.g. "-num_threads 4 -task megablast".
#' @return blast output of 'outfmt 6' and several other column returnd as named data frame.
#'  if all query "No hits found", could not create data.frame
#' @importFrom readr read_delim
#' @importFrom utils tail
#' @examples
#' ## create blast data base
#' # com0 <- c("makeblastdb -in" TAIR10_chr_all.fas "-out" TAIR10 "-dbtype" nucl "-parse_seqids")
#' # system(com0, intern = T)
#'
#' ## sample fasta of rsko package
#' # filepath1 <- system.file("extdata", "AtMlos.fna", package="rskodat")
#' # bndb <-  "~/db/cdna/TAIR10_cdna"
#' # bnout <- rblast(in_f = filepath1, out_f = "-", program = "blastn", db = bndb, "-num_threads 4")
#'
#' ## multiple db
#' # d <- "~/db/cdna"
#' # dbs <- c("TAIR10_cdna", "Vv_ref_rna")
#' # pdbs <- paste0(d, "/", dbs)
#' # bnouts <- rblast(in_f = filepath1, out_f = "-", program = "blastn", db = pdbs, "-num_threads 4")
#' @export
# in_f = filepath1; out_f="-"; program="blastn"; db = pdbs
rblast <- function(in_f, out_f, program, db,  ...){
  # argument check: program PATH ----
  blpath <- suppressWarnings(system(paste("which", program), intern = T))
  if (identical(blpath, character(0))){
    stop(paste0("There is not ", program, ", or the PATH of this does not found."))
  }

  # multiple db or not
  if (out_f != "-" & length(out_f) != length(db)){

    stop("'out_f' must be a vector which has same length of 'db', or '-'. ")

  } else {
    ## data base name ----
    dbname <- sapply(strsplit(db, "/"), function(x)tail(x, 1))

    ## blast search ----
    bout_list <- lapply(seq_along(db), function(i){
      ### command ----
      com <- paste(blpath, "-query", in_f, "-db", db[i], "-out", out_f,
                   "-outfmt", "\"6 std qlen slen sstrand salltitles\"" , ...,  sep = " ")
      cat(paste0(com, " \n"))

      ### execution of system command
      res <- system(com, intern = TRUE)

      ### convert data frame, or read output files ----
      if(out_f=="-"){
        bout <- lapply(res, function(x){unlist(strsplit(x, "\t"))})
        bout <- data.frame(do.call(what = rbind, args = bout),
                           row.names=NULL, stringsAsFactors = F)
        transform(bout, dbname=rep(dbname[i], nrow(bout)))

        } else {
          # in case of read blast output file
          bout <- readr::read_delim(out_f, delim = "\t", col_types = c("ccnnnnnnnnnnnncc"))
        }
    })
  }

  # column names
  bout <- do.call(rbind, bout_list)
  names(bout)[1:16] <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                   "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                   "qlen", "slen", "sstrand","description")

  # return data frame as tbl
  return(bout)
}
