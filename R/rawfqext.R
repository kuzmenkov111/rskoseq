#' rawfqext
#' @description merge NEXTSEQ output fastq, which created by BaseSpace
#' @usage rawfqext(rawfq_dir, prefix)
#' @param rawfq_dir fastq containing directory
#' @param prefix sample name
#' @examples
#' # rawfq_dir <- "~/pub/sampledata/fastq/NEXTSEQ_OUT"
#' # prefix <- c("S1","S2")
#' # rawfqext(rawfq_dir, prefix)
#' @export

rawfqext <- function(rawfq_dir, prefix){

  # collect fastq containing directory -----
  dirs <- list.files(rawfq_dir,  full.names = T)

  # argument check:
  if (length(dirs) != length(prefix)){
    stop("length of fasq directory and prefix must to be same")
  }

  for (i in seq_along(dirs)){
    # collect fasq file path ----
    gzr1 <- list.files(dirs[i], ".*_R1_.*.fastq.gz$", full.names = T)
    gzr2 <- list.files(dirs[i], ".*_R2_.*.fastq.gz$", full.names = T)
    r1 <- list.files(dirs[i], ".*_R1_.*.fastq$", full.names = T)
    r2 <- list.files(dirs[i], ".*_R2_.*.fastq$", full.names = T)

    # still uncompress ----
    if (identical(gzr1, character(0)) ){
      if (!identical(r1, character(0))){ # R1 uncompressed fastq
        com1 <- paste0("cat ",
                       paste(r1, collapse = " "),
                       " | gzip -c > ",
                       rawfq_dir, "/",
                       prefix[i], "_R1.fastq.gz")
        system(com1)

      } else if (!identical(r2, character(0))){ # R2 uncompressed fastq
        com2 <- paste0("cat ",
                       paste(r2, collapse = " "),
                       " | gzip -c > ",
                       rawfq_dir, "/",
                       prefix[i], "_R2.fastq.gz")
        system(com2)

      }

    # compressed fastq gz to merged gz ----
    } else if (!identical(gzr1, character(0))){
      if (!identical(gzr1, character(0))){ # R1

        com1g <- paste0("cat ",
                        paste(gzr1, collapse = " "),
                        " > ",
                        rawfq_dir, "/",
                        prefix[i], "_R1.fastq.gz")
        system(com1g)
        cat(com1g)

      } else if (!identical(gzr2, character(0))){ # R2
        com2g <- paste0("cat ",
                        paste(gzr2, collapse = " "),
                        " > ",
                        rawfq_dir, "/",
                        prefix[i], "_R2.fastq.gz")
        system(com2g)
        cat(paste0(com2g, "\n"))
      }
    }
  }
}


