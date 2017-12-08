#' Counting kmer of nucleotide sequence
#' @description kmer counting from fastq or fasta format file and return kmer table
#' @usage kmer_cnt(filepath, format, k, n, per_site)
#' @return kmer table or kmer count at position
#' @param filepath Input file path
#' @param format Choose "fastq" or "fasta" format
#' @param k Length of kmer
#' @param n Random sampling numbers for ShortRead::FastqSampler
#' @param per_site logical: kmer count per site of all reads of fastq, or not
#' @importFrom ShortRead sread width readFastq FastqSampler yield
#' @importFrom stringr str_sub
#' @importFrom dplyr left_join mutate_at funs %>%
#' @importFrom plyr .
#' @importFrom Biostrings readDNAStringSet subseq
#' @importFrom BiocGenerics table
#' @importFrom stats setNames
#' @examples
#' ## fastq which reads are same
#' # fq <- list.files(system.file("extdata/E-MTAB-1147", package = "ShortRead"), full.names = TRUE)
#' # kmer_dat1 <- kmer_cnt(filepath = fq[1], format = "fastq", k = 7, per_site = TRUE)
#' # kmer_dat2 <- kmer_cnt(filepath = fq[1], format = "fastq", k = 7, per_site = FALSE)
#' @export

kmer_cnt <- function(filepath, format, k, n=NULL, per_site=TRUE){
  # input file format and read sequence ----
  if (format =="fastq" & is.null(n)){
    fas <- ShortRead::sread( ShortRead::readFastq(filepath))

  } else if (format == "fastq" & !is.null(n)){
    sampler <- ShortRead::FastqSampler(filepath, n)
    fas <- ShortRead::sread(ShortRead::yield(sampler))
    close(sampler)

  } else if(format == "fasta"){
    fas <- Biostrings::readDNAStringSet(filepath, format = "fasta")

  }

  # Make all kmer list specified kmer length ----
  nuc <-c("A","T","G","C")
  d <- as.matrix(expand.grid(rep(list(nuc), k)))
  kn <-data.frame(kmer=sapply(1:nrow(d), function(i)paste(d[i,], collapse = "")), stringsAsFactors = F)


  # In case of different read length ----
  kmer = NULL; count = NULL;
  if(length(unique(ShortRead::width(fas))) > 1 && format=="fastq"){
    # create kmer table
    ## A list of kmers per reads ----
    kl <- lapply(fas, function(x){
      st <- 1:(length(x)-k+1)
      ed <- st+k-1
      stringr::str_sub(toString(x), start = st, end = ed)
      })

    ## kmer table ----
    if (per_site == F){
      kmer_dat <- setNames(data.frame(table(unlist(kl)), stringsAsFactors = F), c("kmer", "count")) %>%
        dplyr::mutate(kmer=as.character(kmer)) %>%
        dplyr::left_join(kn, ., by="kmer") %>%
        dplyr::mutate(count= ifelse(is.na(count), 0, count))

    } else if (per_site ==T){
      ## column: read, row:start position kmer matrix ----
      kmer_mat <- sapply(kl, "[", 1:as.numeric(max(names(table(sapply(kl, length))))))

      ## A list of kmer table per site----
      l_ktab <- apply(kmer_mat, 1, function(x){
        data.frame(kmer=names(table(x)), count=as.integer(table(x)), stringsAsFactors = F)
      })
      ## merge all kmer table per site by kmer ----
      kmer_dat <- Reduce(function(x,y){dplyr::left_join(x,y, by="kmer")}, c(list(kn), l_ktab)) %>%
        stats::setNames(., c("kmer", 1:as.numeric(max(names(table(sapply(kl, length))))))) %>%
        dplyr::mutate_at(., -1, dplyr::funs(ifelse(is.na(.), 0, .)))

    }

  # In case of fasta format ----
  }else if(length(unique(ShortRead::width(fas))) >1 && format=="fasta"){
    ## kmer list per seq ----
    kl <- lapply(fas, function(x){
                   st <- 1:(length(x)-k+1)
                   ed <- st+k-1
                   stringr::str_sub(toString(x), start = st, end = ed)
                 })
    ## kmer table ----
    kmer_table <- data.frame(kmer=names(table(unlist(kl))),
                             count=as.integer(table(unlist(kl))), stringsAsFactors = F) %>%
      dplyr::left_join(kn, ., by="kmer") %>%
      mutate_at(., -1, dplyr::funs(ifelse(is.na(.), 0, .)))

  # all read length is same -----
  }else if(length(unique(ShortRead::width(fas)))==1 & format == "fastq"){
    ## start position from read width ----
    starts <- 1:(max(ShortRead::width(fas))-k+1)

    ## A list of kmer table per site -----
    kmer_cnts <- lapply(seq_along(starts), function(i){
        setNames(as.data.frame(BiocGenerics::table(Biostrings::subseq(fas, start =starts[i], width = k)), stringsAsFactors = F),
                 c("kmer", i))
        })

    ## kmer table per site ----
    kmer_dat <- Reduce(function(x,y){dplyr::left_join(x, y, by="kmer")},
                       c(list(kn), kmer_cnts)) %>%
      dplyr::mutate_at(., -1, dplyr::funs(ifelse(is.na(.), 0, .)))
    if (per_site == F){
      kmer_dat <- data.frame(kmer=kmer_dat$kmer, count=apply(kmer_dat[-1], 1, sum))
    }
  }
  # return kmer table
  return(kmer_dat)
}
