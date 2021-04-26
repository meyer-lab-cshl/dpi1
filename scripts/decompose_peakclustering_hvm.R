#################
## libraries ####
#################
library(optparse)
source("scripts/myFastICA.R")



#################
## functions ####
#################
peakClustersDecomposedCtss <- function(ica, gaussian_window_size_half = 20,
                                       bedLine = matrix(NA, ncol = 0,
                                                        nrow = 0)) {
  buf <- lapply(
    ica$signal_components,
    function(component) {
      res <- peakClustersFromCtssVec(
        ica$rescaleS[, component],
        gaussian_window_size_half = gaussian_window_size_half,
        bedLine = bedLine
      )
      if (is.null(dim(res))) {
        return(NA)
      }
      res[, "name"] <- paste(res[, "name"], component, sep = ";comp:")
      res
    }
  )

  idx <- which(sapply(buf, function(res) is.null(dim(res))))
  if (length(idx) == length(buf)) {
    return(NA)
  }
  if (length(idx) > 0) {
    buf <- buf[-idx]
  }

  chrom <- unlist(lapply(buf, function(b) b[, "chrom"]))
  start <- unlist(lapply(buf, function(b) b[, "start"]))
  stop <- unlist(lapply(buf, function(b) b[, "stop"]))
  name <- unlist(lapply(buf, function(b) b[, "name"]))
  score <- unlist(lapply(buf, function(b) b[, "score"]))
  strand <- unlist(lapply(buf, function(b) b[, "strand"]))
  bedTable <- cbind(chrom, start, stop, name, score, strand)
  bedTable
}

peakClustersFromCtssVec <- function(ctssVec, gaussian_window_size_half = 20,
                                    bedLine = matrix(NA, ncol = 0, nrow = 0)) {
  bedTable <- ctssVec2bedTable(ctssVec)
  res <- peakClusters(bedTable,
    gaussian_window_size_half = gaussian_window_size_half,
    bedLine = bedLine
  )
  res
}

ctssVec2bedTable <- function(ctssVec) {
  coords <- t(sapply(
    names(ctssVec),
    function(str) {
      tmp <- strsplit(str, ":")[[1]]
      chrom <- tmp[1]
      tmp <- strsplit(tmp[2], "\\.")[[1]]
      start <- tmp[1]
      tmp <- strsplit(tmp[3], ",")[[1]]
      stop <- tmp[1]
      strand <- tmp[2]
      c(chrom, start, stop, strand)
    }
  ))
  colnames(coords) <- c("chrom", "start", "stop", "strand")
  coords <- as.data.frame(coords, stringsAsFactors = FALSE)
  coords$start <- as.numeric(coords$start)
  coords$stop <- as.numeric(coords$stop)
  score <- ctssVec
  bedTable <- cbind(coords, score)
  bedTable
}

peakClusters <- function(bedTable, gaussian_window_size_half = 20,
                         bedLine = matrix(NA, ncol = 0, nrow = 0)) {
  ### prep
  if (nrow(bedTable) == 0) {
    return(NA)
  }
  pos <- range(bedTable$start)
  vec <- rep(0, pos[2] - pos[1] + 1)
  vec[bedTable$start - pos[1] + 1] <- bedTable$score
  vec[is.na(vec)] <- 0

  ### smoothing & cluster detection
  vec.padding <- c(
    rep(0, gaussian_window_size_half), vec,
    rep(0, gaussian_window_size_half)
  )
  vec.padding.smooth <- filter(
    vec.padding,
    dnorm(seq(-2, 2, by = (2 / gaussian_window_size_half)))
  )
  cutoff <- max(median(vec.padding.smooth, na.rm = TRUE), 3)

  ### assigning cluster IDs
  clst <- rep(0, length(vec.padding.smooth))
  clst[vec.padding.smooth > cutoff] <- 1
  clst_id <- 0
  for (i in 2:length(vec.padding.smooth)) {
    if (clst[i] == 0) {
      next
    }
    if (clst[i - 1] == 0) {
      clst_id <- clst_id + 1
    }
    clst[i] <- clst_id
  }
  clst <- clst[(gaussian_window_size_half + 1):
  (gaussian_window_size_half + length(vec))]
  cluster_ids <- unique(clst)
  cluster_ids <- cluster_ids[cluster_ids != 0]
  if (length(cluster_ids) == 0) {
    return(NA)
  }

  ### prep for output
  chrom <- unique(bedTable$chrom)
  strand <- unique(bedTable$strand)
  chrom <- chrom[!is.na(chrom)]
  strand <- strand[!is.na(strand)]
  if (length(chrom) < 1) {
    return(NA)
  }
  if (length(strand) < 1) {
    return(NA)
  }
  if (length(chrom) > 1) {
    print(sprintf("multiple chroms:%s", paste(chrom, sep = ",", collapse = ",")),
      file = stderr()
    )
    return(NA)
  }
  if (length(strand) > 1) {
    print(sprintf("multiple strands:%s", paste(strand, sep = ",",
                                               collapse = ",")),
      file = stderr()
    )
    return(NA)
  }

  ### cluster output as bed table
  clusterBedTable <- t(sapply(cluster_ids, function(i) {
      idx <- which(clst == i)
      highestPeak <- idx[which.max(vec[idx])[1]] + pos[1] - 1
      cluster_pos <- range(idx) + pos[1] - 1

      origin_cluster <- sprintf(
        "%s:%s..%s,%s", chrom, as.integer(pos[1]),
        as.integer(pos[2] + 1), strand
      )
      if (nrow(bedLine) != 0) {
        origin_cluster <- bedLine$name
      }
      name <- sprintf(
        "%s;peak:%s:%s..%s,%s;peakCounts:%s",
        origin_cluster,
        chrom, as.integer(highestPeak), as.integer(highestPeak + 1), strand,
        max(vec[idx])
      )
      tmp <- c(chrom, cluster_pos[1], cluster_pos[2] + 1, name, 1000, strand)
      names(tmp) <- c("chrom", "start", "stop", "name", "score", "strand")
      tmp
    }
  ))
  clusterBedTable
}

getCtssCounts <- function(bedLine, infile_ctss, force_strandedness = TRUE,
                          noise_subtraction_ratio = 0) {
  if (nrow(bedLine) == 0) {
    return(c())
  }
  strand <- NA
  if (length(grep(".fwd.bw$", infile_ctss)) > 0) {
    strand <- "+"
  }
  if (length(grep(".rev.bw$", infile_ctss)) > 0) {
    strand <- "-"
  }
  if (is.na(strand)) {
    return(c())
  }
  if ((force_strandedness) && (strand != bedLine$strand)) {
    return(c())
  }
  command <- sprintf(
    "bigWigToBedGraph %s /dev/stdout -chrom=%s -start=%s -end=%s",
    infile_ctss, bedLine$chrom, bedLine$start, bedLine$stop
  )
  res <- system(command, intern = TRUE)

  if (length(res) == 0) {
    return(c())
  }
  res <- t(sapply(res, function(str) strsplit(str, "\t")[[1]]))
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  for (i in c(2, 3, 4)) {
    res[, i] <- as.numeric(res[, i])
  }
  colnames(res) <- c("chrom", "start", "stop", "score")

  idx <- which((res$stop - res$start) >= 2)
  for (i in idx) {
    tmp <- as.data.frame(cbind(
      res$chrom[i],
      res$start[i]:(res$stop[i] - 1),
      (res$start[i] + 1):res$stop[i], res$score[i]
    ))
    for (j in c(2, 3, 4)) {
      tmp[, j] <- as.numeric(as.character(tmp[, j]))
    }
    colnames(tmp) <- c("chrom", "start", "stop", "score")
    res <- rbind(res, tmp)
  }
  name <- paste(
    res$chrom, ":",
    as.integer(res$start), "..", as.integer(res$stop),
    ",", strand,
    sep = ""
  )
  res <- as.data.frame(cbind(res$chrom, res$start, res$stop, name, res$score,
                             strand), stringsAsFactors = FALSE)
  if (length(idx) > 0) {
    res <- res[-idx, ]
  }

  colnames(res) <- c("chrom", "start", "stop", "name", "score", "strand")
  res[, "score"] <- as.numeric(res[, "score"])
  res[, "start"] <- as.integer(as.numeric(res[, "start"]))
  res[, "stop"] <- as.integer(as.numeric(res[, "stop"]))
  res <- res[order(res$start), 1:ncol(res), drop = FALSE]

  if (noise_subtraction_ratio > 0) {
    res$score <- res$score - (max(res$score) * noise_subtraction_ratio)
  }
  res <- res[res$score > 0, ]
  res
}

getCtssCountsTable <- function(bedLine, infiles, noise_subtraction_ratio = 0) {
  # prep
  tbl <- c()
  tmpdir <- system("mktemp -d", intern = TRUE)
  tmplist <- sprintf("%s/list.txt", tmpdir)
  tmpmat <- sprintf("%s/mat.txt", tmpdir)
  command <- sprintf(
    "scripts/getCtssCountsTable.sh -c %s -s %s -e %s -i %s > %s",
    bedLine$chrom, bedLine$start, bedLine$stop, tmplist, tmpmat, tmpmat
  )
  if (bedLine$strand == "+") {
    infiles <- infiles[grep(".fwd.bw$", infiles)]
  } else if (bedLine$strand == "-") {
    infiles <- infiles[grep(".rev.bw$", infiles)]
  } else {
    stop(sprintf("ERROR: wrong strand:%s", bedLine$strand))
  }

  # main
  tryCatch({
      write.table(infiles, quote = FALSE, row.names = FALSE, col.names = FALSE,
                  file = tmplist)
      system(command, intern = TRUE)
      tbl <- read.table(tmpmat, row.names = 1, header = TRUE,
                        check.names = FALSE)
      rownames(tbl) <- paste(rownames(tbl), bedLine$strand, sep = ",")
      system(sprintf("rm -rf %s", tmpdir))
    },
    #finally = system(sprintf("rm -rf %s", tmpdir))
    finally=TRUE
  )

  # post processing
  sbtrct <- apply(tbl, 2, max) * noise_subtraction_ratio
  tbl <- t(t(tbl) - sbtrct)
  tbl[tbl < 0] <- 0
  idx1 <- which(apply(tbl, 1, max) > 0)
  idx2 <- which(apply(tbl, 2, max) > 0)
  tbl <- tbl[idx1, idx2, drop = FALSE]
  tbl
}

peakClustersFromCtssVec_print <- function(ctss, gaussian_window_size_half,
                                          bedLine, outfile) {
  ctss <- rowSums(ctss)
  res <- peakClustersFromCtssVec(ctss, gaussian_window_size_half, bedLine)
  if (is.null(dim(res))) {
    write.table(bedLine, file=outfile, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE, append=TRUE)
  } else {
    write.table(res, file=outfile, sep = "\t", quote = FALSE, row.names = FALSE,
                col.names = FALSE, append=TRUE)
  }
}

############
## data ####
############

## command line arguments ####
option_list <- list(
  make_option(c("-a", "--analysis"), action="store",
              dest="analysis",
              type="character", help="Type of analysis: spi | dpi
              [default: %default].", default='spi'),
  make_option(c("-o", "--outfile"), action="store",
              dest="outfile",
              type="character", help="Path the output file [default: %default].",
              default=NULL),
  make_option(c("--icafile"), action="store",
              dest="icafile",
              type="character", help="Path the ica output file [default: %default].",
              default=NULL),
  make_option(c("-c", "--tagclusters"), action="store",
              dest="tagclusters",
              type="character", help="Path to file with tag clusters to be
              decomposed; in `bed.gz` format [default: %default].",
              default=NULL),
  make_option(c("-p", "--ctssprefix"), action="store",
              dest="ctssprefix",
              type="character", help="Prefix of individual bed.gz files with
              transcription start sites [default: %default].",
              default=NULL),
  make_option(c("--path"), action="store",
              dest="path",
              type="character", help=" [default: %default].",
              default=NULL),
  make_option(c("--exclude"), action="store",
              dest="exclude_prefix",
              type="character", help=" [default: %default].",
              default=NULL),
  make_option(c("--pattern"), action="store",
              dest="pattern",
              type="character", help=" [default: %default].",
              default=NULL),
  make_option(c("-w", "--window"), action="store",
              dest="window", type="integer",
              help="Gaussian window size [default: %default].",
              default=5),
  make_option(c( "--bound"), action="store",
              dest="bound", type="integer",
              help=" [default: %default].",
              default=5),
  make_option(c("-l", "--length"), action="store",
              dest="length", type="integer",
              help="Cluster length to decompose [default: %default].",
              default=50),
  make_option(c("-r", "--ratio"), action="store",
              dest="ratio", type="double",
              help="Noise substraction ratio [default: %default].",
              default=0.1),
  make_option(c("--showProgress"), action="store_true",
              dest="verbose",
              default=FALSE, type="logical", help="If set, progress messages
                about analyses are printed to standard out ",
              "[default: %default]."),
  make_option(c("--debug"), action="store_true",
              dest="debug", default=FALSE, type="logical",
              help="If set, predefined arguments are used to test the script",
              "[default: %default].")
)

args <- parse_args(OptionParser(option_list=option_list))

if (args$debug) {
    args <- list()
    #args$outfile <- "/Users/hannah/software/dpi1/_test/snakemake/outPooled/tc.long.spi.bed"
    #args$tagclusters <- "/Users/hannah/software/dpi1/_test/snakemake/outPooled/tc.long.bed.gz"
    #args$ctssprefix <- "/Users/hannah/software/dpi1/_test/snakemake/outPooled/ctssTotalCounts"
    args$outfile <- "/Users/hannah/software/dpi1/_test/snakemake/outPooled/tc.long.decompose/aaaaa.decompose_smoothing.bed"
    args$icafile <- "/Users/hannah/software/dpi1/_test/snakemake/outPooled/tc.long.decompose/aaaaa.ica.txt"
    args$tagclusters <- "/Users/hannah/software/dpi1/_test/snakemake/outPooled/tc.long.files/aaaaa"
    args$path <- "/Users/hannah/software/dpi1/_test/snakemake/outCounts"
    #args$window <- 5
    #args$length <- 50
    #args$ratio <- 0.1
    args$window <- 20
    args$length <- 100
    args$ratio <- 0
    args$verbose <- TRUE
    args$bound <- 5
    args$pattern <- ".bw$"
    args$exclude_prefix <- NA
    verbose = TRUE
    #args$analysis <- 'spi'
    args$analysis <- 'dpi'

}
################
## analysis ####
################

formated <- sapply(seq_along(args), function(x) {
  paste(names(args)[x], ": ", args[x], sep="")
})
if (args$verbose) message(paste(formated, collapse="\n"))

if (file.exists(args$outfile)) file.remove(args$outfile)

if (args$analysis == "spi") {
  base <- read.table(args$tagclusters, sep = "\t", as.is = TRUE, nrow = -1)
  colnames(base) <- c("chrom", "start", "stop", "name", "score", "strand")

  infile_ctss <- c(
    sprintf("%s.fwd.bw", args$ctssprefix),
    sprintf("%s.rev.bw", args$ctssprefix)
  )

  for (i in 1:nrow(base)) {
    bedLine <- base[i, ]
    if (args$verbose) message("#", paste(bedLine, collapse = "\t"))
    if ((bedLine$stop - bedLine$start) < args$length) {
      write.table(bedLine, file=args$outfile, sep = "\t", quote = FALSE,
                  row.names = FALSE, col.names = FALSE, append=TRUE)
      next
    }
    ctss <- getCtssCountsTable(bedLine, infile_ctss, args$ratio)
    peakClustersFromCtssVec_print(ctss, args$window, bedLine, args$outfile)
  }
} else if (args$analysis == "dpi") {
  if (file.exists(args$icafile)) file.remove(args$icafile)
  infile_ctss <- dir(path = args$path, pattern = args$pattern,
                     full.names = TRUE)

  ### special file exclusion
  if (!is.na(args$exclude_prefix)) {
    infile_ctss <- grep(args$exclude_prefix, infile_ctss, value = TRUE,
                        invert = TRUE)
  }

  base <- read.table(args$tagclusters, sep = "\t", as.is = TRUE, nrow = -1)
  colnames(base) <- c("chrom", "start", "stop", "name", "score", "strand")

  for (i in 1:nrow(base)) {
    bedLine <- base[i, ]
    if (args$verbose) message("#", paste(bedLine, collapse = "\t"))
    if ((bedLine$stop - bedLine$start) < args$length) {
      write.table(bedLine, file=args$outfile, sep = "\t", quote = FALSE,
                  row.names = FALSE, col.names = FALSE, append=TRUE)
      next
    }
    ctss <- getCtssCountsTable(bedLine, infile_ctss,
                               noise_subtraction_ratio = args$ratio)
    ica <- myFastICA(ctss, verbose = args$verbose,
                     n.comp.upper_bound = args$bound)

    if (length(ica) == 0) {
      peakClustersFromCtssVec_print(ctss, args$window, bedLine, args$outfile)
    } else {
      res <- peakClustersDecomposedCtss(ica,
                                        gaussian_window_size_half = args$window,
                                        bedLine)
      if (is.null(dim(res))) {
        peakClustersFromCtssVec_print(ctss, args$window, bedLine, args$outfile)
      } else {
        write.table(res, file=args$outfile, sep = "\t", quote = FALSE,
                    row.names = FALSE, col.names = FALSE, append=TRUE)
      }
      write.table(ica$rescaleS, file=args$icafile, sep = "\t", quote = FALSE,
                  row.names = TRUE, col.names = FALSE, append=TRUE)
    }
  }
} else {
  stop ("Analysis type", args$analysis, "not found, has to be spi or dpi")
}

