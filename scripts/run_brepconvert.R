#!/usr/bin/env Rscript
# Run BrepConvert (Mallaby et al. 2023) on one locus so its calls can be
# compared against the topology-aware detector.
#
#   run_brepconvert.R <functional.fasta> <pseudogene.fasta> <transcripts.fasta> <out.csv>
#
# The functional/pseudogene split is supplied by the caller.  We pass the
# expression-confirmed "used" genes as functional rather than the annotation's
# Productive flag, since the annotation is a prediction and demonstrably wrong
# for several genes here.
#
# BrepConvert expects IMGT-gapped sequences.  Our germline is ungapped, which is
# fine: gaps only exist so positions line up with IMGT numbering, and the code
# strips them with gsub(".", "", ...) before doing any sequence work.  Reported
# coordinates are therefore raw offsets into the V gene, not IMGT positions.

suppressPackageStartupMessages({
  library(Biostrings)
  library(BrepConvert)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("usage: run_brepconvert.R <functional.fasta> <pseudogene.fasta> <transcripts.fasta> <out.csv>")
}
functional_fa <- args[1]
pseudogene_fa <- args[2]
transcripts_fa <- args[3]
out_csv <- args[4]

blat <- system.file("exe/blat", package = "BrepConvert")
if (!file.exists(blat)) stop("BLAT executable not found inside BrepConvert")
Sys.chmod(blat, "0755")

seqs <- readDNAStringSet(transcripts_fa)
repertoire <- setNames(as.character(seqs), names(seqs))
cat(sprintf("BrepConvert: %d transcripts, functional=%s, pseudogenes=%s\n",
            length(repertoire), functional_fa, pseudogene_fa))

# Mallaby et al. ran in batches of 3000; keep that behaviour.
batch_size <- 3000
batches <- split(seq_along(repertoire),
                 ceiling(seq_along(repertoire) / batch_size))

results <- list()
for (i in seq_along(batches)) {
  idx <- batches[[i]]
  cat(sprintf("  batch %d/%d (%d sequences)\n", i, length(batches), length(idx)))
  res <- try(batchConvertAnalysis(
    functional = functional_fa,
    pseudogene = pseudogene_fa,
    repertoire = repertoire[idx],
    blat_exec  = blat
  ), silent = TRUE)
  if (inherits(res, "try-error")) {
    cat("    batch failed:", attr(res, "condition")$message, "\n")
    next
  }
  if (!is.null(res) && is.data.frame(res) && nrow(res) > 0) {
    results[[length(results) + 1]] <- res
  }
}

if (length(results) == 0) {
  cat("BrepConvert returned no gene conversion calls\n")
  # still emit an empty file so downstream steps have something to read
  write.csv(data.frame(), out_csv, row.names = FALSE)
} else {
  final <- do.call(rbind, results)
  write.csv(final, out_csv, row.names = FALSE)
  cat(sprintf("wrote %d rows to %s\n", nrow(final), out_csv))
}
