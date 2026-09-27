args <- commandArgs(TRUE)
src <- if (length(args) >= 1) args[1] else "raw/rat/MotrpacRatTraining6moData/data"
out <- if (length(args) >= 2) args[2] else "processed/rat_pkg"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
log <- file(file.path(out,"_summary.tsv"),"w")
writeLines("object\tclass\tnrow\tncol\tcolumns", log)
wt <- function(d, path){ for (c in names(d)) if (is.list(d[[c]])) d[[c]] <- vapply(d[[c]], function(x) paste(x, collapse=";"), "")
  con <- gzfile(path, "w"); write.table(d, con, sep="\t", quote=FALSE, row.names=FALSE, na=""); close(con) }
for (f in list.files(src, pattern="\\.rda$", full.names=TRUE)) {
  e <- new.env(); nm <- load(f, envir=e)
  for (k in nm) {
    v <- get(k, envir=e)
    if (is.data.frame(v)) {
      wt(v, file.path(out, paste0(k, ".tsv.gz")))
      writeLines(paste(k, "data.frame", nrow(v), ncol(v), paste(head(names(v),15), collapse=","), sep="\t"), log)
    } else if (is.list(v) && length(v)>0 && all(vapply(v, is.data.frame, TRUE))) {
      for (s in names(v)) { d <- v[[s]]; wt(d, file.path(out, paste0(k, "__", s, ".tsv.gz")))
        writeLines(paste(paste0(k,"__",s), "data.frame", nrow(d), ncol(d), paste(head(names(d),15), collapse=","), sep="\t"), log) }
    } else {
      saveRDS(v, file.path(out, paste0(k, ".rds")))
      writeLines(paste(k, paste(class(v),collapse="/"), length(v), NA, "", sep="\t"), log)
    }
  }
  flush(log)
}
close(log)
cat("DONE\n")
