src <- "/home/claude/ekg/raw/human/MotrpacHumanPreSuspensionAnalysis/data"
out <- "/home/claude/ekg/processed/human_pkg"
log <- file(file.path(out,"_summary.tsv"),"w")
writeLines("object\tclass\tnrow\tncol\tcolumns", log)
wt <- function(d, path){ for (c in names(d)) if (is.list(d[[c]])) d[[c]] <- vapply(d[[c]], function(x) paste(x, collapse=";"), "")
  con <- gzfile(path, "w"); write.table(d, con, sep="\t", quote=FALSE, row.names=FALSE, na=""); close(con) }
for (f in setdiff(list.files(src, pattern="\\.rda$", full.names=TRUE), file.path(src,"PTMSEA_INPUT.rda"))) {
  e <- new.env(); nm <- tryCatch(load(f, envir=e), error=function(err){writeLines(paste(basename(f),"LOAD_ERROR",NA,NA,conditionMessage(err),sep="\t"),log); character(0)})
  for (k in nm) {
    v <- tryCatch(get(k, envir=e), error=function(err) NULL); if (is.null(v)) next
    if (is.data.frame(v)) {
      wt(v, file.path(out, paste0(k, ".tsv.gz")))
      writeLines(paste(k, "data.frame", nrow(v), ncol(v), paste(head(names(v),15), collapse=","), sep="\t"), log)
    } else if (is.list(v) && length(v)>0 && all(vapply(v, is.data.frame, TRUE))) {
      for (s in names(v)) { d <- v[[s]]; wt(d, file.path(out, paste0(k, "__", s, ".tsv.gz")))
        writeLines(paste(paste0(k,"__",s), "data.frame", nrow(d), ncol(d), paste(head(names(d),15), collapse=","), sep="\t"), log) }
    } else {
      tryCatch(saveRDS(v, file.path(out, paste0(k, ".rds"))), error=function(err) NULL)
      writeLines(paste(k, paste(class(v),collapse="/"), length(v), NA, "", sep="\t"), log)
    }
  }
  flush(log)
}
close(log)
cat("DONE\n")
