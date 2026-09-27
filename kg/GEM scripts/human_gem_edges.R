#!/usr/bin/env Rscript
# How the knowledge graph's gene <-> metabolite edges are built from Human-GEM.
#
# Human-GEM (SysBioChalmers, https://github.com/SysBioChalmers/Human-GEM) is a curated genome-scale model of
# human metabolism. It has no gene-metabolite table; the link runs through reactions:
#
#     gene --(reaction gene rule)--> reaction --(stoichiometric coefficient)--> metabolite
#
#   * Each reaction carries a Boolean gene rule over Ensembl ids, e.g. "ENSG_A or (ENSG_B and ENSG_C)":
#     "or" joins isozymes, "and" joins subunits of one complex.
#   * Each reaction lists its metabolites with a stoichiometric coefficient: negative = consumed (substrate),
#     positive = produced (product). A reaction with lower bound < 0 is reversible.
#   * Metabolite ids are compartment-specific (MAM02403c = L-lactate in the cytosol, MAM02403e = the same
#     molecule outside the cell). The graph collapses compartments: MAM02403.
#
# Edge rule. A gene G and a metabolite M are linked when at least one reaction
#   (1) has G anywhere in its gene rule (isozyme or complex subunit alike), and
#   (2) has M, in any compartment, among its metabolites, and
#   (3) M is not a "currency" metabolite: it appears in <= 60 reactions (ATP, H2O, H+, NAD(H), CoA, ... would
#       otherwise link to almost every gene).
# Edge attributes, aggregated over all reactions that satisfy (1)-(3):
#   n_reactions   number of such reactions
#   role          substrate (M consumed), product (M produced), or either (reversible reaction, or different
#                 roles in different reactions)
#   is_transport  M crosses a compartment boundary in at least one of the reactions (it appears in two
#                 compartments of the same reaction), or the reaction is an exchange with the extracellular
#                 space and its subsystem is a transport subsystem
#   reactions     the reaction ids (first 20, sorted)
#   subsystems    the reactions' subsystems (first 6, sorted)
#
# Usage
#   Rscript human_gem_edges.R download                  # fetch Human-GEM v2.0.1 (~12 MB) into ./data/human_gem
#   Rscript human_gem_edges.R pairs --genes LDHA,SLC16A1,HK2 --metabolites "L-lactate;pyruvate;ATP"
#   Rscript human_gem_edges.R pairs --genes genes.txt --metabolites mets.txt --out-prefix results/query --verbose
#   Rscript human_gem_edges.R build                     # every edge -> ./data/human_gem/gene_metabolite_edges.csv
#
# --genes / --metabolites take either a file (one item per line; a tab-separated file uses its first column;
# blank lines and lines starting with # are skipped) or an inline list. An inline list is split on ";" when it
# contains one, otherwise on ",": use ";" or a file for metabolite names that contain commas
# ("2,3-bisphospho-D-glycerate"). Genes: HGNC symbol or Ensembl id. Metabolites: Human-GEM name (exact,
# case-insensitive), MAM id (with or without compartment letter), KEGG (C00186), HMDB (HMDB0000190 or
# HMDB00190), ChEBI (CHEBI:422) or PubChem CID (bare digits).
#
# `pairs` tests every gene x metabolite combination and writes
#   <prefix>_pairs.tsv     one row per pair: edge TRUE/FALSE, its attributes, and the reason when there is none
#   <prefix>_evidence.tsv  one row per pair x shared reaction: equation, coefficients, reversibility, whether the
#                          gene is an isozyme or complex subunit, and what the reaction contributes to the edge
#   <prefix>_matrix.tsv    genes x metabolites, n_reactions of the edge (0 = no edge, empty = input not resolved)
#
# `pairs` and `build` download the model first if it is missing; the parsed model is cached as
# human_gem_parsed.rds next to it (the first run parses the 8 MB YAML, ~30 s). Requires R >= 4.1 and the
# `yaml` package. Functions can also be used from an R session: source("human_gem_edges.R");
# gem <- load_gem(); res <- query_pairs(gem, c("LDHA", "HK2"), c("L-lactate", "ATP")).

GEM_VERSION <- "v2.0.1"
RAW_URL <- "https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/%s/model/%s"
FILES <- c("Human-GEM.yml", "metabolites.tsv", "genes.tsv")
CURRENCY_DEGREE <- 60        # metabolites in more reactions than this get no gene edges
MAX_LISTED_REACTIONS <- 20
MAX_LISTED_SUBSYSTEMS <- 6

`%||%` <- function(a, b) if (is.null(a)) b else a

script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
}
DEFAULT_DIR <- file.path(script_dir(), "data", "human_gem")


# ---------------------------------------------------------------- step 0: download

download_gem <- function(gem_dir = DEFAULT_DIR, version = GEM_VERSION, force = FALSE) {
  # the model YAML (reactions, metabolites, gene rules) and the metabolite / gene annotation tables
  dir.create(gem_dir, recursive = TRUE, showWarnings = FALSE)
  old <- options(timeout = max(600, getOption("timeout"))); on.exit(options(old))
  for (f in FILES) {
    dest <- file.path(gem_dir, f)
    if (file.exists(dest) && !force) { message("  have ", dest); next }
    url <- sprintf(RAW_URL, version, f)
    message("  get  ", url)
    tmp <- paste0(dest, ".part")
    download.file(url, tmp, mode = "wb", quiet = TRUE)
    file.rename(tmp, dest)
  }
}

ensure_downloaded <- function(gem_dir) {
  if (!all(file.exists(file.path(gem_dir, FILES)))) {
    message("Human-GEM not found in ", gem_dir, "; downloading ", GEM_VERSION)
    download_gem(gem_dir)
  }
}


# ---------------------------------------------------------------- step 1: load the model

base_id <- function(met) sub("[a-z]$", "", met)   # MAM02403c -> MAM02403 (one-letter compartment suffix)

first_by <- function(values, keys) {                # named vector key -> value, first occurrence wins
  keep <- !duplicated(keys) & !is.na(keys) & keys != ""
  setNames(values[keep], keys[keep])
}

load_gem <- function(gem_dir = DEFAULT_DIR) {
  ensure_downloaded(gem_dir)
  cache <- file.path(gem_dir, "human_gem_parsed.rds")
  if (file.exists(cache) && file.mtime(cache) >= max(file.mtime(file.path(gem_dir, FILES))))
    return(readRDS(cache))

  yml <- file.path(gem_dir, "Human-GEM.yml")
  message("parsing ", yml, " (first run only)")
  # YAML 1.1 reads a bare NO as boolean FALSE, and NO is nitric oxide's name in Human-GEM: keep such words as text
  keep_text <- function(x) x
  model <- yaml::read_yaml(yml, handlers = list("bool#yes" = keep_text, "bool#no" = keep_text))

  # metabolites: compartment-specific id -> name, compartment
  mets <- model$metabolites
  met_id <- vapply(mets, function(m) m$id, "")
  met_name <- setNames(vapply(mets, function(m) as.character(m$name %||% ""), ""), met_id)

  # reactions: id, name, gene rule, subsystems, reversibility; and the stoichiometry as a long table
  rxns <- model$reactions
  rid <- vapply(rxns, function(r) r$id, "")
  rx <- data.frame(
    rxn = rid,
    name = vapply(rxns, function(r) as.character(r$name %||% ""), ""),
    rule = vapply(rxns, function(r) trimws(as.character(r$gene_reaction_rule %||% "")), ""),
    subsystem = vapply(rxns, function(r) paste(unlist(r$subsystem), collapse = ";"), ""),
    reversible = vapply(rxns, function(r) as.numeric(r$lower_bound %||% 0) < 0, TRUE),
    stringsAsFactors = FALSE)
  st <- data.frame(
    rxn = rep(rid, vapply(rxns, function(r) length(r$metabolites), 0L)),
    met = unlist(lapply(rxns, function(r) names(r$metabolites)), use.names = FALSE),
    coef = as.numeric(unlist(lapply(rxns, function(r) unlist(r$metabolites)), use.names = FALSE)),
    stringsAsFactors = FALSE)
  st$base <- base_id(st$met)
  st$comp <- substring(st$met, nchar(st$met))

  # reaction -> genes named in its rule
  rule_genes <- lapply(regmatches(rx$rule, gregexpr("ENSG[0-9]+", rx$rule)), unique)
  rg <- data.frame(rxn = rep(rx$rxn, lengths(rule_genes)), ensembl = unlist(rule_genes), stringsAsFactors = FALSE)

  # genes: Ensembl id <-> symbol
  gt <- read.delim(file.path(gem_dir, "genes.tsv"), colClasses = "character", quote = "\"", na.strings = character())

  # metabolite cross-references (compartment-free id), used to resolve user input
  mt <- read.delim(file.path(gem_dir, "metabolites.tsv"), colClasses = "character", quote = "", na.strings = character())
  xr <- do.call(rbind, lapply(c("metKEGGID", "metHMDBID", "metChEBIID", "metPubChemID"), function(col) {
    v <- strsplit(mt[[col]], ";", fixed = TRUE)
    data.frame(key = toupper(trimws(unlist(v))), base = rep(mt$metsNoComp, lengths(v)), stringsAsFactors = FALSE)
  }))

  # reaction degree of each compartment-free metabolite: the currency test
  ub <- unique(st[c("base", "rxn")])
  degree <- c(table(ub$base))

  base_name <- first_by(unname(met_name), base_id(met_id))
  gem <- list(
    version = model$metaData$version,
    rx = rx, st = st, rg = rg,
    st_by_rxn = split(seq_len(nrow(st)), st$rxn),
    met_name = met_name, base_name = base_name,
    base_of_name = first_by(names(base_name), tolower(base_name)),
    comps_of = tapply(st$comp, st$base, function(x) paste(sort(unique(x)), collapse = ",")),
    symbol = setNames(gt$geneSymbols, gt$genes),
    ens_of_symbol = first_by(gt$genes, toupper(gt$geneSymbols)),
    xref = first_by(xr$base, xr$key),
    degree = degree)
  saveRDS(gem, cache)
  gem
}


# ---------------------------------------------------------------- step 2: reactions -> (metabolite, gene) contributions

contributions <- function(gem, ensembl = NULL) {
  # one row per (reaction, compartment-free metabolite, rule gene) for non-currency metabolites,
  # with the metabolite's role in that reaction and whether the reaction moves it between compartments
  rg <- gem$rg
  if (!is.null(ensembl)) rg <- rg[rg$ensembl %in% ensembl, ]
  st <- gem$st[gem$st$rxn %in% rg$rxn, ]
  st <- st[gem$degree[st$base] <= CURRENCY_DEGREE, ]
  if (!nrow(st)) return(data.frame(rxn = character(), base = character(), role = character(),
                                   transport = logical(), ensembl = character()))
  rxi <- match(st$rxn, gem$rx$rxn)
  st$role <- ifelse(gem$rx$reversible[rxi], "either", ifelse(st$coef < 0, "substrate", "product"))
  n_comp_rxn <- tapply(st$comp, st$rxn, function(x) length(unique(x)))
  has_e <- tapply(st$comp, st$rxn, function(x) any(x == "e"))
  key <- paste(st$rxn, st$base)
  crosses <- tapply(st$comp, key, function(x) length(unique(x)) > 1)[key]
  exchange <- has_e[st$rxn] & n_comp_rxn[st$rxn] > 1 & grepl("ransport", gem$rx$subsystem[rxi])
  st$transport <- as.vector(crosses | exchange)

  first <- !duplicated(key)
  rb <- data.frame(rxn = st$rxn[first], base = st$base[first], stringsAsFactors = FALSE)
  rk <- key[first]
  rb$role <- as.vector(tapply(st$role, key, function(x) if (length(unique(x)) > 1) "either" else x[1])[rk])
  rb$transport <- as.vector(tapply(st$transport, key, any)[rk])
  merge(rb, rg, by = "rxn")
}


# ---------------------------------------------------------------- step 3: aggregate over reactions -> edges

aggregate_edges <- function(gem, cb) {
  if (!nrow(cb)) return(data.frame(met_id = character(), met_name = character(), ensembl = character(),
                                   symbol = character(), n_reactions = integer(), role = character(),
                                   is_transport = logical(), reactions = character(), subsystems = character()))
  cb <- cb[order(cb$base, cb$ensembl, method = "radix"), ]
  k <- paste(cb$base, cb$ensembl)
  grp <- split(seq_len(nrow(cb)), factor(k, levels = unique(k)))
  first <- vapply(grp, function(i) i[1], 0L)
  sub_of <- setNames(gem$rx$subsystem, gem$rx$rxn)
  sym <- unname(gem$symbol[cb$ensembl[first]])
  data.frame(
    met_id = cb$base[first],
    met_name = unname(gem$base_name[cb$base[first]]),
    ensembl = cb$ensembl[first],
    symbol = ifelse(is.na(sym), "", sym),
    n_reactions = vapply(grp, function(i) length(unique(cb$rxn[i])), 0L),
    role = vapply(grp, function(i) { r <- unique(cb$role[i]); if (length(r) > 1) "either" else r }, ""),
    is_transport = vapply(grp, function(i) any(cb$transport[i]), TRUE),
    reactions = vapply(grp, function(i)
      paste(head(sort(unique(cb$rxn[i]), method = "radix"), MAX_LISTED_REACTIONS), collapse = ";"), ""),
    subsystems = vapply(grp, function(i) {
      s <- unlist(strsplit(sub_of[unique(cb$rxn[i])], ";", fixed = TRUE))
      paste(head(sort(unique(s[s != ""]), method = "radix"), MAX_LISTED_SUBSYSTEMS), collapse = ";")
    }, ""),
    row.names = NULL, stringsAsFactors = FALSE)
}

build_edges <- function(gem) aggregate_edges(gem, contributions(gem))


# ---------------------------------------------------------------- input resolution

resolve_genes <- function(gem, q) {
  up <- toupper(trimws(q))
  out <- ifelse(grepl("^ENSG[0-9]+$", up), up, NA_character_)
  i <- is.na(out); out[i] <- gem$ens_of_symbol[up[i]]
  unname(out)
}

resolve_metabolites <- function(gem, q) {
  s <- trimws(q)
  out <- ifelse(grepl("^MAM[0-9]+[a-z]?$", s), base_id(s), NA_character_)
  up <- sub("^HMDB([0-9]{5})$", "HMDB00\\1", toupper(s))          # old 5-digit HMDB ids
  i <- is.na(out); out[i] <- gem$xref[up[i]]
  i <- is.na(out); out[i] <- gem$base_of_name[tolower(s[i])]
  out[!is.na(out) & !(out %in% names(gem$base_name))] <- NA_character_
  unname(out)
}

name_suggestions <- function(gem, q, k = 3) {
  n <- unique(unname(gem$base_name[grepl(tolower(trimws(q)), tolower(gem$base_name), fixed = TRUE)]))
  head(n[order(nchar(n))], k)
}


# ---------------------------------------------------------------- per-reaction evidence

equation <- function(gem, rxn) {
  s <- gem$st[gem$st_by_rxn[[rxn]], ]
  term <- function(j) {
    n <- paste0(gem$met_name[s$met[j]], "[", s$comp[j], "]")
    ifelse(abs(s$coef[j]) == 1, n, paste(abs(s$coef[j]), n))
  }
  arrow <- if (gem$rx$reversible[match(rxn, gem$rx$rxn)]) "<=>" else "=>"
  paste(paste(term(which(s$coef < 0)), collapse = " + "), arrow, paste(term(which(s$coef > 0)), collapse = " + "))
}

gene_context <- function(rule, g) {
  # complex subunit if an "and" joins the terms at any bracket level enclosing the gene; else isozyme
  toks <- regmatches(rule, gregexpr("\\(|\\)|\\band\\b|\\bor\\b|ENSG[0-9]+", rule))[[1]]
  if (length(unique(grep("^ENSG", toks, value = TRUE))) == 1) return("sole gene")
  d <- cumsum(toks == "(") - cumsum(toks == ")")
  i <- match(g, toks)
  for (L in 0:d[i]) {
    lo <- i; while (lo > 1 && d[lo - 1] >= L) lo <- lo - 1
    hi <- i; while (hi < length(toks) && d[hi + 1] >= L) hi <- hi + 1
    span <- lo:hi
    if (any(toks[span] == "and" & d[span] == L)) return("complex subunit (AND)")
  }
  "isozyme (OR)"
}


# ---------------------------------------------------------------- gene list x metabolite list

query_pairs <- function(gem, genes, metabolites) {
  genes <- unique(trimws(genes)); metabolites <- unique(trimws(metabolites))
  g_ens <- resolve_genes(gem, genes)
  m_id <- resolve_metabolites(gem, metabolites)

  # edges: the same contribution + aggregation code as `build`, restricted to the queried genes and metabolites
  cb <- contributions(gem, unique(na.omit(g_ens)))
  cb <- cb[cb$base %in% m_id, ]
  edges <- aggregate_edges(gem, cb)

  # all reactions shared by a queried gene and metabolite, including those that do not count (currency)
  rg_q <- gem$rg[gem$rg$ensembl %in% g_ens, ]
  st_q <- unique(gem$st[gem$st$base %in% m_id, c("rxn", "base")])
  sh <- merge(rg_q, st_q, by = "rxn")
  n_shared <- table(factor(paste(sh$ensembl, sh$base)))

  gi <- rep(seq_along(genes), each = length(metabolites))
  mi <- rep(seq_along(metabolites), times = length(genes))
  p <- data.frame(gene_query = genes[gi], ensembl = g_ens[gi], metabolite_query = metabolites[mi], met_id = m_id[mi],
                  stringsAsFactors = FALSE)
  sym <- unname(gem$symbol[p$ensembl]); p$symbol <- ifelse(is.na(sym), "", sym)
  p$met_name <- ifelse(is.na(p$met_id), "", unname(gem$base_name[p$met_id]))
  ei <- match(paste(p$met_id, p$ensembl), paste(edges$met_id, edges$ensembl))
  p$edge <- !is.na(ei)
  for (col in c("n_reactions", "role", "is_transport", "reactions", "subsystems")) p[[col]] <- edges[[col]][ei]
  p$n_reactions[!p$edge & !is.na(p$ensembl) & !is.na(p$met_id)] <- 0L
  ns <- n_shared[paste(p$ensembl, p$met_id)]
  p$n_shared_reactions <- ifelse(is.na(ns), 0L, as.integer(ns))
  p$met_degree <- as.integer(gem$degree[p$met_id])

  suggest <- vapply(metabolites, function(q) {
    s <- name_suggestions(gem, q)
    if (length(s)) paste0("metabolite not matched; names containing it: ", paste(s, collapse = "; ")) else "metabolite not matched"
  }, "", USE.NAMES = FALSE)
  p$reason <- ifelse(p$edge, "edge",
    ifelse(is.na(p$ensembl), "gene not recognised",
    ifelse(is.na(p$met_id), suggest[mi],
    ifelse(!p$ensembl %in% gem$rg$ensembl, "gene in no Human-GEM reaction",
    ifelse(p$met_degree > CURRENCY_DEGREE, sprintf("currency metabolite (%d reactions > %d)", p$met_degree, CURRENCY_DEGREE),
    "no shared reaction")))))

  # display labels: the resolved name, with the query in brackets when it was an id or differed from the name
  label <- function(q, name) ifelse(is.na(name) | name == "", q, ifelse(tolower(q) == tolower(name), name, paste0(name, " (", q, ")")))
  p$gene_label <- label(p$gene_query, p$symbol)
  p$metabolite_label <- label(p$metabolite_query, p$met_name)

  # evidence: one row per resolved (gene, metabolite) x shared reaction
  ck <- match(paste(sh$rxn, sh$base, sh$ensembl), paste(cb$rxn, cb$base, cb$ensembl))
  rxi <- match(sh$rxn, gem$rx$rxn)
  ev <- data.frame(
    symbol = unname(gem$symbol[sh$ensembl]), ensembl = sh$ensembl,
    met_id = sh$base, met_name = unname(gem$base_name[sh$base]),
    rxn = sh$rxn, rxn_name = gem$rx$name[rxi], subsystem = gem$rx$subsystem[rxi],
    equation = vapply(sh$rxn, function(r) equation(gem, r), "", USE.NAMES = FALSE),
    coefficients = mapply(function(r, b) {
      s <- gem$st[gem$st_by_rxn[[r]], ]; s <- s[s$base == b, ]
      paste(sprintf("%s %+g", s$met, s$coef), collapse = ", ")
    }, sh$rxn, sh$base, USE.NAMES = FALSE),
    reversible = gem$rx$reversible[rxi],
    gene_context = mapply(gene_context, gem$rx$rule[rxi], sh$ensembl, USE.NAMES = FALSE),
    contributes = !is.na(ck),
    role = cb$role[ck], is_transport = cb$transport[ck],
    stringsAsFactors = FALSE)
  ev <- ev[order(ev$symbol, ev$met_name, ev$rxn, method = "radix"), ]

  # matrix: n_reactions of the edge; 0 = no edge; NA = gene or metabolite not resolved
  mat <- matrix(p$n_reactions, nrow = length(genes), byrow = TRUE,
                dimnames = list(p$gene_label[mi == 1], p$metabolite_label[gi == 1]))

  list(pairs = p, evidence = ev, matrix = mat)
}


# ---------------------------------------------------------------- CLI

read_list <- function(x) {
  if (file.exists(x)) {
    v <- sub("\t.*$", "", readLines(x, warn = FALSE))
    v <- trimws(v)
    v <- v[v != "" & !startsWith(v, "#")]
  } else {
    v <- trimws(strsplit(x, if (grepl(";", x, fixed = TRUE)) ";" else ",", fixed = TRUE)[[1]])
    v <- v[v != ""]
  }
  unique(v)
}

parse_opts <- function(args) {
  opt <- list(); i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (!startsWith(a, "--")) stop("unexpected argument: ", a)
    k <- substring(a, 3)
    if (i < length(args) && !startsWith(args[i + 1], "--")) { opt[[k]] <- args[i + 1]; i <- i + 2 }
    else { opt[[k]] <- TRUE; i <- i + 1 }
  }
  opt
}

write_tsv <- function(df, path, row_names = FALSE) {
  write.table(df, path, sep = "\t", quote = FALSE, row.names = row_names, col.names = if (row_names) NA else TRUE, na = "")
}

print_verbose <- function(res) {
  p <- res$pairs; ev <- res$evidence
  for (j in which(p$n_shared_reactions > 0)) {
    e <- ev[ev$ensembl == p$ensembl[j] & ev$met_id == p$met_id[j], ]
    cat(sprintf("\n%s -- %s (%s): %s\n", p$symbol[j], p$met_name[j], p$met_id[j],
                if (p$edge[j]) sprintf("EDGE, %d reactions, role %s, transport %s", p$n_reactions[j], p$role[j], p$is_transport[j]) else p$reason[j]))
    for (r in seq_len(nrow(e))) {
      cat(sprintf("  %s  %s  [%s]\n      %s\n      %s; %s; gene is %s -> %s\n", e$rxn[r], e$rxn_name[r], e$subsystem[r],
                  e$equation[r], e$coefficients[r], if (e$reversible[r]) "reversible" else "irreversible", e$gene_context[r],
                  if (e$contributes[r]) sprintf("role %s, transport %s", e$role[r], e$is_transport[r]) else "nothing (currency metabolite)"))
    }
  }
}

USAGE <- "usage:
  Rscript human_gem_edges.R download [--version v2.0.1] [--force]
  Rscript human_gem_edges.R pairs --genes <file|list> --metabolites <file|list> [--out-prefix gem_query] [--verbose]
  Rscript human_gem_edges.R build [--out <gem-dir>/gene_metabolite_edges.csv]
common: --gem-dir <dir> (default: data/human_gem next to this script)
See the header of this script for input formats and outputs.
"

main <- function(args = commandArgs(TRUE)) {
  if (!length(args) || args[1] %in% c("-h", "--help", "help")) { cat(USAGE); return(0L) }
  cmd <- args[1]; opt <- parse_opts(args[-1])
  gem_dir <- opt[["gem-dir"]] %||% DEFAULT_DIR

  if (cmd == "download") {
    download_gem(gem_dir, opt$version %||% GEM_VERSION, isTRUE(opt$force)); return(0L)
  }
  if (!cmd %in% c("pairs", "build")) { cat(USAGE); return(1L) }
  gem <- load_gem(gem_dir)

  if (cmd == "build") {
    edges <- build_edges(gem)
    out <- opt$out %||% file.path(gem_dir, "gene_metabolite_edges.csv")
    write.csv(edges, out, row.names = FALSE)
    cat(sprintf("Human-GEM %s: %d gene-metabolite edges (%d metabolites x %d genes) -> %s\n", gem$version, nrow(edges),
                length(unique(edges$met_id)), length(unique(edges$ensembl)), out))
    return(0L)
  }

  if (is.null(opt$genes) || is.null(opt$metabolites) || isTRUE(opt$genes) || isTRUE(opt$metabolites)) {
    cat("pairs needs --genes and --metabolites\n\n", USAGE, sep = ""); return(1L)
  }
  res <- query_pairs(gem, read_list(opt$genes), read_list(opt$metabolites))
  p <- res$pairs
  cat(sprintf("Human-GEM %s: %d genes x %d metabolites = %d pairs, %d edges\n\n", gem$version,
              nrow(res$matrix), ncol(res$matrix), nrow(p), sum(p$edge)))
  op <- options(width = 200); on.exit(options(op))
  unres <- unique(p$reason[p$reason %in% c("gene not recognised", "gene in no Human-GEM reaction") | startsWith(p$reason, "metabolite not matched")])
  for (r in unres) {
    who <- if (startsWith(r, "metabolite")) unique(p$metabolite_query[p$reason == r]) else unique(p$gene_label[p$reason == r])
    cat(sprintf("unresolved: %s: %s\n", paste(who, collapse = ", "), r))
  }
  e <- p[p$edge, ]
  if (nrow(e)) {
    cat("\nedges\n")
    print(data.frame(gene = e$gene_label, metabolite = e$metabolite_label, n_rxn = e$n_reactions, role = e$role,
                     transport = e$is_transport, subsystems = e$subsystems), row.names = FALSE, right = FALSE)
  }
  ne <- p[!p$edge & !p$reason %in% unres, ]
  if (nrow(ne)) {
    cat("\nno edge (resolved pairs)\n")
    tab <- table(ne$reason); for (r in names(tab)) cat(sprintf("  %4d  %s\n", tab[[r]], r))
  }
  if (nrow(res$matrix) <= 30 && ncol(res$matrix) <= 12) {
    cat("\nn_reactions per edge (0 = no edge, NA = not resolved)\n"); print(res$matrix)
  }
  if (isTRUE(opt$verbose)) print_verbose(res)

  prefix <- opt[["out-prefix"]] %||% "gem_query"
  if (dirname(prefix) != ".") dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  write_tsv(p, paste0(prefix, "_pairs.tsv"))
  write_tsv(res$evidence, paste0(prefix, "_evidence.tsv"))
  write_tsv(res$matrix, paste0(prefix, "_matrix.tsv"), row_names = TRUE)
  cat(sprintf("\nwrote %s_pairs.tsv, %s_evidence.tsv, %s_matrix.tsv\n", prefix, prefix, prefix))
  0L
}

if (sys.nframe() == 0L) quit(status = main())
