# src/harmonization/drift/drift_audit.R — acquisition audit: settings that can make tubes differ
# Base R. Used by scripts/drift/11b_acquisition_audit.R. Collects, per tube, what the cytometer and FlowJo
# did to it (detector voltages, compensation, CS&T, FlowJo display transformation, sample quality) and how
# these settings overlap with the design factors, so the variance partition only models sources that can be
# told apart and reports the rest as explanations.

#' Acquisition settings per tube from a long keyword table (`kw`: filename, key, value; keys upper case) and
#' the panel map (`pm`: filename, channel, canonical). Voltages are reported per canonical marker
#' ("V_<marker>"), so cohorts with different channel names line up. Returns one row per file.
acquisition_settings <- function(kw, pm) {
  get <- function(k) { x <- kw[kw$key == k, ]; setNames(x$value, x$filename) }
  files <- unique(kw$filename)
  pick <- function(k) { v <- get(k); unname(v[files]) }
  out <- data.frame(filename = files, cytometer = pick("$CYT"), serial = pick("CYTNUM"), software = pick("CREATOR"),
                    cst_setup = pick("CST SETUP DATE"), cst_beads_lot = pick("CST BEADS LOT ID"),
                    cst_expired = !is.na(pick("CST PERFORMANCE EXPIRED")), config = pick("CYTOMETER CONFIG NAME"),
                    experiment = pick("EXPERIMENT NAME"), stringsAsFactors = FALSE)
  sp <- pick("SPILL"); out$compensation_id <- ifelse(is.na(sp), NA, as.integer(factor(sp)))
  pn <- kw[grepl("^\\$P[0-9]+N$", kw$key), ]; pv <- kw[grepl("^\\$P[0-9]+V$", kw$key), ]
  pn$n <- sub("^\\$P([0-9]+)N$", "\\1", pn$key); pv$n <- sub("^\\$P([0-9]+)V$", "\\1", pv$key)
  vv <- merge(pn[, c("filename", "n", "value")], pv[, c("filename", "n", "value")], by = c("filename", "n"), suffixes = c("_channel", "_voltage"))
  vv <- merge(vv, pm[!is.na(pm$canonical), c("filename", "channel", "canonical")], by.x = c("filename", "value_channel"), by.y = c("filename", "channel"))
  if (nrow(vv)) {
    vw <- reshape(data.frame(filename = vv$filename, marker = vv$canonical, v = suppressWarnings(as.numeric(vv$value_voltage))),
                  idvar = "filename", timevar = "marker", direction = "wide")
    names(vw) <- sub("^v\\.", "V_", names(vw))
    out <- merge(out, vw, by = "filename", all.x = TRUE)
    vc <- grep("^V_", names(out), value = TRUE)
    out$voltage_set <- as.integer(factor(do.call(paste, c(out[vc], sep = "/"))))
  }
  out
}

#' Sample quality per tube from the FlowJo gate counts (all_gate_counts.csv: filename, Population, Parent,
#' Count, ParentCount): viability = cells in the live-cell population / cells in its parent gate, and the
#' number of live cells. `live_population` is matched as the last path element (cfg$neg_population).
sample_quality <- function(counts, live_population) {
  leaf <- sub(".*/", "", counts$Population)
  x <- counts[tolower(leaf) == tolower(sub(".*/", "", live_population)), ]
  x <- x[!duplicated(x$filename), ]
  data.frame(filename = x$filename, viability = x$Count / x$ParentCount, live_cells = x$Count, stringsAsFactors = FALSE)
}

#' FlowJo display transformations per sample and parameter, read from the workspace XML (FlowJo 10:
#' <transforms:biex|linear|log|arcsinh ...> inside each <Sample>, with the parameter in data-type:name).
#' Returns workspace, filename, parameter, transform and its attributes as "key=value;..." (one row each).
flowjo_transformations <- function(wsp_path) {
  s <- paste(readLines(wsp_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  samples <- regmatches(s, gregexpr("(?s)<Sample[ >].*?</Sample>", s, perl = TRUE))[[1]]
  rows <- lapply(samples, function(blk) {
    uri <- regmatches(blk, regexpr('<DataSet[^>]*uri="[^"]+"', blk))
    fn <- if (length(uri)) utils::URLdecode(basename(sub('.*uri="([^"]+)"', "\\1", uri))) else NA_character_
    tr <- regmatches(blk, gregexpr("(?s)<transforms:(biex|linear|log|fasinh|arcsinh|logicle|hyperlog)[^>]*>.*?</transforms:\\1>", blk, perl = TRUE))[[1]]
    if (!length(tr)) return(NULL)
    do.call(rbind, lapply(tr, function(t) {
      type <- sub("^<transforms:([a-z]+).*$", "\\1", t)
      head <- sub("(?s)>.*$", "", t, perl = TRUE)
      at <- regmatches(head, gregexpr('transforms:[A-Za-z]+="[^"]*"', head))[[1]]
      at <- paste(sub('^transforms:', "", gsub('"', "", at)), collapse = ";")
      par <- regmatches(t, regexpr('data-type:name="[^"]+"', t)); par <- if (length(par)) sub('.*="([^"]+)"', "\\1", par) else NA_character_
      data.frame(workspace = basename(wsp_path), filename = fn, parameter = sub("^Comp-", "", par), transform = type,
                 settings = at, stringsAsFactors = FALSE)
    }))
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) return(data.frame(workspace = character(0), filename = character(0), parameter = character(0), transform = character(0), settings = character(0)))
  unique(out)
}

#' How much knowing factor A determines factor B, for every pair of columns `cols` of `d`:
#' 1 - H(B | A) / H(B) (entropy based; 1 = A fully determines B, 0 = independent). Rows are A, columns B.
determination_matrix <- function(d, cols) {
  H <- function(x) { p <- table(x) / length(x); -sum(p * log(p)) }
  Hc <- function(y, x) { s <- split(y, x); sum(vapply(s, function(v) length(v) * H(v), numeric(1))) / length(y) }
  M <- matrix(NA_real_, length(cols), length(cols), dimnames = list(names(cols), names(cols)))
  for (a in names(cols)) for (b in names(cols)) {
    x <- as.character(d[[cols[[a]]]]); y <- as.character(d[[cols[[b]]]]); ok <- !is.na(x) & !is.na(y)
    hy <- H(y[ok]); M[a, b] <- if (a == b) 1 else if (hy == 0) NA_real_ else 1 - Hc(y[ok], x[ok]) / hy
  }
  M
}

#' Share of a numeric variable explained by a factor alone (eta squared).
eta_squared <- function(y, g) {
  ok <- is.finite(y) & !is.na(g); y <- y[ok]; g <- as.character(g[ok])
  1 - sum(tapply(y, g, function(v) sum((v - mean(v))^2))) / sum((y - mean(y))^2)
}
