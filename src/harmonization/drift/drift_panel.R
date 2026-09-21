# =============================================================================
# src/harmonization/drift/drift_panel.R — reconcile 12/13/14-parameter panels
# =============================================================================
# Base R. Maps every file's $PnN/$PnS pairs onto the canonical marker names in
# config/mtb_drift.yaml so all later steps address channels by marker. Files
# whose panel lacks a required marker are flagged rather than dropped here.
# =============================================================================

#' Parameter table of one file from its keywords: index, channel ($PnN), marker ($PnS).
fcs_parameter_table <- function(kw) {
  idx <- as.integer(sub("^\\$P(\\d+)N$", "\\1",
                        grep("^\\$P\\d+N$", names(kw), value = TRUE)))
  idx <- sort(idx)
  data.frame(
    index   = idx,
    channel = unname(kw[paste0("$P", idx, "N")]),
    marker  = unname(ifelse(is.na(kw[paste0("$P", idx, "S")]), "",
                            kw[paste0("$P", idx, "S")])),
    stringsAsFactors = FALSE)
}

#' Assign canonical names. `marker_patterns` is the `markers` block of the config:
#' a named list of regex vectors. Match on $PnS first, then on $PnN.
assign_canonical <- function(ptab, marker_patterns) {
  ptab$canonical <- NA_character_
  for (canon in names(marker_patterns)) {
    pats <- marker_patterns[[canon]]
    hit <- rep(FALSE, nrow(ptab))
    for (p in pats) {
      hit <- hit | grepl(p, ptab$marker, ignore.case = TRUE, perl = TRUE) |
                   grepl(p, ptab$channel, ignore.case = TRUE, perl = TRUE)
    }
    hit <- hit & is.na(ptab$canonical) & !grepl("-H$|-W$", ptab$channel)
    if (sum(hit) > 1L) {
      # prefer an area channel, then the first
      a <- hit & grepl("-A$", ptab$channel)
      hit <- if (any(a)) a & !duplicated(a) & a else hit & !duplicated(hit) & hit
      hit[which(hit)[-1]] <- FALSE
    }
    ptab$canonical[hit] <- canon
  }
  ptab
}

#' One long table over all files, plus a per-file panel signature.
build_panel_map <- function(tube_table, marker_patterns) {
  maps <- lapply(seq_len(nrow(tube_table)), function(i) {
    kw <- tryCatch(fcs_keywords(tube_table$path[i]), error = function(e) character(0))
    pt <- assign_canonical(fcs_parameter_table(kw), marker_patterns)
    pt$filename <- tube_table$filename[i]
    pt
  })
  pm <- do.call(rbind, maps)
  sig <- vapply(split(pm, pm$filename), function(d)
    paste(d$channel, collapse = ";"), character(1))
  list(map = pm,
       signature = data.frame(filename = names(sig), panel_signature = unname(sig),
                              stringsAsFactors = FALSE))
}

#' Files missing any of `required` canonical markers.
panel_gaps <- function(panel_map, required) {
  have <- tapply(panel_map$canonical, panel_map$filename,
                 function(x) required[!required %in% x])
  gaps <- Filter(length, have)
  data.frame(filename = names(gaps),
             missing = vapply(gaps, paste, character(1), collapse = ","),
             row.names = NULL, stringsAsFactors = FALSE)
}
