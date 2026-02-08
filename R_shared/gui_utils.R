# =============================================================================
# R_shared/gui_utils.R — GUI Utility Functions (tcltk with CLI fallback)
# =============================================================================
# Deduplicated from: Functions.R, 01_split_clusters.R, GLMM, Fine_clustering
# =============================================================================

#' Check if session supports interactive GUI
#' @return Logical TRUE if tcltk GUI is available
is_gui_available <- function() {
  interactive() && requireNamespace("tcltk", quietly = TRUE) &&
    tryCatch({ tcltk::tclRequire("Tk"); TRUE }, error = function(e) FALSE)
}

#' Multi-selection GUI dialog with CLI fallback
#' @param choices Character vector of items to choose from
#' @param title Dialog title
#' @param multiple Allow multiple selections
#' @param cli_fallback Use CLI when GUI unavailable (default TRUE)
#' @return Character vector of selected items, or NULL if cancelled
tk_multiselect <- function(choices, title = "Select items",
                           multiple = TRUE, cli_fallback = TRUE) {
  if (is_gui_available()) {
    tryCatch({
      return(.tk_multiselect_gui(choices, title, multiple))
    }, error = function(e) {
      message("GUI failed: ", e$message)
    })
  }

  if (cli_fallback) {
    return(.cli_multiselect(choices, title, multiple))
  }

  message("No interactive selection available, returning all choices")
  choices
}

# Internal: tcltk GUI multi-select
.tk_multiselect_gui <- function(choices, title, multiple) {
  require(tcltk)
  tt <- tktoplevel()
  tkwm.title(tt, title)
  tkwm.geometry(tt, "400x300")

  main_frame <- tkframe(tt)
  tkpack(main_frame, expand = TRUE, fill = "both")

  select_mode <- if (multiple) "extended" else "single"
  listbox <- tklistbox(main_frame, selectmode = select_mode, height = 10, width = 30)
  scrollbar <- tkscrollbar(main_frame, orient = "vertical")
  tkconfigure(listbox, yscrollcommand = function(...) tkset(scrollbar, ...))
  tkconfigure(scrollbar, command = function(...) tkyview(listbox, ...))
  tkpack(scrollbar, side = "right", fill = "y")
  tkpack(listbox, side = "left", expand = TRUE, fill = "both")
  for (item in choices) tkinsert(listbox, "end", item)

  selected_items <- NULL
  ok_function <- function() {
    idx <- as.integer(tkcurselection(listbox)) + 1
    selected_items <<- if (length(idx) > 0) choices[idx] else character(0)
    tkdestroy(tt)
  }
  cancel_function <- function() { selected_items <<- NULL; tkdestroy(tt) }

  btn_frame <- tkframe(main_frame)
  tkpack(btn_frame, side = "bottom", fill = "x")
  ok_btn <- tkbutton(btn_frame, text = "OK", command = ok_function)
  cancel_btn <- tkbutton(btn_frame, text = "Cancel", command = cancel_function)
  tkpack(ok_btn, side = "left", padx = 5, pady = 5)
  tkpack(cancel_btn, side = "left", padx = 5, pady = 5)

  tkbind(listbox, "<Double-Button-1>", function() ok_function())
  tkbind(tt, "<Return>", ok_function)
  tkbind(tt, "<Escape>", cancel_function)
  tkfocus(listbox)
  tkgrab.set(tt)
  tkwait.window(tt)
  selected_items
}

# Internal: CLI fallback for multi-select
.cli_multiselect <- function(choices, title, multiple) {
  cat("\n===", title, "===\n")
  for (i in seq_along(choices)) {
    cat(sprintf("  %d. %s\n", i, choices[i]))
  }
  if (multiple) {
    cat("Enter numbers (comma-separated), 'all', or 'q' to cancel: ")
  } else {
    cat("Enter number, or 'q' to cancel: ")
  }
  input <- readline()
  if (tolower(input) == "q") return(NULL)
  if (tolower(input) == "all") return(choices)
  nums <- as.integer(strsplit(trimws(input), "[,\\s]+")[[1]])
  nums <- nums[!is.na(nums) & nums >= 1 & nums <= length(choices)]
  if (length(nums) == 0) return(character(0))
  choices[nums]
}

#' Marker selection popup with checkbox scrollbar
#' @param marker_info List with 'indices' and 'names' elements
#' @param title Dialog title
#' @return Integer vector of selected marker indices, or NULL if cancelled
select_markers <- function(marker_info, title = "Select Markers") {
  if (!is_gui_available()) {
    return(.cli_select_markers(marker_info, title))
  }

  tryCatch({
    require(tcltk)
    tt <- tktoplevel()
    tkwm.title(tt, title)
    tkwm.geometry(tt, "400x350")

    main_frame <- tkframe(tt)
    tkpack(main_frame, expand = TRUE, fill = "both", padx = 10, pady = 10)

    label <- tklabel(main_frame, text = title,
                     font = tkfont.create(size = 12, weight = "bold"))
    tkpack(label, side = "top", pady = 5)

    cb_frame <- tkframe(main_frame)
    tkpack(cb_frame, expand = TRUE, fill = "both")

    # Create scrollable checkboxes
    canvas <- tkcanvas(cb_frame, bg = "white")
    scrollbar <- tkscrollbar(cb_frame, orient = "vertical",
                             command = function(...) tkyview(canvas, ...))
    tkconfigure(canvas, yscrollcommand = function(...) tkset(scrollbar, ...))
    tkpack(scrollbar, side = "right", fill = "y")
    tkpack(canvas, side = "left", expand = TRUE, fill = "both")

    frame <- tkframe(canvas, bg = "white")
    tkcreate(canvas, "window", 0, 0, anchor = "nw", window = frame)

    vars <- vector("list", length(marker_info$indices))
    for (i in seq_along(marker_info$indices)) {
      vars[[i]] <- tclVar("1")
      txt <- paste0("Column ", marker_info$indices[i], ": ", marker_info$names[i])
      cb <- tkcheckbutton(frame, text = txt, variable = vars[[i]],
                          anchor = "w", bg = "white")
      tkpack(cb, side = "top", fill = "x")
    }
    tkbind(frame, "<Configure>", function() {
      bbox <- tkbbox(canvas, "all")
      tkconfigure(canvas, scrollregion = bbox)
    })

    selected_markers <- NULL
    ok_function <- function() {
      selected_markers <<- marker_info$indices[
        which(sapply(vars, tclvalue) == "1")
      ]
      tkdestroy(tt)
    }
    cancel_function <- function() { selected_markers <<- NULL; tkdestroy(tt) }

    ok_btn <- tkbutton(main_frame, text = "OK", command = ok_function)
    cancel_btn <- tkbutton(main_frame, text = "Cancel", command = cancel_function)
    tkpack(ok_btn, side = "left", padx = 5, pady = 5)
    tkpack(cancel_btn, side = "left", padx = 5, pady = 5)
    tkbind(tt, "<Return>", ok_function)
    tkbind(tt, "<Escape>", cancel_function)
    tkfocus(tt)
    tkgrab.set(tt)
    tkwait.window(tt)
    selected_markers
  }, error = function(e) {
    warning("Marker selection popup failed: ", e$message)
    .cli_select_markers(marker_info, title)
  })
}

# Internal: CLI fallback for marker selection
.cli_select_markers <- function(marker_info, title) {
  cat("\n===", title, "===\n")
  for (i in seq_along(marker_info$indices)) {
    cat(sprintf("  %d. Column %d: %s\n", i, marker_info$indices[i], marker_info$names[i]))
  }
  cat("Enter numbers (comma-separated), 'all', or 'q' to cancel: ")
  input <- readline()
  if (tolower(input) == "q") return(NULL)
  if (tolower(input) == "all") return(marker_info$indices)
  nums <- as.integer(strsplit(trimws(input), "[,\\s]+")[[1]])
  nums <- nums[!is.na(nums) & nums >= 1 & nums <= length(marker_info$indices)]
  if (length(nums) == 0) return(marker_info$indices)
  marker_info$indices[nums]
}

#' Simplified marker selection by mode
#' @param mode One of "all", "popup", "terminal"
#' @param filtered_markers Character vector of marker names
#' @param filtered_indices Integer vector of marker column indices
#' @return Integer vector of selected marker indices
select_markers_by_mode <- function(mode = c("all", "popup", "terminal"),
                                   filtered_markers, filtered_indices) {
  mode <- match.arg(mode)
  switch(mode,
    "all" = filtered_indices,
    "popup" = {
      cat("\nFiltered markers available:\n")
      purrr::iwalk(filtered_markers,
                   ~ cat(sprintf("%d. Column %d: %s\n", .y, filtered_indices[.y], .x)))
      tryCatch({
        selected <- select_markers(list(indices = filtered_indices, names = filtered_markers))
        if (is.null(selected) || length(selected) == 0) filtered_indices else selected
      }, error = function(e) { cat("Popup failed, using all markers\n"); filtered_indices })
    },
    "terminal" = {
      cat("\nAvailable markers:\n")
      purrr::iwalk(filtered_markers,
                   ~ cat(sprintf("%d. Column %d: %s\n", .y, filtered_indices[.y], .x)))
      cat("Enter marker numbers (comma-separated, or 'all'): ")
      input <- readline()
      if (input == "all") filtered_indices else {
        nums <- as.numeric(strsplit(input, ",")[[1]])
        filtered_indices[nums[nums %in% seq_along(filtered_indices)]]
      }
    }
  )
}

#' GUI column selection with CLI fallback
#' @param column_names Character vector of available columns
#' @param title Dialog title
#' @param multiple Allow multiple selections
#' @return Character vector of selected column names
select_columns_gui <- function(column_names, title = "Select Columns", multiple = TRUE) {
  tk_multiselect(column_names, title = title, multiple = multiple, cli_fallback = TRUE)
}
