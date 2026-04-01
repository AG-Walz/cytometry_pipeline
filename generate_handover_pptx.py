#!/usr/bin/env python3
"""Generate a handover PowerPoint presentation for the Cytometry Pipeline."""

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE

# ── Colour palette ──────────────────────────────────────────────────────────
WHITE       = RGBColor(0xFF, 0xFF, 0xFF)
BLACK       = RGBColor(0x00, 0x00, 0x00)
DARK_BG     = RGBColor(0x1B, 0x2A, 0x4A)  # deep navy
LIGHT_BG    = RGBColor(0xF5, 0xF7, 0xFA)  # very light grey
ACCENT1     = RGBColor(0x3B, 0x82, 0xF6)  # blue
ACCENT2     = RGBColor(0x10, 0xB9, 0x81)  # green
ACCENT3     = RGBColor(0xF5, 0x9E, 0x0B)  # amber
ACCENT4     = RGBColor(0xEF, 0x44, 0x44)  # red
ACCENT5     = RGBColor(0x8B, 0x5C, 0xF6)  # purple
ACCENT6     = RGBColor(0x06, 0xB6, 0xD4)  # cyan
GRAY        = RGBColor(0x64, 0x74, 0x8B)  # slate
LIGHT_GRAY  = RGBColor(0xCB, 0xD5, 0xE1)
MID_GRAY    = RGBColor(0x94, 0xA3, 0xB8)

# Module colours
MOD_SHARED  = RGBColor(0x64, 0x74, 0x8B)
MOD_LOAD    = RGBColor(0x3B, 0x82, 0xF6)
MOD_PRE     = RGBColor(0x10, 0xB9, 0x81)
MOD_HARM    = RGBColor(0xF5, 0x9E, 0x0B)
MOD_CLUST   = RGBColor(0xEF, 0x44, 0x44)
MOD_ANAL    = RGBColor(0x8B, 0x5C, 0xF6)
MOD_VIZ     = RGBColor(0x06, 0xB6, 0xD4)
MOD_META    = RGBColor(0xF9, 0x73, 0x16)

prs = Presentation()
prs.slide_width  = Inches(13.333)
prs.slide_height = Inches(7.5)

# ── Helper functions ────────────────────────────────────────────────────────

def add_bg(slide, color=LIGHT_BG):
    bg = slide.background
    fill = bg.fill
    fill.solid()
    fill.fore_color.rgb = color

def add_shape_bg(slide, left, top, width, height, color, alpha=None):
    shape = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, left, top, width, height)
    shape.fill.solid()
    shape.fill.fore_color.rgb = color
    shape.line.fill.background()
    if alpha is not None:
        shape.fill.fore_color.brightness = 0
    return shape

def add_textbox(slide, left, top, width, height, text, font_size=14,
                bold=False, color=BLACK, alignment=PP_ALIGN.LEFT,
                font_name="Calibri", line_spacing=1.15):
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(font_size)
    p.font.bold = bold
    p.font.color.rgb = color
    p.font.name = font_name
    p.alignment = alignment
    p.space_after = Pt(2)
    return txBox

def add_bullet_list(slide, left, top, width, height, items, font_size=12,
                    color=BLACK, bullet_color=None, font_name="Calibri",
                    bold_prefix=True, line_spacing=1.1):
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = True
    for i, item in enumerate(items):
        if i == 0:
            p = tf.paragraphs[0]
        else:
            p = tf.add_paragraph()
        p.space_after = Pt(2)
        p.space_before = Pt(1)
        p.level = 0

        # Split on " -- " to bold the prefix
        if " -- " in item and bold_prefix:
            parts = item.split(" -- ", 1)
            run1 = p.add_run()
            run1.text = parts[0]
            run1.font.size = Pt(font_size)
            run1.font.bold = True
            run1.font.color.rgb = color
            run1.font.name = font_name
            run2 = p.add_run()
            run2.text = " -- " + parts[1]
            run2.font.size = Pt(font_size)
            run2.font.bold = False
            run2.font.color.rgb = color
            run2.font.name = font_name
        else:
            run = p.add_run()
            run.text = item
            run.font.size = Pt(font_size)
            run.font.bold = False
            run.font.color.rgb = color
            run.font.name = font_name
    return txBox

def add_colored_box(slide, left, top, width, height, title, items,
                    box_color=ACCENT1, title_color=WHITE, text_color=BLACK,
                    title_size=13, text_size=10):
    # Header bar
    hdr = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, left, top, width, Inches(0.38))
    hdr.fill.solid()
    hdr.fill.fore_color.rgb = box_color
    hdr.line.fill.background()
    tf = hdr.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = title
    p.font.size = Pt(title_size)
    p.font.bold = True
    p.font.color.rgb = title_color
    p.font.name = "Calibri"
    p.alignment = PP_ALIGN.CENTER

    # Body
    body = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                                  left, top + Inches(0.36), width, height - Inches(0.36))
    body.fill.solid()
    body.fill.fore_color.rgb = WHITE
    body.line.color.rgb = box_color
    body.line.width = Pt(1.5)
    tf = body.text_frame
    tf.word_wrap = True
    tf.margin_left = Inches(0.1)
    tf.margin_right = Inches(0.05)
    tf.margin_top = Inches(0.05)
    tf.margin_bottom = Inches(0.05)
    for i, item in enumerate(items):
        if i == 0:
            p = tf.paragraphs[0]
        else:
            p = tf.add_paragraph()
        p.space_after = Pt(1)
        p.space_before = Pt(0)
        # bold prefix
        if " -- " in item:
            parts = item.split(" -- ", 1)
            r1 = p.add_run()
            r1.text = parts[0]
            r1.font.size = Pt(text_size)
            r1.font.bold = True
            r1.font.color.rgb = text_color
            r1.font.name = "Calibri"
            r2 = p.add_run()
            r2.text = " -- " + parts[1]
            r2.font.size = Pt(text_size)
            r2.font.color.rgb = text_color
            r2.font.name = "Calibri"
        else:
            r = p.add_run()
            r.text = item
            r.font.size = Pt(text_size)
            r.font.color.rgb = text_color
            r.font.name = "Calibri"


def add_table(slide, left, top, width, height, headers, rows,
              header_color=ACCENT1, header_text_color=WHITE,
              alt_row_color=RGBColor(0xEF, 0xF6, 0xFF), font_size=10):
    n_rows = len(rows) + 1
    n_cols = len(headers)
    table_shape = slide.shapes.add_table(n_rows, n_cols, left, top, width, height)
    table = table_shape.table

    # Set column widths proportionally
    col_total = sum(1 for _ in headers)

    # Header row
    for j, h in enumerate(headers):
        cell = table.cell(0, j)
        cell.text = h
        cell.fill.solid()
        cell.fill.fore_color.rgb = header_color
        for paragraph in cell.text_frame.paragraphs:
            paragraph.font.size = Pt(font_size)
            paragraph.font.bold = True
            paragraph.font.color.rgb = header_text_color
            paragraph.font.name = "Calibri"
            paragraph.alignment = PP_ALIGN.CENTER

    # Data rows
    for i, row in enumerate(rows):
        for j, val in enumerate(row):
            cell = table.cell(i + 1, j)
            cell.text = str(val)
            if i % 2 == 1:
                cell.fill.solid()
                cell.fill.fore_color.rgb = alt_row_color
            else:
                cell.fill.solid()
                cell.fill.fore_color.rgb = WHITE
            for paragraph in cell.text_frame.paragraphs:
                paragraph.font.size = Pt(font_size)
                paragraph.font.color.rgb = BLACK
                paragraph.font.name = "Calibri"
                paragraph.alignment = PP_ALIGN.LEFT
    return table_shape


def add_arrow(slide, left, top, width, height, color=ACCENT1):
    shape = slide.shapes.add_shape(MSO_SHAPE.DOWN_ARROW, left, top, width, height)
    shape.fill.solid()
    shape.fill.fore_color.rgb = color
    shape.line.fill.background()
    return shape


def add_flow_box(slide, left, top, width, height, text, color=ACCENT1, text_color=WHITE, font_size=11):
    shape = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, left, top, width, height)
    shape.fill.solid()
    shape.fill.fore_color.rgb = color
    shape.line.fill.background()
    tf = shape.text_frame
    tf.word_wrap = True
    tf.margin_left = Inches(0.05)
    tf.margin_right = Inches(0.05)
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(font_size)
    p.font.bold = True
    p.font.color.rgb = text_color
    p.font.name = "Calibri"
    p.alignment = PP_ALIGN.CENTER
    tf.paragraphs[0].space_before = Pt(0)
    tf.paragraphs[0].space_after = Pt(0)
    return shape


# ═══════════════════════════════════════════════════════════════════════════
# SLIDE 1: TITLE
# ═══════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])  # blank
add_bg(slide, DARK_BG)

add_textbox(slide, Inches(1), Inches(1.5), Inches(11), Inches(1.2),
            "Cytometry Data Analysis Pipeline", font_size=44, bold=True,
            color=WHITE, alignment=PP_ALIGN.CENTER)
add_textbox(slide, Inches(1), Inches(2.8), Inches(11), Inches(0.6),
            "Handover Meeting", font_size=28,
            color=ACCENT6, alignment=PP_ALIGN.CENTER)

# Subtitle details
add_textbox(slide, Inches(1), Inches(4.0), Inches(11), Inches(0.5),
            "Modular flow cytometry framework  |  67 R scripts  |  5 pipelines  |  11 shared modules",
            font_size=16, color=MID_GRAY, alignment=PP_ALIGN.CENTER)

# Date
add_textbox(slide, Inches(1), Inches(5.5), Inches(11), Inches(0.4),
            "February 2026", font_size=14, color=MID_GRAY, alignment=PP_ALIGN.CENTER)


# ═══════════════════════════════════════════════════════════════════════════
# SLIDE 2: TABLE OF CONTENTS
# ═══════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "Agenda", font_size=32, bold=True, color=DARK_BG)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), ACCENT1)

toc_items = [
    "1.  Project Overview & Architecture",
    "2.  Shared Utilities (R_shared) -- 11 modules",
    "3.  Data Loading -- 6 scripts",
    "4.  Preprocessing & QC -- 13 scripts",
    "5.  Panel Harmonization & Batch Correction -- 5 scripts",
    "6.  Clustering, Validation & Annotation -- 9 scripts",
    "7.  Statistical Analysis -- 12 scripts",
    "8.  Visualization -- 13 scripts",
    "9.  Metadata QC -- 4 scripts",
    "10. Pipeline A: Raw FCS Full Analysis",
    "11. Pipeline B: CSV Expression Quick Stats",
    "12. Pipeline C: CSV Gating Abundance",
    "13. Pipeline D: FlowJo Comparison",
    "14. Pipeline E: Preprocessed FCS Fast Path",
    "15. Pipeline Comparison Matrix",
    "16. Data Formats, Bridge Functions & Dispatchers",
    "17. Decision Tree: Which Pipeline to Use",
    "18. Key Takeaways & Next Steps",
]
add_bullet_list(slide, Inches(1.2), Inches(1.2), Inches(10), Inches(5.8),
                toc_items, font_size=15, color=DARK_BG, bold_prefix=False)


# ═══════════════════════════════════════════════════════════════════════════
# SLIDE 3: PROJECT OVERVIEW
# ═══════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "1. Project Overview", font_size=32, bold=True, color=DARK_BG)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), ACCENT1)

# Left column: description
add_textbox(slide, Inches(0.5), Inches(1.2), Inches(6), Inches(0.5),
            "What is this?", font_size=18, bold=True, color=DARK_BG)
desc = [
    "Modular flow cytometry analysis framework for immune cell phenotyping",
    "Supports 5 input types: raw FCS, preprocessed FCS, CSV expression, CSV gating, FlowJo WSP",
    "67 standalone R scripts organized in 8 functional modules",
    "5 pipeline orchestrators that chain scripts into complete workflows",
    "Config-driven: swap clustering, dim-reduction, batch correction via one parameter",
    "GUI dialogs (tcltk) with CLI fallback for headless/server execution",
    "All I/O via data.table; semicolon separator by default",
]
add_bullet_list(slide, Inches(0.7), Inches(1.7), Inches(5.8), Inches(3.5),
                desc, font_size=13, color=DARK_BG)

# Right column: directory tree
add_textbox(slide, Inches(7), Inches(1.2), Inches(5.5), Inches(0.5),
            "Directory Structure", font_size=18, bold=True, color=DARK_BG)
tree = [
    "src/R_shared/ -- 11 shared utility modules",
    "src/loading/ -- 6 data loaders",
    "src/preprocessing/ -- 13 QC & transformation scripts",
    "src/harmonization/ -- 5 batch correction scripts",
    "src/clustering/ -- 9 clustering & annotation scripts",
    "src/analysis/ -- 12 statistical analysis scripts",
    "src/visualization/ -- 13 plotting scripts",
    "src/metadata_qc/ -- 4 metadata QC scripts",
    "pipelines/ -- 5 pipeline orchestrators (A-E)",
    "config/ -- Pipeline + standalone configs",
    "docs/ -- Composability guide, API ref, methods",
    "test_data/ -- Synthetic test data",
]
add_bullet_list(slide, Inches(7.2), Inches(1.7), Inches(5.5), Inches(4),
                tree, font_size=12, color=DARK_BG)

# Philosophy box
add_colored_box(slide, Inches(0.5), Inches(5.5), Inches(12.3), Inches(1.5),
                "Design Philosophy", [
                    "Scripts are building blocks -- each has one purpose with documented inputs/outputs",
                    "Pipelines are thin orchestrators -- example workflows, not the only way to use scripts",
                    "Bridge functions convert between data formats -- enabling free mix-and-match across pipelines",
                    "Method dispatchers allow swapping algorithms via config -- no code changes needed",
                ], box_color=ACCENT1, text_size=12)


# ═══════════════════════════════════════════════════════════════════════════
# SLIDE 4: ARCHITECTURE -- HIGH LEVEL DATA FLOW
# ═══════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "Architecture: High-Level Data Flow", font_size=32, bold=True, color=DARK_BG)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), ACCENT1)

# Flow boxes
y = 1.3
box_h = 0.55
arrow_h = 0.28
gap = 0.08
center_x = 4.5
box_w = 4.0
arrow_w = 0.35

stages = [
    ("Loading", MOD_LOAD, "FCS / CSV / FlowJo WSP files"),
    ("Preprocessing", MOD_PRE, "Compensation, QC, transformation, gating"),
    ("Harmonization", MOD_HARM, "Panel unification, batch correction"),
    ("Clustering", MOD_CLUST, "Leiden / FlowSOM / PhenoGraph / DBSCAN"),
    ("Analysis", MOD_ANAL, "Abundance, enrichment, diff. expression, survival"),
    ("Visualization", MOD_VIZ, "UMAP, t-SNE, heatmaps, barplots, etc."),
]

for i, (label, color, desc) in enumerate(stages):
    cur_y = y + i * (box_h + arrow_h + gap)
    add_flow_box(slide, Inches(center_x), Inches(cur_y), Inches(box_w), Inches(box_h),
                 label, color=color, font_size=14)
    # description to the right
    add_textbox(slide, Inches(center_x + box_w + 0.3), Inches(cur_y + 0.05),
                Inches(4), Inches(box_h),
                desc, font_size=12, color=GRAY)
    # arrow down
    if i < len(stages) - 1:
        add_arrow(slide, Inches(center_x + box_w / 2 - arrow_w / 2),
                  Inches(cur_y + box_h), Inches(arrow_w), Inches(arrow_h),
                  color=LIGHT_GRAY)

# Left side: data format annotations
add_textbox(slide, Inches(0.5), Inches(1.3), Inches(4), Inches(0.4),
            "Data Formats at Each Stage:", font_size=14, bold=True, color=DARK_BG)
formats = [
    "flowFrame list -- after Loading, Preprocessing",
    "Numeric matrix -- after bridge function",
    "Gating CSV -- after clustering (per-sample files)",
    "Tracking matrix -- sample x cluster cell counts",
    "Long data.frame -- tidy format for ggplot2",
]
add_bullet_list(slide, Inches(0.5), Inches(1.7), Inches(3.8), Inches(3.5),
                formats, font_size=11, color=DARK_BG)

# Side modules
add_flow_box(slide, Inches(0.5), Inches(5.5), Inches(2.5), Inches(0.45),
             "Metadata QC", color=MOD_META, font_size=12)
add_textbox(slide, Inches(3.1), Inches(5.55), Inches(3), Inches(0.4),
            "Acquisition QC, gating comparison", font_size=11, color=GRAY)

add_flow_box(slide, Inches(0.5), Inches(6.1), Inches(2.5), Inches(0.45),
             "R_shared Utilities", color=MOD_SHARED, font_size=12)
add_textbox(slide, Inches(3.1), Inches(6.15), Inches(4), Inches(0.4),
            "I/O, GUI, memory, stats, dispatchers (used by all modules)",
            font_size=11, color=GRAY)


# ═══════════════════════════════════════════════════════════════════════════
# SLIDE 5: SHARED UTILITIES
# ═══════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "2. Shared Utilities -- src/R_shared/ (11 modules)", font_size=28, bold=True, color=DARK_BG)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), MOD_SHARED)

# Left column
modules_left = [
    ("io_utils.R", "File I/O, CSV detection, column classification, die(), make_outdir()"),
    ("memory_utils.R", "Memory monitoring, chunked CSV reading for large files"),
    ("fcs_utils.R", "FCS reading, event counting, compensation extraction, load_fcs_list()"),
    ("gui_utils.R", "tcltk dialogs with CLI readline fallback, marker selection"),
    ("flowframe_utils.R", "FlowFrame manipulation, build_combined_expr_matrix(), cellID generation"),
    ("cluster_utils.R", "Cluster naming, label conversion, binary-to-single-cluster"),
]
modules_right = [
    ("plot_utils.R", "Multi-format saving (PDF/PNG/SVG), filename sanitization, barplots"),
    ("stats_utils.R", "Fisher exact, mixed-effects, NB-GLMM, PCA, UMAP helpers"),
    ("clustering_dispatch.R", "run_clustering() -- routes to Leiden/FlowSOM/PhenoGraph/DBSCAN"),
    ("dimred_dispatch.R", "run_dimred() -- routes to UMAP/t-SNE/diffusion map"),
    ("batch_correction_dispatch.R", "run_batch_correction_dispatch() -- Harmony/CytoNorm/ComBat/none"),
]

y_start = 1.2
for i, (name, desc) in enumerate(modules_left):
    add_colored_box(slide, Inches(0.4), Inches(y_start + i * 1.0), Inches(6.2), Inches(0.9),
                    name, [desc], box_color=MOD_SHARED, text_size=11, title_size=12)

for i, (name, desc) in enumerate(modules_right):
    add_colored_box(slide, Inches(6.8), Inches(y_start + i * 1.0), Inches(6.2), Inches(0.9),
                    name, [desc], box_color=MOD_SHARED, text_size=11, title_size=12)

add_textbox(slide, Inches(6.8), Inches(y_start + 5 * 1.0 + 0.1), Inches(6), Inches(0.4),
            "All scripts source these via:  shared_dir <- file.path(dirname(getwd()), \"src\", \"R_shared\")",
            font_size=11, color=GRAY)


# ═══════════════════════════════════════════════════════════════════════════
# SLIDE 6: LOADING
# ═══════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "3. Data Loading -- src/loading/ (6 scripts)", font_size=28, bold=True, color=DARK_BG)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), MOD_LOAD)

headers = ["Script", "Main Function(s)", "Input", "Output", "Used By"]
rows = [
    ["load_fcs_raw.R", "load_fcs_raw()", "Directory of .fcs files",
     "List of flowFrames + sample names", "Pipeline A, E"],
    ["load_csv_expression.R", "load_csv_expression()", "Directory of CSV expression files",
     "Named list of data.tables", "Pipeline B"],
    ["load_csv_gating.R", "load_csv_gating()\nprocess_single_csv_gating()", "CSVs with gating columns + FCS files",
     "Per-sample gating CSVs or long-format CSVs", "Pipeline C"],
    ["load_manual_gating.R", "load_manual_gating()", "CSVs with manual gating columns",
     "Named list of label vectors", "Custom workflows"],
    ["load_flowjo_wsp.R", "load_flowjo_wsp()\nload_multiple_workspaces()", ".wsp files + FCS directory",
     "GatingSet objects, population stats/strings", "Pipeline D"],
    ["create_flowframes.R", "create_flowframes_from_matrices()\ncreate_merged_flowset()", "Named list of expression matrices",
     "List of flowFrames, flowSet", "Bridge function"],
]
add_table(slide, Inches(0.4), Inches(1.2), Inches(12.5), Inches(5.5),
          headers, rows, header_color=MOD_LOAD, font_size=11)


# ═══════════════════════════════════════════════════════════════════════════
# SLIDE 7-8: PREPROCESSING (split across 2 slides)
# ═══════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "4. Preprocessing & QC -- src/preprocessing/ (1/2)", font_size=28, bold=True, color=DARK_BG)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), MOD_PRE)

headers = ["Script", "Main Function(s)", "Input", "Output"]
rows_p1 = [
    ["compensation.R", "apply_compensation()\nextract_compensation_matrix()", "flowFrame + comp. matrix",
     "Compensated flowFrame"],
    ["transformation.R", "transform_flowframe()\nfind_optimal_b_cofactor()", "flowFrame, markers, method",
     "Transformed flowFrame (arcsinh/logicle/biexp/log)"],
    ["spectral_unmixing.R", "unmix_spectral()\nbuild_unmixing_matrix()", "Single-stain controls + spectral data",
     "Unmixed flowFrame with virtual channels"],
    ["bead_normalization.R", "normalize_beads()\nidentify_beads()", "flowFrame + bead channel names",
     "Normalized flowFrame + bead indices + factors"],
    ["doublet_detection.R", "detect_doublets()", "flowFrame with scatter params",
     "Filtered flowFrame + doublet indices + singlet rate"],
    ["live_dead_gating.R", "gate_live_cells()\nfind_density_valley()", "flowFrame + viability marker",
     "Live flowFrame + dead indices + viability rate"],
    ["quality_control.R", "quality_control()", "flowFrame",
     "QC-filtered flowFrame (flowAI/PeacoQC)"],
]
add_table(slide, Inches(0.4), Inches(1.2), Inches(12.5), Inches(5.8),
          headers, rows_p1, header_color=MOD_PRE, font_size=11)

# Page 2
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "4. Preprocessing & QC -- src/preprocessing/ (2/2)", font_size=28, bold=True, color=DARK_BG)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), MOD_PRE)

rows_p2 = [
    ["downsampling.R", "downsample_flowframes()", "List of flowFrames + target count",
     "Downsampled flowFrame list"],
    ["normalization.R", "z_score_scaling()\nquantile_norm()\nharmony_normalize()", "flowFrame or matrix + batch vector",
     "Normalized flowFrame or matrix"],
    ["supercell.R", "run_supercell()\nrun_supercell_per_cluster()", "Expr. matrix + sample labels + markers",
     "Supercell matrix + labels + cell map + sizes"],
    ["automated_gating.R", "run_auto_gating()", "flowFrame + channels + optional template",
     "Gated flowFrame + gate coords + population stats"],
    ["gate_utils.R", "apply_gate_sequence()\nboolean_gate()", "flowFrame + gate list + expression",
     "Filtered flowFrame + gate results"],
    ["preprocess_pipeline.R", "run_preprocessing()", "flowFrame list + config params",
     "Preprocessed flowFrame list (QC+transform+norm)"],
]
add_table(slide, Inches(0.4), Inches(1.2), Inches(12.5), Inches(5),
          headers, rows_p2, header_color=MOD_PRE, font_size=11)

add_textbox(slide, Inches(0.4), Inches(5.7), Inches(12), Inches(0.5),
            "Typical preprocessing chain:  compensation -> transformation -> QC -> doublet removal -> live/dead gating -> downsampling -> normalization",
            font_size=13, bold=True, color=MOD_PRE)


# ═══════════════════════════════════════════════════════════════════════════
# SLIDE 9: HARMONIZATION
# ═══════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "5. Panel Harmonization & Batch Correction -- src/harmonization/ (5 scripts)", font_size=26, bold=True, color=DARK_BG)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), MOD_HARM)

headers = ["Script", "Main Function(s)", "Input", "Output"]
rows = [
    ["panel_comparison.R", "compare_panels_detailed()", "flowFrame list",
     "Common/all markers, presence matrix, panel groups"],
    ["panel_homogenization.R", "homogenize_flowframes()\nreplace_markernames()", "flowFrame list + target order + replacement map",
     "Homogenized flowFrame list with unified panel"],
    ["batch_correction.R", "run_batch_correction()", "flowFrame list + batch labels + markers",
     "Corrected matrix + sample/batch labels (Harmony)"],
    ["cytonorm.R", "run_cytonorm()", "Train/validation flowFrame lists + batch labels",
     "Corrected flowFrames + CytoNorm model"],
    ["combat.R", "run_combat()", "Expression matrix + batch vector + optional covariates",
     "Corrected matrix (empirical Bayes)"],
]
add_table(slide, Inches(0.4), Inches(1.2), Inches(12.5), Inches(4.2),
          headers, rows, header_color=MOD_HARM, font_size=11)

# Add note about dispatcher
add_colored_box(slide, Inches(0.4), Inches(5.5), Inches(12.5), Inches(1.5),
                "Batch Correction Dispatcher", [
                    "run_batch_correction_dispatch(data, method, batch_labels, ...) in R_shared/batch_correction_dispatch.R",
                    "Routes to: harmony | cytonorm | combat | none -- based on config param batch_method",
                    "Accepts BOTH flowFrame list and numeric matrix -- auto-converts as needed",
                    "Returns: list(corrected_matrix, sample_labels, batch_labels, markers, corrected_frames)",
                ], box_color=MOD_HARM, text_size=11)


# ═══════════════════════════════════════════════════════════════════════════
# SLIDE 10-11: CLUSTERING
# ═══════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "6. Clustering, Validation & Annotation -- src/clustering/ (9 scripts)", font_size=26, bold=True, color=DARK_BG)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), MOD_CLUST)

headers = ["Script", "Main Function(s)", "Input", "Output"]
rows = [
    ["leiden.R", "run_leiden()\nclusters_to_gating_csv()", "Expr. matrix, n_neighbors, resolution",
     "Integer cluster IDs + UMAP coords + per-sample gating CSVs"],
    ["flowsom.R", "run_flowsom()\nrun_flowsom_flowframes()", "Expr. matrix or flowFrame list, n_clusters",
     "Integer cluster IDs + FlowSOM object"],
    ["phenograph.R", "run_phenograph()", "Expr. matrix, k (nearest neighbors)",
     "Integer cluster IDs"],
    ["dbscan.R", "run_dbscan()\nestimate_eps()", "Expr. matrix, eps, min_pts",
     "Cluster IDs + noise indices + noise fraction"],
    ["meta_clustering.R", "run_resolution_sweep()\nconsensus_clusters()", "Expr. matrix + method + resolution range",
     "Best clusters + scores + all clusterings"],
    ["cluster_validation.R", "validate_clusters()", "Expr. matrix + cluster IDs + optional true labels",
     "Metrics table (silhouette, CH, DB, ARI)"],
    ["cell_annotation.R", "annotate_clusters()\nannotate_by_correlation()", "Expr. matrix + cluster IDs + reference profiles",
     "Annotation table (cluster, type, confidence, medians)"],
    ["flowjo_import.R", "import_flowjo_gates()", ".wsp file + FCS directory",
     "GatingSet + assignments + populations"],
    ["cluster_assignment.R", "get_cluster_assignments()", "Source type + file paths/objects",
     "Data frame (cell_id, cluster, cluster_name)"],
]
add_table(slide, Inches(0.4), Inches(1.15), Inches(12.5), Inches(5.2),
          headers, rows, header_color=MOD_CLUST, font_size=10)

add_colored_box(slide, Inches(0.4), Inches(6.3), Inches(12.5), Inches(0.9),
                "Clustering Dispatcher", [
                    "run_clustering(expr_matrix, method, params) -- all methods return integer vector of cluster IDs",
                    "Dispatcher coerces Leiden output from numeric to integer; DBSCAN noise is cluster 0",
                ], box_color=MOD_CLUST, text_size=11)


# ═══════════════════════════════════════════════════════════════════════════
# SLIDE 12-13: ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "7. Statistical Analysis -- src/analysis/ (1/2)", font_size=28, bold=True, color=DARK_BG)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), MOD_ANAL)

headers = ["Script", "Main Function(s)", "Input", "Output"]
rows = [
    ["cluster_splitting.R", "split_by_cluster()\ntracking_matrix_to_long()", "Per-cell CSVs + cluster regex",
     "Per-cluster CSVs + tracking matrix + summary"],
    ["expression_matrices.R", "create_expression_matrices()\ncombine_expression_matrices()", "Cluster-split CSVs + markers",
     "Named list of expression matrices per cluster"],
    ["differential_abundance.R", "build_abundance_matrix()\nrun_differential_abundance()", "Per-cell CSVs + conditions mapping",
     "Abundance stats + volcano plots (GLMM)"],
    ["cluster_enrichment.R", "run_cluster_enrichment()", "Tracking matrix + conditions",
     "Enrichment results (Fisher exact / mixed-effects logistic)"],
    ["marker_statistics.R", "run_marker_statistics()", "Expression matrices + group labels",
     "Per-marker stats (Wilcoxon/KW/lmer) + volcano plots"],
    ["differential_expression.R", "run_differential_expression()", "Expression matrix + conditions",
     "DE results per marker (limma/lmer/Wilcoxon)"],
]
add_table(slide, Inches(0.4), Inches(1.2), Inches(12.5), Inches(5),
          headers, rows, header_color=MOD_ANAL, font_size=11)

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "7. Statistical Analysis -- src/analysis/ (2/2)", font_size=28, bold=True, color=DARK_BG)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), MOD_ANAL)

rows2 = [
    ["effect_sizes.R", "compute_effect_sizes()", "Expression matrix + group labels",
     "Cohen's d, Cliff's delta, Hedges' g per marker"],
    ["correlation_analysis.R", "run_correlation_analysis()", "Expression matrix + cluster IDs",
     "Spearman/Pearson correlation matrices per cluster"],
    ["diffcyt_analysis.R", "run_diffcyt()", "flowFrame list + experiment info + design",
     "diffcyt DA/DS results (edgeR/limma-based)"],
    ["survival_analysis.R", "run_survival_analysis()", "Cluster abundances + survival data",
     "Cox PH results + Kaplan-Meier plots"],
    ["abundance_plotting.R", "plot_abundance()", "Tracking matrix + sample metadata",
     "Stacked, proportional, heatmap, totals plots"],
    ["cell_ratio.R", "compute_cell_ratios()", "Cluster abundances + ratio definitions",
     "Cell type ratio table + comparisons"],
]
add_table(slide, Inches(0.4), Inches(1.2), Inches(12.5), Inches(5),
          headers, rows2, header_color=MOD_ANAL, font_size=11)

add_textbox(slide, Inches(0.4), Inches(5.5), Inches(12), Inches(1),
            "Statistical methods available: Fisher's exact, mixed-effects logistic regression, "
            "negative binomial GLMM (lme4), Wilcoxon, Kruskal-Wallis, limma, diffcyt (edgeR), "
            "Cox proportional hazards, Cohen's d, Cliff's delta, Hedges' g, Spearman/Pearson correlations",
            font_size=12, color=MOD_ANAL)


# ═══════════════════════════════════════════════════════════════════════════
# SLIDE 14-15: VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "8. Visualization -- src/visualization/ (1/2)", font_size=28, bold=True, color=DARK_BG)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), MOD_VIZ)

headers = ["Script", "Main Function(s)", "Input", "Output"]
rows = [
    ["umap.R", "plot_umap()\nplot_umap_from_matrix()", "flowFrame list OR expr. matrix + labels",
     "UMAP plots (PNG/PDF/SVG) + coordinates"],
    ["tsne.R", "plot_tsne()\nplot_tsne_from_matrix()", "flowFrame list OR expr. matrix",
     "t-SNE plots + coordinates"],
    ["diffusion_map.R", "plot_diffusion_map()", "Expr. matrix + labels",
     "Diffusion map plots + coordinates (destiny pkg)"],
    ["heatmaps.R", "plot_heatmap()", "Expression matrix + cluster/sample labels",
     "Hierarchical cluster heatmaps"],
    ["density_plots.R", "plot_density()\nplot_density_from_matrix()", "flowFrame list OR matrix + markers",
     "1D/2D density plots per marker"],
    ["biaxial_plots.R", "plot_biaxial()\nplot_biaxial_from_matrix()", "flowFrame list OR matrix + marker pairs",
     "2D scatter with density/marker/cluster coloring"],
    ["violin_plots.R", "plot_violin()", "Expr. matrix + cluster/condition labels",
     "Violin/box/sina plots per marker"],
]
add_table(slide, Inches(0.4), Inches(1.15), Inches(12.5), Inches(5.2),
          headers, rows, header_color=MOD_VIZ, font_size=11)

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "8. Visualization -- src/visualization/ (2/2)", font_size=28, bold=True, color=DARK_BG)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), MOD_VIZ)

rows2 = [
    ["radar_plots.R", "plot_radar()", "Cluster median expression + markers",
     "Spider/radar charts per cluster"],
    ["contour_plots.R", "plot_contour()", "2D coordinates (UMAP/tSNE) + values",
     "2D density contour overlays"],
    ["barplots.R", "plot_composition_barplot()", "Tracking matrix (long format)",
     "Stacked composition barplots"],
    ["upset_plots.R", "plot_upset()", "Binary membership matrix",
     "UpSet intersection plots"],
    ["venn_diagrams.R", "plot_venn()", "List of sets or binary matrix",
     "Venn/Euler diagrams"],
    ["pca.R", "plot_pca()", "Expr. matrix + labels",
     "PCA plots + variance explained + loadings"],
]
add_table(slide, Inches(0.4), Inches(1.15), Inches(12.5), Inches(4.5),
          headers, rows2, header_color=MOD_VIZ, font_size=11)

add_colored_box(slide, Inches(0.4), Inches(5.2), Inches(12.5), Inches(1.5),
                "Dual-Input Pattern", [
                    "UMAP, t-SNE, density, biaxial all accept EITHER flowFrame list OR numeric matrix",
                    "flowFrame path: plot_umap(flow_frames, ...) -- internal extraction",
                    "Matrix path: plot_umap_from_matrix(expr_matrix, labels, ...) -- direct, faster",
                    "Dimred Dispatcher: run_dimred(expr_matrix, method) -- unified interface for all three",
                ], box_color=MOD_VIZ, text_size=11)


# ═══════════════════════════════════════════════════════════════════════════
# SLIDE 16: METADATA QC
# ═══════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "9. Metadata QC -- src/metadata_qc/ (4 scripts)", font_size=28, bold=True, color=DARK_BG)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), MOD_META)

headers = ["Script", "Main Function(s)", "Input", "Output"]
rows = [
    ["check_metadata.R", "check_metadata()", "Sample metadata table",
     "Validation report (missing, inconsistent values)"],
    ["check_acquisition.R", "check_acquisition()", "FCS files or metadata CSV",
     "Acquisition date/parameter QC plots"],
    ["overlap_statistics.R", "compute_fleiss_kappa()", "Multiple gating assignments",
     "Inter-rater reliability (Fleiss' kappa)"],
    ["compare_gatings.R", "compare_gatings()", "Two or more gating CSVs",
     "Gating comparison visualizations (Venn/UpSet)"],
]
add_table(slide, Inches(0.4), Inches(1.2), Inches(12.5), Inches(3.5),
          headers, rows, header_color=MOD_META, font_size=12)

add_textbox(slide, Inches(0.5), Inches(5.0), Inches(12), Inches(1),
            "These scripts are primarily used in Pipeline D (FlowJo comparison) and as standalone QC checks "
            "before running any analysis pipeline. They help validate data quality and gating consistency.",
            font_size=14, color=DARK_BG)


# ═══════════════════════════════════════════════════════════════════════════
# SLIDE 17: PIPELINE A
# ═══════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "10. Pipeline A: Raw FCS -> Full Analysis", font_size=28, bold=True, color=DARK_BG)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), ACCENT4)

# Flow diagram - left side
stages_a = [
    ("load_fcs_raw()", MOD_LOAD),
    ("run_preprocessing()", MOD_PRE),
    ("homogenize_flowframes()", MOD_HARM),
    ("run_batch_correction_dispatch()", MOD_HARM),
    ("downsample_flowframes()", MOD_PRE),
    ("build_combined_expr_matrix()", MOD_SHARED),
    ("run_clustering()", MOD_CLUST),
    ("clusters_to_gating_csv()", MOD_CLUST),
    ("split_by_cluster()", MOD_ANAL),
    ("run_differential_abundance()", MOD_ANAL),
    ("run_dimred()", MOD_VIZ),
]

box_w_a = 3.0
box_h_a = 0.35
arrow_h_a = 0.15
x_start = 0.5
y_start_a = 1.1

for i, (func, color) in enumerate(stages_a):
    cur_y = y_start_a + i * (box_h_a + arrow_h_a + 0.02)
    add_flow_box(slide, Inches(x_start), Inches(cur_y), Inches(box_w_a), Inches(box_h_a),
                 func, color=color, text_color=WHITE, font_size=10)
    if i < len(stages_a) - 1:
        add_arrow(slide, Inches(x_start + box_w_a / 2 - 0.12),
                  Inches(cur_y + box_h_a), Inches(0.24), Inches(arrow_h_a),
                  color=LIGHT_GRAY)

# Right side: details
add_textbox(slide, Inches(4.0), Inches(1.1), Inches(4.5), Inches(0.4),
            "Input", font_size=16, bold=True, color=DARK_BG)
add_textbox(slide, Inches(4.0), Inches(1.45), Inches(4.5), Inches(0.4),
            "Raw .fcs files in a directory", font_size=13, color=GRAY)

add_textbox(slide, Inches(4.0), Inches(2.0), Inches(4.5), Inches(0.4),
            "Swappable Methods", font_size=16, bold=True, color=DARK_BG)
methods = [
    "Clustering: leiden | flowsom | phenograph | dbscan",
    "Dim. reduction: umap | tsne | diffusion_map",
    "Batch correction: harmony | cytonorm | combat | none",
]
add_bullet_list(slide, Inches(4.0), Inches(2.4), Inches(4.5), Inches(1.5),
                methods, font_size=12, color=DARK_BG)

add_textbox(slide, Inches(4.0), Inches(3.8), Inches(4.5), Inches(0.4),
            "Outputs", font_size=16, bold=True, color=DARK_BG)
outputs = [
    "Preprocessed flowFrame list",
    "Per-sample gating CSVs with cluster columns",
    "Per-cluster expression matrices",
    "Differential abundance results + volcano plots",
    "UMAP/t-SNE/diffusion map coordinates + plots",
]
add_bullet_list(slide, Inches(4.0), Inches(4.2), Inches(4.5), Inches(2.5),
                outputs, font_size=12, color=DARK_BG)

# Far right: summary box
add_colored_box(slide, Inches(8.8), Inches(1.1), Inches(4.2), Inches(3.0),
                "Most Comprehensive Pipeline", [
                    "Covers the entire workflow from raw instrument data to publication-ready figures",
                    "Handles QC, compensation, transformation, panel unification, batch correction, "
                    "clustering, statistical analysis, and visualization",
                    "All method choices configurable via a single config file",
                ], box_color=ACCENT4, text_size=11)

add_colored_box(slide, Inches(8.8), Inches(4.3), Inches(4.2), Inches(2.5),
                "Output Directories", [
                    "01_preprocessing/ -- QC'd flowFrames",
                    "02_harmonization/ -- unified panels",
                    "03_clustering/ -- cluster assignments",
                    "04_analysis/ -- stats + split clusters",
                    "05_visualization/ -- plots + coordinates",
                ], box_color=GRAY, text_size=11)


# ═══════════════════════════════════════════════════════════════════════════
# SLIDE 18: PIPELINE B
# ═══════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "11. Pipeline B: CSV Expression -> Quick Stats", font_size=28, bold=True, color=DARK_BG)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), ACCENT1)

# Flow
stages_b = [
    ("load_csv_expression()", MOD_LOAD),
    ("run_marker_statistics()", MOD_ANAL),
    ("combine_expression_matrices()", MOD_ANAL),
    ("plot_heatmap()", MOD_VIZ),
]
y_s = 1.8
for i, (func, color) in enumerate(stages_b):
    cx = 1.5
    cy = y_s + i * 1.2
    add_flow_box(slide, Inches(cx), Inches(cy), Inches(3.5), Inches(0.5),
                 func, color=color, font_size=13)
    if i < len(stages_b) - 1:
        add_arrow(slide, Inches(cx + 1.5), Inches(cy + 0.5), Inches(0.3), Inches(0.6),
                  color=LIGHT_GRAY)

# Details
add_colored_box(slide, Inches(6.0), Inches(1.5), Inches(6.8), Inches(2.0),
                "Pipeline B Summary", [
                    "Input: CSV expression matrices (one file per sample, rows=cells, cols=markers)",
                    "No preprocessing, clustering, or gating -- purely exploratory",
                    "Quickest pipeline: runs in under 5 minutes",
                    "Best for: initial marker validation, QC checks, expression profiling",
                ], box_color=ACCENT1, text_size=12)

add_colored_box(slide, Inches(6.0), Inches(3.8), Inches(6.8), Inches(1.5),
                "Outputs", [
                    "Marker statistics CSV -- per-marker mean, median, percentiles across samples",
                    "Expression heatmap -- hierarchical clustering of markers x samples",
                ], box_color=GRAY, text_size=12)

add_textbox(slide, Inches(6.0), Inches(5.6), Inches(6.8), Inches(1),
            "Tip: Can be extended with custom scripts. Add clustering via run_clustering() "
            "on the combined matrix, then export via clusters_to_gating_csv() to feed into Pipeline C.",
            font_size=13, color=MOD_ANAL)


# ═══════════════════════════════════════════════════════════════════════════
# SLIDE 19: PIPELINE C
# ═══════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "12. Pipeline C: CSV Gating -> Cluster Abundance", font_size=28, bold=True, color=DARK_BG)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), ACCENT2)

stages_c = [
    ("load_csv_gating()", MOD_LOAD),
    ("split_by_cluster()", MOD_ANAL),
    ("run_cluster_enrichment()\nOR run_differential_abundance()", MOD_ANAL),
    ("plot_composition_barplot()", MOD_VIZ),
    ("plot_abundance()", MOD_VIZ),
]
y_s = 1.5
for i, (func, color) in enumerate(stages_c):
    cx = 0.5
    cy = y_s + i * 1.1
    add_flow_box(slide, Inches(cx), Inches(cy), Inches(4.5), Inches(0.5),
                 func, color=color, font_size=12)
    if i < len(stages_c) - 1:
        add_arrow(slide, Inches(cx + 2.0), Inches(cy + 0.5), Inches(0.3), Inches(0.5),
                  color=LIGHT_GRAY)

add_colored_box(slide, Inches(5.5), Inches(1.3), Inches(7.3), Inches(2.2),
                "Pipeline C Summary", [
                    "Input: CSVs with expression columns + binary cluster/gating columns",
                    "Skips clustering -- works with pre-existing gates (FlowJo, Pipeline A, manual)",
                    "Branching analysis: enrichment (Fisher exact) OR abundance (GLMM)",
                    "Generates comprehensive abundance/composition visualizations",
                    "Configurable analysis_method: 'enrichment' or 'abundance'",
                ], box_color=ACCENT2, text_size=12)

add_colored_box(slide, Inches(5.5), Inches(3.8), Inches(7.3), Inches(2.5),
                "Outputs", [
                    "Per-cluster CSVs split by population",
                    "Statistical results: enrichment (Fisher) or abundance (NB-GLMM) + volcano plots",
                    "Composition barplots (stacked)",
                    "Abundance plots: stacked, proportional, heatmap, totals",
                ], box_color=GRAY, text_size=12)


# ═══════════════════════════════════════════════════════════════════════════
# SLIDE 20: PIPELINE D
# ═══════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "13. Pipeline D: FlowJo Comparison -> Gating Validation", font_size=28, bold=True, color=DARK_BG)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), ACCENT3)

stages_d = [
    ("load_multiple_workspaces()", MOD_LOAD),
    ("compare_gatings()", MOD_META),
    ("compute_fleiss_kappa()", MOD_META),
]
y_s = 2.0
for i, (func, color) in enumerate(stages_d):
    cx = 1.5
    cy = y_s + i * 1.3
    add_flow_box(slide, Inches(cx), Inches(cy), Inches(4.0), Inches(0.55),
                 func, color=color, font_size=14)
    if i < len(stages_d) - 1:
        add_arrow(slide, Inches(cx + 1.7), Inches(cy + 0.55), Inches(0.3), Inches(0.65),
                  color=LIGHT_GRAY)

add_colored_box(slide, Inches(6.5), Inches(1.5), Inches(6.3), Inches(2.2),
                "Pipeline D Summary", [
                    "Input: FlowJo .wsp workspace files + corresponding .fcs files",
                    "Specialized for comparing gating strategies across analysts/methods",
                    "Does NOT perform clustering or statistical abundance testing",
                    "Requires study-specific config (no default config exists)",
                ], box_color=ACCENT3, text_size=12)

add_colored_box(slide, Inches(6.5), Inches(4.0), Inches(6.3), Inches(1.8),
                "Outputs", [
                    "Gating comparison visualizations (Venn/Euler diagrams)",
                    "Fleiss' kappa statistics for inter-rater agreement",
                    "Population assignment concordance tables",
                ], box_color=GRAY, text_size=12)

add_textbox(slide, Inches(0.5), Inches(5.8), Inches(12), Inches(1),
            "Typical use: Validating manual gating consistency between labs, comparing automated "
            "vs. manual gating, assessing gating protocol reproducibility across cytometrists.",
            font_size=14, color=DARK_BG)


# ═══════════════════════════════════════════════════════════════════════════
# SLIDE 21: PIPELINE E
# ═══════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "14. Pipeline E: Preprocessed FCS -> Fast UMAP/Clustering", font_size=28, bold=True, color=DARK_BG)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), ACCENT5)

stages_e = [
    ("load_fcs_raw()", MOD_LOAD),
    ("homogenize_flowframes()", MOD_HARM),
    ("downsample_flowframes()", MOD_PRE),
    ("build_combined_expr_matrix()", MOD_SHARED),
    ("run_clustering()", MOD_CLUST),
    ("run_dimred()", MOD_VIZ),
]
y_s = 1.5
for i, (func, color) in enumerate(stages_e):
    cx = 1.0
    cy = y_s + i * 0.85
    add_flow_box(slide, Inches(cx), Inches(cy), Inches(3.8), Inches(0.45),
                 func, color=color, font_size=12)
    if i < len(stages_e) - 1:
        add_arrow(slide, Inches(cx + 1.6), Inches(cy + 0.45), Inches(0.3), Inches(0.3),
                  color=LIGHT_GRAY)

add_colored_box(slide, Inches(5.5), Inches(1.3), Inches(7.3), Inches(2.5),
                "Pipeline E Summary", [
                    "Input: Already-preprocessed .fcs files (no QC/compensation needed)",
                    "Lightweight version of Pipeline A -- skips QC and preprocessing",
                    "Use when FCS files preprocessed externally (CellRanger, FlowJo, Cytobank)",
                    "Fast path to visualization and clustering",
                    "Same swappable methods as Pipeline A (clustering + dimred dispatchers)",
                ], box_color=ACCENT5, text_size=12)

add_colored_box(slide, Inches(5.5), Inches(4.1), Inches(7.3), Inches(1.5),
                "Outputs", [
                    "Harmonized flowFrame list",
                    "Cluster assignments",
                    "Dimensionality reduction coordinates + plots",
                ], box_color=GRAY, text_size=12)

add_colored_box(slide, Inches(5.5), Inches(5.9), Inches(7.3), Inches(1.0),
                "Pipeline A vs E", [
                    "Pipeline A: QC + Preprocess + Harmonize + Cluster + Analyze + Visualize",
                    "Pipeline E: Harmonize + Cluster + Visualize (skips QC, preprocess, analysis)",
                ], box_color=ACCENT5, text_size=11)


# ═══════════════════════════════════════════════════════════════════════════
# SLIDE 22: PIPELINE COMPARISON MATRIX
# ═══════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "15. Pipeline Comparison Matrix", font_size=28, bold=True, color=DARK_BG)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), ACCENT1)

headers = ["Feature", "Pipeline A", "Pipeline B", "Pipeline C", "Pipeline D", "Pipeline E"]
check = "\u2713"
cross = "\u2717"
rows = [
    ["Input type", "Raw FCS", "CSV expression", "CSV with gates", "FlowJo WSP", "Preprocessed FCS"],
    ["QC / Preprocessing", check, cross, cross, cross, cross],
    ["Panel harmonization", check, cross, cross, cross, check],
    ["Batch correction", check, cross, cross, cross, "Optional"],
    ["Clustering", check, cross, cross, cross, check],
    ["Differential abundance", check, cross, check, cross, cross],
    ["Enrichment stats", cross, cross, check, cross, cross],
    ["Marker statistics", cross, check, cross, cross, cross],
    ["Gating comparison", cross, cross, cross, check, cross],
    ["Visualization type", "UMAP/tSNE/DM", "Heatmap", "Barplots", "Venn/UpSet", "UMAP/tSNE/DM"],
    ["Best use case", "Full analysis", "Quick QC", "Existing gates", "Method validation", "Fast iteration"],
]
add_table(slide, Inches(0.4), Inches(1.2), Inches(12.5), Inches(5.5),
          headers, rows, header_color=ACCENT1, font_size=11)


# ═══════════════════════════════════════════════════════════════════════════
# SLIDE 23: DATA FORMATS & BRIDGE FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "16. Data Formats, Bridge Functions & Dispatchers", font_size=26, bold=True, color=DARK_BG)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), ACCENT1)

# Data formats
add_textbox(slide, Inches(0.5), Inches(1.1), Inches(6), Inches(0.4),
            "5 Core Data Formats", font_size=18, bold=True, color=DARK_BG)
formats = [
    ("flowFrame List", "Named list of flowCore::flowFrame objects (one per sample)", MOD_LOAD),
    ("Numeric Matrix", "Rows=cells, cols=markers. Standard input for stats/clustering", MOD_ANAL),
    ("Gating CSV", "Semicolon-separated; expression + binary cluster columns", MOD_CLUST),
    ("Tracking Matrix", "Wide data.frame: rows=samples, cols=clusters, vals=counts", MOD_VIZ),
    ("Long Data.frame", "Tidy: sample-cluster rows with file, cluster, count cols", MOD_VIZ),
]
for i, (name, desc, color) in enumerate(formats):
    add_flow_box(slide, Inches(0.5), Inches(1.55 + i * 0.65), Inches(2.0), Inches(0.5),
                 name, color=color, font_size=11)
    add_textbox(slide, Inches(2.7), Inches(1.55 + i * 0.65), Inches(4), Inches(0.5),
                desc, font_size=11, color=DARK_BG)

# Bridge functions table
add_textbox(slide, Inches(0.5), Inches(4.9), Inches(6), Inches(0.4),
            "Bridge Functions", font_size=18, bold=True, color=DARK_BG)

bridge_headers = ["Function", "From", "To", "Location"]
bridge_rows = [
    ["build_combined_expr_matrix()", "flowFrame list", "Numeric matrix", "R_shared/flowframe_utils.R"],
    ["create_flowframes_from_matrices()", "Numeric matrix", "flowFrame list", "loading/create_flowframes.R"],
    ["clusters_to_gating_csv()", "Matrix + clusters", "Per-sample CSVs", "clustering/leiden.R"],
    ["tracking_matrix_to_long()", "Wide tracking matrix", "Long data.frame", "analysis/cluster_splitting.R"],
    ["run_supercell()", "Full matrix", "Compressed matrix", "preprocessing/supercell.R"],
]
add_table(slide, Inches(0.4), Inches(5.3), Inches(7), Inches(2.0),
          bridge_headers, bridge_rows, header_color=ACCENT2, font_size=10)

# Dispatchers
add_textbox(slide, Inches(7.5), Inches(1.1), Inches(5.5), Inches(0.4),
            "3 Method Dispatchers", font_size=18, bold=True, color=DARK_BG)

add_colored_box(slide, Inches(7.5), Inches(1.55), Inches(5.5), Inches(1.5),
                "run_clustering(expr_matrix, method, params)", [
                    "Routes to: leiden | flowsom | phenograph | dbscan",
                    "Always returns: integer vector of cluster IDs",
                    "Config param: clustering_method",
                ], box_color=MOD_CLUST, text_size=11)

add_colored_box(slide, Inches(7.5), Inches(3.2), Inches(5.5), Inches(1.5),
                "run_dimred(expr_matrix, method, ...)", [
                    "Routes to: umap | tsne | diffusion_map",
                    "Returns: list(coords, method, plots, extra)",
                    "Config param: dimred_method",
                ], box_color=MOD_VIZ, text_size=11)

add_colored_box(slide, Inches(7.5), Inches(4.85), Inches(5.5), Inches(1.7),
                "run_batch_correction_dispatch(data, method, batch, ...)", [
                    "Routes to: harmony | cytonorm | combat | none",
                    "Accepts flowFrame list OR matrix (auto-converts)",
                    "Returns: list(corrected_matrix, sample_labels, ...)",
                    "Config param: batch_method",
                ], box_color=MOD_HARM, text_size=11)


# ═══════════════════════════════════════════════════════════════════════════
# SLIDE 24: DECISION TREE
# ═══════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "17. Decision Tree: Which Pipeline to Use?", font_size=28, bold=True, color=DARK_BG)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), ACCENT1)

# Root
add_flow_box(slide, Inches(4.5), Inches(1.2), Inches(4.0), Inches(0.5),
             "What is your input data?", color=DARK_BG, font_size=14)

# Level 1 branches
branches = [
    (0.3, 2.5, ".fcs files", ACCENT4),
    (3.3, 2.5, "CSV expression\n(no gates)", ACCENT1),
    (6.3, 2.5, "CSV with\ngating columns", ACCENT2),
    (9.8, 2.5, "FlowJo .wsp\nworkspaces", ACCENT3),
]
for (x, y, text, color) in branches:
    add_flow_box(slide, Inches(x), Inches(y), Inches(2.7), Inches(0.7),
                 text, color=color, font_size=12)

# Level 2 - FCS branch
add_flow_box(slide, Inches(0.3), Inches(3.7), Inches(2.7), Inches(0.5),
             "Need QC?", color=GRAY, font_size=12)
add_flow_box(slide, Inches(0.1), Inches(4.6), Inches(1.4), Inches(0.7),
             "YES\nPipeline A", color=ACCENT4, font_size=11)
add_flow_box(slide, Inches(1.7), Inches(4.6), Inches(1.4), Inches(0.7),
             "NO\nPipeline E", color=ACCENT5, font_size=11)

# Level 2 - CSV expression branch
add_flow_box(slide, Inches(3.3), Inches(3.7), Inches(2.7), Inches(0.5),
             "Need clustering?", color=GRAY, font_size=12)
add_flow_box(slide, Inches(3.1), Inches(4.6), Inches(1.4), Inches(0.7),
             "NO\nPipeline B", color=ACCENT1, font_size=11)
add_flow_box(slide, Inches(4.7), Inches(4.6), Inches(1.4), Inches(0.7),
             "YES\nB + cluster\n-> Pipeline C", color=ACCENT2, font_size=11)

# Level 2 - CSV gating
add_flow_box(slide, Inches(6.6), Inches(3.7), Inches(2.4), Inches(0.7),
             "Pipeline C\nAbundance/Enrichment", color=ACCENT2, font_size=11)

# Level 2 - FlowJo
add_flow_box(slide, Inches(9.8), Inches(3.7), Inches(2.7), Inches(0.7),
             "Pipeline D\nGating Comparison", color=ACCENT3, font_size=11)

# Bottom note
add_colored_box(slide, Inches(0.4), Inches(5.7), Inches(12.5), Inches(1.4),
                "Remember: Pipelines are examples -- you can mix and match scripts freely!", [
                    "The composability guide (docs/composability-guide.md) shows 8 custom workflow examples",
                    "Any script chain is valid if data formats match at each connection point",
                    "Bridge functions convert between formats: flowFrame <-> matrix <-> CSV <-> tracking matrix",
                    "See docs/api-reference.md for complete function signatures and parameter documentation",
                ], box_color=ACCENT1, text_size=12)


# ═══════════════════════════════════════════════════════════════════════════
# SLIDE 25: KEY TAKEAWAYS
# ═══════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, DARK_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "18. Key Takeaways", font_size=32, bold=True, color=WHITE)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), ACCENT6)

takeaways = [
    ("67 standalone scripts", "Each script has one purpose with documented inputs/outputs. "
     "Read the top-of-file annotations for I/O contracts."),
    ("5 pipeline orchestrators", "Thin example workflows. Pipeline A is most comprehensive; "
     "Pipeline B is fastest. All configurable."),
    ("3 method dispatchers", "Swap clustering (Leiden/FlowSOM/PhenoGraph/DBSCAN), "
     "dimensionality reduction (UMAP/t-SNE/diffusion map), and batch correction "
     "(Harmony/CytoNorm/ComBat) via a single config parameter."),
    ("5 data format contracts", "flowFrame list, numeric matrix, gating CSV, tracking matrix, "
     "long data.frame. Bridge functions convert between them."),
    ("Config-driven", "All parameters in config files under config/. "
     "Study-specific configs in config/study_params/. No code changes needed for different analyses."),
    ("Composable by design", "Mix and match scripts freely across pipeline boundaries. "
     "See docs/composability-guide.md for the full data flow diagram and 8 custom workflow examples."),
]

for i, (title, desc) in enumerate(takeaways):
    y = 1.2 + i * 1.0
    add_flow_box(slide, Inches(0.5), Inches(y), Inches(3.0), Inches(0.5),
                 title, color=ACCENT6, font_size=13)
    add_textbox(slide, Inches(3.8), Inches(y), Inches(9), Inches(0.85),
                desc, font_size=13, color=LIGHT_GRAY)


# ═══════════════════════════════════════════════════════════════════════════
# SLIDE 26: REFERENCE & NEXT STEPS
# ═══════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, DARK_BG)
add_textbox(slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
            "Reference & Next Steps", font_size=32, bold=True, color=WHITE)
add_shape_bg(slide, Inches(0.5), Inches(0.85), Inches(2), Inches(0.04), ACCENT6)

# Key documentation
add_textbox(slide, Inches(0.7), Inches(1.3), Inches(5.5), Inches(0.4),
            "Key Documentation", font_size=20, bold=True, color=ACCENT6)
docs = [
    "CLAUDE.md -- Project overview (root of repo)",
    "docs/composability-guide.md -- Data flow, bridge functions, custom workflows",
    "docs/methods-guide.md -- When to use each statistical method (with references)",
    "docs/api-reference.md -- Complete function signatures and parameters",
    "docs/methods/*.md -- Deep-dive per method (FlowSOM, PhenoGraph, etc.)",
    "config/*.R -- Pipeline + standalone config templates",
    "Each src/ subproject has its own CLAUDE.md and tests/ directory",
]
add_bullet_list(slide, Inches(0.9), Inches(1.8), Inches(5.5), Inches(3.5),
                docs, font_size=13, color=LIGHT_GRAY)

# Environment
add_textbox(slide, Inches(7), Inches(1.3), Inches(5.5), Inches(0.4),
            "Environment Setup", font_size=20, bold=True, color=ACCENT6)
env = [
    "conda env create -f environment.yml",
    "conda activate cytometry_env",
    "Rscript install_extra_packages.R",
]
add_bullet_list(slide, Inches(7.2), Inches(1.8), Inches(5.5), Inches(1.5),
                env, font_size=13, color=LIGHT_GRAY)

# Key packages
add_textbox(slide, Inches(7), Inches(3.3), Inches(5.5), Inches(0.4),
            "Key R Packages", font_size=20, bold=True, color=ACCENT6)
pkgs = [
    "flowCore, flowWorkspace -- FCS handling",
    "data.table -- Fast I/O",
    "lme4 -- Mixed-effects models",
    "uwot, Rtsne, destiny -- Dim. reduction",
    "ggplot2, patchwork, viridis -- Plotting",
    "FlowSOM, dbscan -- Clustering",
    "diffcyt, limma, sva -- Differential analysis",
    "survival, survminer -- Survival analysis",
    "SuperCellCyto -- Metacell compression",
]
add_bullet_list(slide, Inches(7.2), Inches(3.8), Inches(5.5), Inches(3.5),
                pkgs, font_size=12, color=LIGHT_GRAY)

# Footer
add_textbox(slide, Inches(1), Inches(6.5), Inches(11), Inches(0.5),
            "Questions?", font_size=28, bold=True, color=ACCENT6, alignment=PP_ALIGN.CENTER)


# ═══════════════════════════════════════════════════════════════════════════
# SAVE
# ═══════════════════════════════════════════════════════════════════════════
output_path = "/mnt/volume/nacho/cytometry/Cytometry_Pipeline_Handover.pptx"
prs.save(output_path)
print(f"Presentation saved to {output_path}")
print(f"Total slides: {len(prs.slides)}")
