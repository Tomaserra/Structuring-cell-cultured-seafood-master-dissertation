# =============================================================================
#  MyHC2 quantification (C2C12, DAPI / phalloidin / MyHC2) - plots and statistics
#   1  MyHC2-positive area (%)                       bars        (image level)
#   2  Fusion index (%)                              bars        (image level)
#   3  Myotube morphology: length, width, aspect ratio  bars + individual myotubes
#   4  Alignment A (myotubes):  % of myotubes within +-15 deg of the reference, mean angular deviation  + rose diagrams
#   5  Alignment B (structure tensor, OrientationJ-like): order parameter S, alignment score (+-10 deg)
#  Every metric gets one figure per platform (culture plate / seeded on hydrogels / encapsulated) and one with all platforms.
#  Style = myogenin script: ink colours of the characterization script (same colour for controls and ES of one ink),
#  ES bars with fluorescent-yellow border + stripes, detailed captions inside the plots, brackets only for significant results.
#
#  INPUT
#   (1) the combined summary CSV (all ALL_images_summary.csv of the macro put together, like the myogenin sheet:
#       one block per platform, each with its own header row, first column = image names).
#   (2) the macro output folders, with one sub-folder per image:  <image>/<image>_MyHC2_objects.csv  (per-myotube values)
#       Leave these files where the macro wrote them - do NOT combine them: the script finds them by itself
#       (give the folder(s) in 'results_dirs'). The *_nuclei_assignment.csv files are not needed for the plots.
#  OUTPUT  <out_dir>/<metric folder>/ figures (PDF + PNG)  and  <out_dir>/data/  (tables, captions)
# =============================================================================

# -----------------------------------------------------------------------------
# 0. OPTIONS
# -----------------------------------------------------------------------------
mhc2_root <- "/Users/tomas_serra/Library/CloudStorage/GoogleDrive-tomasfserra1@gmail.com/My Drive/Tese Drive/Current Thesis Files/C2C12/Immunostainings/Quantifications/MHC2"

summary_csv  <- file.path(mhc2_root, "Cell quantifications - Immuno DAPI + Phalloidin + MHC2 .csv")   # combined summary (NA = choose in a window)
#  folders that contain the macro results (one sub-folder per image, with <image>_MyHC2_objects.csv inside).
#  The script searches these folders recursively, so the name of the results folder inside each one does not matter.
results_dirs <- file.path(mhc2_root, c("SeedPlate-MHC2", "SeedHydro-MHC2", "MHC2-EN30k-c13.3MC0Gel"))
out_dir      <- file.path(mhc2_root, "MHC2_QuantPlots")      # all figures and tables are saved here

density_plate    <- "2,000 cells/cm\u00B2"
density_hydrogel <- "200,000 cells/cm\u00B2"
es_unit <- "V"
myog_medium <- "PM3"      # medium of the hydrogels (colour looked up as ink + medium)

# --- ink colours (copied from the characterization script) -------------------------------------------
ink_table <- data.frame(
  medium = c("PM3", "PM3", "PM3", "PM3", "PM3", "PM3", "PM3",
             "L15-PM", "L15-PM", "L15-PM", "L15-PM", "DMEM-PM", "DMEM-PM", "MACK-PM"),
  ink    = c("30k-c13.3MC0Gel", "30k-c13.3MC100Gel", "72k-c40MC0Gel", "90k-c40MC0Gel", "90k-c40MC100Gel", "90k-c40MC200Gel", "90k-c40MC375Gel",
             "90k-c", "72k-c40MC0Gel", "72k-c40MC100Gel", "72k-c40MC375Gel", "72k-c40MC0Gel", "72k-c40MC100Gel", "72k-c40MC0Gel"),
  color  = c("#ffcad8", "#ffae9c", "#fa7b8e", "#9f0019", "#e43954", "#ff742f", "#ffc14a",
             "#003b74", "#43b6ff", "#009692", "#b7e8d8", "#a600a5", "#eb67ea", "#afc40b"),
  stringsAsFactors = FALSE)
plate_color <- "#cdcdcd"
ink_rx <- "[0-9]+k-c[0-9._]+MC[0-9]+Gel"

# --- electrical stimulation (ES): fluorescent-yellow border + stripes -------------------------------------
es_yellow   <- "#fffd3c"
es_border_w <- 1
es_stripe_w <- 0.6
es_stripe_d <- 0.4

# --- statistics ----------------------------------------------------------------------------------------
t_var_equal <- TRUE            # Student's t-test (image-level metrics: n = 2 images per condition)
morph_test  <- "t-test"      # test for the individual-myotube metrics (morphology): "wilcoxon" (Mann-Whitney) or "t-test"
min_myotubes_test <- 3         # a condition needs at least this many myotubes to be tested
#  Each ES condition is compared with the expansion control and with the differentiation control of its own panel.
#  Significance: ns p >= 0.05; * 0.01-0.05; ** 0.001-0.01; *** 0.0001-0.001; **** < 0.0001
#  Non-significant comparisons get NO bracket on the plot; they are listed in the caption.
run_tests <- TRUE

fig_height <- 9
base_size  <- 11
sd_alpha   <- 0.4           # translucency of the range lines / SD bars
save_pdf   <- FALSE
make_platform_figures <- TRUE  # (only used when make_bar_figures = TRUE)

#  WHAT IS PLOTTED (the numbers of the other metrics go into the tables of the thesis: lengths, widths, aspect ratios, % parallel...)
make_bar_figures       <- FALSE  # the old bar figure of every metric (length, width, aspect, % parallel, S, ...). OFF: values are in tables
make_summary_figures   <- TRUE   # NEW: 1_Area_and_fusion (all platforms, one figure) + 2_Alignment_summary (one figure per platform and ink)
plot_rose_diagrams     <- TRUE   # rose diagrams: ONE figure per platform, all inks and conditions in a grid (angles are pooled, not averaged)
make_morphology_glyphs <- FALSE  # to-scale rectangles of the myotubes. OFF
make_histograms        <- FALSE  # orientation histograms. OFF
rose_bin   <- 10               # bin width (deg) of the rose diagrams

#  DIRECTION OF THE ELECTRIC FIELD in the images, in degrees counter-clockwise from the horizontal axis of the image.
#  Your C-shaped electrodes "( cells )" face each other left/right, so the field runs between them horizontally -> 0.
#  Only valid if the images were saved with the electrodes on the left and right (check this on one image of an ES plate).
#  NA = use the reference angle that the macro stored (0 unless you changed it in the macro dialog).
field_direction_deg <- NA

# -----------------------------------------------------------------------------
# 1. PACKAGES AND READING
# -----------------------------------------------------------------------------
need <- function(p) if (!requireNamespace(p, quietly = TRUE)) stop(sprintf('Package "%s" is missing: install.packages("%s")', p, p), call. = FALSE)
invisible(lapply(c("readr", "dplyr", "tidyr", "stringr", "ggplot2", "ggprism", "patchwork", "tibble"), need))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(stringr); library(ggplot2); library(ggprism); library(patchwork) })

if (is.na(summary_csv) || !file.exists(summary_csv)) summary_csv <- file.choose()
if (is.na(out_dir)) out_dir <- file.path(dirname(normalizePath(summary_csv)), "MyHC2 results")
dir.create(file.path(out_dir, "data"), recursive = TRUE, showWarnings = FALSE)
message("Summary CSV: ", summary_csv, "\nOutput:      ", out_dir)

first_line <- readLines(summary_csv, n = 1, warn = FALSE, encoding = "UTF-8")
delim <- if (str_count(first_line, ";") > str_count(first_line, ",")) ";" else ","
raw0 <- readr::read_delim(summary_csv, delim = delim, col_types = readr::cols(.default = readr::col_character()),
                          trim_ws = TRUE, show_col_types = FALSE, name_repair = "minimal")
to_num <- function(x) {
  x <- str_replace_all(as.character(x), "\\s", "")
  x <- ifelse(str_detect(x, ",") & !str_detect(x, "\\."), str_replace(x, ",", "."), x)
  suppressWarnings(as.numeric(x))
}
col_total <- "Total_DAPI_nuclei"
if (!col_total %in% names(raw0)) stop('Column "', col_total, '" not found. Columns are: ', paste(names(raw0), collapse = ", "))

# first column = image names; its header cell is the title of the FIRST block; every further block has its own title + header row
block_first <- names(raw0)[1]

names(raw0)[1] <- "Image"
raw0$row_id <- seq_len(nrow(raw0))
is_hdr <- dplyr::coalesce(raw0[[col_total]] == col_total, FALSE)
block <- character(nrow(raw0)); cur <- block_first
for (i in seq_len(nrow(raw0))) { if (is_hdr[i]) cur <- raw0$Image[i]; block[i] <- cur }
raw0$block <- block
raw <- raw0[!is_hdr & !is.na(raw0$Image) & nzchar(str_trim(raw0$Image)), ]
message("Blocks found: ", paste(unique(raw$block), collapse = " | "), "\nImages read: ", nrow(raw))
if (anyDuplicated(raw$Image)) message("!! Duplicated image names in the summary: ", paste(unique(raw$Image[duplicated(raw$Image)]), collapse = ", "))

# -----------------------------------------------------------------------------
# 2. PLATFORM, INK, CELL DENSITY AND CONDITION (block title + image name)
# -----------------------------------------------------------------------------
rx_encaps <- "(^|[_ .-])(10|25)M([_ .-]|$)"
ink_names <- ink_table$ink[ink_table$medium == myog_medium]
ink_cols  <- c(setNames(ink_table$color[ink_table$medium == myog_medium], ink_names), "Culture plate (no hydrogel)" = plate_color)
get_ink <- function(x) {
  tok <- str_replace_all(str_extract(x, regex(ink_rx, ignore_case = TRUE)), "_", ".")
  hit <- ink_names[match(tolower(tok), tolower(ink_names))]
  ifelse(is.na(hit), NA_character_, hit)
}
parse_meta <- function(x, blk) {
  volt <- str_match(x, "(?<![0-9A-Za-z.])([0-9]+(?:[._][0-9]+)?)\\s*V(?![A-Za-z])")[, 2]
  is_enc <- str_detect(x, rx_encaps)
  plat <- case_when(str_detect(blk, regex("encaps", ignore_case = TRUE)) ~ "Encapsulated",
                    str_detect(blk, regex("hydro", ignore_case = TRUE)) ~ "Hydrogel-seeded",
                    str_detect(blk, regex("plate", ignore_case = TRUE)) ~ "Plate",
                    is_enc ~ "Encapsulated", str_detect(x, regex("SeedPlate", ignore_case = TRUE)) ~ "Plate", TRUE ~ "Hydrogel-seeded")
  tibble(platform = plat,
         cells = ifelse(plat == "Encapsulated", paste0(str_match(x, rx_encaps)[, 3], "M cells/mL"), NA_character_),
         condition = case_when(
           str_detect(x, regex("Control[_ .-]*Exp", ignore_case = TRUE)) ~ "Control (expansion)",
           str_detect(x, regex("Control[_ .-]*Diff", ignore_case = TRUE)) ~ "Control (differentiation)",
           !is.na(volt) ~ paste0("ES ", str_replace(volt, "_", "."), " ", es_unit),
           TRUE ~ NA_character_),
         ink_found = str_extract(x, regex(ink_rx, ignore_case = TRUE)))
}
meta <- parse_meta(raw$Image, raw$block)
meta$ink <- ifelse(meta$platform == "Plate", "Culture plate (no hydrogel)", get_ink(raw$Image))
chk_tab <- bind_cols(tibble(row_id = raw$row_id, Image = raw$Image), meta)
problems <- chk_tab %>% filter(is.na(condition) | (platform != "Plate" & is.na(ink)) | (platform == "Encapsulated" & is.na(cells)))
if (nrow(problems)) {
  print(as.data.frame(problems %>% select(row_id, Image, platform, ink_found, ink, cells, condition)))
  stop("The rows above could not be fully understood (ink not in 'ink_table', or no ControlExpan / ControlDiff / '<n>V' / 10M / 25M in the name).")
}

# numeric columns of the summary that the plots use
metric_cols <- c("MyHC2_area_fraction_pct", "Fusion_index_pct", "Pct_within_15deg", "Pct_within_10deg", "Mean_abs_angle_dev_deg",
                 "Order_parameter_S_objects", "ST_order_parameter_S", "ST_alignment_score_10deg_pct", "ST_alignment_score_15deg_pct",
                 "N_myotubes", "N_MyHC2_objects", "Total_DAPI_nuclei", "Reference_angle_deg", "Myotube_min_nuclei", "Alignment_min_aspect_ratio",
                 "Mean_myotube_length_um", "Mean_myotube_width_um", "Mean_myotube_aspect_ratio", "ST_nematic_magnitude_R", "ST_dominant_angle_deg")
for (cc in setdiff(metric_cols, names(raw))) message("Column not in the summary (that metric will be skipped): ", cc)
img <- tibble(row_id = raw$row_id, Image = raw$Image) %>% bind_cols(meta %>% select(platform, cells, condition, ink))
for (cc in intersect(metric_cols, names(raw))) img[[cc]] <- to_num(raw[[cc]])
img <- img %>% mutate(platform = factor(platform, levels = c("Plate", "Hydrogel-seeded", "Encapsulated")))
volt_levels <- unique(img$condition[str_detect(img$condition, "^ES ")])
volt_levels <- volt_levels[order(as.numeric(str_match(volt_levels, "^ES ([0-9.]+)")[, 2]))]
cond_levels <- c("Control (expansion)", "Control (differentiation)", volt_levels)
cond_short <- setNames(c("Expansion control", "Differentiation control", volt_levels), cond_levels)
img <- img %>% mutate(condition = factor(condition, levels = cond_levels),
                      cond_group = case_when(condition == "Control (expansion)" ~ "expansion", condition == "Control (differentiation)" ~ "diff", TRUE ~ "ES"),
                      panel = case_when(platform == "Plate" ~ "Culture plate", platform == "Hydrogel-seeded" ~ sub("MC", "\nMC", ink),
                                        TRUE ~ paste0(cells, "\n", sub("MC", "\nMC", ink))),
                      panel_plain = case_when(platform == "Plate" ~ "Culture plate", platform == "Hydrogel-seeded" ~ ink, TRUE ~ paste0(cells, " \u00B7 ", ink)),
                      group_id = paste(platform, ink, coalesce(cells, "-"), sep = "|"))
message("\nWhat the script understood (check the ink column):")
print(as.data.frame(img %>% select(row_id, platform, ink, cells, condition, Image) %>% arrange(platform, ink, cells, condition)), row.names = FALSE)
nt <- img %>% count(group_id, condition) %>% filter(n != 2)
if (nrow(nt)) { message("\n!! Conditions with n != 2 images:"); print(as.data.frame(nt)) } else message("\nEvery condition has n = 2 images.")

# data check: do the two images of a condition really differ? (the structure tensor depends only on the MyHC2 channel:
# identical values in both images mean the same image was analysed twice)
if (all(c("ST_order_parameter_S", "Total_DAPI_nuclei") %in% names(img))) {
  dchk <- img %>% filter(!is.na(ST_order_parameter_S)) %>% group_by(group_id, condition) %>%
    filter(n() == 2) %>% summarise(Image_1 = Image[1], Image_2 = Image[2], same_structure_tensor = n_distinct(round(ST_order_parameter_S, 6)) == 1,
                                   nuclei_1 = Total_DAPI_nuclei[1], nuclei_2 = Total_DAPI_nuclei[2], .groups = "drop")
  if (nrow(dchk) && any(dchk$same_structure_tensor)) {
    message("\n!! ", sum(dchk$same_structure_tensor), " of ", nrow(dchk), " conditions have two images with IDENTICAL structure-tensor values. ",
            "That is only possible if both images share the same MyHC2 channel (same field analysed twice, probably with different thresholds):")
    print(as.data.frame(dchk %>% filter(same_structure_tensor) %>% select(Image_1, Image_2, nuclei_1, nuclei_2)), row.names = FALSE)
  }
  readr::write_csv(dchk %>% select(-group_id), file.path(out_dir, "data", "MyHC2_data_check_duplicate_images.csv"))
}

first_val <- function(col, default) { v <- unique(na.omit(img[[col]])); if (length(v)) v[1] else default }
ref_deg   <- first_val("Reference_angle_deg", 0)
min_nuc   <- first_val("Myotube_min_nuclei", 2)
min_ar    <- first_val("Alignment_min_aspect_ratio", 2)
if (length(unique(na.omit(img$Reference_angle_deg))) > 1) message("!! Different reference angles in the summary: ", paste(unique(na.omit(img$Reference_angle_deg)), collapse = ", "))
field_deg <- if (is.na(field_direction_deg)) ref_deg else field_direction_deg
ref_txt <- if (field_deg == 0) "the horizontal axis of the image" else paste0("a direction ", field_deg, "\u00B0 counter-clockwise from the horizontal axis of the image")
if (field_deg != ref_deg) message("Field direction ", field_deg, " deg differs from the macro's reference (", ref_deg, " deg): angle-based metrics are recomputed from the per-myotube angles; the alignment score +-10 deg (macro) cannot be recomputed and is skipped.")

# -----------------------------------------------------------------------------
# 3. PER-MYOTUBE FILES (found automatically in the macro output folders)
# -----------------------------------------------------------------------------
search_dirs <- if (all(is.na(results_dirs))) dirname(normalizePath(summary_csv)) else results_dirs
obj_files <- unique(unlist(lapply(search_dirs, function(dd) list.files(dd, pattern = "_MyHC2_objects\\.csv$", recursive = TRUE, full.names = TRUE))))
obj <- NULL
if (length(obj_files)) {
  obj <- bind_rows(lapply(obj_files, function(f) {
    o <- readr::read_csv(f, col_types = readr::cols(.default = readr::col_character()), show_col_types = FALSE, name_repair = "minimal")
    if (!"Image" %in% names(o)) o$Image <- sub("_MyHC2_objects\\.csv$", "", basename(f))
    o
  }))
  for (cc in c("Area_um2", "Feret_length_um", "MinFeret_width_um", "AspectRatio_Feret", "Angle_deg", "N_nuclei", "Is_myotube", "Used_for_alignment"))
    if (cc %in% names(obj)) obj[[cc]] <- to_num(obj[[cc]])
  miss <- setdiff(unique(obj$Image), img$Image)
  if (length(miss)) message("!! Per-myotube files of images that are not in the summary (ignored): ", paste(miss, collapse = ", "))
  nofile <- setdiff(img$Image, unique(obj$Image))
  if (length(nofile)) message("Images of the summary without any myotube row (no myotubes, or file not found): ", paste(nofile, collapse = ", "))
  dimg <- img %>% distinct(Image, .keep_all = TRUE) %>% select(Image, platform, cells, condition, cond_group, ink, panel, panel_plain, group_id, Reference_angle_deg)
  obj <- obj %>% inner_join(dimg, by = "Image")
  message("Per-myotube files read: ", length(obj_files), " (", nrow(obj), " objects, ", sum(obj$Is_myotube == 1, na.rm = TRUE), " myotubes)")
} else message("!! No *_MyHC2_objects.csv found in: ", paste(search_dirs, collapse = " | "), "\n   Set 'results_dirs' (section 0). Morphology and the rose diagrams are skipped; the other plots still work.")
obj_myo <- if (!is.null(obj)) obj %>% filter(Is_myotube == 1) else NULL

# ---- alignment relative to the ELECTRIC FIELD: parallel (0 deg) and perpendicular (90 deg) --------------------------------
rel_to_field <- function(a) ((a - field_deg + 90) %% 180) - 90      # -90..90: 0 = parallel to the field, +-90 = perpendicular to it
if (!is.null(obj)) {
  al <- obj %>% filter(Used_for_alignment == 1, !is.na(Angle_deg)) %>% mutate(rel = rel_to_field(Angle_deg)) %>% group_by(Image) %>%
    summarise(Pct_parallel15 = 100 * mean(abs(rel) <= 15), Pct_perp15 = 100 * mean(abs(rel) >= 75), Mean_angle_to_field = mean(abs(rel)), .groups = "drop")
  img <- img %>% left_join(al, by = "Image")
} else {
  img$Pct_parallel15 <- if ("Pct_within_15deg" %in% names(img)) img$Pct_within_15deg else NA_real_
  img$Pct_perp15 <- NA_real_
  img$Mean_angle_to_field <- if ("Mean_abs_angle_dev_deg" %in% names(img)) img$Mean_abs_angle_dev_deg else NA_real_
}
# structure-tensor order parameter relative to the field: S = R cos 2(theta_dominant - theta_field)  (= the macro's S when field = macro reference)
img$ST_S_field <- if (all(c("ST_nematic_magnitude_R", "ST_dominant_angle_deg") %in% names(img)))
  img$ST_nematic_magnitude_R * cos(2 * (img$ST_dominant_angle_deg - field_deg) * pi / 180) else img$ST_order_parameter_S
if (field_deg != ref_deg && "ST_alignment_score_10deg_pct" %in% names(img)) img$ST_alignment_score_10deg_pct <- NA_real_

# -----------------------------------------------------------------------------
# 4. GENERIC MACHINERY (statistics, brackets, bars, captions)
# -----------------------------------------------------------------------------
sig_symbol <- function(p) {
  if (is.na(p)) "NA" else if (p < 0.0001) "****" else if (p < 0.001) "***" else
    if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}
ctrl_short <- c("Control (differentiation)" = "differentiation control", "Control (expansion)" = "expansion control")
base_theme <- function() {
  ggprism::theme_prism(base_size = base_size, base_fontface = "plain") +
    theme(plot.title = element_text(face = "bold"),
          plot.caption = element_text(colour = "grey40", size = base_size - 2, hjust = 0, face = "plain", margin = margin(t = 12)),
          plot.caption.position = "plot", strip.background = element_blank(), strip.text = element_text(face = "bold"))
}
wrap_parts <- function(parts, w_cm, chars_per_cm = 5) {
  paste(vapply(parts, function(p) paste(strwrap(p, width = max(25, floor(w_cm * chars_per_cm))), collapse = "\n"), ""), collapse = "\n")
}
cap_theme <- function() theme(plot.caption = element_text(colour = "grey40", size = base_size - 2, hjust = 0, face = "plain", margin = margin(t = 12)), plot.caption.position = "plot")
# FIGURE CAPTIONS: set draw_caption <- TRUE to print captions inside the figures (default: not drawn; the LaTeX caption carries the text,
# and the caption texts are still written to the *_figure_captions.txt file).
draw_caption <- FALSE
n_lines <- function(txt) length(strsplit(txt, "\n", fixed = TRUE)[[1]])
cap_height <- function(cap) if (draw_caption) 0.38 * n_lines(cap) + 0.6 else 0
cap_draw <- function(cap) if (draw_caption) cap else NULL
data_path <- function(name) file.path(out_dir, "data", name)
save_fig <- function(p, folder, name, w, h) {
  dd <- file.path(out_dir, folder); dir.create(dd, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(dd, name)
  ggsave(paste0(f, ".png"), p, width = w, height = h, units = "cm", dpi = 300, bg = "white")
  if (save_pdf) ggsave(paste0(f, ".pdf"), p, width = w, height = h, units = "cm", bg = "white")
}

# metric specifications ---------------------------------------------------------------------------------
#  level "image": value = column of the summary (n = 2 images per condition; bars = mean, circles = images, line = range)
#  level "myotube": value = column of the per-myotube files (bars = mean +- SD of the myotubes, small dots = individual myotubes)
specs <- list(
  list(id = "area", folder = "1_MyHC2_area", level = "image", col = "MyHC2_area_fraction_pct", ylab = "MyHC2-positive area (%)",
       what = "MyHC2-positive area (%) = total area of the MyHC2-positive objects (thresholded, cleaned and outlined MyHC2 mask) / image area \u00D7 100."),
  list(id = "fusion", folder = "2_Fusion_index", level = "image", col = "Fusion_index_pct", ylab = "Fusion index (%)",
       what = paste0("Fusion index (%) = DAPI nuclei inside MyHC2-positive myotubes (MyHC2-positive objects containing at least ", min_nuc, " nuclei) / total DAPI nuclei \u00D7 100.")),
  list(id = "length", folder = "3_Morphology", level = "myotube", col = "Feret_length_um", ylab = "Myotube length (\u00B5m)",
       what = paste0("Myotube length = maximum Feret diameter of each myotube (MyHC2-positive object with at least ", min_nuc, " nuclei).")),
  list(id = "width", folder = "3_Morphology", level = "myotube", col = "MinFeret_width_um", ylab = "Myotube width (\u00B5m)",
       what = paste0("Myotube width = minimum Feret diameter of each myotube (MyHC2-positive object with at least ", min_nuc, " nuclei).")),
  list(id = "aspect", folder = "3_Morphology", level = "myotube", col = "AspectRatio_Feret", ylab = "Myotube aspect ratio",
       what = paste0("Aspect ratio = Feret length / minimum Feret width of each myotube (MyHC2-positive object with at least ", min_nuc, " nuclei).")),
  list(id = "within15", folder = "4_Alignment_A_myotubes", level = "image", col = "Pct_parallel15", ylab = "Myotubes parallel to the field (\u00B115\u00B0) (%)",
       what = paste0("Alignment A: percentage of myotubes (aspect ratio \u2265 ", min_ar, ") whose long axis lies within \u00B115\u00B0 of the direction of the electric field (", ref_txt, "), i.e. oriented along the field and perpendicular to the electrodes.")),
  list(id = "perp15", folder = "4_Alignment_A_myotubes", level = "image", col = "Pct_perp15", ylab = "Myotubes perpendicular to the field (\u00B115\u00B0) (%)",
       what = paste0("Alignment A: percentage of myotubes (aspect ratio \u2265 ", min_ar, ") whose long axis lies within \u00B115\u00B0 of the direction perpendicular to the electric field (field direction: ", ref_txt, "), i.e. parallel to the electrodes.")),
  list(id = "meandev", folder = "4_Alignment_A_myotubes", level = "image", col = "Mean_angle_to_field", ylab = "Mean angle to the field (\u00B0)",
       what = paste0("Alignment A: mean absolute angle between the long axis of the myotubes (aspect ratio \u2265 ", min_ar, ") and the electric-field direction (", ref_txt, "); 0\u00B0 = parallel to the field, 90\u00B0 = perpendicular, about 45\u00B0 = random orientation.")),
  list(id = "S", folder = "5_Alignment_B_structure_tensor", level = "image", col = "ST_S_field", ylab = "Order parameter S (+1 parallel, \u22121 perpendicular)",
       what = paste0("Alignment B (structure tensor of the MyHC2 channel, as computed by OrientationJ): order parameter S = mean cos 2(\u03B8 \u2212 \u03B8field) of the oriented signal, with \u03B8field = the electric-field direction (", ref_txt, "); +1 = all signal parallel to the field, 0 = random, \u22121 = all perpendicular to the field.")),
  list(id = "score10", folder = "5_Alignment_B_structure_tensor", level = "image", col = "ST_alignment_score_10deg_pct", ylab = "Alignment score \u00B110\u00B0 of the field (%)",
       what = paste0("Alignment B (structure tensor of the MyHC2 channel): percentage of the oriented signal within \u00B110\u00B0 of the electric-field direction (", ref_txt, ").")))

# one spec -> values, summaries, tests, brackets ---------------------------------------------------------
prep_spec <- function(sp) {
  if (sp$level == "image") {
    if (!sp$col %in% names(img)) { message("Skipped (column missing): ", sp$id); return(NULL) }
    v <- img %>% filter(!is.na(.data[[sp$col]])) %>% mutate(value = .data[[sp$col]])
  } else {
    if (is.null(obj_myo) || !sp$col %in% names(obj_myo)) { message("Skipped (no per-myotube data): ", sp$id); return(NULL) }
    v <- obj_myo %>% filter(!is.na(.data[[sp$col]])) %>% mutate(value = .data[[sp$col]])
  }
  if (!nrow(v)) { message("Skipped (no values): ", sp$id); return(NULL) }
  v <- v %>% mutate(condition = factor(as.character(condition), levels = cond_levels))
  sm <- v %>% group_by(platform, ink, cells, panel, panel_plain, group_id, condition, cond_group) %>%
    summarise(n = n(), mean = mean(value), sd = ifelse(n() > 1, sd(value), 0), lo = min(value), hi = max(value), .groups = "drop") %>%
    arrange(platform, ink, cells, condition) %>%
    group_by(group_id) %>% mutate(pos = match(as.character(condition), cond_levels[cond_levels %in% as.character(condition)])) %>% ungroup()
  sm$bar_lo <- if (sp$level == "image") sm$lo else sm$mean - sm$sd
  sm$bar_hi <- if (sp$level == "image") sm$hi else sm$mean + sm$sd
  # tests: each ES condition vs each control of the same panel
  tests <- NULL
  if (run_tests) {
    rows <- list(); min_n <- if (sp$level == "image") 2 else min_myotubes_test
    use_w <- sp$level == "myotube" && morph_test == "wilcoxon"
    for (g in unique(v$group_id)) {
      vg <- v[v$group_id == g, ]
      for (e in unique(as.character(vg$condition[vg$cond_group == "ES"]))) for (cn in c("Control (differentiation)", "Control (expansion)")) {
        a <- vg$value[as.character(vg$condition) == e]; b <- vg$value[as.character(vg$condition) == cn]
        if (length(a) < min_n || length(b) < min_n) next
        tt <- tryCatch(if (use_w) wilcox.test(a, b, exact = FALSE) else t.test(a, b, var.equal = t_var_equal), error = function(err) NULL)
        if (is.null(tt) || is.na(tt$p.value)) next
        rows[[length(rows) + 1]] <- data.frame(metric = sp$id, platform = as.character(vg$platform[1]), panel = vg$panel[1], panel_plain = vg$panel_plain[1], group_id = g,
                                               es = e, control = cn, n_es = length(a), n_control = length(b), mean_es = mean(a), mean_control = mean(b),
                                               difference = mean(a) - mean(b), test = if (use_w) "Wilcoxon rank-sum" else if (t_var_equal) "Student t" else "Welch t",
                                               p_value = tt$p.value, significance = sig_symbol(tt$p.value), stringsAsFactors = FALSE)
      }
    }
    if (length(rows)) tests <- do.call(rbind, rows)
  }
  scale_y <- max(c(sm$bar_hi, sm$hi, 0), na.rm = TRUE) - min(c(sm$bar_lo, sm$lo, 0), na.rm = TRUE)
  if (!is.finite(scale_y) || scale_y <= 0) scale_y <- 1
  brk <- NULL
  if (!is.null(tests) && any(tests$significance != "ns")) {
    sg <- tests[tests$significance != "ns", ]; out <- list()
    for (g in unique(sg$group_id)) {
      s <- sg[sg$group_id == g, ]
      present <- cond_levels[cond_levels %in% as.character(v$condition[v$group_id == g])]
      top <- max(c(v$value[v$group_id == g], sm$bar_hi[sm$group_id == g]))
      s$x1 <- match(s$es, present); s$x2 <- match(s$control, present); s$span <- abs(s$x1 - s$x2)
      s <- s[order(s$span), ]
      s$y <- top + 0.05 * scale_y + (seq_len(nrow(s)) - 1) * 0.085 * scale_y
      out[[length(out) + 1]] <- s
    }
    brk <- do.call(rbind, out)
  }
  y_max <- max(c(sm$bar_hi, sm$hi, if (!is.null(brk)) brk$y + 0.09 * scale_y), na.rm = TRUE) * 1.04
  y_min <- min(0, min(c(sm$bar_lo, sm$lo), na.rm = TRUE) * 1.08)
  if (sp$id %in% c("area", "fusion", "within15", "perp15", "score10")) y_max <- max(y_max, 0)
  list(spec = sp, v = v, sm = sm, tests = tests, brk = brk, y_min = y_min, y_max = y_max, scale_y = scale_y)
}

ns_caption <- function(tb) {
  if (is.null(tb) || !nrow(tb)) return(NULL)
  if (!any(tb$significance == "ns")) return(NULL)
  "Other comparisons were not significant (p \u2265 0.05)."
}
stat_legend <- function(sp) {
  if (!run_tests) return(NULL)
  paste0("Each ES condition vs each control: Student's t-test; all conditions of a panel: one-way ANOVA (data/MyHC2_anova_by_panel.csv). ",
         "* p < 0.05; ** p < 0.01; *** p < 0.001; **** p < 0.0001.")
}
bars_text <- function(sp) {
  if (sp$level == "image") "Mean of two images (circles: images; line: range)."
  else "Mean \u00B1 SD of myotubes (dots: myotubes; n under each label)."
}
es_text <- "Yellow border and stripes: ES; plain bars: controls."
rep_desc <- c(Plate = "culture plates (2D, ", `Hydrogel-seeded` = "the top of the hydrogels (", Encapsulated = "")
cap_where <- function(pf) switch(pf,
                                 Plate = paste0("culture plates (", density_plate, ")."),
                                 `Hydrogel-seeded` = paste0("cells seeded on hydrogels (", density_hydrogel, ")."),
                                 Encapsulated = "cells encapsulated (10 or 25 million cells/mL).")
caption_for <- function(S, pf, w) {
  tb <- if (!is.null(S$tests)) S$tests[S$tests$platform == pf, ] else NULL
  wrap_parts(c(paste0(S$spec$ylab, ": ", cap_where(pf)), bars_text(S$spec), es_text, stat_legend(S$spec), ns_caption(tb)), w)
}

# plotting -----------------------------------------------------------------------------------------------
stripe_segments <- function(sg, k) {
  out <- list()
  for (r in seq_len(nrow(sg))) {
    h <- sg$mean[r]; x0 <- sg$pos[r] - 0.36; x1 <- sg$pos[r] + 0.36
    if (!is.finite(h) || abs(h) < 1e-9) next
    ylo <- min(0, h); yhi <- max(0, h)
    for (cc in seq(x0 - yhi / k, x1 - ylo / k, by = es_stripe_d)) {
      xs <- max(x0, cc + ylo / k); xe <- min(x1, cc + yhi / k)
      if (xs < xe) out[[length(out) + 1]] <- data.frame(panel = sg$panel[r], xs = xs, xe = xe, ys = k * (xs - cc), ye = k * (xe - cc))
    }
  }
  if (length(out)) do.call(rbind, out) else NULL
}
metric_bar_plot <- function(S, pf, fig_w, fig_h, show_y = TRUE, x_units_total = NULL, cap_h = 0, y_lim = NULL) {
  sp <- S$spec
  s <- S$sm %>% filter(platform == pf); i <- S$v %>% filter(platform == pf)
  xlab_of <- function(df) {
    lab <- cond_short[as.character(df$condition)]
    if (sp$level == "myotube") lab <- paste0(lab, "\n(n = ", df$n, ")")
    lab
  }
  s$x_lab <- xlab_of(s)
  lev <- unique(s$x_lab[order(as.integer(s$condition), s$x_lab)])
  s$x <- factor(s$x_lab, levels = lev)
  if (sp$level == "myotube") {
    key <- s %>% select(group_id, condition, x_lab)
    i <- i %>% left_join(key, by = c("group_id", "condition")); i$x <- factor(i$x_lab, levels = lev)
  } else i$x <- factor(cond_short[as.character(i$condition)], levels = lev)
  s$panel <- factor(s$panel, levels = unique(s$panel)); i$panel <- factor(i$panel, levels = levels(s$panel))
  n_pan <- nlevels(s$panel)
  if (is.null(x_units_total)) x_units_total <- nrow(s) + 1.2 * n_pan
  y_lim <- if (is.null(y_lim)) c(S$y_min, S$y_max) else y_lim
  panel_w <- max(3, fig_w - 2.2 - 0.5 * n_pan); panel_h <- max(3, fig_h - cap_h - 5.9)
  k <- (diff(y_lim) / panel_h) * (panel_w / x_units_total)
  s_es <- s %>% filter(cond_group == "ES"); s_ct <- s %>% filter(cond_group != "ES")
  p <- ggplot(s, aes(x = x)) + geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_col(data = s_ct, aes(y = mean, fill = ink), colour = "black", linewidth = 0.4, width = 0.72) +
    geom_col(data = s_es, aes(y = mean, fill = ink), colour = es_yellow, linewidth = es_border_w, width = 0.72)
  st <- stripe_segments(as.data.frame(s_es), k)
  if (!is.null(st)) {
    st$panel <- factor(st$panel, levels = levels(s$panel))
    p <- p + geom_segment(data = st, aes(x = xs, xend = xe, y = ys, yend = ye), colour = es_yellow, linewidth = es_stripe_w, inherit.aes = FALSE, lineend = "butt")
  }
  if (sp$level == "image") {
    p <- p + geom_linerange(aes(ymin = bar_lo, ymax = bar_hi), linewidth = 0.5, colour = "black", alpha = sd_alpha) +
      geom_point(data = i, aes(y = value), shape = 21, fill = "white", colour = "black", size = 2.4, stroke = 0.6, position = position_jitter(width = 0.07, height = 0, seed = 7))
  } else {
    p <- p + geom_point(data = i, aes(y = value), shape = 16, colour = "grey25", size = 0.9, alpha = 0.35, position = position_jitter(width = 0.2, height = 0, seed = 7)) +
      geom_errorbar(aes(ymin = bar_lo, ymax = bar_hi), width = 0.25, linewidth = 0.5, colour = "black", alpha = sd_alpha)
  }
  if (!is.null(S$brk)) {
    b <- S$brk %>% filter(platform == pf); b$panel <- factor(b$panel, levels = levels(s$panel))
    tick <- 0.02 * S$scale_y
    if (nrow(b)) p <- p +
      geom_segment(data = b, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
      geom_segment(data = b, aes(x = x1, xend = x1, y = y, yend = y - tick), inherit.aes = FALSE, linewidth = 0.4) +
      geom_segment(data = b, aes(x = x2, xend = x2, y = y, yend = y - tick), inherit.aes = FALSE, linewidth = 0.4) +
      geom_text(data = b, aes(x = (x1 + x2) / 2, y = y, label = significance), inherit.aes = FALSE, vjust = -0.2, size = 3.5)
  }
  p + scale_fill_manual(values = ink_cols, guide = "none") +
    scale_y_continuous(expand = c(0, 0), breaks = scales::pretty_breaks(5), guide = "prism_offset_minor") +
    coord_cartesian(ylim = y_lim) + facet_grid(. ~ panel, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = if (show_y) sp$ylab else NULL) + base_theme() +
    theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1, lineheight = 0.9), strip.text = element_text(face = "bold", size = 9, lineheight = 0.9),
          panel.spacing.x = unit(0.9, "lines"))
}
banner <- function(line1, line2, w_cm) {
  l1 <- paste(strwrap(line1, width = max(14, floor(w_cm * 3.3))), collapse = "\n")
  l2 <- paste(strwrap(line2, width = max(18, floor(w_cm * 4.4))), collapse = "\n")
  ggplot() + annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = "grey86", colour = NA) +
    annotate("text", x = 0.5, y = 0.94, vjust = 1, label = l1, fontface = "bold", size = 3.7, lineheight = 0.95) +
    annotate("text", x = 0.5, y = 0.06, vjust = 0, label = l2, size = 3.1, colour = "grey20", lineheight = 0.95) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) + theme_void()
}
banner_text <- function(pf) {
  if (pf == "Plate") return(c("CELLS ON CULTURE PLATES (2D)", paste0("seeded at ", density_plate)))
  if (pf == "Hydrogel-seeded") return(c("CELLS SEEDED ON TOP OF HYDROGELS", paste0("seeded at ", density_hydrogel)))
  cl <- sort(unique(na.omit(img$cells[img$platform == "Encapsulated"])))
  c("CELLS ENCAPSULATED IN HYDROGELS (3D)", paste0("encapsulated at ", paste(str_replace(cl, " cells/mL", ""), collapse = " and "), " cells/mL of hydrogel"))
}
platform_column <- function(S, pf, fig_w, fig_h, show_y = TRUE, x_units_total = NULL, cap_h = 0, y_lim = NULL) {
  bt <- banner_text(pf)
  nl <- n_lines(paste(strwrap(bt[1], width = max(14, floor(fig_w * 3.3))), collapse = "\n")) + n_lines(paste(strwrap(bt[2], width = max(18, floor(fig_w * 4.4))), collapse = "\n"))
  wrap_plots(banner(bt[1], bt[2], fig_w), metric_bar_plot(S, pf, fig_w, fig_h, show_y, x_units_total, cap_h, y_lim), ncol = 1, heights = c(0.065 * nl + 0.06, 1))
}
n_bars <- function(S, pf) { s <- S$sm %>% filter(platform == pf); nrow(s) + 1.2 * dplyr::n_distinct(s$panel) }


run_metric <- function(S) {
  sp <- S$spec; caps <- character(0)
  if (!make_bar_figures) return(caps)
  pfs <- levels(img$platform)[levels(img$platform) %in% as.character(S$sm$platform)]
  fig_prefix <- c(Plate = "01_Cells_on_culture_plates", `Hydrogel-seeded` = "02_Cells_seeded_on_hydrogels", Encapsulated = "03_Cells_encapsulated_in_hydrogels")
  if (make_platform_figures) for (pf in pfs) {
    w <- max(12, 0.95 * n_bars(S, pf) + 2.8)
    cap <- caption_for(S, pf, w); capH <- cap_height(cap)
    h <- fig_height + 1.2 + capH - 0.6
    p <- platform_column(S, pf, w, h, cap_h = capH) + plot_annotation(caption = cap_draw(cap), theme = cap_theme())
    save_fig(p, sp$folder, paste0(sp$id, "_", fig_prefix[[pf]]), w, h)
    caps <- c(caps, paste0(sp$id, "_", fig_prefix[[pf]], ": ", gsub("\n", " ", cap)))
  }
  if (length(pfs) > 1) {
    bars <- sapply(pfs, function(pf) n_bars(S, pf)); w <- min(34, 0.95 * sum(bars) + 3 * length(pfs))
    tb <- S$tests
    ns_all <- ns_caption(tb)
    cap <- wrap_parts(c(paste0(sp$ylab, ": plates, seeded and encapsulated cells."), bars_text(sp), es_text, stat_legend(sp), ns_all), w)
    capH <- cap_height(cap); h <- fig_height + 1.5 + capH - 0.6
    cols <- lapply(seq_along(pfs), function(j) platform_column(S, pfs[j], w * bars[j] / sum(bars), h, show_y = j == 1, x_units_total = bars[j], cap_h = capH))
    allp <- wrap_plots(cols, nrow = 1, widths = bars) + plot_annotation(caption = cap_draw(cap), theme = cap_theme())
    save_fig(allp, sp$folder, paste0(sp$id, "_04_All_platforms"), w, h)
    caps <- c(caps, paste0(sp$id, "_04_All_platforms: ", gsub("\n", " ", cap)))
  }
  caps
}

# -----------------------------------------------------------------------------
# 5. RUN ALL METRICS
# -----------------------------------------------------------------------------
all_caps <- character(0); all_tests <- list(); all_sum <- list()
for (sp in specs) {
  S <- prep_spec(sp); if (is.null(S)) next
  message("\n[", sp$id, "] ", sp$ylab)
  all_caps <- c(all_caps, run_metric(S))
  if (!is.null(S$tests)) { all_tests[[sp$id]] <- S$tests; message("   tests: ", nrow(S$tests), " comparisons, ", sum(S$tests$significance != "ns"), " significant") }
  all_sum[[sp$id]] <- S$sm %>% mutate(metric = sp$id, level = sp$level) %>% select(metric, level, platform, panel_plain, condition, n, mean, sd, lo, hi)
}
if (length(all_sum)) readr::write_csv(bind_rows(all_sum) %>% mutate(across(where(is.numeric), ~ round(.x, 4))), data_path("MyHC2_summary_by_condition.csv"))
if (length(all_tests)) readr::write_csv(bind_rows(all_tests) %>% select(-panel) %>% mutate(across(where(is.numeric), ~ round(.x, 5))), data_path("MyHC2_statistical_tests_ES_vs_controls.csv"))
readr::write_csv(img %>% select(-panel, -group_id) %>% mutate(across(where(is.numeric), ~ round(.x, 4))), data_path("MyHC2_image_level_data.csv"))
if (!is.null(obj_myo)) readr::write_csv(obj_myo %>% select(Image, platform, ink, cells, condition, Object, Area_um2, Feret_length_um, MinFeret_width_um, AspectRatio_Feret, Angle_deg, N_nuclei, Used_for_alignment),
                                        data_path("MyHC2_myotube_level_data.csv"))

# -----------------------------------------------------------------------------
# 6. ROSE DIAGRAMS: orientation of the individual myotubes relative to the reference direction
# -----------------------------------------------------------------------------
#  Axial data (a myotube at +20 deg and one at -160 deg have the same orientation) -> shown mirrored. The dashed line is the
#  reference direction (0 deg); a bar pointing along the line = a myotube parallel to the reference. Bars: % of the myotubes
#  of the condition (aspect ratio >= the minimum used in the macro) in each 10-degree bin.
rose_data <- NULL
if (!is.null(obj)) {
  ro <- obj %>% filter(Used_for_alignment == 1, !is.na(Angle_deg)) %>%
    mutate(rel = rel_to_field(Angle_deg),
           bin = floor((rel + 90) / rose_bin) * rose_bin - 90 + rose_bin / 2)
  if (nrow(ro)) {
    tot <- ro %>% count(platform, panel, panel_plain, group_id, ink, condition, name = "n_tot")
    rose_data <- ro %>% count(platform, panel, panel_plain, group_id, ink, condition, bin, name = "cnt") %>%
      left_join(tot, by = c("platform", "panel", "panel_plain", "group_id", "ink", "condition")) %>% mutate(pct = 100 * cnt / n_tot)
  }
}
make_rose <- function(pf) {
  rd <- rose_data %>% filter(platform == pf)
  if (!nrow(rd)) return(NULL)
  rd <- rd %>% mutate(es = str_detect(as.character(condition), "^ES "))
  full <- bind_rows(rd %>% mutate(th = (bin + 360) %% 360), rd %>% mutate(th = (bin + 180 + 360) %% 360))
  full$condition <- factor(as.character(full$condition), levels = cond_levels)
  full$panel <- factor(full$panel, levels = unique(rd$panel))
  nlab <- rd %>% distinct(panel, condition, n_tot) %>% mutate(th = 0, pct = 0)
  nlab$condition <- factor(as.character(nlab$condition), levels = cond_levels); nlab$panel <- factor(nlab$panel, levels = levels(full$panel))
  ymax <- max(full$pct) * 1.05
  p <- ggplot(full, aes(x = th, y = pct)) +
    geom_rect(data = data.frame(xmin = c(0, 345, 165), xmax = c(15, 360, 195)), aes(xmin = xmin, xmax = xmax), ymin = 0, ymax = ymax, fill = "grey80", alpha = 0.6, inherit.aes = FALSE) +
    geom_rect(data = data.frame(xmin = c(75, 255), xmax = c(105, 285)), aes(xmin = xmin, xmax = xmax), ymin = 0, ymax = ymax, fill = "#cfe3f5", alpha = 0.7, inherit.aes = FALSE) +
    geom_col(aes(fill = ink, colour = es, linewidth = es), width = rose_bin) +
    geom_vline(xintercept = c(0, 180), linetype = "dashed", colour = "grey15", linewidth = 0.4) +
    geom_text(data = nlab, aes(label = paste0("n = ", n_tot)), x = 315, y = ymax * 0.93, size = 2.5, colour = "grey20", inherit.aes = FALSE) +
    scale_fill_manual(values = ink_cols, guide = "none") +
    scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = es_yellow), guide = "none") + scale_linewidth_manual(values = c(`FALSE` = 0.25, `TRUE` = 0.8), guide = "none") +
    scale_x_continuous(limits = c(0, 360), breaks = c(0, 90, 180, 270), labels = c("0\u00B0", "90\u00B0", "", "")) +
    scale_y_continuous(limits = c(0, ymax)) +
    coord_polar(theta = "x", start = -pi / 2, direction = -1) +
    facet_grid(panel ~ condition, labeller = labeller(condition = function(x) cond_short[x])) +
    labs(x = NULL, y = NULL) + theme_minimal(base_size = base_size - 2) +
    theme(axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_line(colour = "grey85", linewidth = 0.3), panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold", size = 8.5, lineheight = 0.9), strip.text.y = element_text(angle = 0), axis.text.x = element_text(size = 7), panel.spacing = unit(0.4, "lines"))
  list(plot = p, ncol = length(unique(rd$condition)), nrow = length(unique(rd$panel)), tb = rd)
}
if (plot_rose_diagrams && !is.null(rose_data)) for (pf in levels(img$platform)[levels(img$platform) %in% as.character(rose_data$platform)]) {
  R <- make_rose(pf); if (is.null(R)) next
  w <- R$ncol * 3.3 + 3.2
  cap <- wrap_parts(c(paste0("Orientation of the individual myotubes: ", cap_where(pf)),
                      paste0("Rose diagrams (axial data, shown mirrored): bars give the percentage of the myotubes of each condition (aspect ratio \u2265 ", min_ar, ") whose axis lies in each ", rose_bin,
                             "\u00B0 bin relative to the electric-field direction; the dashed line is the field direction (", ref_txt, "), so bars along it are myotubes parallel to the field (perpendicular to the electrodes) and bars along the perpendicular axis are myotubes parallel to the electrodes; grey sectors: \u00B115\u00B0 around the field direction, blue sectors: \u00B115\u00B0 around the perpendicular direction. The same axis is used for the unstimulated controls. n = number of myotubes."),
                      "Bar colour identifies the ink (culture plates: grey; same colour for controls and ES); bars with a fluorescent-yellow border: electrical stimulation (ES) at the indicated voltage; plain bars: controls (expansion control and differentiation control, no ES)."), w)
  capH <- cap_height(cap); h <- R$nrow * 3.4 + 1.6 + capH
  save_fig(R$plot + plot_annotation(caption = cap_draw(cap), theme = cap_theme()), "4_Alignment_A_myotubes", paste0("rose_", c(Plate = "01_Cells_on_culture_plates", `Hydrogel-seeded` = "02_Cells_seeded_on_hydrogels", Encapsulated = "03_Cells_encapsulated_in_hydrogels")[[pf]]), w, h)
  all_caps <- c(all_caps, paste0("rose_", pf, ": ", gsub("\n", " ", cap)))
}

# -----------------------------------------------------------------------------
# 7. EASY-TO-READ FIGURES: what length / width / angle mean
# -----------------------------------------------------------------------------
ell_poly <- function(cx, cy, a, b, ang = 0, n = 90) {
  t <- seq(0, 2 * pi, length.out = n); x <- a * cos(t); y <- b * sin(t)
  data.frame(x = cx + x * cos(ang) - y * sin(ang), y = cy + x * sin(ang) + y * cos(ang))
}
# 7a. schematic: how to read the morphology and alignment figures ------------------------------------------
make_schematic <- function() {
  myo_col <- "#e43954"; nuc_col <- "#2748c9"
  arr <- arrow(ends = "both", length = unit(0.14, "cm"), type = "closed")
  xs <- c(-2.4, -0.8, 0.8, 2.4)
  nuc <- do.call(rbind, lapply(seq_along(xs), function(i) cbind(ell_poly(xs[i], 0, 0.55, 0.38), id = i)))
  pA <- ggplot() +
    geom_polygon(data = ell_poly(0, 0, 4, 0.95), aes(x, y), fill = myo_col, alpha = 0.55, colour = "black", linewidth = 0.5) +
    geom_polygon(data = nuc, aes(x, y, group = id), fill = nuc_col, colour = NA) +
    annotate("segment", x = -4, xend = 4, y = -1.75, yend = -1.75, arrow = arr, linewidth = 0.5) +
    annotate("text", x = 0, y = -2.35, label = "Length = maximum Feret diameter\n(longest distance across the myotube)", size = 3.2, lineheight = 0.95) +
    annotate("segment", x = 4.7, xend = 4.7, y = -0.95, yend = 0.95, arrow = arr, linewidth = 0.5) +
    annotate("text", x = 5.0, y = 0, label = "Width =\nminimum\nFeret\ndiameter", hjust = 0, size = 3.2, lineheight = 0.95) +
    annotate("text", x = 0, y = 2.25, label = "Aspect ratio = length / width", size = 3.5, fontface = "bold") +
    annotate("text", x = 0, y = 1.65, label = paste0("Myotube = MyHC2-positive object with at least ", min_nuc, " nuclei (blue)"), size = 2.9, colour = "grey25") +
    coord_fixed(xlim = c(-4.6, 7.4), ylim = c(-3.1, 2.7), clip = "off") + theme_void() +
    labs(title = "A   Morphology") + theme(plot.title = element_text(face = "bold", size = 12))
  th <- 28 * pi / 180; a15 <- seq(-15, 15, length.out = 20) * pi / 180
  wedge <- function(base, rad = 4.7) { a <- base + a15; data.frame(x = c(0, rad * cos(a)), y = c(0, rad * sin(a))) }
  arc <- data.frame(x = 1.9 * cos(seq(0, th, length.out = 30)), y = 1.9 * sin(seq(0, th, length.out = 30)))
  tt <- seq(-1, 1, length.out = 40)
  el_l <- data.frame(x = -6.2 - 0.55 * (1 - tt^2), y = 3.3 * tt); el_r <- data.frame(x = 6.2 + 0.55 * (1 - tt^2), y = 3.3 * tt)
  fld <- function(y) annotate("segment", x = -5.3, xend = 5.3, y = y, yend = y, colour = "grey55", linewidth = 0.35, arrow = arrow(length = unit(0.12, "cm"), type = "closed"))
  pB <- ggplot() +
    geom_polygon(data = wedge(0), aes(x, y), fill = "grey82", alpha = 0.8) + geom_polygon(data = wedge(pi), aes(x, y), fill = "grey82", alpha = 0.8) +
    geom_polygon(data = wedge(pi / 2), aes(x, y), fill = "#cfe3f5", alpha = 0.8) + geom_polygon(data = wedge(-pi / 2), aes(x, y), fill = "#cfe3f5", alpha = 0.8) +
    geom_path(data = el_l, aes(x, y), linewidth = 2.4, lineend = "round") + geom_path(data = el_r, aes(x, y), linewidth = 2.4, lineend = "round") +
    fld(-2.4) + fld(2.4) +
    annotate("segment", x = -5.3, xend = 5.3, y = 0, yend = 0, linetype = "dashed", linewidth = 0.5, arrow = arrow(length = unit(0.14, "cm"), type = "closed")) +
    geom_polygon(data = ell_poly(0, 0, 3.3, 0.7, th), aes(x, y), fill = myo_col, alpha = 0.55, colour = "black", linewidth = 0.5) +
    geom_path(data = arc, aes(x, y), linewidth = 0.5) +
    annotate("text", x = 2.35 * cos(th / 2), y = 2.35 * sin(th / 2) - 0.1, label = "\u03B8", size = 4.2, fontface = "italic") +
    annotate("text", x = -6.7, y = 4.0, label = "electrode", size = 3, colour = "grey25") + annotate("text", x = 6.7, y = 4.0, label = "electrode", size = 3, colour = "grey25") +
    annotate("text", x = 0, y = 4.0, label = "Electric field E (reference direction, 0\u00B0)", size = 3.3, fontface = "bold") +
    annotate("text", x = 0, y = -3.9, label = "\u03B8 = angle between the myotube axis and the electric field", size = 3, fontface = "bold") +
    annotate("text", x = 0, y = -4.5, label = "0\u00B0 = parallel to the field (perpendicular to the electrodes)  \u00B7  90\u00B0 = perpendicular to the field (parallel to the electrodes)", size = 2.7) +
    annotate("text", x = 0, y = -5.1, label = "grey sectors: within \u00B115\u00B0 of parallel  \u00B7  blue sectors: within \u00B115\u00B0 of perpendicular", size = 2.7, colour = "grey30") +
    annotate("text", x = 0, y = -5.7, label = "Order parameter S = mean cos 2\u03B8:  +1 = all parallel to the field \u00B7 0 = random \u00B7 \u22121 = all perpendicular", size = 2.7) +
    coord_fixed(xlim = c(-8.2, 8.2), ylim = c(-6.1, 4.4), clip = "off") + theme_void() +
    labs(title = "B   Alignment relative to the electric field") + theme(plot.title = element_text(face = "bold", size = 12))
  (pA | pB) + plot_layout(widths = c(1, 1.4)) + plot_annotation(title = "How to read the morphology and alignment figures",
                                                                theme = theme(plot.title = element_text(face = "bold", size = 13)))
}
dir.create(file.path(out_dir, "0_How_to_read"), showWarnings = FALSE, recursive = TRUE)
save_fig(make_schematic(), "0_How_to_read", "How_to_read_morphology_and_alignment", 30, 11)

# 7b. myotubes drawn to scale: one rectangle per condition = mean length x mean width -------------------------
rect_stripes <- function(x0, x1, y0, y1, d) {
  h <- y1 - y0; out <- list()
  for (cc in seq(x0 - h, x1, by = d)) {
    xs <- max(x0, cc); xe <- min(x1, cc + h)
    if (xs < xe) out[[length(out) + 1]] <- data.frame(xs = xs, xe = xe, ys = y0 + (xs - cc), ye = y0 + (xe - cc))
  }
  if (length(out)) do.call(rbind, out) else NULL
}
glyph_figure <- function(pf) {
  g <- obj_myo %>% filter(platform == pf, !is.na(Feret_length_um), !is.na(MinFeret_width_um)) %>%
    group_by(panel_plain, ink, cells, group_id, condition, cond_group) %>%
    summarise(n = n(), L = mean(Feret_length_um), W = mean(MinFeret_width_um), AR = mean(AspectRatio_Feret), .groups = "drop") %>%
    arrange(ink, cells, condition)
  if (!nrow(g)) return(NULL)
  Lmax <- max(g$L); S <- 44; gap <- 52; y <- 0; rows <- list(); hdr <- list()
  for (gid in unique(g$group_id)) {
    gg <- g[g$group_id == gid, ]
    y <- y - gap; hdr[[length(hdr) + 1]] <- data.frame(y = y, label = gg$panel_plain[1])
    for (r in seq_len(nrow(gg))) { y <- y - S; rows[[length(rows) + 1]] <- cbind(as.data.frame(gg[r, ]), y = y) }
  }
  rows <- do.call(rbind, rows); hdr <- do.call(rbind, hdr)
  rows$es <- rows$cond_group == "ES"; rows$lab <- cond_short[as.character(rows$condition)]
  rows$txt <- sprintf("L %.0f \u00B5m \u00B7 W %.0f \u00B5m \u00B7 AR %.1f \u00B7 n = %d", rows$L, rows$W, rows$AR, rows$n)
  xr <- Lmax * 1.9; y_end <- min(rows$y) - S; y_lo <- y_end - 45
  st <- do.call(rbind, lapply(which(rows$es), function(i) rect_stripes(0, rows$L[i], rows$y[i] - rows$W[i] / 2, rows$y[i] + rows$W[i] / 2, d = min(12, max(4, rows$W[i] * 0.55)))))
  p <- ggplot() +
    geom_rect(data = rows[!rows$es, ], aes(xmin = 0, xmax = L, ymin = y - W / 2, ymax = y + W / 2, fill = ink), colour = "black", linewidth = 0.4) +
    geom_rect(data = rows[rows$es, ], aes(xmin = 0, xmax = L, ymin = y - W / 2, ymax = y + W / 2, fill = ink), colour = es_yellow, linewidth = es_border_w)
  if (!is.null(st)) p <- p + geom_segment(data = st, aes(x = xs, xend = xe, y = ys, yend = ye), colour = es_yellow, linewidth = es_stripe_w, lineend = "butt")
  bt <- banner_text(pf)
  p <- p + geom_text(data = rows, aes(x = -Lmax * 0.03, y = y, label = lab), hjust = 1, size = 3.1) +
    geom_text(data = rows, aes(x = Lmax * 1.04, y = y, label = txt), hjust = 0, size = 2.9, colour = "grey15") +
    geom_text(data = hdr, aes(x = -Lmax * 0.03, y = y, label = label), hjust = 1, fontface = "bold", size = 3.3) +
    annotate("segment", x = 0, xend = 100, y = y_end, yend = y_end, linewidth = 1) +
    annotate("text", x = 50, y = y_end - 20, label = "100 \u00B5m", size = 3) +
    scale_fill_manual(values = ink_cols, guide = "none") +
    coord_fixed(xlim = c(0, xr), ylim = c(y_lo, 5), expand = FALSE, clip = "off") + theme_void() +
    labs(title = bt[1], subtitle = bt[2]) +
    theme(plot.margin = margin(0.4, 0.4, 0.4, 4.6, "cm"), plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 9, colour = "grey25"),
          plot.title.position = "plot")
  panel_w <- 12.5
  list(plot = p, w = panel_w + 5.0, panel_h = (5 - y_lo) * panel_w / xr)
}
if (make_morphology_glyphs && !is.null(obj_myo) && nrow(obj_myo)) for (pf in levels(img$platform)[levels(img$platform) %in% as.character(obj_myo$platform)]) {
  G <- glyph_figure(pf); if (is.null(G)) next
  cap <- wrap_parts(c(paste0("Myotube morphology: ", cap_where(pf)),
                      paste0("Each rectangle is drawn to scale (the same \u00B5m scale horizontally and vertically): its length is the mean maximum Feret diameter and its height the mean minimum Feret diameter (width) of the myotubes of the condition (MyHC2-positive objects with at least ", min_nuc, " nuclei); L = length, W = width, AR = mean aspect ratio (length / width), n = number of myotubes of the two images. See 0_How_to_read for the definitions."),
                      es_text, "The means with their SD, the individual myotubes and the statistics are in the bar figures of this folder."), G$w)
  capH <- cap_height(cap)
  save_fig(G$plot + plot_annotation(caption = cap_draw(cap), theme = cap_theme()), "3_Morphology",
           paste0("morphology_summary_", c(Plate = "01_Cells_on_culture_plates", `Hydrogel-seeded` = "02_Cells_seeded_on_hydrogels", Encapsulated = "03_Cells_encapsulated_in_hydrogels")[[pf]]),
           G$w, G$panel_h + 2.4 + capH)
  all_caps <- c(all_caps, paste0("morphology_summary_", pf, ": ", gsub("\n", " ", cap)))
}

# 7c. orientation histograms (angle relative to the reference, +-15 deg band marked, numbers written in each panel) -----------------
if (make_histograms && !is.null(rose_data)) {
  hstat <- ro %>% group_by(platform, panel, condition) %>% summarise(n = n(), within15 = 100 * mean(abs(rel) <= 15), perp15 = 100 * mean(abs(rel) >= 75), meandev = mean(abs(rel)), .groups = "drop")
  for (pf in levels(img$platform)[levels(img$platform) %in% as.character(rose_data$platform)]) {
    rd <- rose_data %>% filter(platform == pf) %>% mutate(es = str_detect(as.character(condition), "^ES "))
    st <- hstat %>% filter(platform == pf)
    rd$condition <- factor(as.character(rd$condition), levels = cond_levels); st$condition <- factor(as.character(st$condition), levels = cond_levels)
    rd$panel <- factor(rd$panel, levels = unique(rd$panel)); st$panel <- factor(st$panel, levels = levels(rd$panel))
    st$lab <- sprintf("n = %d\n%.0f%% parallel (\u00B115\u00B0)\n%.0f%% perpendicular (\u00B115\u00B0)\nmean angle %.0f\u00B0", st$n, st$within15, st$perp15, st$meandev)
    ymax <- max(rd$pct) * 2.2
    p <- ggplot(rd, aes(x = bin, y = pct)) +
      annotate("rect", xmin = -15, xmax = 15, ymin = 0, ymax = Inf, fill = "grey85", alpha = 0.7) +
      annotate("rect", xmin = c(-92, 75), xmax = c(-75, 92), ymin = 0, ymax = Inf, fill = "#cfe3f5", alpha = 0.7) +
      geom_col(aes(fill = ink, colour = es, linewidth = es), width = rose_bin) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
      geom_label(data = st, aes(label = lab), x = -90, y = ymax, hjust = 0, vjust = 1, size = 2.4, lineheight = 0.95, fill = "white", alpha = 0.85, label.size = 0, label.padding = unit(0.12, "lines"), inherit.aes = FALSE) +
      scale_fill_manual(values = ink_cols, guide = "none") +
      scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = es_yellow), guide = "none") + scale_linewidth_manual(values = c(`FALSE` = 0.25, `TRUE` = 0.8), guide = "none") +
      scale_x_continuous(breaks = c(-90, -45, 0, 45, 90), limits = c(-92, 92)) + scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
      facet_grid(panel ~ condition, labeller = labeller(condition = function(x) cond_short[x])) +
      labs(x = "Myotube angle relative to the electric field (\u00B0)     0\u00B0 = parallel to the field \u00B7 \u00B190\u00B0 = perpendicular to the field", y = "Myotubes (%)") + base_theme() +
      theme(strip.text = element_text(face = "bold", size = 8.5, lineheight = 0.9), strip.text.y = element_text(angle = 0), axis.text = element_text(size = 7.5), panel.spacing = unit(0.5, "lines"))
    ncl <- length(unique(rd$condition)); nrw <- length(unique(rd$panel)); w <- ncl * 3.7 + 3.2
    cap <- wrap_parts(c(paste0("Orientation of the individual myotubes: ", cap_where(pf)),
                        paste0("Bars: percentage of the myotubes of each condition (aspect ratio \u2265 ", min_ar, ") in each ", rose_bin, "\u00B0 class of angle relative to the electric-field direction (",
                               ref_txt, "). The dashed line is the field direction (0\u00B0 = parallel to the field, perpendicular to the electrodes); the grey band marks \u00B115\u00B0 around it and the blue bands \u00B115\u00B0 around the perpendicular direction (parallel to the electrodes). The same axis is used for the unstimulated controls. n = number of myotubes; mean angle = mean absolute angle to the field."),
                        "Bar colour identifies the ink (culture plates: grey; same colour for controls and ES); bars with a fluorescent-yellow border: electrical stimulation (ES); plain bars: controls (expansion control and differentiation control, no ES)."), w)
    h <- nrw * 3.1 + 2.8 + 0.38 * n_lines(cap)
    save_fig(p + plot_annotation(caption = cap_draw(cap), theme = cap_theme()), "4_Alignment_A_myotubes",
             paste0("orientation_histogram_", c(Plate = "01_Cells_on_culture_plates", `Hydrogel-seeded` = "02_Cells_seeded_on_hydrogels", Encapsulated = "03_Cells_encapsulated_in_hydrogels")[[pf]]), w, h)
    all_caps <- c(all_caps, paste0("orientation_histogram_", pf, ": ", gsub("\n", " ", cap)))
  }
}

if (make_summary_figures) {
  # =============================================================================
  # 8. COMPACT FIGURES (replace the many per-metric bar figures; the values are in the table of the thesis)
  #    1_Area_and_fusion/area_fusion_all_platforms     MyHC2 area (solid bar) and fusion index (hatched bar) of every condition,
  #                                                    all platforms in ONE figure, with the platform headers
  #    2_Alignment_summary/alignment_summary_<panel>   ONE figure per platform and ink: the alignment parameters of both methods,
  #                                                    each averaged over the images of the condition
  #    (the rose diagrams of section 6 stay: one figure per platform, all inks and conditions)
  #  Needs (defined above): out_dir, save_fig(), ink_cols, plate_color, es_yellow, banner(), banner_text()
  #  Both quantities are percentages, so they share ONE y axis (area: % of the image; fusion index: % of the nuclei). They are NOT
  #  parts of one another (different denominators), so the bars stand side by side and are not stacked.
  # =============================================================================
  local({
    cl2 <- c("Control (expansion)", "Control (differentiation)", "ES 0.5 V", "ES 1 V", "ES 1.5 V", "ES 2 V")
    cs2 <- c("Exp.", "Diff.", "0.5 V", "1 V", "1.5 V", "2 V")
    I2 <- readr::read_csv(file.path(out_dir, "data", "MyHC2_image_level_data.csv"), show_col_types = FALSE)
    T2 <- readr::read_csv(file.path(out_dir, "data", "MyHC2_statistical_tests_ES_vs_controls.csv"), show_col_types = FALSE)
    I2 <- I2 %>% mutate(panel2 = dplyr::case_when(platform == "Plate" ~ "Culture plate", platform == "Hydrogel-seeded" ~ ink,
                                                  TRUE ~ paste0(cells, " \u00B7 ", ink)),
                        cond2 = factor(condition, levels = cl2))
    fill_for <- function(i) { v <- ink_cols[i]; if (length(v) == 0 || is.na(v)) plate_color else unname(v) }
    pn <- I2 %>% distinct(platform, panel2, ink) %>% mutate(o = match(platform, c("Plate", "Hydrogel-seeded", "Encapsulated"))) %>% arrange(o, panel2)
    ttl <- function(pf, p) ifelse(pf == "Plate", "Culture plate", ifelse(pf == "Hydrogel-seeded", sub("MC", "\nMC", p), sub(" cells/mL \u00B7 ", "\n", p)))
    star_of <- function(tm, panel, es_c) {            # strongest significance of one ES condition against any control
      s <- T2$significance[T2$metric == tm & T2$panel_plain == panel & T2$es == es_c & T2$significance != "ns"]
      if (!length(s)) "" else s[which.max(nchar(s))]
    }
    # ---------------- (1) area + fusion index in one panel per ink ----------------------------------------------
    ymax_all <- max(I2$Fusion_index_pct, na.rm = TRUE) * 1.14
    n_pan <- nrow(pn); n_cond <- sapply(seq_len(n_pan), function(i) length(unique(I2$cond2[I2$panel2 == pn$panel2[i] & I2$platform == pn$platform[i]])))
    fig_w <- 28; fig_h <- 11.5
    unit_w <- (fig_w - 2) / (sum(n_cond) + 0.7 * n_pan); panel_h <- fig_h - 4.2
    k <- ymax_all * unit_w / panel_h                  # slope of the hatching lines (data units) that looks like 45 degrees
    hatch <- function(x0, x1, h, d = 0.11) {
      out <- list(); if (!is.finite(h) || h <= 0) return(NULL)
      for (cc in seq(x0 - h / k, x1, by = d)) { xs <- max(x0, cc); xe <- min(x1, cc + h / k)
      if (xs < xe) out[[length(out) + 1]] <- data.frame(xs = xs, xe = xe, ys = k * (xs - cc), ye = k * (xe - cc)) }
      if (length(out)) do.call(rbind, out) else NULL
    }
    panel_af <- function(i, show_y) {      # fusion index only (the MyHC2 area is no longer shown: it was used to segment the objects)
      d <- I2 %>% filter(panel2 == pn$panel2[i], platform == pn$platform[i])
      present <- cl2[cl2 %in% as.character(d$cond2)]
      m <- do.call(rbind, lapply(seq_along(present), function(j) { x <- d[as.character(d$cond2) == present[j], ]
      data.frame(j = j, cond = present[j], es = grepl("^ES", present[j]),
                 f = mean(x$Fusion_index_pct), f_lo = min(x$Fusion_index_pct), f_hi = max(x$Fusion_index_pct)) }))
      m$sf <- sapply(m$cond, function(cc) if (grepl("^ES", cc)) star_of("fusion", pn$panel2[i], cc) else "")
      bw <- 0.34; col <- fill_for(pn$ink[i])
      pts <- do.call(rbind, lapply(seq_along(present), function(j) { x <- d[as.character(d$cond2) == present[j], ]
      data.frame(x = j, v = x$Fusion_index_pct) }))
      ggplot() +
        geom_rect(data = m, aes(xmin = j - bw, xmax = j + bw, ymin = 0, ymax = f, colour = es, linewidth = es), fill = col) +
        scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = es_yellow), guide = "none") + scale_linewidth_manual(values = c(`FALSE` = 0.4, `TRUE` = 1.1), guide = "none") +
        geom_linerange(data = m, aes(x = j, ymin = f_lo, ymax = f_hi), linewidth = 0.3, alpha = sd_alpha) +
        geom_point(data = pts, aes(x, v), shape = 21, fill = "white", size = 1.1, stroke = 0.3) +
        geom_text(data = m, aes(x = j, y = pmax(f_hi, f) + ymax_all * 0.015, label = sf), vjust = 0, size = 3.2) +
        scale_x_continuous(breaks = seq_along(present), labels = cs2[match(present, cl2)], limits = c(0.4, length(present) + 0.6), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, ymax_all), expand = c(0, 0)) +
        labs(x = NULL, y = if (show_y) "Fusion index (% of the nuclei)" else NULL, title = ttl(pn$platform[i], pn$panel2[i])) + theme_classic(base_size = 8) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), plot.title = element_text(size = 7, face = "bold"), axis.line = element_line(linewidth = 0.4),
              axis.text.y = if (show_y) element_text() else element_blank(), axis.ticks.y = if (show_y) element_line() else element_blank())
    }
    key <- ggplot() + theme_void() + coord_cartesian(xlim = c(0, 14), ylim = c(0, 1), expand = FALSE) +
      annotate("text", x = 1.2, y = 0.5, label = "Fusion index: share of the nuclei inside myotubes (%)", hjust = 0, size = 2.8) +
      annotate("text", x = 9.6, y = 0.5, label = "yellow border: ES;  * differs from a control", hjust = 0, size = 2.4, colour = "grey30")
    groups <- lapply(c("Plate", "Hydrogel-seeded", "Encapsulated"), function(pf) {
      idx <- which(pn$platform == pf); if (!length(idx)) return(NULL)
      wd <- n_cond[idx] + 0.7; wtot <- sum(wd) * unit_w
      row <- wrap_plots(lapply(seq_along(idx), function(q) panel_af(idx[q], show_y = (pf == "Plate" && q == 1))), nrow = 1, widths = wd)
      bt <- banner_text(pf)
      list(plot = wrap_plots(banner(bt[1], bt[2], wtot), row, ncol = 1, heights = c(0.36, 1)), w = sum(wd))
    })
    groups <- Filter(Negate(is.null), groups)
    fig1 <- wrap_plots(wrap_plots(lapply(groups, `[[`, "plot"), nrow = 1, widths = sapply(groups, `[[`, "w")), key, ncol = 1, heights = c(1, 0.06))
    save_fig(fig1, "1_Fusion_index", "fusion_index_all_platforms", fig_w, fig_h)
    # ---------------- (2) alignment summary: one figure per platform and ink ------------------------------------
    bar_panel <- function(d, col, tm, ylab, title, ref = NULL) {
      d <- d %>% mutate(v = .data[[col]]) %>% filter(!is.na(v))
      if (!nrow(d)) return(ggplot() + theme_void())
      present <- cl2[cl2 %in% as.character(d$cond2)]; d$cond2 <- factor(as.character(d$cond2), levels = present)
      s <- d %>% group_by(cond2) %>% summarise(m = mean(v), sd = sd(v), .groups = "drop") %>% mutate(es = grepl("^ES", cond2), sdv = ifelse(is.na(sd), 0, sd))
      p <- ggplot(s, aes(cond2, m)) +
        geom_col(fill = fill_for(d$ink[1]), aes(colour = es, linewidth = es), width = 0.72) +
        scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = es_yellow), guide = "none") + scale_linewidth_manual(values = c(`FALSE` = 0.4, `TRUE` = 1.1), guide = "none") +
        geom_errorbar(aes(ymin = m - sdv, ymax = m + sdv), width = 0.22, linewidth = 0.35, alpha = sd_alpha) +
        geom_point(data = d, aes(cond2, v), shape = 21, fill = "white", size = 1.3, stroke = 0.3, inherit.aes = FALSE)
      if (!is.null(ref)) p <- p + geom_hline(yintercept = ref, linetype = "dashed", colour = "grey35", linewidth = 0.35)
      tt <- T2 %>% filter(metric == tm, panel_plain == d$panel2[1], significance != "ns", es %in% present, control %in% present)
      top <- max(c(s$m + s$sdv, d$v), na.rm = TRUE); rng <- diff(range(c(0, top)))
      if (nrow(tt)) {
        tt <- tt %>% mutate(x1 = match(control, present), x2 = match(es, present), span = abs(x1 - x2)) %>% arrange(span) %>% mutate(y = top + rng * (0.06 + 0.11 * (row_number() - 1)))
        p <- p + geom_segment(data = tt, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.3) +
          geom_segment(data = tt, aes(x = x1, xend = x1, y = y, yend = y - rng * 0.025), inherit.aes = FALSE, linewidth = 0.3) +
          geom_segment(data = tt, aes(x = x2, xend = x2, y = y, yend = y - rng * 0.025), inherit.aes = FALSE, linewidth = 0.3) +
          geom_text(data = tt, aes(x = (x1 + x2) / 2, y = y, label = significance), inherit.aes = FALSE, vjust = -0.15, size = 3)
      }
      p + scale_x_discrete(labels = cs2[match(present, cl2)]) + scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
        labs(x = NULL, y = ylab, title = title) + theme_classic(base_size = 8) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), plot.title = element_text(size = 7, face = "bold"), axis.line = element_line(linewidth = 0.4))
    }
    for (i in seq_len(n_pan)) {
      d <- I2 %>% filter(panel2 == pn$panel2[i], platform == pn$platform[i])
      if (!any(grepl("^ES", d$condition))) next
      f <- bar_panel(d, "Pct_parallel15", "within15", "Parallel \u00B115\u00B0 (%)", "A  Myotubes: parallel", ref = 100 * 30 / 180) +
        bar_panel(d, "Pct_perp15", "perp15", "Perpendicular \u00B115\u00B0 (%)", "A  Myotubes: perpendicular", ref = 100 * 30 / 180) +
        bar_panel(d, "Mean_angle_to_field", "meandev", "Mean angle to field (\u00B0)", "A  Myotubes: mean angle", ref = 45) +
        bar_panel(d, "ST_S_field", "S", "Order parameter S", "B  Structure tensor: S", ref = 0) +
        bar_panel(d, "ST_alignment_score_10deg_pct", "score10", "Score \u00B110\u00B0 (%)", "B  Structure tensor: score", ref = 100 * 20 / 180) + plot_layout(nrow = 1)
      bt <- banner_text(pn$platform[i])
      save_fig(wrap_plots(banner(bt[1], paste0(bt[2], "  \u00B7  ", gsub("\n", " ", ttl(pn$platform[i], pn$panel2[i]))), 27), f, ncol = 1, heights = c(0.2, 1)),
               "2_Alignment_summary", paste0("alignment_summary_", gsub("[^A-Za-z0-9]+", "_", paste(pn$platform[i], pn$panel2[i]))), 27, 7)
    }
  })
}
writeLines(all_caps, data_path("MyHC2_figure_captions.txt"))
message("\nDone. Figures in: ", normalizePath(out_dir), "\n      tables and captions in: ", normalizePath(file.path(out_dir, "data")))



# ---- one-way ANOVA across the conditions of each panel (three or more conditions) -------------------
local({
  fi <- file.path(out_dir, "data", "MyHC2_image_level_data.csv"); fm <- file.path(out_dir, "data", "MyHC2_myotube_level_data.csv")
  res <- list()
  if (file.exists(fi)) { d <- readr::read_csv(fi, show_col_types = FALSE)
  for (v in intersect(c("MyHC2_area_fraction_pct","Fusion_index_pct","Pct_within_15deg","Mean_abs_angle_dev_deg","ST_order_parameter_S","ST_alignment_score_10deg_pct"), names(d)))
    for (g in unique(paste(d$platform, d$ink, d$cells, sep = " | "))) {
      x <- d[paste(d$platform, d$ink, d$cells, sep = " | ") == g, ]
      if (length(unique(x$condition)) >= 3 && sum(!is.na(x[[v]])) > length(unique(x$condition)))
        res[[length(res) + 1]] <- data.frame(level = "image", panel = g, metric = v, p_anova = summary(aov(x[[v]] ~ factor(x$condition)))[[1]][["Pr(>F)"]][1]) } }
  if (file.exists(fm)) { d <- readr::read_csv(fm, show_col_types = FALSE)
  for (v in intersect(c("Feret_length_um","MinFeret_width_um","AspectRatio_Feret"), names(d)))
    for (g in unique(paste(d$platform, d$ink, d$cells, sep = " | "))) {
      x <- d[paste(d$platform, d$ink, d$cells, sep = " | ") == g, ]
      if (length(unique(x$condition)) >= 3 && sum(!is.na(x[[v]])) > length(unique(x$condition)))
        res[[length(res) + 1]] <- data.frame(level = "myotube", panel = g, metric = v, p_anova = summary(aov(x[[v]] ~ factor(x$condition)))[[1]][["Pr(>F)"]][1]) } }
  if (length(res)) readr::write_csv(do.call(rbind, res), file.path(out_dir, "data", "MyHC2_anova_by_panel.csv"))
})

