# =============================================================================
#  Myogenin-positive fraction (%) - C2C12 immunostaining (DAPI / phalloidin / myogenin)
#   Plots: cells on culture plates | cells seeded on hydrogels | cells encapsulated in hydrogels | all together
#   Statistics: Student's t-test, every electrical-stimulation (ES) condition vs the expansion control and vs the
#               differentiation control (same panel); only significant comparisons get a bracket, the others are
#               listed in the caption (same convention as the characterization script).
#  Plot style: GraphPad Prism-like ('ggprism'); ink colours = those of the characterization script.
#  Input : CSV exported from your sheet (blocks SeedPlate / SeedHydrogel / Encapsulation, each with its own header row).
#          Used: Total_DAPI_nuclei and Myog_positive_nuclei (background-ROI method: threshold = mean + 2 SD).
#  Output: <out_dir>/  figures (PDF + PNG)  and  data/  (csv tables, captions)
# =============================================================================

# -----------------------------------------------------------------------------
# 0. OPTIONS
# -----------------------------------------------------------------------------
csv_path <- '/Users/tomas_serra/Library/CloudStorage/GoogleDrive-tomasfserra1@gmail.com/My Drive/Tese Drive/Current Thesis Files/C2C12/Immunostainings/Quantifications/Myogenin/Cell quantifications - Immuno DAPI + Phalloidin + Myogenin .csv'            # path to the CSV, or NA to choose the file in a window
out_dir  <- '/Users/tomas_serra/Library/CloudStorage/GoogleDrive-tomasfserra1@gmail.com/My Drive/Tese Drive/Current Thesis Files/C2C12/Immunostainings/Quantifications/Myogenin/Myog_Quant_Plots'            # NA = folder "Myogenin results" next to the CSV

col_total <- "Total_DAPI_nuclei"
col_pos   <- "Myog_positive_nuclei"
col_frac  <- "Myog_positive_fraction_pct"      # only used to cross-check positives / total

density_plate    <- "2,000 cells/cm\u00B2"
density_hydrogel <- "200,000 cells/cm\u00B2"   # <-- CONFIRM this number
es_unit <- "V"
myog_medium <- "PM3"      # medium of the hydrogels in the immunostainings (the colour is looked up as ink + medium)

# --- ink colours: copied from the characterization script (ink + medium -> colour) ---------------------
ink_table <- data.frame(
  medium = c("PM3", "PM3", "PM3", "PM3", "PM3", "PM3", "PM3",
             "L15-PM", "L15-PM", "L15-PM", "L15-PM", "DMEM-PM", "DMEM-PM", "MACK-PM"),
  ink    = c("30k-c13.3MC0Gel", "30k-c13.3MC100Gel", "72k-c40MC0Gel", "90k-c40MC0Gel", "90k-c40MC100Gel", "90k-c40MC200Gel", "90k-c40MC375Gel",
             "90k-c", "72k-c40MC0Gel", "72k-c40MC100Gel", "72k-c40MC375Gel", "72k-c40MC0Gel", "72k-c40MC100Gel", "72k-c40MC0Gel"),
  color  = c("#ffcad8", "#ffae9c", "#fa7b8e", "#9f0019", "#e43954", "#ff742f", "#ffc14a",
             "#003b74", "#43b6ff", "#009692", "#b7e8d8", "#a600a5", "#eb67ea", "#afc40b"),
  stringsAsFactors = FALSE)
plate_color <- "#cdcdcd"                       # culture plate: no ink -> neutral grey
ink_rx <- "[0-9]+k-c[0-9._]+MC[0-9]+Gel"        # how an ink name looks inside an image name (e.g. 90k-c40MC375Gel)

# --- electrical stimulation (ES) bars: fluorescent-yellow border + stripes ----------------------------
es_yellow   <- "#fffd3c"
es_border_w <- 1          # border thickness
es_stripe_w <- 0.6          # stripe thickness
es_stripe_d <- 0.4         # distance between stripes (x units; 1 = distance between two bars)

# --- statistics (same convention as the characterization script) ---------------------------------------
t_var_equal <- TRUE         # TRUE = Student's t-test; FALSE = Welch's
#  Significance: ns p >= 0.05; * 0.01-0.05; ** 0.001-0.01; *** 0.0001-0.001; **** < 0.0001
#  Non-significant comparisons get NO bracket on the plot; they are listed in the caption.
run_tests <- TRUE

fig_height <- 9             # cm (single platform plot; the caption adds height)
base_size  <- 11            # font size
bracket_gap  <- 3           # distance of the first bracket above the tallest point (% units)
bracket_step <- 5.5         # vertical spacing between stacked brackets
bracket_tick <- 1.5         # length of the bracket end ticks
save_pdf <- FALSE

# -----------------------------------------------------------------------------
# 1. PACKAGES AND READING
# -----------------------------------------------------------------------------
need <- function(p) if (!requireNamespace(p, quietly = TRUE)) stop(sprintf('Package "%s" is missing: install.packages("%s")', p, p), call. = FALSE)
invisible(lapply(c("readr", "dplyr", "tidyr", "stringr", "ggplot2", "ggprism", "patchwork", "tibble"), need))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(stringr); library(ggplot2); library(ggprism); library(patchwork) })

if (is.na(csv_path) || !file.exists(csv_path)) csv_path <- file.choose()
if (is.na(out_dir)) out_dir <- file.path(dirname(normalizePath(csv_path)), "Myogenin results")
dir.create(file.path(out_dir, "data"), recursive = TRUE, showWarnings = FALSE)
message("CSV:    ", csv_path, "\nOutput: ", out_dir)

first_line <- readLines(csv_path, n = 1, warn = FALSE, encoding = "UTF-8")
delim <- if (str_count(first_line, ";") > str_count(first_line, ",")) ";" else ","
raw0 <- readr::read_delim(csv_path, delim = delim, col_types = readr::cols(.default = readr::col_character()),
                          trim_ws = TRUE, show_col_types = FALSE, name_repair = "minimal")
to_num <- function(x) {
  x <- str_replace_all(as.character(x), "\\s", "")
  x <- ifelse(str_detect(x, ",") & !str_detect(x, "\\."), str_replace(x, ",", "."), x)
  suppressWarnings(as.numeric(x))
}
for (cc in c(col_total, col_pos)) if (!cc %in% names(raw0)) stop('Column "', cc, '" not found. Columns are: ', paste(names(raw0), collapse = ", "))

# first column = image names; its header cell is the title of the FIRST block, each further block has its own title + header row
block_first <- names(raw0)[1]
names(raw0)[1] <- "Image"
raw0$row_id <- seq_len(nrow(raw0))
is_hdr <- dplyr::coalesce(raw0[[col_total]] == col_total, FALSE)
block <- character(nrow(raw0)); cur <- block_first
for (i in seq_len(nrow(raw0))) { if (is_hdr[i]) cur <- raw0$Image[i]; block[i] <- cur }
raw0$block <- block
raw <- raw0[!is_hdr & !is.na(raw0$Image) & nzchar(str_trim(raw0$Image)), ]
message("Blocks found: ", paste(unique(raw$block), collapse = " | "), "\nImages read: ", nrow(raw))

img <- tibble(row_id = raw$row_id, Image = raw$Image, block = raw$block, total = to_num(raw[[col_total]]), pos = to_num(raw[[col_pos]]),
              frac_stated = if (col_frac %in% names(raw)) to_num(raw[[col_frac]]) else NA_real_) %>%
  mutate(frac = 100 * pos / total)
if (anyNA(img$frac)) stop("Missing nuclei counts for rows: ", paste(img$row_id[is.na(img$frac)], collapse = ", "))
bad <- which(!is.na(img$frac_stated) & abs(img$frac_stated - img$frac) > 0.05)
if (length(bad)) message("!! Stated fraction differs from positives/total for rows: ", paste(img$row_id[bad], collapse = ", ")) else
  message("Per-image check: positives / total x 100 equals the fraction in your sheet for every image.")

# -----------------------------------------------------------------------------
# 2. PLATFORM, INK, CELL DENSITY AND CONDITION  (read from the block title and the image name)
# -----------------------------------------------------------------------------
rx_encaps <- "(^|[_ .-])(10|25)M([_ .-]|$)"          # "..._10M_..." / "..._25M_..." = cells per mL of hydrogel
ink_names <- ink_table$ink[ink_table$medium == myog_medium]
ink_cols  <- setNames(ink_table$color[ink_table$medium == myog_medium], ink_names)
ink_cols  <- c(ink_cols, "Culture plate (no hydrogel)" = plate_color)

get_ink <- function(x) {
  tok <- str_extract(x, regex(ink_rx, ignore_case = TRUE))
  tok <- str_replace_all(tok, "_", ".")
  hit <- ink_names[match(tolower(tok), tolower(ink_names))]
  ifelse(is.na(hit), NA_character_, hit)
}
parse_meta <- function(x, blk) {
  volt <- str_match(x, "(?<![0-9.])([0-9]+(?:[._][0-9]+)?)\\s*V(?![A-Za-z])")[, 2]
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
meta <- parse_meta(img$Image, img$block)
meta$ink <- ifelse(meta$platform == "Plate", "Culture plate (no hydrogel)", get_ink(img$Image))
chk_tab <- bind_cols(img %>% select(row_id, Image), meta)
problems <- chk_tab %>% filter(is.na(condition) | (platform != "Plate" & is.na(ink)) | (platform == "Encapsulated" & is.na(cells)))
if (nrow(problems)) {
  print(as.data.frame(problems %>% select(row_id, Image, platform, ink_found, ink, cells, condition)))
  stop("The rows above could not be fully understood.\n",
       "  - ink NA: the ink in the image name (column ink_found) is not in 'ink_table' for medium '", myog_medium, "' (section 0).\n",
       "  - condition / cells NA: the image name has no ControlExpan / ControlDiff / '<n>V' / 10M / 25M.")
}
d <- bind_cols(img %>% select(-Image), chk_tab %>% select(Image, platform, ink, cells, condition)) %>%
  mutate(platform = factor(platform, levels = c("Plate", "Hydrogel-seeded", "Encapsulated")))
volt_levels <- unique(d$condition[str_detect(d$condition, "^ES ")])
volt_levels <- volt_levels[order(as.numeric(str_match(volt_levels, "^ES ([0-9.]+)")[, 2]))]
cond_levels <- c("Control (expansion)", "Control (differentiation)", volt_levels)
cond_short <- setNames(c("Expansion control", "Differentiation control", volt_levels), cond_levels)
d <- d %>% mutate(condition = factor(condition, levels = cond_levels),
                  cond_group = case_when(condition == "Control (expansion)" ~ "expansion", condition == "Control (differentiation)" ~ "diff", TRUE ~ "ES"),
                  panel = case_when(platform == "Plate" ~ "Culture plate", platform == "Hydrogel-seeded" ~ sub("MC", "\nMC", ink),
                                    TRUE ~ paste0(cells, "\n", sub("MC", "\nMC", ink))),
                  panel_plain = case_when(platform == "Plate" ~ "Culture plate", platform == "Hydrogel-seeded" ~ ink, TRUE ~ paste0(cells, " \u00B7 ", ink)),
                  group_id = paste(platform, ink, coalesce(cells, "-"), sep = "|"))
message("\nWhat the script understood (check this table - the ink column must be filled for every hydrogel image):")
print(as.data.frame(d %>% select(row_id, platform, ink, cells, condition, Image) %>% arrange(platform, ink, cells, condition)), row.names = FALSE)
d <- d %>% group_by(group_id, condition) %>% mutate(rep_n = row_number(), n_img = n()) %>% ungroup()
nt <- d %>% distinct(group_id, condition, n_img) %>% filter(n_img != 2)
if (nrow(nt)) { message("\n!! Conditions with n != 2 images:"); print(as.data.frame(nt)) } else message("\nEvery condition has n = 2 images.")

# -----------------------------------------------------------------------------
# 3. SUMMARY AND CHECK OF YOUR MEANS
# -----------------------------------------------------------------------------
data_path <- function(name) file.path(out_dir, "data", name)
summ <- d %>% group_by(platform, ink, cells, panel, group_id, condition, cond_group) %>%
  summarise(panel_plain = first(panel_plain), n_images = n(), mean_frac = mean(frac), min_frac = min(frac), max_frac = max(frac), sd_frac = sd(frac), .groups = "drop") %>%
  arrange(platform, ink, cells, condition) %>%
  group_by(group_id) %>% mutate(pos = match(as.character(condition), cond_levels[cond_levels %in% as.character(condition)])) %>% ungroup()
readr::write_csv(summ %>% mutate(across(where(is.numeric), ~ round(.x, 3))), data_path("myogenin_summary_by_condition.csv"))

cand <- names(raw)[str_detect(names(raw), regex("(mean|avg|average|m[e\u00E9]dia).*(frac|pos|%)|(frac|pos|%).*(mean|avg|average|m[e\u00E9]dia)", ignore_case = TRUE)) &
                     !str_detect(names(raw), regex("raw|BGsub|intens|MFI|median|bleed|otsu", ignore_case = TRUE)) & names(raw) != col_frac]
for (cn in cand) {
  chk <- tibble(row_id = raw$row_id, in_sheet = to_num(raw[[cn]])) %>% filter(!is.na(in_sheet)) %>%
    left_join(d %>% select(row_id, Image, platform, ink, cells, condition, group_id), by = "row_id") %>%
    left_join(summ %>% select(group_id, condition, recalculated_mean = mean_frac), by = c("group_id", "condition")) %>%
    mutate(difference = in_sheet - recalculated_mean, ok = abs(difference) < 0.05)
  readr::write_csv(chk %>% select(row_id, Image, platform, ink, cells, condition, in_sheet, recalculated_mean, difference, ok) %>%
                     mutate(across(where(is.numeric), ~ round(.x, 3))), data_path("myogenin_mean_check.csv"))
  message("\nMean check, column '", cn, "': ", sum(chk$ok), " of ", nrow(chk), " averages in your sheet match the mean of the two images of that condition",
          if (any(!chk$ok)) paste0("  (!! ", sum(!chk$ok), " differ - see myogenin_mean_check.csv)") else "  (all correct).")
}

# -----------------------------------------------------------------------------
# 4. STATISTICS: each ES condition vs the expansion control and vs the differentiation control
# -----------------------------------------------------------------------------
sig_symbol <- function(p) {
  if (is.na(p)) "NA" else if (p < 0.0001) "****" else if (p < 0.001) "***" else
    if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}
tests <- NULL
if (run_tests) {
  rows <- list()
  for (g in unique(d$group_id)) {
    dg <- d[d$group_id == g, ]
    for (e in unique(as.character(dg$condition[dg$cond_group == "ES"]))) {
      for (cn in c("Control (differentiation)", "Control (expansion)")) {
        a <- dg$frac[as.character(dg$condition) == e]; b <- dg$frac[as.character(dg$condition) == cn]
        if (length(a) < 2 || length(b) < 2) next
        tt <- tryCatch(t.test(a, b, var.equal = t_var_equal), error = function(err) NULL)
        if (is.null(tt) || is.na(tt$p.value)) { message("Test not possible (identical values): ", e, " vs ", cn, " [", dg$panel_plain[1], "]"); next }
        rows[[length(rows) + 1]] <- data.frame(platform = as.character(dg$platform[1]), panel = dg$panel[1], panel_plain = dg$panel_plain[1], group_id = g, es = e, control = cn,
                                               mean_es = mean(a), mean_control = mean(b), difference_pp = mean(a) - mean(b),
                                               t = unname(tt$statistic), df = unname(tt$parameter), p_value = tt$p.value,
                                               significance = sig_symbol(tt$p.value), stringsAsFactors = FALSE)
      }
    }
  }
  if (length(rows)) {
    tests <- do.call(rbind, rows)
    readr::write_csv(tests %>% select(-panel) %>% rename(panel = panel_plain) %>% mutate(across(where(is.numeric), ~ round(.x, 4))), data_path("myogenin_ttests_ES_vs_controls.csv"))
    message("\nStudent's t-tests (ES vs controls): ", nrow(tests), " comparisons, ", sum(tests$significance != "ns"), " significant.")
  } else message("No ES condition with a control in the same panel: no tests.")
}
ctrl_short <- c("Control (differentiation)" = "differentiation control", "Control (expansion)" = "expansion control")
ns_caption <- function(tb, with_panel = FALSE) {
  if (is.null(tb) || !nrow(tb)) return(NULL)
  if (!any(tb$significance == "ns")) return(NULL)
  "Other comparisons were not significant (p \u2265 0.05)."
}
stat_legend <- function() paste0(
  "Each ES condition vs each control: Student's t-test; all conditions of a panel: one-way ANOVA (see the tables). ",
  "* p < 0.05; ** p < 0.01; *** p < 0.001; **** p < 0.0001.")

# brackets (only significant comparisons), stacked per panel
brk <- NULL
if (!is.null(tests)) {
  sg <- tests[tests$significance != "ns", , drop = FALSE]
  if (nrow(sg)) {
    out <- list()
    for (g in unique(sg$group_id)) {
      s <- sg[sg$group_id == g, ]
      present <- cond_levels[cond_levels %in% as.character(d$condition[d$group_id == g])]
      top <- max(d$frac[d$group_id == g])
      s$x1 <- match(s$es, present); s$x2 <- match(s$control, present); s$span <- abs(s$x1 - s$x2)
      s <- s[order(s$span), ]
      s$y <- top + bracket_gap + (seq_len(nrow(s)) - 1) * bracket_step
      out[[length(out) + 1]] <- s
    }
    brk <- do.call(rbind, out)
  }
}
y_top <- max(106, if (!is.null(brk)) max(brk$y) + 9 else 0)

# -----------------------------------------------------------------------------
# 5. PLOTTING
# -----------------------------------------------------------------------------
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
y_lab <- "Myogenin-positive fraction (%)"
# FIGURE CAPTIONS: set draw_caption <- TRUE to print captions inside the figures (default: not drawn; the LaTeX caption carries the text,
# and the caption texts are still written to the *_figure_captions.txt file).
draw_caption <- FALSE
n_lines <- function(txt) length(strsplit(txt, "\n", fixed = TRUE)[[1]])
cap_height <- function(cap) if (draw_caption) 0.38 * n_lines(cap) + 0.6 else 0
cap_draw <- function(cap) if (draw_caption) cap else NULL

# diagonal stripes inside the ES bars (drawn with plain ggplot2 segments; k = slope in data units for ~45 degrees on the page)
stripe_segments <- function(sg, k) {
  out <- list()
  for (r in seq_len(nrow(sg))) {
    h <- sg$mean_frac[r]; x0 <- sg$pos[r] - 0.36; x1 <- sg$pos[r] + 0.36
    if (!is.finite(h) || h <= 0) next
    for (cc in seq(x0 - h / k, x1, by = es_stripe_d)) {
      xs <- max(x0, cc); xe <- min(x1, cc + h / k)
      if (xs < xe) out[[length(out) + 1]] <- data.frame(panel = sg$panel[r], xs = xs, xe = xe, ys = k * (xs - cc), ye = k * (xe - cc))
    }
  }
  if (length(out)) do.call(rbind, out) else NULL
}

bar_plot <- function(pf, fig_w, fig_h, show_y = TRUE, x_units_total = NULL, cap_h = 0, key = TRUE) {
  s <- summ %>% filter(platform == pf); i <- d %>% filter(platform == pf)
  s$x <- factor(cond_short[as.character(s$condition)], levels = unname(cond_short))
  i$x <- factor(cond_short[as.character(i$condition)], levels = unname(cond_short))
  s$panel <- factor(s$panel, levels = unique(s$panel)); i$panel <- factor(i$panel, levels = levels(s$panel))
  n_pan <- nlevels(s$panel)
  if (is.null(x_units_total)) x_units_total <- nrow(s) + 1.2 * n_pan
  panel_w <- max(3, fig_w - 2.2 - 0.5 * n_pan); panel_h <- max(3, fig_h - cap_h - 5.9)
  k <- (y_top / panel_h) * (panel_w / x_units_total)               # data slope that looks like 45 degrees
  s_es <- s %>% filter(cond_group == "ES"); s_ct <- s %>% filter(cond_group != "ES")
  p <- ggplot(s, aes(x = x)) +
    geom_col(data = s_ct, aes(y = mean_frac, fill = ink), colour = "black", linewidth = 0.4, width = 0.72) +
    geom_col(data = s_es, aes(y = mean_frac, fill = ink), colour = es_yellow, linewidth = es_border_w, width = 0.72)
  st <- stripe_segments(as.data.frame(s_es), k)
  if (!is.null(st)) {
    st$panel <- factor(st$panel, levels = levels(s$panel))
    p <- p + geom_segment(data = st, aes(x = xs, xend = xe, y = ys, yend = ye), colour = es_yellow, linewidth = es_stripe_w,
                          inherit.aes = FALSE, lineend = "butt")
  }
  p <- p +
    geom_linerange(aes(ymin = min_frac, ymax = max_frac), linewidth = 0.5, colour = "black", alpha = sd_alpha) +
    geom_point(data = i, aes(y = frac), shape = 21, fill = "white", colour = "black", size = 2.4, stroke = 0.6,
               position = position_jitter(width = 0.07, height = 0, seed = 7))
  if (!is.null(brk)) {
    b <- brk %>% filter(platform == pf); b$panel <- factor(b$panel, levels = levels(s$panel))
    if (nrow(b)) p <- p +
        geom_segment(data = b, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
        geom_segment(data = b, aes(x = x1, xend = x1, y = y, yend = y - bracket_tick), inherit.aes = FALSE, linewidth = 0.4) +
        geom_segment(data = b, aes(x = x2, xend = x2, y = y, yend = y - bracket_tick), inherit.aes = FALSE, linewidth = 0.4) +
        geom_text(data = b, aes(x = (x1 + x2) / 2, y = y, label = significance), inherit.aes = FALSE, vjust = -0.2, size = 3.5)
  }
  p +
    scale_fill_manual(values = ink_cols, name = NULL,
                      guide = if (key) guide_legend(keywidth = unit(0.45, "cm"), keyheight = unit(0.45, "cm")) else "none") +
    scale_y_continuous(breaks = seq(0, 100, 20), expand = c(0, 0), guide = "prism_offset_minor") +
    coord_cartesian(ylim = c(0, y_top)) +
    facet_grid(. ~ panel, scales = "free_x", space = "free_x", labeller = labeller(panel = strip_text)) +
    labs(x = NULL, y = if (show_y) y_lab else NULL) +
    base_theme() +
    theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1), strip.text = element_text(face = "bold", size = 9, lineheight = 0.9),
          panel.spacing.x = unit(0.9, "lines"))
}
# ink names are NOT written above the panels any more: the ink is identified by the colour key (beside the plot or beside the whole figure);
# the strips only keep 'Culture plate' and the cell density of the encapsulated cultures
strip_text <- function(x) ifelse(grepl("^[0-9]+M cells/mL", x), sub("\n.*$", "", x), ifelse(grepl("^[0-9]+k-c", x), "", x))
sd_alpha <- 0.4                      # translucency of the range lines (the 'SD' bars of the other scripts)
key_w    <- 4.4                      # cm added for the colour key beside a figure
# colour key (swatch + ink name) of the inks of the given platforms
ink_key_plot <- function(pfs) {
  kk <- unique(d[d$platform %in% pfs, c("platform", "ink")]); kk$lab <- ifelse(kk$platform == "Plate", "Culture plate", kk$ink)
  kk <- kk[!duplicated(kk$lab), ]; kk$y <- rev(seq_len(nrow(kk))); kk$col <- ifelse(kk$platform == "Plate", plate_color, unname(ink_cols[kk$ink]))
  g <- ggplot(kk) + geom_tile(aes(x = 0.5, y = y, fill = col), width = 0.7, height = 0.78, colour = "black", linewidth = 0.4) +
    geom_text(aes(x = 1.0, y = y, label = lab), hjust = 0, size = 3.1) + scale_fill_identity() +
    coord_cartesian(xlim = c(0, 7), ylim = c(0.4, nrow(kk) + 0.6), expand = FALSE) + theme_void()
  wrap_plots(plot_spacer(), g, plot_spacer(), ncol = 1, heights = c(1.2, max(1.2, 0.45 * nrow(kk)), 4))
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
  cl <- sort(unique(na.omit(d$cells[d$platform == "Encapsulated"])))
  c("CELLS ENCAPSULATED IN HYDROGELS (3D)", paste0("encapsulated at ", paste(str_replace(cl, " cells/mL", ""), collapse = " and "), " cells/mL of hydrogel"))
}
platform_column <- function(pf, fig_w, fig_h, show_y = TRUE, x_units_total = NULL, cap_h = 0, key = TRUE) {
  bt <- banner_text(pf)
  nl <- n_lines(paste(strwrap(bt[1], width = max(14, floor(fig_w * 3.3))), collapse = "\n")) + n_lines(paste(strwrap(bt[2], width = max(18, floor(fig_w * 4.4))), collapse = "\n"))
  wrap_plots(banner(bt[1], bt[2], fig_w), bar_plot(pf, fig_w, fig_h, show_y, x_units_total, cap_h, key), ncol = 1, heights = c(0.065 * nl + 0.06, 1))
}
n_bars <- function(pf) { s <- summ %>% filter(platform == pf); nrow(s) + 1.2 * dplyr::n_distinct(s$panel) }

# captions -------------------------------------------------------------------------------------------
rep_desc <- c(Plate = "two fields of view (Center and Corner)", `Hydrogel-seeded` = "two fields of view (Center and Corner)",
              Encapsulated = "two z-planes of the same construct")
cap_what <- c(
  Plate = paste0("Myogenin-positive fraction (%), culture plates (", density_plate, ")."),
  `Hydrogel-seeded` = paste0("Myogenin-positive fraction (%), cells seeded on hydrogels (", density_hydrogel, ")."),
  Encapsulated = "Myogenin-positive fraction (%), cells encapsulated (10 or 25 million cells/mL).")
cap_method <- "Mean of two images (circles: images; line: range). Yellow border and stripes: ES; plain bars: controls."
cap_bars <- function(pf) NULL
caption_for <- function(pf, w) {
  tb <- if (!is.null(tests)) tests[tests$platform == pf, ] else NULL
  wrap_parts(c(cap_what[[pf]], cap_method, if (!is.null(tests)) stat_legend(), ns_caption(tb, with_panel = TRUE)), w)
}
save_fig <- function(p, name, w, h) {
  f <- file.path(out_dir, name)
  ggsave(paste0(f, ".png"), p, width = w, height = h, units = "cm", dpi = 300, bg = "white")
  if (save_pdf) ggsave(paste0(f, ".pdf"), p, width = w, height = h, units = "cm", bg = "white")
}

fig_name <- c(Plate = "01_Cells_on_culture_plates", `Hydrogel-seeded` = "02_Cells_seeded_on_hydrogels", Encapsulated = "03_Cells_encapsulated_in_hydrogels")
all_captions <- character(0)
present_pf <- levels(d$platform)[levels(d$platform) %in% as.character(d$platform)]
for (pf in present_pf) {
  w <- max(12, 0.95 * n_bars(pf) + 2.8)
  cap <- caption_for(pf, w)
  h <- fig_height + 1.2 + 0.38 * n_lines(cap)
  p <- platform_column(pf, w, h, cap_h = cap_height(cap)) + plot_annotation(caption = cap_draw(cap), theme = cap_theme())
  save_fig(p, fig_name[[pf]], w + key_w, h)
  all_captions <- c(all_captions, paste0(fig_name[[pf]], ": ", gsub("\n", " ", cap)))
}
if (length(present_pf) > 1) {
  bars <- sapply(present_pf, n_bars); w <- min(34, 0.95 * sum(bars) + 3 * length(present_pf))
  ns_all <- if (!is.null(tests)) ns_caption(tests, with_panel = TRUE)
  cap <- wrap_parts(c("Myogenin-positive fraction (%) of C2C12 cells cultured on culture plates, seeded on top of hydrogels and encapsulated in hydrogels. ",
                      cap_method,
                      "Bars show the mean of the two images of each condition (n = 2), open circles the individual images and the line their range. Bar colour identifies the ink (same colour for controls and ES); bars with a fluorescent-yellow border and stripes: electrical stimulation (ES) at the indicated voltage; plain bars: controls (expansion control and differentiation control, no ES).",
                      "The platforms differ in seeding density and dimensionality: each condition should be compared with the controls of its own panel, not across panels.",
                      if (!is.null(tests)) stat_legend(), ns_all), w)
  h <- fig_height + 1.5 + 0.38 * n_lines(cap)
  cols <- lapply(seq_along(present_pf), function(j) platform_column(present_pf[j], w * bars[j] / sum(bars), h, show_y = j == 1, x_units_total = bars[j], cap_h = cap_height(cap), key = FALSE))
  allp <- wrap_plots(wrap_plots(cols, nrow = 1, widths = bars), ink_key_plot(present_pf), nrow = 1, widths = c(w, key_w)) + plot_annotation(caption = cap_draw(cap), theme = cap_theme())
  save_fig(allp, "04_All_platforms", w + key_w, h)
  all_captions <- c(all_captions, paste0("04_All_platforms: ", gsub("\n", " ", cap)))
}
writeLines(all_captions, data_path("myogenin_figure_captions.txt"))
message("\nDone. Figures in: ", normalizePath(out_dir), "\n      tables and captions in: ", normalizePath(file.path(out_dir, "data")))

