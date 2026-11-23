# =============================================================================
#  Hydrogel ink characterization
#   A) Swelling degree over time          -> mean ± SD (n = 3)
#   B) Water content                      -> mean ± SD, Student's t-test
#   C) Compression tests                  -> stress-strain curves, linear fits,
#                                            Young's modulus: mean ± SD, Student's t-test
#   D) Printability factor                -> Pr = L^2 / (16 A), dot plots with the adequate range
#   F) Rheology                           -> viscosity, thixotropy (recovery), gelation (mean ± SD)
#   E) Comparison groups for Results & Discussion (section 7)
#  Plot style: GraphPad Prism-like (package 'ggprism'; install once with install.packages("ggprism"))
#  Output: <results_root>/<technique>/  (plots)  and  <technique>/data/  (csv tables)
# =============================================================================

library(ggplot2)
if (!requireNamespace("ggprism", quietly = TRUE))
  message("Package 'ggprism' not installed -> using a Prism-like fallback theme. ",
          "For the real Prism style run once: install.packages(\"ggprism\")")

# -----------------------------------------------------------------------------
# 1. DATA ENTRY — SWELLING DEGREE (edit weights/colours only)
# -----------------------------------------------------------------------------
#  'medium' groups the inks into separate plots (PM3, L15-PM, DMEM-PM).
#  'color' is any R colour name or hex code; it is used for that ink in EVERY
#  plot (swelling, water content, compression, comparison groups).
#  Replicates are NOT matched across days -> swelling uses the day means
#  (see sd_method in the options to understand the size of the error bars).
#  Measurements not yet taken -> NA (that day is simply skipped).

inks <- list(
  # ---------------- PM3 ----------------
  list(medium = "PM3", ink = "30k-c13.3MC0Gel",   color = "#ffcad8",
       day0 = c(494.5333333, 476.4, 348.6666667),
       day1 = c(424, 277.7, 356.6),
       day3 = c(364.4, 274, 185.3),
       day7 = c(294.5, 144.9, 125.8),
       day10 = c(161.6, 135, 68.2)),
  
  list(medium = "PM3", ink = "30k-c13.3MC100Gel", color = "#ffae9c",
       day0 = c(440.5, 291.2333333, 334.8333333),
       day1 = c(212.2, 280.6, 180.4),
       day3 = c(148.2, 198.4, 110.9),
       day7 = c(204.6, 147.9, 114.7),
       day10 = c(196.8, 138.6, 113)),
  
  list(medium = "PM3", ink = "72k-c40MC0Gel",     color = "#fa7b8e",
       day0 = c(549.3, 589.3, 503.4),
       day1 = c(496.3, 550.1, 522.9),
       day3 = c(497.6, 541.5, 475.3),
       day7 = c(515.2, 510.8, 476.8),
       day10 = c(407.3, 490.5, 432.4)),
  
  list(medium = "PM3", ink = "90k-c40MC0Gel",     color = "#9f0019",
       day0 = c(635.5, 582.2666667, 473.2),
       day1 = c(540.8333333, 632.9666667, 458.6),
       day3 = c(444.5, 542.2, 596.2),
       day7 = c(537.4, 455.5, 511.6),
       day10 = c(509, 439.4, 544.7)),
  
  list(medium = "PM3", ink = "90k-c40MC100Gel",   color = "#e43954",
       day0 = c(583.7, 618, 572.8666667),
       day1 = c(542.3666667, 589.7, 527.3333333),
       day3 = c(512.2, 593.1, 520.5),
       day7 = c(508.8, 560.9, 457.8),
       day10 = c(510.3, 460, 559)),
  
  list(medium = "PM3", ink = "90k-c40MC200Gel",   color = "#ff742f",
       day0 = c(440.7666667, 556.6, 503.2333333),
       day1 = c(398.1, 456.2, 530.7),
       day3 = c(444.9, 509.2, 381.3),
       day7 = c(511.1, 362.8, 441.5),
       day10 = c(442, 387.9, 500.3)),
  
  list(medium = "PM3", ink = "90k-c40MC375Gel",   color = "#ffc14a",
       day0 = c(543.5, 571.1, 520.7),
       day1 = c(516.2, 488.7, 491.1),
       day3 = c(483.6, 487, 483.8),
       day7 = c(481.1, 498.6, 460.6),
       day10 = c(451.1, 447.8, 428.3)),
  
  # ---------------- L15-PM ----------------
  list(medium = "L15-PM", ink = "90k-c",           color = "#003b74",
       day0 = c(526, 657.6, 536.5666667),
       day1 = c(514.7, 642.1, 511.3),
       day3 = c(609.4, 471.2, 481.3),
       day7 = c(512.5, 495.1, 610.7),
       day10 = c(482.3, 621.5, 452.6)),
  
  list(medium = "L15-PM", ink = "72k-c40MC0Gel",   color = "#43b6ff",
       day0 = c(522.6666667, 527.0333333, 497.4666667),
       day1 = c(516.4, 528.1, 492),
       day3 = c(498.5, 494.3, 466.7),
       day7 = c(501.5, 496.5, 469.3),
       day10 = c(498.5, 469, 469.2)),
  
  list(medium = "L15-PM", ink = "72k-c40MC100Gel", color = "#009692",
       day0 = c(492.4666667, 519, 646.1333333),
       day1 = c(489.4, 622.4, 525.4),
       day3 = c(464.8, 560.9, 416.5),
       day7 = c(589, 473.2, 419.6),
       day10 = c(569.5, 423.5, 448.4)),
  
  list(medium = "L15-PM", ink = "72k-c40MC375Gel", color = "#b7e8d8",
       day0 = c(521.9666667, 520.3666667, 558.2),
       day1 = c(496.9, 471.5, 544.1),
       day3 = c(430.2, 424.8, 446),
       day7 = c(426.8, 238.7, 434.6),
       day10 = c(412.6, 404.8, 410.4)),
  
  # ---------------- DMEM-PM ----------------
  list(medium = "DMEM-PM", ink = "72k-c40MC0Gel",   color = "#a600a5",
       day0 = c(638.5, 613.6, 455.6),
       day1 = c(454, 630.5, 600.3),
       day3 = c(491.1, 454.1, 450.1),
       day7 = c(490.6, 425.6, 445.7),
       day10 = c(434.2, 474.2, 415.3)),
  
  list(medium = "DMEM-PM", ink = "72k-c40MC100Gel", color = "#eb67ea",
       day0 = c(515.8, 537.6, 508.2),
       day1 = c(528.5, 505.5, 498.8),
       day3 = c(431.2, 592.4, 550),
       day7 = c(584.8, 528.5, 396.5),
       day10 = c(412.5, 548.3, 568.8)),
  
  # ---------------- MACK-PM (only printability so far; swelling not measured) ----------------
  list(medium = "MACK-PM", ink = "72k-c40MC0Gel",   color = "#afc40b",
       day0 = c(NA, NA, NA), day1 = c(NA, NA, NA), day3 = c(NA, NA, NA),
       day7 = c(NA, NA, NA), day10 = c(NA, NA, NA))
)

# -----------------------------------------------------------------------------
# 1b. DATA ENTRY — WATER CONTENT (edit weights only; colours come from above)
# -----------------------------------------------------------------------------
#  wo = initial (wet) mass, wd = dry mass, same units, n = 3 each.
#  medium + ink must match an entry of 'inks' above (that is how the colour is found).
#  Not measured yet -> NA (ink is skipped).

water <- list(
  # ---------------- PM3 ----------------
  list(medium = "PM3", ink = "30k-c13.3MC0Gel",   wo = c(252, 263.2, 198), wd = c(3.4, 2.8, 2)),
  list(medium = "PM3", ink = "30k-c13.3MC100Gel", wo = c(248.1, 166.6, 239.1), wd = c(7.7, 4.7, 7)),
  list(medium = "PM3", ink = "72k-c40MC0Gel",     wo = c(533, 432, 529.3), wd = c(12.9, 10.2, 11.7)),
  list(medium = "PM3", ink = "90k-c40MC0Gel",     wo = c(572.8, 560.9, 422.8), wd = c(18.4, 12.7, 15.1)),
  list(medium = "PM3", ink = "90k-c40MC100Gel",   wo = c(534.5, 534.2, 510.3), wd = c(13.6, 16.8, 12.9)),
  list(medium = "PM3", ink = "90k-c40MC200Gel",   wo = c(509.2, 592.7, 594.2), wd = c(13.4, 15.9, 9.8)),
  list(medium = "PM3", ink = "90k-c40MC375Gel",   wo = c(547.3, 459.4, 432.8), wd = c(27.2, 20.9, 11.2)),
  
  # ---------------- L15-PM ----------------
  list(medium = "L15-PM", ink = "90k-c",           wo = c(462.9, 550.4, 484.2), wd = c(14.3, 15.1, 18.5)),
  list(medium = "L15-PM", ink = "72k-c40MC0Gel",   wo = c(469.2, 461.6, 469.6), wd = c(13, 13.7, 13.6)),
  list(medium = "L15-PM", ink = "72k-c40MC100Gel", wo = c(479.5, 498.1, 476.5), wd = c(17.5, 17.5, 16)),
  list(medium = "L15-PM", ink = "72k-c40MC375Gel", wo = c(543.5, 815.6, 532), wd = c(43.1, 27.6, 23.1)),
  
  # ---------------- DMEM-PM ----------------
  list(medium = "DMEM-PM", ink = "72k-c40MC0Gel",   wo = c(506.5, 514.2, 507.9), wd = c(16.8, 16.8, 17.1)),
  list(medium = "DMEM-PM", ink = "72k-c40MC100Gel", wo = c(449.6, 440.4, 558.2), wd = c(10.2, 12.6, 13.3))
)

# -----------------------------------------------------------------------------
# 1c. DATA ENTRY — COMPRESSION TESTS: specimen dimensions (mm)
# -----------------------------------------------------------------------------
#  Force/displacement are read automatically from 'comp_dir' (section 2).
#  Values are in the order of the replicate files _01, _02, _03.
#   height   = measured height (mm); only used for the strain if comp_L0_source = "measured"
#              (default is now "size_contact", the machine height at the start of loading)
#   diameter = diameter (mm)           -> area A = pi * d^2 / 4 -> stress = F / A
#  Missing diameter -> that replicate is skipped. Missing height -> the first
#  'Size_mm' value of the csv is used instead (with a message).

comp_dims <- list(
  # ---------------- PM3 ----------------
  list(medium = "PM3", ink = "30k-c13.3MC0Gel",   height = c(7.6, 5.7, 5.2),  diameter = c(9, 9, 9)),
  list(medium = "PM3", ink = "30k-c13.3MC100Gel", height = c(5.9, 6.38, 6.2), diameter = c(9.5, 8.8, 9)),
  list(medium = "PM3", ink = "72k-c40MC0Gel",     height = c(7.3, 7, 5.5),    diameter = c(8.5, 9, 9)),
  list(medium = "PM3", ink = "90k-c40MC0Gel",     height = c(7.4, 7.5, 7),    diameter = c(9, 9, 9)),
  list(medium = "PM3", ink = "90k-c40MC100Gel",   height = c(7.2, 7.6, 7.5),  diameter = c(9, 9, 9)),
  list(medium = "PM3", ink = "90k-c40MC200Gel",   height = c(7.2, 6.5, 7.1),  diameter = c(9, 9, 8.5)),
  list(medium = "PM3", ink = "90k-c40MC375Gel",   height = c(8, 6.7, 7.3),    diameter = c(9, 9, 9)),
  # ---------------- L15-PM ----------------
  list(medium = "L15-PM", ink = "90k-c",           height = c(9.3, 7, 10.1), diameter = c(9.3, 9.1, 9.2)),
  list(medium = "L15-PM", ink = "72k-c40MC0Gel",   height = c(8, 9.4, 8.1),  diameter = c(9.5, 10, 9.6)),
  list(medium = "L15-PM", ink = "72k-c40MC100Gel", height = c(8, 9, 9.5),    diameter = c(9.5, 9.2, 9)),
  list(medium = "L15-PM", ink = "72k-c40MC375Gel", height = c(7.1, 6, 6.5),  diameter = c(9.5, 9.5, 9)),
  # ---------------- DMEM-PM ----------------
  list(medium = "DMEM-PM", ink = "72k-c40MC0Gel",   height = c(7.8, 7.3, 7.8), diameter = c(9, 9.2, 8.9)),
  # !! CHECK: 4 values here (all other inks have 3). Only files _01.._04 that exist are used.
  list(medium = "DMEM-PM", ink = "72k-c40MC100Gel", height = c(7.8, 6.5, 6.5, 7.3), diameter = c(9.5, 9, 9, 9))
)

# -----------------------------------------------------------------------------
# 2. OPTIONS
# -----------------------------------------------------------------------------
days        <- c(1, 3, 7, 10)          # swelling time points
comp_dir <- '/Users/tomas_serra/Library/CloudStorage/GoogleDrive-tomasfserra1@gmail.com/My Drive/Tese Drive/Current Thesis Files/Ink Characterizations/Compression Tests'
#  ALL results go here, in one sub-folder per technique ("Swelling degree", "Water content",
#  "Compression tests", "Printability", "Viscosity", "Thixotropy", "Gelation",
#  "Comparisons/<group>"); the csv tables are in the 'data' sub-folder of each.
results_root <- '/Users/tomas_serra/Library/CloudStorage/GoogleDrive-tomasfserra1@gmail.com/My Drive/Tese Drive/Current Thesis Files/Ink Characterizations/Characterization results'
out_dir  <- results_root
fig_width   <- 12                  # cm (single plot)
fig_height  <- 9                   # cm
#  INK KEY (colour + ink name). Ink names are NOT written on the x axis any more:
#   - plots with several inks (per medium, all media): the key (colour + name) is drawn beside the plot;
#   - comparison plots (Comparisons/<group>/): no key and no names at all; the key is saved as a separate small image,
#     <group>/cmp_<id>_key.png (and Comparisons/Keys/<set>.png for the sets below) to be placed next to the plots in the thesis.
key_w       <- 5.6                 # cm added to the width of bar plots for the key at the side
cmp_skip    <- c("stress_strain")   # comparison plots NOT made (add "thixotropy" to skip the thixotropy comparison curves)
#  ONE key image with ALL inks (colour + name): Comparisons/Keys/cmp_key_all_inks.png -> use it with every comparison plot, in Chapter 4 and in the appendix
base_size   <- 11                  # font size

# --- FIGURE PANELS FOR THE THESIS (section 8) ---------------------------------
#  The readability problem came from saving the all-media figures 36 cm wide and
#  including them at ~16 cm: LaTeX shrank them to 44 %, so 11 pt text printed at ~5 pt.
#  The panels below are saved at the width they are USED at, so nothing is shrunk.
make_panels      <- TRUE
panel_dir        <- "Panels"       # sub-folder of out_dir where the panels are written
panel_width_cm   <- 17             # = \linewidth of the thesis; do NOT scale in LaTeX
panel_height_cm  <- 7.4            # height of the PLOT row; the key row is added on top of this
panel_base_size  <- 9              # font size inside the panels
panel_dpi        <- 600
panel_widths     <- c(2, 1)        # (a) : (b) in the first row; (a) takes two thirds of the width
panel_point_size <- 0.9            # dots of the curve/dot plots inside the panels
key_ncol         <- 7              # ink key: 7 entries per row, so the 14 inks fit in 2 rows
key_text_size    <- 2.3
E_band_3D        <- c(1, 3)        # kPa, optimum for 3D culture   (Costantini et al. 2017)
E_band_2D        <- c(12, 15)      # kPa, optimum for 2D substrates (Brock et al. 2022; Serena et al. 2010)
E_band_3D_fill   <- "#FFD43B"      # yellow, 3D window
E_band_2D_fill   <- "#B0B0B0"      # grey, 2D window
visc_ref_rate    <- 100            # shear rate (1/s) of the viscosity summary panel; Fonseca et al. 2020 recommend comparing extrusion bioinks at 100 1/s
gel_ref_temps    <- c(30, 15)      # temperatures (C) of the G' panel: printing temp and end of cooling
ph_csv           <- '/Users/tomas_serra/Library/CloudStorage/GoogleDrive-tomasfserra1@gmail.com/My Drive/Tese Drive/Current Thesis Files/Ink Characterizations/Experiments - pH.csv'

# --- comparison plots ---------------------------------------------------------
#  The comparison figures were removed from the thesis: the panels above show every
#  ink in one plot, so the per-group comparisons repeated information. Set to TRUE
#  to produce them again.
run_comparison_plots <- FALSE

# --- degradation ---------------------------------------------------------------
#  The assay is reported as mass loss (degradation), not swelling:
#  Mass loss (%) = (W0 - Wt)/W0 x 100 = -(swelling). The computation is unchanged;
#  only the sign and the axis label of the plot change.
report_mass_loss <- TRUE

# --- alternative ink names --------------------------------------------------------
#  Some files/folders use a different name for an ink than the 'inks' list above.
#  Format: c("name used in the file" = "name used in 'inks'"). Applied to the
#  compression, printability and rheology (incl. gelation) file/folder names.
ink_aliases <- c("90k-c0MC0Gel" = "90k-c")
#  name shown in the keys and legends (the thesis name of the ink); the internal label is kept for the files
ink_display <- c(
  "90k-c (L15-PM)"             = "FD[hL]+G0",
  "90k-c0MC0Gel (L15-PM)"      = "FD[hL]+G0",
  "72k-c40MC0Gel (L15-PM)"     = "FD[mM]+G0",
  "72k-c40MC100Gel (L15-PM)"   = "FD[mM]+G1",
  "72k-c40MC375Gel (L15-PM)"   = "FD[mM]+G3.75",
  "72k-c40MC0Gel (MACK-PM)"    = "FM[mM]+G0",
  "72k-c40MC0Gel (DMEM-PM)"    = "C[mM]+G0",
  "72k-c40MC100Gel (DMEM-PM)"  = "C[mM]+G1",
  "72k-c40MC0Gel (PM3)"        = "C(-)[mM]+G0",
  "90k-c40MC0Gel (PM3)"        = "C(-)[hM]+G0",
  "90k-c40MC100Gel (PM3)"      = "C(-)[hM]+G1",
  "90k-c40MC200Gel (PM3)"      = "C(-)[hM]+G2",
  "90k-c40MC375Gel (PM3)"      = "C(-)[hM]+G3.75",
  "30k-c13.3MC0Gel (PM3)"      = "C(-)[lL]+G0",
  "30k-c13.3MC100Gel (PM3)"    = "C(-)[lL]+G1")

# --- statistics (all t-tests) ---------------------------------------------------
t_var_equal <- TRUE                # TRUE = Student's t-test; FALSE = Welch's
#  Significance: ns p >= 0.05; * 0.01-0.05; ** 0.001-0.01; *** 0.0001-0.001; **** < 0.0001
#  Non-significant comparisons get NO bracket on the plot; they are listed in the caption.

# --- swelling: how the mean ± SD of the swelling degree is obtained ----------------
#  Swelling (%) = (W_t - W_0) / W_0 * 100.  W_0 and W_t are not matched per
#  hydrogel, so the swelling is calculated from the day means and the SD has to
#  carry the spread of BOTH days (that is why the bars are long: the hydrogels
#  differ a lot in size/weight, e.g. CV of W_0 ~ 5-20 %):
#   "propagated": SD_S = 100 * (mean_t/mean_0) * sqrt((sd_t/mean_t)^2 + (sd_0/mean_0)^2)
#   "day_t_only": SD_S = 100 * sd_t / mean_0            (ignores the spread of day 0)
#   "paired"    : swelling of each hydrogel i = (W_t,i - W_0,i)/W_0,i * 100, then
#                 mean ± SD of the n values. ONLY valid if value i of day 0 and of
#                 day t are the same hydrogel (usually gives much smaller SD).
#  All three SDs are written to swelling_sd_breakdown.csv for comparison.
sd_method <- "propagated"

# --- water content ------------------------------------------------------------
#  "means" : WC = (mean Wo - mean Wd) / mean Wo * 100, SD by error propagation.
#            Use this when Wo and Wd are NOT matched per hydrogel.
#  "paired": WC per hydrogel (Wo_i, Wd_i), then mean ± SD of the 3 values.
wc_method   <- "means"
show_points <- TRUE                # individual hydrogels as dots (only possible in "paired")

#  t-tests for the water-content plots (computed from mean, SD and n of each ink).
#  Each pair as c("ink (medium)", "ink (medium)"), both in the same medium.
wc_comparisons <- list(
  c("90k-c40MC0Gel (PM3)",     "30k-c13.3MC0Gel (PM3)"),
  c("30k-c13.3MC0Gel (PM3)",   "30k-c13.3MC100Gel (PM3)"),
  c("90k-c0MC0Gel (L15-PM)",   "72k-c40MC0Gel (L15-PM)"),
  c("90k-c40MC0Gel (PM3)",     "90k-c40MC100Gel (PM3)"),
  c("90k-c40MC0Gel (PM3)",     "90k-c40MC200Gel (PM3)"),
  c("90k-c40MC0Gel (PM3)",     "90k-c40MC375Gel (PM3)"),
  c("72k-c40MC0Gel (L15-PM)",  "72k-c40MC100Gel (L15-PM)"),
  c("72k-c40MC0Gel (L15-PM)",  "72k-c40MC375Gel (L15-PM)"),
  c("72k-c40MC0Gel (DMEM-PM)", "72k-c40MC100Gel (DMEM-PM)")
)
# set to list() for plots without significance brackets

# --- water content plot: y axis ------------------------------------------------
wc_y_min        <- 90              # y axis starts at 90 %
wc_bracket_gap  <- 0.6             # distance of first bracket above the tallest bar (% units)
wc_bracket_step <- 0.9             # vertical spacing between stacked brackets
wc_bracket_tick <- 0.25            # length of the bracket end ticks
wc_label_room   <- 0.9             # room above the top bracket for the asterisks
wc_value_room   <- 1.6             # room for the printed values between the bars and the first bracket (% units)

# --- compression tests ---------------------------------------------------------
comp_dir       <- '/Users/tomas_serra/Library/CloudStorage/GoogleDrive-tomasfserra1@gmail.com/My Drive/Tese Drive/Current Thesis Files/Ink Characterizations/Compression Tests'
lin_max_strain <- 0.15             # linear fit region: 0 to 15 % strain
r2_min         <- 0.9              # fits with R² below this are flagged
curve_step     <- 0.005            # strain step (mm/mm) of the averaged curve (~ acquisition spacing)
curve_point_size <- 1.1            # size of the dots in the stress-strain curves
curve_show_sd    <- TRUE           # thin SD bars on each dot of the mean curve

#  t-tests for the Young's modulus plots (computed on the replicate moduli, n = 3).
#  The gelatin pairs of 'wc_comparisons' PLUS the pairs that the text argues from:
#  dilution, gelatin in the diluted ink and methylcellulose. A bracket can only be
#  drawn between two inks of the SAME printing medium, because the plot is faceted
#  by medium; the comparisons across media are in the appendix table.
E_comparisons <- wc_comparisons   # the same pre-specified pairs; see the list above

#  Start-of-loading detection (initial noise / stabilisation removal).
#  The noise level is estimated from the first 'onset_n_base' points. Loading
#  starts at the first point where the force stays above
#  baseline + onset_k * noise for 'onset_n_consec' consecutive points.
#  All points before it are removed; the displacement at that point is
#  subtracted (new zero) and the force baseline (median of removed points) is
#  subtracted from the force.
onset_n_base   <- 10
onset_k        <- 2
onset_n_consec <- 5
#  Initial height L0 used for the strain (strain = corrected displacement / L0):
#   "size_contact"           : DEFAULT. Size_mm at the detected start of loading (gap between the
#                              plates when the plate touches the gel, after noise removal)
#   "measured"               : the height you typed in 'comp_dims' 
#   "measured_minus_removed" : measured height minus the displacement removed before the start of loading
#   "size_first"             : Size_mm of the first row of the csv (gap between the plates at the start)
#  The 3 candidate heights are written side by side to compression_onsets.csv.
comp_L0_source <- "size_contact"

#  Manual override of the start point (use after checking the
#  'compression_onset_check_*' plots). onset_mm = raw displacement (mm) where
#  loading really starts. Example:
#   list(medium = "PM3", ink = "30k-c13.3MC0Gel", rep = 1, onset_mm = 1.05),
onset_manual <- list(
)

# --- printability --------------------------------------------------------------
#  pr_dir/<ink>_<medium>/Results_1.txt, Results_2.txt, ...  (one txt per printed mesh;
#  columns: 1 = index, 2 = pore area, 3 = pore perimeter).
pr_dir        <- '/Users/tomas_serra/Library/CloudStorage/GoogleDrive-tomasfserra1@gmail.com/My Drive/Tese Drive/Current Thesis Files/Ink Characterizations/Printability'
pr_n_pores    <- 20                # per mesh: keep the N pores with Pr closest to 1 (15 or 20; same for all)
pr_n_meshes   <- 3                 # meshes expected per ink; fewer -> flagged
pr_adequate   <- c(0.9, 1.1)       # grey region: adequate printability (Ouyang et al. 2016)
pr_seed       <- 1                 # random pairing of pores between replicate meshes (fixed -> reproducible)

# --- comparison groups (section 7) ------------------------------------------------
#  labels : inks compared in the plot, written "ink (medium)"; the first one is
#           the REFERENCE for the t-tests (override with ref = "...").
#  folder : name of the sub-folder (inside "Comparisons") where this group's plots are saved.
#  tests  : "compression" (Young's modulus + stress-strain curves),
#           "water" (water content), "swelling" (swelling degree),
#           "printability" (Pr dot plot). Rheology (viscosity, thixotropy, gelation) is added
#           to every group automatically (see 'rheo_in_groups').
#  Inks that have no data yet (e.g. MACK-PM, 72k-c40MC50Gel) are skipped with a message.
comparison_groups <- list(
  list(id = "1_MC_printability", folder = "MC addition effect",
       title = "MC improves printability (L15-PM)",
       labels = c("90k-c (L15-PM)", "72k-c40MC0Gel (L15-PM)"),
       tests = c("compression", "printability")),
  
  list(id = "2_solvent_composition", folder = "Solvent composition effect",
       title = "Effect of solvent composition (72k-c40MC0Gel)",
       labels = c("72k-c40MC0Gel (L15-PM)", "72k-c40MC0Gel (MACK-PM)",
                  "72k-c40MC0Gel (DMEM-PM)", "72k-c40MC0Gel (PM3)"),
       tests = c("compression", "water", "swelling", "printability")),
  
  list(id = "3_DMEM_to_Media3", folder = "DMEM vs Media 3",
       title = "DMEM vs Media 3 as solvent",
       labels = c("72k-c40MC0Gel (DMEM-PM)", "72k-c40MC100Gel (DMEM-PM)",
                  "72k-c40MC0Gel (PM3)", "72k-c40MC100Gel (PM3)", "90k-c40MC0Gel (PM3)"),
       tests = c("compression", "water", "swelling", "printability")),
  
  list(id = "4a_gelatin_DLEC", folder = "Gelatin effect - DLEC inks",
       title = "Effect of gelatin — DLEC inks (L15-PM)",
       labels = c("72k-c40MC0Gel (L15-PM)", "72k-c40MC50Gel (L15-PM)",
                  "72k-c40MC100Gel (L15-PM)", "72k-c40MC375Gel (L15-PM)"),
       tests = c("compression", "water", "swelling", "printability")),
  
  list(id = "4b_gelatin_C2C12", folder = "Gelatin effect - C2C12 inks",
       title = "Effect of gelatin — C2C12 inks (PM3)",
       labels = c("90k-c40MC0Gel (PM3)", "90k-c40MC100Gel (PM3)", "90k-c40MC200Gel (PM3)",
                  "90k-c40MC375Gel (PM3)", "72k-c40MC0Gel (PM3)", "72k-c40MC100Gel (PM3)"),
       tests = c("compression", "water", "swelling", "printability")),
  
  list(id = "5_material_decrease", folder = "Decrease of k-c and MC",
       title = "Decrease of k-c and MC content (PM3)",
       labels = c("90k-c40MC0Gel (PM3)", "30k-c13.3MC0Gel (PM3)"),
       tests = c("compression", "water", "swelling", "printability"))
)

# =============================================================================
# 3. SETUP  (no need to edit below this line)
# =============================================================================
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ink_label <- function(x) paste0(x$ink, " (", x$medium, ")")
all_labels  <- sapply(inks, ink_label)
all_media   <- unique(sapply(inks, `[[`, "medium"))
ink_colors  <- setNames(sapply(inks, `[[`, "color"), all_labels)
safe_name   <- function(m) gsub("[^A-Za-z0-9_-]", "_", m)
disp        <- function(l) { i <- match(l, names(ink_display)); ifelse(is.na(i), l, unname(ink_display[i])) }
strip_lab   <- function(strip) function(l) { l <- disp(l); if (strip) sub(" \\(.*\\)$", "", l) else l }

# --- where files go: <out_dir>/<technique folder>/ (plots) and .../data/ (csv) -------------
cmp_folders <- setNames(vapply(comparison_groups, function(g) if (!is.null(g$folder)) g$folder else g$id, ""),
                        vapply(comparison_groups, function(g) g$id, ""))
result_folder <- function(name) {
  n <- sub("\\.(png|csv)$", "", name)
  if (grepl("^cmp_", n)) {
    ids <- names(cmp_folders)[startsWith(n, paste0("cmp_", names(cmp_folders), "_"))]
    return(file.path("Comparisons", if (length(ids)) cmp_folders[[ids[1]]] else "Other"))
  }
  rules <- c("^comparison_groups" = "Comparisons", "^swelling" = "Swelling degree",
             "^water_content" = "Water content", "^compression|^youngs_modulus" = "Compression tests",
             "^printability" = "Printability", "^pH" = "pH", "^rheology_viscosity" = "Viscosity",
             "^rheology_(thixotropy|recovery)" = "Thixotropy", "^rheology_gelation" = "Gelation")
  for (r in names(rules)) if (grepl(r, n)) return(rules[[r]])
  "Other"
}
out_path <- function(name) {                     # csv tables
  d <- file.path(out_dir, result_folder(name), "data")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.path(d, name)
}
save_plot <- function(p, name, w, h) {
  d <- file.path(out_dir, result_folder(name))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(d, paste0(name, ".png")), p, width = w, height = h, units = "cm", dpi = 300)
  if (!is.null(.cap_env$last)) cat(name, ": ", gsub("\n", " ", .cap_env$last), "\n", file = file.path(out_dir, "figure_captions.txt"), append = TRUE, sep = "")
}

# Prism-like theme: ggprism::theme_prism (thick axes, outward ticks, no grid, no top/right box)
base_theme <- function() {
  th <- if (requireNamespace("ggprism", quietly = TRUE)) {
    ggprism::theme_prism(base_size = base_size, base_fontface = "plain")
  } else {
    theme_classic(base_size = base_size) +
      theme(axis.line = element_line(linewidth = 0.8), axis.ticks = element_line(linewidth = 0.8),
            axis.ticks.length = unit(0.2, "cm"), axis.text = element_text(colour = "black"))
  }
  th + theme(plot.title = element_text(face = "bold"),
             plot.caption = element_text(colour = "grey40", size = base_size - 2, hjust = 0, face = "plain", margin = margin(t = 12)),
             plot.caption.position = "plot",
             strip.background = element_blank(), strip.text = element_text(face = "bold"))
}

# FIGURE CAPTIONS: set draw_caption <- TRUE to print captions inside the figures.
# By default they are NOT drawn (the LaTeX caption carries the text); the text is written to figure_captions.txt.
draw_caption <- FALSE
#  every plot title starts with the characterization technique; SD bars/ribbons are translucent (alpha below)
sd_alpha   <- 0.4
tech_lab   <- c(viscosity = "Flow curve", thixotropy = "Thixotropy (three-interval test)", gelation = "Gelation")
tm         <- function(tech, m) paste0(tech, " \u2013 ", m)                                    # per-medium title: "<technique> \u2013 <medium>"
tt         <- function(tech, title, w) wrap_parts(paste0(tech, ": ", title), w, chars_per_cm = 4)   # comparison title

# ---- ink key (colour + ink name) -------------------------------------------------------------
.key_mode <- "side"       # "side" = key beside the plot; "none" = no key (comparison plots; the key is a separate image)
key_guides <- function() if (.key_mode == "none") guides(colour = "none", fill = "none") else NULL
ink_key_gg <- function(labs, text_size = 4.4) {
  labs <- labs[labs %in% names(ink_colors)]
  d <- data.frame(y = rev(seq_along(labs)), lab = disp(labs), col = unname(ink_colors[labs]), stringsAsFactors = FALSE)
  W <- 1.3 + 0.24 * max(nchar(d$lab))                # cm: swatch + longest name
  list(plot = ggplot(d) +
         geom_tile(aes(x = 0.55, y = y, fill = col), width = 0.8, height = 0.78, colour = "black", linewidth = 0.4) +
         geom_text(aes(x = 1.1, y = y, label = lab), hjust = 0, size = text_size) +
         scale_fill_identity() + coord_cartesian(xlim = c(0, W), ylim = c(0.45, length(labs) + 0.55), expand = FALSE) + theme_void(),
       w = W, h = 0.9 * length(labs) + 0.1)
}
save_key <- function(labs, name, folder = NULL) {
  k <- ink_key_gg(labs); if (!length(labs[labs %in% names(ink_colors)])) return(invisible())
  d <- file.path(out_dir, if (is.null(folder)) result_folder(name) else folder); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(d, paste0(name, ".png")), k$plot, width = k$w, height = k$h, units = "cm", dpi = 300, bg = "white")
  invisible()
}
.cap_env <- new.env()
cap_out <- function(parts, w_cm, ...) { txt <- wrap_parts(parts, w_cm, ...); .cap_env$last <- txt; if (draw_caption) txt else NULL }
# wrap caption/title text to the width of the saved figure (so it is not cut off)
wrap_parts <- function(parts, w_cm, chars_per_cm = 5) {
  paste(vapply(parts, function(p) paste(strwrap(p, width = max(25, floor(w_cm * chars_per_cm))), collapse = "\n"), ""),
        collapse = "\n")
}
stat_legend <- function() paste0(
  "Student's t-test between two inks; one-way ANOVA among three or more. ",
  "* p < 0.05; ** p < 0.01; *** p < 0.001; **** p < 0.0001.")

# one-way ANOVA from mean, SD and n of each group (three or more groups)
anova_p_summary <- function(m, s, n) {
  k <- length(m); if (k < 3 || any(is.na(m)) || any(is.na(s))) return(NA_real_)
  N <- sum(n); gm <- sum(n * m) / N
  ssb <- sum(n * (m - gm)^2); ssw <- sum((n - 1) * s^2)
  pf((ssb / (k - 1)) / (ssw / (N - k)), k - 1, N - k, lower.tail = FALSE)
}
fmt_p <- function(p) if (is.na(p)) NA_character_ else if (p < 0.0001) "< 0.0001" else paste0("= ", signif(p, 2))
# caption sentence with the ANOVA of every medium (only media with three or more inks)
anova_by_medium <- function(df, mcol, scol, ncol = NULL) {
  if (!"medium" %in% names(df)) df$medium <- "all"
  out <- character(0)
  for (m in unique(as.character(df$medium))) {
    d <- df[as.character(df$medium) == m, , drop = FALSE]
    n <- if (!is.null(ncol) && ncol %in% names(d)) d[[ncol]] else rep(3, nrow(d))
    p <- anova_p_summary(d[[mcol]], d[[scol]], n)
    if (!is.na(p)) out <- c(out, paste0(if (length(unique(df$medium)) > 1) paste0(m, " ") else "", "p ", fmt_p(p)))
  }
  if (!length(out)) return(NULL)
  paste0("One-way ANOVA across inks: ", paste(out, collapse = "; "), ".")
}

# caption sentence for comparisons that are not significant (they get no bracket)
ns_caption <- function(tab, strip = TRUE) {
  if (is.null(tab) || !nrow(tab)) return(NULL)
  if (!any(tab$significance == "ns")) return(NULL)
  "Other comparisons were not significant (p \u2265 0.05)."
}
tab_in <- function(tab, col, v) if (is.null(tab)) NULL else tab[as.character(tab[[col]]) == v, , drop = FALSE]

# map alternative ink names (see 'ink_aliases' in the options) to the name used in 'inks'
canon_ink <- function(ink) {
  i <- match(ink, names(ink_aliases))
  ifelse(is.na(i), ink, unname(ink_aliases[i]))
}

sig_symbol <- function(p) {
  if (is.na(p)) "NA" else if (p < 0.0001) "****" else if (p < 0.001) "***" else
    if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}
# two-sample t-test from mean, SD and n
t_from_summary <- function(m1, s1, n1, m2, s2, n2, var_equal = TRUE) {
  if (var_equal) {
    df <- n1 + n2 - 2
    se <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / df * (1 / n1 + 1 / n2))
  } else {
    se <- sqrt(s1^2 / n1 + s2^2 / n2)
    df <- se^4 / ((s1^2 / n1)^2 / (n1 - 1) + (s2^2 / n2)^2 / (n2 - 1))
  }
  tval <- (m1 - m2) / se
  c(t = tval, df = df, p = 2 * pt(-abs(tval), df))
}

#  Significance brackets.
#   pairs    : list of c(label_a, label_b)
#   summ     : data frame with columns label (FACTOR, levels = plotting order),
#              <group_col> (panel) and top (mean + SD of the bar)
#   test_fun : function(label_a, label_b) -> c(t=, df=, p=) or NULL
#   rel      : TRUE -> gap/step/tick are fractions of the tallest bar of the panel
#  audit: test EVERY pair of inks within a printing medium, write the full table and
#  warn about significant pairs that are not in the pre-specified bracket list
audit_all_pairs <- function(summ, test_fun, group_col, pairs, file) {
  lab <- as.character(summ$label); grp <- as.character(summ[[group_col]])
  out <- list()
  for (g in unique(grp)) {
    ls <- unique(lab[grp == g]); if (length(ls) < 2) next
    cb <- utils::combn(ls, 2)
    for (k in seq_len(ncol(cb))) {
      tt <- test_fun(cb[1, k], cb[2, k]); if (is.null(tt) || is.na(tt["p"])) next
      out[[length(out) + 1]] <- data.frame(group = g, a = cb[1, k], b = cb[2, k],
                                           t = unname(tt["t"]), df = unname(tt["df"]), p_value = unname(tt["p"]),
                                           significance = sig_symbol(tt["p"]), stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) return(invisible(NULL))
  tab <- do.call(rbind, out)
  write.csv(tab, out_path(file), row.names = FALSE)
  key <- vapply(pairs, function(p) paste(sort(p), collapse = " | "), character(1))
  mis <- tab[tab$p_value < 0.05 &
               !(paste(pmin(tab$a, tab$b), pmax(tab$a, tab$b), sep = " | ") %in% key), ]
  if (nrow(mis)) {
    message("Significant pairs NOT drawn on the plot (see ", file, "):")
    for (i in seq_len(nrow(mis)))
      message("   ", mis$a[i], " vs ", mis$b[i], "  p = ", signif(mis$p_value[i], 3))
  } else message("All significant within-medium pairs are on the plot (", file, ").")
  invisible(tab)
}

build_brackets <- function(pairs, summ, test_fun, group_col, gap, step, tick, rel = FALSE) {
  br <- list(); tb <- list()
  lab <- as.character(summ$label)
  grp <- as.character(summ[[group_col]])
  for (pr in pairs) {
    if (!all(pr %in% lab)) { message("Comparison skipped (no data): ", pr[1], " vs ", pr[2]); next }
    ga <- grp[lab == pr[1]][1]; gb <- grp[lab == pr[2]][1]
    if (ga != gb) { message("Skipped (different panels): ", pr[1], " vs ", pr[2]); next }
    tt <- test_fun(pr[1], pr[2])
    if (is.null(tt) || is.na(tt["p"])) { message("Test not possible (n < 2): ", pr[1], " vs ", pr[2]); next }
    present <- levels(droplevels(summ$label[grp == ga]))
    row <- data.frame(g = ga, x1 = match(pr[1], present), x2 = match(pr[2], present),
                      p = unname(tt["p"]), sig = sig_symbol(tt["p"]), stringsAsFactors = FALSE)
    tab <- data.frame(g = ga, a = pr[1], b = pr[2], t = unname(tt["t"]), df = unname(tt["df"]),
                      p_value = unname(tt["p"]), significance = sig_symbol(tt["p"]), stringsAsFactors = FALSE)
    names(row)[1] <- group_col; names(tab)[1] <- group_col
    br[[length(br) + 1]] <- row; tb[[length(tb) + 1]] <- tab
  }
  if (!length(br)) return(list(br = NULL, tab = NULL))
  br <- do.call(rbind, br); tab <- do.call(rbind, tb)
  br <- br[br$sig != "ns", , drop = FALSE]                 # non-significant: no bracket (reported in the caption)
  if (!nrow(br)) return(list(br = NULL, tab = tab))
  br$span <- abs(br$x2 - br$x1); br$y <- NA_real_; br$tick <- NA_real_
  for (g in unique(br[[group_col]])) {                    # narrow brackets lowest
    idx <- which(br[[group_col]] == g); idx <- idx[order(br$span[idx])]
    top <- max(summ$top[grp == g], na.rm = TRUE)
    f <- if (rel) abs(top) else 1
    br$y[idx]    <- top + gap * f + (seq_along(idx) - 1) * step * f
    br$tick[idx] <- tick * f
  }
  list(br = br, tab = tab)
}

#  Bar plot (mean ± SD, optional replicate dots, optional significance brackets)
make_bar_plot <- function(df, mean_col, sd_col, ylab, caption, strip_medium = TRUE,
                          rep_df = NULL, rep_col = NULL, br_df = NULL,
                          ylim_zoom = NULL, ybreaks = waiver(), room = 0,
                          show_values = FALSE, value_fmt = "%.1f", value_size = 2.3, value_gap = 0.06) {
  df$.m <- df[[mean_col]]; df$.s <- df[[sd_col]]
  has_br <- !is.null(br_df) && nrow(br_df) > 0
  p <- ggplot(df, aes(x = label, y = .m, fill = label)) +
    geom_col(colour = "black", linewidth = 0.4, width = 0.7) +
    geom_errorbar(aes(ymin = .m - .s, ymax = .m + .s), width = 0.25, linewidth = 0.5, alpha = sd_alpha)
  if (!is.null(rep_df) && nrow(rep_df) > 0) {
    rep_df$.v <- rep_df[[rep_col]]
    p <- p + geom_point(data = rep_df, aes(x = label, y = .v), inherit.aes = FALSE,
                        shape = 21, fill = "white", colour = "black", size = 1.6,
                        position = position_jitter(width = 0.1, height = 0))
  }
  if (has_br) {
    p <- p +
      geom_segment(data = br_df, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
      geom_segment(data = br_df, aes(x = x1, xend = x1, y = y, yend = y - tick), inherit.aes = FALSE, linewidth = 0.4) +
      geom_segment(data = br_df, aes(x = x2, xend = x2, y = y, yend = y - tick), inherit.aes = FALSE, linewidth = 0.4) +
      geom_text(data = br_df, aes(x = (x1 + x2) / 2, y = y, label = sig),
                inherit.aes = FALSE, vjust = -0.2, size = 3.5)
  }
  if (isTRUE(show_values)) {
    #  the value sits above the SD bar AND above every replicate dot, with a margin
    vt <- df$.m + ifelse(is.finite(df$.s), df$.s, 0)
    if (!is.null(rep_df) && nrow(rep_df) > 0) {
      rmax <- tapply(rep_df[[rep_col]], as.character(rep_df$label), max, na.rm = TRUE)
      vt <- pmax(vt, unname(rmax[as.character(df$label)]), na.rm = TRUE)
    }
    base0 <- if (!is.null(ylim_zoom)) ylim_zoom else 0          # gap relative to the visible axis range
    df$.vy <- vt + value_gap * (max(vt, na.rm = TRUE) - base0)
    p <- p + geom_text(data = df, aes(x = label, y = .vy, label = sprintf(value_fmt, .m)),
                       inherit.aes = FALSE, angle = 90, hjust = 0, vjust = 0.5,
                       size = value_size, colour = "black")
  }
  p <- p + scale_fill_manual(values = ink_colors, name = NULL, labels = disp,
                             guide = if (.key_mode == "none") "none" else guide_legend(keywidth = unit(0.45, "cm"), keyheight = unit(0.45, "cm"))) +
    scale_x_discrete()
  if (!is.null(ylim_zoom)) {
    ytop <- max(c(df$.m + df$.s, if (has_br) br_df$y), na.rm = TRUE) + room
    if (isTRUE(show_values)) ytop <- max(ytop, max(df$.vy, na.rm = TRUE) + 0.18 * (ytop - ylim_zoom))
    p <- p + scale_y_continuous(breaks = ybreaks) +
      coord_cartesian(ylim = c(ylim_zoom, ytop), expand = FALSE)
  } else {
    p <- p + scale_y_continuous(expand = expansion(mult = c(0, if (isTRUE(show_values)) 0.34 else 0.12)))
  }
  p + labs(x = NULL, y = ylab, caption = caption) + base_theme() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "right")
}

# =============================================================================
# 4. A) SWELLING DEGREE
# =============================================================================
rows <- list()
for (x in inks) {
  w0_all <- x$day0
  w0 <- w0_all[!is.na(w0_all)]
  m0 <- mean(w0); s0 <- sd(w0)
  for (d in days) {
    wt_all <- x[[paste0("day", d)]]
    wt <- wt_all[!is.na(wt_all)]
    if (length(wt) == 0) next
    mt <- mean(wt); st <- sd(wt)
    sd_prop <- 100 * (mt / m0) * sqrt((st / mt)^2 + (s0 / m0)^2)
    sd_dayt <- 100 * st / m0
    if (length(w0_all) == length(wt_all)) {
      v <- (wt_all - w0_all) / w0_all * 100
      sw_pair <- mean(v, na.rm = TRUE); sd_pair <- sd(v, na.rm = TRUE)
    } else { sw_pair <- NA; sd_pair <- NA }
    if (sd_method == "paired") {
      swell <- sw_pair; swell_sd <- sd_pair
    } else {
      swell <- (mt - m0) / m0 * 100
      swell_sd <- if (sd_method == "propagated") sd_prop else sd_dayt
    }
    rows[[length(rows) + 1]] <- data.frame(
      medium = x$medium, ink = x$ink, label = ink_label(x), day = d,
      n_day0 = length(w0), n_day = length(wt),
      mean_w0 = m0, sd_w0 = s0, cv_w0_pct = s0 / m0 * 100,
      mean_wt = mt, sd_wt = st, cv_wt_pct = st / mt * 100,
      swelling_pct = swell, swelling_sd = swell_sd,
      sd_propagated = sd_prop, sd_day_t_only = sd_dayt, sd_paired = sd_pair)
  }
}
res <- do.call(rbind, rows)
res$medium <- factor(res$medium, levels = all_media)
res$label  <- factor(res$label,  levels = all_labels)

write.csv(res, out_path("swelling_degree_summary.csv"), row.names = FALSE)
write.csv(res[, c("medium", "ink", "day", "mean_w0", "sd_w0", "cv_w0_pct", "mean_wt", "sd_wt", "cv_wt_pct",
                  "swelling_pct", "sd_propagated", "sd_day_t_only", "sd_paired")],
          out_path("swelling_sd_breakdown.csv"), row.names = FALSE)
message("Swelling SD breakdown (sd_method = '", sd_method, "'): CV = SD / mean of the weights")
print(res[, c("label", "day", "cv_w0_pct", "cv_wt_pct", "swelling_pct", "sd_propagated", "sd_day_t_only", "sd_paired")],
      digits = 3, row.names = FALSE)

make_swelling_plot <- function(df, strip_medium = TRUE, cap_w = 12) {
  if (isTRUE(report_mass_loss)) df$swelling_pct <- -df$swelling_pct   # mass loss = -(swelling)
  y_lab <- if (isTRUE(report_mass_loss)) "Mass loss (%)" else "Swelling degree (%)"
  ggplot(df, aes(x = day, y = swelling_pct, colour = label, group = label)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_line(linewidth = 0.6) +
    geom_errorbar(aes(ymin = swelling_pct - swelling_sd, ymax = swelling_pct + swelling_sd),
                  width = 0.2, linewidth = 0.5, alpha = sd_alpha) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = ink_colors, name = "Ink", labels = strip_lab(strip_medium)) +
    scale_x_continuous(breaks = days, labels = days, expand = expansion(add = 0.6)) +
    labs(x = "Time (days)", y = y_lab,
         caption = cap_out("Data shown as mean ± SD (n = 3).", cap_w)) +
    base_theme() + theme(legend.position = "right") + key_guides()
}

for (m in levels(res$medium)) {
  df <- droplevels(subset(res, medium == m)); if (nrow(df) == 0) next
  p  <- make_swelling_plot(df, cap_w = fig_width) + ggtitle(tm("Swelling degree", m))
  save_plot(p, paste0("swelling_", safe_name(m)), fig_width, fig_height)
  print(p)
}
p_all <- make_swelling_plot(res, strip_medium = FALSE, cap_w = 3 * fig_width) + facet_wrap(~ medium, nrow = 1) + ggtitle("Swelling degree")
save_plot(p_all, "swelling_all_media", 3 * fig_width, fig_height + 2)
print(p_all)

# =============================================================================
# 5. B) WATER CONTENT
# =============================================================================
wc_rows <- list(); wc_rep <- list()
for (x in water) {
  lab <- paste0(x$ink, " (", x$medium, ")")
  if (!lab %in% all_labels) { warning("Water-content ink not found in 'inks' (no colour): ", lab); next }
  if (wc_method == "paired") {
    ok <- !is.na(x$wo) & !is.na(x$wd)
    if (!any(ok)) next
    vals <- (x$wo[ok] - x$wd[ok]) / x$wo[ok] * 100
    m <- mean(vals); s <- sd(vals); n <- length(vals)
    wc_rep[[length(wc_rep) + 1]] <- data.frame(medium = x$medium, label = lab, value = vals)
  } else {
    wo <- x$wo[!is.na(x$wo)]; wd <- x$wd[!is.na(x$wd)]
    if (length(wo) == 0 || length(wd) == 0) next
    mo <- mean(wo); md <- mean(wd); so <- sd(wo); sdd <- sd(wd)
    m <- (mo - md) / mo * 100
    s <- 100 * (md / mo) * sqrt((sdd / md)^2 + (so / mo)^2)
    n <- min(length(wo), length(wd))
  }
  wc_rows[[length(wc_rows) + 1]] <- data.frame(
    medium = x$medium, ink = x$ink, label = lab, n = n,
    water_content_pct = m, water_content_sd = s)
}

wc_test <- function(a, b) {
  A <- wc[as.character(wc$label) == a, ]; B <- wc[as.character(wc$label) == b, ]
  if (nrow(A) == 0 || nrow(B) == 0 || A$n < 2 || B$n < 2) return(NULL)
  t_from_summary(A$water_content_pct, A$water_content_sd, A$n,
                 B$water_content_pct, B$water_content_sd, B$n, t_var_equal)
}

if (length(wc_rows) == 0) {
  message("No water-content data entered yet — fill in 'water' (section 1b).")
} else {
  wc <- do.call(rbind, wc_rows)
  wc$medium <- factor(wc$medium, levels = all_media)
  wc$label  <- factor(wc$label,  levels = all_labels)
  wc <- wc[order(wc$label), ]
  write.csv(wc, out_path("water_content_summary.csv"), row.names = FALSE)
  print(wc[, c("medium", "ink", "water_content_pct", "water_content_sd")], digits = 4)
  
  wc_rep <- if (length(wc_rep)) do.call(rbind, wc_rep) else NULL
  if (!is.null(wc_rep)) {
    wc_rep$medium <- factor(wc_rep$medium, levels = all_media)
    wc_rep$label  <- factor(wc_rep$label,  levels = all_labels)
  }
  if (!show_points) wc_rep <- NULL
  
  #  brackets start above the printed values (panel): top = mean + SD (and replicates) + label room
  wc_hi <- wc$water_content_pct + ifelse(is.finite(wc$water_content_sd), wc$water_content_sd, 0)
  if (!is.null(wc_rep) && nrow(wc_rep)) {
    wc_rmx <- tapply(wc_rep$value, as.character(wc_rep$label), max, na.rm = TRUE)
    wc_hi  <- pmax(wc_hi, unname(wc_rmx[as.character(wc$label)]), na.rm = TRUE)
  }
  wc_sum <- data.frame(medium = wc$medium, label = wc$label, top = wc_hi + wc_value_room)
  bw <- build_brackets(wc_comparisons, wc_sum, wc_test, "medium",
                       gap = wc_bracket_gap, step = wc_bracket_step, tick = wc_bracket_tick, rel = FALSE)
  if (!is.null(bw$br)) bw$br$medium <- factor(bw$br$medium, levels = all_media)
  audit_all_pairs(wc_sum, wc_test, "medium", wc_comparisons, "water_content_all_pairwise.csv")
  if (!is.null(bw$tab)) {
    write.csv(bw$tab, out_path("water_content_ttests.csv"), row.names = FALSE)
    print(bw$tab, digits = 3)
  }
  
  wc_plot <- function(df, rp, br, cap_w, nsn = NULL) {
    make_bar_plot(df, "water_content_pct", "water_content_sd", ylab = "Water content (%)",
                  caption = cap_out(c("Mean ± SD (n = 3).", stat_legend(), anova_by_medium(df, "water_content_pct", "water_content_sd", "n"), nsn), cap_w),
                  rep_df = rp, rep_col = "value", br_df = br,
                  ylim_zoom = wc_y_min, ybreaks = seq(wc_y_min, 100, 2), room = wc_label_room)
  }
  for (m in levels(wc$medium)) {
    df <- droplevels(subset(wc, medium == m)); if (nrow(df) == 0) next
    rp <- if (!is.null(wc_rep)) droplevels(subset(wc_rep, medium == m)) else NULL
    bp <- if (!is.null(bw$br)) bw$br[bw$br$medium == m, ] else NULL
    w  <- max(fig_width * 0.6, 3 + 1.8 * nrow(df))
    p  <- wc_plot(df, rp, bp, w, ns_caption(tab_in(bw$tab, "medium", m))) + ggtitle(tm("Water content", m))
    save_plot(p, paste0("water_content_", safe_name(m)), w + key_w, fig_height + 1.5)
    print(p)
  }
  p_all_wc <- wc_plot(wc, wc_rep, bw$br, 3 * fig_width, ns_caption(bw$tab, strip = FALSE)) + facet_grid(. ~ medium, scales = "free_x", space = "free_x") + ggtitle("Water content")
  save_plot(p_all_wc, "water_content_all_media", 3 * fig_width + key_w, fig_height + 1.5)
  print(p_all_wc)
}

# =============================================================================
# 6. C) COMPRESSION TESTS
# =============================================================================
#  csv columns: SetName, Cycle, Time_S, Size_mm, Displacement_mm, Force_N
#  Only the first "Compress" cycle is used ("Recover" rows are ignored).
#  stress (kPa) = F (N) / A (mm²) * 1000        strain = dL (mm) / L0 (mm)
#  Young's modulus E = slope of stress vs strain in 0-15 % strain (kPa)
#
#  Data frames created (also saved as csv):
#   comp_raw    : every compression point of every file, uncorrected
#   comp_onsets : detected start of loading per replicate (what was removed)
#   comp_corr   : corrected data (noise removed, displacement re-zeroed,
#                 force offset subtracted) + strain and stress
#   comp_fits   : linear fit 0-15 % per replicate (E, intercept, R²)
#   comp_avg    : mean ± SD stress-strain curve per ink (n replicates)
#   E_sum       : Young's modulus per ink, mean ± SD of the replicates

find_onset <- function(force) {
  n  <- length(force)
  nb <- max(3, min(onset_n_base, floor(n / 3)))
  base  <- median(force[1:nb])
  noise <- mad(force[1:nb])
  if (!is.finite(noise) || noise == 0) noise <- sd(force[1:nb])
  if (!is.finite(noise) || noise == 0) noise <- 1e-6
  thr   <- base + onset_k * noise
  r     <- rle(force > thr)
  ends  <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  ok    <- which(r$values & r$lengths >= onset_n_consec)
  if (length(ok) == 0) return(list(idx = 1L, thr = thr, found = FALSE))
  list(idx = starts[ok[1]], thr = thr, found = TRUE)
}

fmt_eq <- function(a, b) sprintf("y = %.2fx %s %.3f", a, ifelse(b < 0, "\u2212", "+"), abs(b))

n_cap <- function(df) paste0("Mean ± SD (n = ", paste(unique(df$n), collapse = "/"), ")")

# mean stress-strain curve as dots (one dot per acquired/averaged strain step), N/mm² vs mm/mm
make_curve_plot <- function(df, strip_medium = TRUE, cap_w = 12) {
  df$s_mean <- df$stress_mean / 1000                    # kPa -> N/mm²
  df$s_sd   <- df$stress_sd / 1000
  p <- ggplot(df, aes(x = strain, y = s_mean, colour = label)) +
    annotate("rect", xmin = 0, xmax = lin_max_strain, ymin = -Inf, ymax = Inf,
             fill = "grey85", alpha = 0.4)
  if (curve_show_sd)
    p <- p + geom_errorbar(aes(ymin = s_mean - s_sd, ymax = s_mean + s_sd),
                           width = 0, alpha = 0.35, linewidth = 0.3)
  p + geom_point(size = curve_point_size) +
    scale_colour_manual(values = ink_colors, name = "Ink", labels = strip_lab(strip_medium)) +
    labs(x = "Strain (mm/mm)", y = expression("Stress (N/mm"^2*")"),
         caption = cap_out(paste0(n_cap(df), ". Grey area: linear region (0–", lin_max_strain * 100, " % strain)."), cap_w)) +
    base_theme() + theme(legend.position = "right") + key_guides()
}

run_compression <- function() {
  if (!dir.exists(comp_dir)) {
    message("Compression folder not found (is the drive mounted?): ", comp_dir); return(invisible())
  }
  files <- list.files(comp_dir, pattern = "_[0-9]+Data\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) { message("No *_0XData.csv files found in: ", comp_dir); return(invisible()) }
  
  # file name "<ink>_<medium>_<rep>Data.csv"; '_' inside the ink name stands for '.'
  # (30k-c13_3MC0Gel_PM3_01Data.csv -> ink 30k-c13.3MC0Gel, medium PM3, replicate 1)
  parsed <- do.call(rbind, lapply(files, function(f) {
    b <- sub("Data\\.csv$", "", basename(f))
    m <- regmatches(b, regexec("^(.*)_([^_]+)_([0-9]+)$", b))[[1]]
    if (length(m) == 0) { message("File name not understood, skipped: ", basename(f)); return(NULL) }
    data.frame(file = f, ink = canon_ink(gsub("_", ".", m[2])), medium = m[3],
               rep_n = as.integer(m[4]), stringsAsFactors = FALSE)
  }))
  if (is.null(parsed)) { message("Could not parse any file names."); return(invisible()) }
  
  dims <- list()
  for (x in comp_dims) dims[[paste0(x$ink, " (", x$medium, ")")]] <- x
  man <- if (length(onset_manual)) do.call(rbind, lapply(onset_manual, as.data.frame)) else NULL
  
  # ---------------- (a) read everything -> comp_raw, comp_onsets, comp_corr ----------------
  raw_l <- list(); ons_l <- list(); corr_l <- list()
  for (i in seq_len(nrow(parsed))) {
    pr  <- parsed[i, ]
    lab <- paste0(pr$ink, " (", pr$medium, ")")
    if (!lab %in% all_labels) { warning("Not in 'inks' (no colour), skipped: ", lab); next }
    dm   <- dims[[lab]]
    diam <- if (is.null(dm) || pr$rep_n > length(dm$diameter)) NA else dm$diameter[pr$rep_n]
    hgt  <- if (is.null(dm) || pr$rep_n > length(dm$height))   NA else dm$height[pr$rep_n]
    if (is.na(diam)) { message("No diameter entered for ", lab, " replicate ", pr$rep_n, " -> skipped"); next }
    
    d <- read.csv(pr$file, stringsAsFactors = FALSE, strip.white = TRUE)
    d <- d[grepl("Compress", d$Cycle), ]
    d <- d[d$Cycle == d$Cycle[1], ]                       # first compression cycle only
    d <- d[order(d$Time_S), ]
    if (nrow(d) < 10) { message("Too few points in ", basename(pr$file)); next }
    if (is.na(hgt)) {
      hgt <- d$Size_mm[1]
      message("No height entered for ", lab, " replicate ", pr$rep_n, " -> using first Size_mm = ", hgt)
    }
    
    # start of loading: automatic, or manual override
    on <- find_onset(d$Force_N); method <- if (on$found) "auto" else "none found (kept all)"
    if (!is.null(man)) {
      mm <- man[man$medium == pr$medium & man$ink == pr$ink & man$rep == pr$rep_n, ]
      if (nrow(mm)) { on$idx <- which.min(abs(d$Displacement_mm - mm$onset_mm[1])); method <- "manual" }
    }
    k        <- on$idx
    disp0    <- d$Displacement_mm[k]
    f_offset <- if (k > 1) median(d$Force_N[seq_len(k - 1)]) else 0
    L0       <- switch(comp_L0_source,
                       measured = hgt, measured_minus_removed = hgt - disp0,
                       size_first = d$Size_mm[1], size_contact = d$Size_mm[k], hgt)
    area     <- pi * diam^2 / 4                           # mm²
    
    raw_l[[length(raw_l) + 1]] <- data.frame(
      medium = pr$medium, ink = pr$ink, label = lab, rep_n = pr$rep_n, point = seq_len(nrow(d)),
      Time_S = d$Time_S, Size_mm = d$Size_mm, Displacement_mm = d$Displacement_mm,
      Force_N = d$Force_N, removed = seq_len(nrow(d)) < k, file = basename(pr$file))
    ons_l[[length(ons_l) + 1]] <- data.frame(
      medium = pr$medium, ink = pr$ink, label = lab, rep_n = pr$rep_n, method = method,
      onset_point = k, n_removed = k - 1, onset_disp_mm = disp0, force_offset_N = f_offset,
      height_mm = hgt, size_first_mm = d$Size_mm[1], size_contact_mm = d$Size_mm[k],
      L0_used_mm = L0, L0_source = comp_L0_source, diameter_mm = diam, area_mm2 = area)
    
    keep <- k:nrow(d)
    dc <- d$Displacement_mm[keep] - disp0
    fc <- d$Force_N[keep] - f_offset
    corr_l[[length(corr_l) + 1]] <- data.frame(
      medium = pr$medium, ink = pr$ink, label = lab, rep_n = pr$rep_n,
      Time_S = d$Time_S[keep], Displacement_corr_mm = dc, Force_corr_N = fc,
      strain = dc / L0, strain_pct = dc / L0 * 100,
      stress_kPa = fc / area * 1000)                      # N/mm² -> kPa
  }
  if (length(corr_l) == 0) {
    message("No compression data could be processed — check 'comp_dims' (section 1c)."); return(invisible())
  }
  
  fix <- function(df) {
    df$medium <- factor(df$medium, levels = all_media)
    df$label  <- factor(df$label,  levels = all_labels)
    if (!is.null(df$rep_n)) df$rep_lab <- paste("Replicate", df$rep_n)
    df
  }
  comp_raw    <- fix(do.call(rbind, raw_l))
  comp_onsets <- fix(do.call(rbind, ons_l))
  comp_corr   <- fix(do.call(rbind, corr_l))
  comp_onsets <- comp_onsets[order(comp_onsets$label, comp_onsets$rep_n), ]
  write.csv(comp_raw,    out_path("compression_raw_data.csv"), row.names = FALSE)
  write.csv(comp_onsets, out_path("compression_onsets.csv"), row.names = FALSE)
  write.csv(comp_corr,   out_path("compression_corrected_data.csv"), row.names = FALSE)
  print(comp_onsets[, c("medium", "ink", "rep_n", "method", "n_removed", "onset_disp_mm", "force_offset_N")], digits = 3)
  
  # ---------------- (b) onset check plots (raw force, removed points in grey) ----------------
  ons_pts <- merge(comp_onsets[, c("label", "rep_n", "onset_disp_mm")],
                   unique(comp_raw[, c("label", "rep_n", "medium", "rep_lab")]))
  for (m in levels(comp_raw$medium)) {
    df <- droplevels(subset(comp_raw, medium == m)); if (nrow(df) == 0) next
    op <- droplevels(subset(ons_pts, medium == m))
    p <- ggplot(df, aes(Displacement_mm, Force_N)) +
      geom_point(aes(colour = removed), size = 0.9) +
      geom_vline(data = op, aes(xintercept = onset_disp_mm), linetype = "dashed", colour = "red") +
      scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "grey70"),
                          labels = c("kept", "removed"), name = NULL) +
      facet_grid(label ~ rep_lab, scales = "free",
                 labeller = labeller(label = function(l) sub(" \\(.*\\)$", "", l))) +
      labs(x = "Raw displacement (mm)", y = "Raw force (N)",
           caption = "Red dashed line: detected start of loading") +
      ggtitle(paste0(m, " — start-of-loading check")) + base_theme() +
      theme(strip.text.y = element_text(angle = 0, hjust = 0), legend.position = "top")
    save_plot(p, paste0("compression_onset_check_", safe_name(m)), 22, 3.5 * nlevels(df$label) + 3)
  }
  
  # ---------------- (c) linear fit 0-15 % per replicate -> E ----------------
  fit_l <- list()
  for (g in split(comp_corr, list(comp_corr$label, comp_corr$rep_n), drop = TRUE)) {
    w <- g$strain >= 0 & g$strain <= lin_max_strain
    if (sum(w) < 3) { message("Too few points in the linear region for ", g$label[1], " rep ", g$rep_n[1]); next }
    fit <- lm(stress_kPa ~ strain, data = g[w, ])
    fit_l[[length(fit_l) + 1]] <- data.frame(
      medium = as.character(g$medium[1]), ink = g$ink[1], label = as.character(g$label[1]),
      rep_n = g$rep_n[1], n_points = sum(w),
      E_kPa = unname(coef(fit)[2]), intercept_kPa = unname(coef(fit)[1]),
      r2 = summary(fit)$r.squared)
  }
  comp_fits <- fix(do.call(rbind, fit_l))
  comp_fits$r2_ok <- comp_fits$r2 > r2_min
  comp_fits <- comp_fits[order(comp_fits$label, comp_fits$rep_n), ]
  write.csv(comp_fits, out_path("compression_fits_per_replicate.csv"), row.names = FALSE)
  print(comp_fits[, c("medium", "ink", "rep_n", "E_kPa", "intercept_kPa", "r2")], digits = 4)
  if (any(!comp_fits$r2_ok)) {
    message("Fits with R² <= ", r2_min, " (check the onset plots / consider a manual onset):")
    print(comp_fits[!comp_fits$r2_ok, c("label", "rep_n", "E_kPa", "r2")], digits = 3)
  }
  
  E_sum <- do.call(rbind, lapply(split(comp_fits, droplevels(comp_fits$label)), function(g)
    data.frame(medium = as.character(g$medium[1]), ink = g$ink[1], label = as.character(g$label[1]),
               n = nrow(g), E_mean = mean(g$E_kPa), E_sd = sd(g$E_kPa))))
  E_sum <- fix(E_sum)
  write.csv(E_sum, out_path("youngs_modulus_summary.csv"), row.names = FALSE)
  print(E_sum[, c("medium", "ink", "n", "E_mean", "E_sd")], digits = 4)
  
  # ---------------- (d) mean ± SD curve per ink (replicates interpolated on a common strain grid) ----------------
  avg_l <- list()
  for (g in split(comp_corr, droplevels(comp_corr$label))) {
    reps <- split(g, g$rep_n)
    smax <- min(sapply(reps, function(r) max(r$strain)))   # common strain range of all replicates
    grid <- seq(0, smax, by = curve_step)
    mat  <- sapply(reps, function(r) approx(r$strain, r$stress_kPa, xout = grid, ties = mean, rule = 2)$y)
    mat  <- matrix(mat, nrow = length(grid))
    avg_l[[length(avg_l) + 1]] <- data.frame(
      medium = as.character(g$medium[1]), label = as.character(g$label[1]),
      strain = grid, strain_pct = grid * 100, n = ncol(mat),
      stress_mean = rowMeans(mat), stress_sd = if (ncol(mat) > 1) apply(mat, 1, sd) else 0)
  }
  comp_avg <- fix(do.call(rbind, avg_l))
  write.csv(comp_avg, out_path("compression_mean_curves.csv"), row.names = FALSE)
  
  # ---------------- plot 1: mean stress-strain curves (dots; strain mm/mm, stress N/mm²) ----------------
  for (m in levels(comp_avg$medium)) {
    df <- droplevels(subset(comp_avg, medium == m)); if (nrow(df) == 0) next
    p <- make_curve_plot(df, cap_w = fig_width + 2) + ggtitle(tm("Compression curves", m))
    save_plot(p, paste0("compression_curves_", safe_name(m)), fig_width + 2, fig_height + 1)
    print(p)
  }
  p_c_all <- make_curve_plot(comp_avg, strip_medium = FALSE, cap_w = 3 * fig_width) +
    facet_wrap(~ medium, nrow = 1) + ggtitle("Compression curves")
  save_plot(p_c_all, "compression_curves_all_media", 3 * fig_width, fig_height + 2)
  print(p_c_all)
  
  # ---------------- plot 2: linear regression 0-15 % on the mean curve, with equation ----------------
  lin <- droplevels(subset(comp_avg, strain <= lin_max_strain + 1e-9))
  lab_l <- list()
  for (g in split(lin, lin$label, drop = TRUE)) {
    f  <- lm(stress_mean ~ strain, data = g)
    es <- E_sum[as.character(E_sum$label) == as.character(g$label[1]), ]
    r2 <- summary(f)$r.squared
    lab_l[[length(lab_l) + 1]] <- data.frame(
      medium = g$medium[1], label = g$label[1],
      txt = paste0(fmt_eq(coef(f)[2], coef(f)[1]), "\n",
                   sprintf("R\u00b2 = %.3f%s", r2, ifelse(r2 > r2_min, "", paste0(" (< ", r2_min, ")"))), "\n",
                   sprintf("E = %.2f \u00b1 %.2f kPa (n = %d)", es$E_mean, es$E_sd, es$n)))
  }
  lin_lab <- fix(do.call(rbind, lab_l))
  
  make_fit_plot <- function(pts, lb) {
    ggplot(pts, aes(x = strain, y = stress_mean)) +
      geom_errorbar(aes(ymin = stress_mean - stress_sd, ymax = stress_mean + stress_sd),
                    width = 0, colour = "grey55", linewidth = 0.3, alpha = sd_alpha) +
      geom_point(aes(fill = label), shape = 21, colour = "black", size = 1.8, stroke = 0.3) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.5) +
      geom_text(data = lb, aes(x = -Inf, y = Inf, label = txt), inherit.aes = FALSE,
                hjust = -0.05, vjust = 1.15, size = 2.8, lineheight = 0.95) +
      facet_wrap(~ label, scales = "free", ncol = 2,
                 labeller = labeller(label = function(l) sub(" \\(.*\\)$", "", l))) +
      scale_fill_manual(values = ink_colors, guide = "none") +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.45))) +
      labs(x = "Strain (mm/mm)", y = "Stress (kPa)",
           caption = cap_out(paste0("Mean ± SD curve; line: linear fit 0–", lin_max_strain * 100,
                                    " % strain. E = mean ± SD of the replicate slopes."), 18)) +
      base_theme()
  }
  for (m in levels(lin$medium)) {
    pts <- droplevels(subset(lin, medium == m)); if (nrow(pts) == 0) next
    lb  <- droplevels(subset(lin_lab, medium == m))
    p <- make_fit_plot(pts, lb) + ggtitle(m)
    save_plot(p, paste0("compression_linear_fits_", safe_name(m)),
              18, 5.5 * ceiling(nlevels(pts$label) / 2) + 2.5)
    print(p)
  }
  
  # ---------------- plot 2b (supplementary): per-replicate fits with equations ----------------
  fitpts <- subset(comp_corr, strain <= lin_max_strain)
  comp_fits$txt <- paste0(fmt_eq(comp_fits$E_kPa, comp_fits$intercept_kPa), "\n",
                          sprintf("R\u00b2 = %.3f", comp_fits$r2))
  for (m in levels(fitpts$medium)) {
    pts <- droplevels(subset(fitpts, medium == m)); if (nrow(pts) == 0) next
    ft  <- droplevels(subset(comp_fits, medium == m))
    p <- ggplot(pts, aes(strain, stress_kPa)) +
      geom_point(aes(fill = label), shape = 21, colour = "black", size = 1.5, stroke = 0.3) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.5) +
      geom_text(data = ft, aes(x = -Inf, y = Inf, label = txt), inherit.aes = FALSE,
                hjust = -0.05, vjust = 1.15, size = 2.6) +
      facet_grid(label ~ rep_lab, scales = "free",
                 labeller = labeller(label = function(l) sub(" \\(.*\\)$", "", l))) +
      scale_fill_manual(values = ink_colors, guide = "none") +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
      labs(x = "Strain (mm/mm)", y = "Stress (kPa)") + ggtitle(paste0(m, " — individual replicates")) +
      base_theme() + theme(strip.text.y = element_text(angle = 0, hjust = 0))
    save_plot(p, paste0("compression_linear_fits_replicates_", safe_name(m)),
              24, 3.8 * nlevels(pts$label) + 2)
  }
  
  # ---------------- plot 3: Young's modulus (mean ± SD, Student's t-test) ----------------
  E_test <- function(a, b) {
    ea <- comp_fits$E_kPa[as.character(comp_fits$label) == a]
    eb <- comp_fits$E_kPa[as.character(comp_fits$label) == b]
    if (length(ea) < 2 || length(eb) < 2) return(NULL)
    tt <- t.test(ea, eb, var.equal = t_var_equal)
    c(t = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value)
  }
  #  brackets start above the printed values: top = highest of (mean + SD, replicates) + label room
  E_hi  <- E_sum$E_mean + ifelse(is.finite(E_sum$E_sd), E_sum$E_sd, 0)
  E_rmx <- tapply(comp_fits$E_kPa, as.character(comp_fits$label), max, na.rm = TRUE)
  E_hi  <- pmax(E_hi, unname(E_rmx[as.character(E_sum$label)]), na.rm = TRUE)
  E_top <- data.frame(medium = E_sum$medium, label = E_sum$label,
                      top = E_hi + 0.25 * max(E_hi, na.rm = TRUE))
  bE <- build_brackets(E_comparisons, E_top, E_test, "medium", gap = 0.04, step = 0.09, tick = 0.02, rel = TRUE)
  if (!is.null(bE$br)) bE$br$medium <- factor(bE$br$medium, levels = all_media)
  audit_all_pairs(E_top, E_test, "medium", E_comparisons, "youngs_modulus_all_pairwise.csv")
  if (!is.null(bE$tab)) {
    write.csv(bE$tab, out_path("youngs_modulus_ttests.csv"), row.names = FALSE)
    print(bE$tab, digits = 3)
  }
  
  E_plot <- function(df, rp, br, cap_w, nsn = NULL) {
    make_bar_plot(df, "E_mean", "E_sd", ylab = "Young's modulus (kPa)",
                  caption = cap_out(c("Mean ± SD (n = 3).", stat_legend(), anova_by_medium(df, "E_mean", "E_sd"), nsn), cap_w),
                  rep_df = rp, rep_col = "E_kPa", br_df = br)
  }
  for (m in levels(E_sum$medium)) {
    df <- droplevels(subset(E_sum, medium == m)); if (nrow(df) == 0) next
    rp <- droplevels(subset(comp_fits, medium == m))
    bp <- if (!is.null(bE$br)) bE$br[bE$br$medium == m, ] else NULL
    w  <- max(fig_width * 0.6, 3 + 1.8 * nrow(df))
    p  <- E_plot(df, rp, bp, w, ns_caption(tab_in(bE$tab, "medium", m))) + ggtitle(tm("Young's modulus", m))
    save_plot(p, paste0("youngs_modulus_", safe_name(m)), w + key_w, fig_height + 1.5)
    print(p)
  }
  p_E_all <- E_plot(E_sum, comp_fits, bE$br, 3 * fig_width, ns_caption(bE$tab, strip = FALSE)) +
    facet_grid(. ~ medium, scales = "free_x", space = "free_x") + ggtitle("Young's modulus")
  save_plot(p_E_all, "youngs_modulus_all_media", 3 * fig_width + key_w, fig_height + 1.5)
  print(p_E_all)
  
  # keep the data frames and tests in the R session (used by section 7)
  list2env(list(comp_raw = comp_raw, comp_onsets = comp_onsets, comp_corr = comp_corr,
                comp_fits = comp_fits, comp_avg = comp_avg, E_sum = E_sum, E_test = E_test),
           envir = globalenv())
  invisible()
}
run_compression()

# =============================================================================
# 6b. D) PRINTABILITY FACTOR
# =============================================================================
#  Per pore:  C = 4*pi*A / L^2 (circularity);  Pr = (pi/4) / C = L^2 / (16 * A)
#   Pr = 1 square pore | Pr < 1 circular pore (under-gelled) | Pr > 1 over-gelled
#  Processing:
#   1. pr_raw    : every pore of every mesh (area, perimeter, Pr)
#   2. per mesh keep the 'pr_n_pores' pores closest to Pr = 1 (meshes with fewer
#      pores keep all of them and are flagged)
#   3. pr_points : the points that are plotted. The pores of the replicate meshes
#      are NOT the same pores, so they are randomly paired (shuffled, fixed seed)
#      and every point is the MEAN of the replicate meshes.
#   4. pr_status : per ink: meshes obtained, pores per mesh, flags
#      0 meshes -> shown as an empty "not printable" slot; < 3 meshes -> plotted, flagged †;
#      a mesh with < N measurable pores -> flagged ‡.

pr_flag_mesh <- "\u2020"
pr_flag_pore <- "\u2021"
norm_key <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))

read_results <- function(f) {
  l1 <- readLines(f, n = 1, warn = FALSE)
  if (length(l1) == 0) return(NULL)
  sep  <- if (grepl("\t", l1)) "\t" else if (grepl(",", l1)) "," else ""
  skip <- if (grepl("[A-Za-z]", l1)) 1 else 0            # header line
  d <- tryCatch(read.table(f, header = FALSE, sep = sep, skip = skip, fill = TRUE, quote = "",
                           comment.char = "", stringsAsFactors = FALSE, strip.white = TRUE),
                error = function(e) NULL)
  if (is.null(d) || ncol(d) < 3 || nrow(d) == 0) return(NULL)
  num <- function(v) suppressWarnings(as.numeric(if (sep == ",") as.character(v) else gsub(",", ".", as.character(v))))
  out <- data.frame(area = num(d[[2]]), perimeter = num(d[[3]]))
  out[is.finite(out$area) & is.finite(out$perimeter) & out$area > 0 & out$perimeter > 0, ]
}

make_pr_plot <- function(pts, lev, np = NULL, strip_medium = TRUE, cap_w = 12) {
  #  pts : processed points (label, medium, Pr_mean, n_meshes, suffix)
  #  lev : x-axis order (labels, printed + not printable inks)
  #  np  : data frame (label, medium) of inks with 0 meshes -> empty slot "not printable"
  pts$status <- factor(ifelse(pts$n_meshes < pr_n_meshes, "fewer meshes", "all meshes"),
                       levels = c("all meshes", "fewer meshes"))
  any_mesh <- any(pts$n_meshes < pr_n_meshes)
  any_pore <- any(grepl(pr_flag_pore, pts$suffix, fixed = TRUE))
  sfx <- unique(data.frame(label = as.character(pts$label), suffix = pts$suffix, stringsAsFactors = FALSE))
  sfx <- setNames(sfx$suffix, sfx$label)
  lab_fun <- function(l) { s_ <- unname(sfx[l]); s_[is.na(s_)] <- ""; s_ }       # only the flags (e.g. the mesh/pore marks) stay under the dots; the ink names are in the key
  d <- pts[, c("label", "medium", "Pr_mean", "status")]
  has_np <- !is.null(np) && nrow(np) > 0
  if (has_np) {
    np <- data.frame(label = as.character(np$label), medium = as.character(np$medium), stringsAsFactors = FALSE)
    d <- rbind(d, data.frame(label = np$label, medium = np$medium, Pr_mean = NA_real_,
                             status = factor("all meshes", levels = levels(d$status))))
  }
  d$label  <- factor(as.character(d$label), levels = lev)
  d$medium <- factor(as.character(d$medium), levels = all_media)
  cap <- c(paste0("Each dot is the mean of the replicate meshes of the ink for one pore (pores of different meshes randomly paired); ",
                  "per mesh, the ", pr_n_pores, " pores closest to Pr = 1 were used."),
           paste0("Grey region: range of adequate ink printability (", pr_adequate[1], "\u2013", pr_adequate[2],
                  ") (Ouyang et al. 2016). Dotted line: Pr = 1 (perfectly square pore); Pr < 1: circular pore (under-gelled); ",
                  "Pr > 1: over-gelled pore."))
  if (any_mesh) cap <- c(cap, paste0(pr_flag_mesh, " fewer than ", pr_n_meshes, " meshes."))
  if (any_pore) cap <- c(cap, paste0(pr_flag_pore, " a mesh with fewer than ", pr_n_pores, " pores."))
  if (has_np) cap <- c(cap, "Not printable: no mesh.")
  p <- ggplot(d, aes(x = label, y = Pr_mean, fill = label, shape = status)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = pr_adequate[1], ymax = pr_adequate[2],
             fill = "grey85", alpha = 0.6) +
    geom_hline(yintercept = 1, linetype = "dotted", colour = "grey30") +
    geom_point(colour = "black", size = 2, stroke = 0.3, na.rm = TRUE,
               position = position_jitter(width = 0.18, height = 0, seed = pr_seed))
  if (has_np)
    p <- p + geom_text(data = d[is.na(d$Pr_mean), ], aes(x = label, y = 1, label = "not printable"),
                       inherit.aes = FALSE, angle = 90, colour = "grey35", fontface = "italic", size = 3.3)
  p + scale_fill_manual(values = ink_colors, name = NULL, labels = disp,
                        guide = if (.key_mode == "none") "none" else guide_legend(order = 1, keywidth = unit(0.45, "cm"), keyheight = unit(0.45, "cm"),
                                                                                  override.aes = list(shape = 21, size = 3.2, colour = "black"))) +
    scale_shape_manual(values = c("all meshes" = 21, "fewer meshes" = 24), name = NULL,
                       labels = c("all meshes" = paste0(pr_n_meshes, " meshes"),
                                  "fewer meshes" = paste0("fewer than ", pr_n_meshes, " meshes")),
                       guide = if (any_mesh) "legend" else "none") +
    scale_x_discrete(labels = lab_fun) +
    labs(x = NULL, y = "Printability factor (Pr)", caption = cap_out(cap, cap_w)) +
    base_theme() + theme(axis.text.x = element_text(size = base_size - 1, margin = margin(t = 1)), axis.ticks.x = element_blank(),
                         legend.position = if (.key_mode == "none") "top" else "right")
}

run_printability <- function() {
  if (!dir.exists(pr_dir)) {
    message("Printability folder not found (is the drive mounted?): ", pr_dir); return(invisible())
  }
  ink_keys <- setNames(all_labels, norm_key(sapply(inks, function(x) paste0(x$ink, x$medium))))
  for (a in names(ink_aliases)) for (x in inks)            # folders named with an alternative ink name
    if (x$ink == ink_aliases[[a]]) ink_keys[norm_key(paste0(a, x$medium))] <- ink_label(x)
  dirs <- list.dirs(pr_dir, recursive = FALSE, full.names = TRUE)
  raw_l <- list(); fold_l <- list()
  for (d in dirs) {
    b   <- trimws(basename(d))
    lab <- unname(ink_keys[norm_key(b)])
    if (is.na(lab)) { message("Printability folder not matched to an ink (skipped): ", b); next }
    x   <- inks[[match(lab, all_labels)]]
    fl  <- list.files(d, pattern = "^Results[_ ]?[0-9]+\\.txt$", ignore.case = TRUE,
                      recursive = TRUE, full.names = TRUE)
    n_ok <- 0
    for (f in fl) {
      rep_n <- as.integer(sub("^Results[_ ]?([0-9]+)\\.txt$", "\\1", basename(f), ignore.case = TRUE))
      dat <- read_results(f)
      if (is.null(dat) || nrow(dat) == 0) { message("No valid pores in ", lab, " / ", basename(f)); next }
      n_ok <- n_ok + 1
      raw_l[[length(raw_l) + 1]] <- data.frame(
        medium = x$medium, ink = x$ink, label = lab, rep_n = rep_n, pore = seq_len(nrow(dat)),
        area = dat$area, perimeter = dat$perimeter,
        Pr = dat$perimeter^2 / (16 * dat$area), file = basename(f))
    }
    fold_l[[length(fold_l) + 1]] <- data.frame(medium = x$medium, ink = x$ink, label = lab,
                                               n_files = length(fl), n_meshes = n_ok)
  }
  if (length(raw_l) == 0) { message("No printability data could be read."); return(invisible()) }
  
  fix <- function(df) {
    df$medium <- factor(df$medium, levels = all_media)
    df$label  <- factor(df$label,  levels = all_labels)
    df
  }
  pr_raw <- fix(do.call(rbind, raw_l))
  
  # ---- keep the pores closest to Pr = 1 in each mesh ----
  pr_raw$abs_dev  <- abs(pr_raw$Pr - 1)
  pr_raw$rank     <- ave(pr_raw$abs_dev, pr_raw$label, pr_raw$rep_n,
                         FUN = function(v) rank(v, ties.method = "first"))
  pr_raw$selected <- pr_raw$rank <= pr_n_pores
  
  # ---- mesh table ----
  mesh <- do.call(rbind, lapply(split(pr_raw, list(pr_raw$label, pr_raw$rep_n), drop = TRUE), function(g)
    data.frame(label = as.character(g$label[1]), rep_n = g$rep_n[1], n_pores = nrow(g),
               n_used = sum(g$selected), mean_Pr_all = mean(g$Pr), mean_Pr_used = mean(g$Pr[g$selected]))))
  mesh$label <- factor(mesh$label, levels = all_labels)
  mesh <- mesh[order(mesh$label, mesh$rep_n), ]
  
  # ---- random pairing of the selected pores between replicate meshes -> mean ----
  set.seed(pr_seed)
  pts_l <- list()
  for (lab in unique(as.character(pr_raw$label))) {
    g    <- pr_raw[as.character(pr_raw$label) == lab & pr_raw$selected, ]
    vals <- lapply(split(g$Pr, g$rep_n), function(v) v[sample.int(length(v))])
    J    <- max(lengths(vals))
    mat  <- do.call(cbind, lapply(vals, function(v) c(v, rep(NA_real_, J - length(v)))))
    pts_l[[length(pts_l) + 1]] <- data.frame(
      medium = g$medium[1], ink = g$ink[1], label = lab, point = seq_len(J),
      Pr_mean = rowMeans(mat, na.rm = TRUE), n_replicates = rowSums(!is.na(mat)))
  }
  pr_points <- fix(do.call(rbind, pts_l))
  
  # ---- status / flags per ink ----
  pr_status <- fix(do.call(rbind, fold_l))
  pr_status$min_pores <- sapply(as.character(pr_status$label), function(l)
    if (any(mesh$label == l)) min(mesh$n_pores[mesh$label == l]) else NA)
  pr_status$n_points  <- sapply(as.character(pr_status$label), function(l) sum(as.character(pr_points$label) == l))
  pr_status$flag_meshes <- pr_status$n_meshes > 0 & pr_status$n_meshes < pr_n_meshes
  pr_status$flag_pores  <- pr_status$n_meshes > 0 & !is.na(pr_status$min_pores) & pr_status$min_pores < pr_n_pores
  pr_status$suffix <- paste0(ifelse(pr_status$flag_meshes, paste0(" ", pr_flag_mesh), ""),
                             ifelse(pr_status$flag_pores,  paste0(" ", pr_flag_pore), ""))
  pr_status$printable <- ifelse(pr_status$n_meshes == 0, "not printable (no mesh)",
                                ifelse(pr_status$flag_meshes | pr_status$flag_pores, "partially printable", "printable"))
  pr_status <- pr_status[order(pr_status$label), ]
  
  m_idx <- match(as.character(pr_points$label), as.character(pr_status$label))
  pr_points$n_meshes <- pr_status$n_meshes[m_idx]
  pr_points$suffix   <- pr_status$suffix[m_idx]
  pr_points <- pr_points[order(pr_points$label, pr_points$point), ]
  
  write.csv(pr_raw,    out_path("printability_raw_pores.csv"), row.names = FALSE)
  write.csv(mesh,      out_path("printability_meshes.csv"), row.names = FALSE)
  write.csv(pr_points, out_path("printability_points.csv"), row.names = FALSE)
  write.csv(pr_status, out_path("printability_status.csv"), row.names = FALSE)
  message("Printability status per ink:")
  print(pr_status[, c("label", "n_meshes", "min_pores", "n_points", "printable")], row.names = FALSE)
  
  # ---- plots (inks without any mesh are shown as empty "not printable" slots) ----
  for (m in levels(pr_status$medium)) {
    df  <- droplevels(subset(pr_points, medium == m))
    nop <- subset(pr_status, medium == m & n_meshes == 0)
    if (nrow(df) == 0 && nrow(nop) == 0) next
    lev <- all_labels[all_labels %in% c(as.character(df$label), as.character(nop$label))]
    w   <- max(fig_width * 0.6, 3 + 1.8 * length(lev))
    p   <- make_pr_plot(df, lev, nop, cap_w = w) + ggtitle(tm("Printability", m))
    save_plot(p, paste0("printability_", safe_name(m)), w + key_w, fig_height + 3)
    print(p)
  }
  nop_all <- subset(pr_status, n_meshes == 0)
  lev_all <- all_labels[all_labels %in% c(as.character(pr_points$label), as.character(nop_all$label))]
  p_all <- make_pr_plot(pr_points, lev_all, nop_all, strip_medium = TRUE, cap_w = 3 * fig_width) +
    facet_grid(. ~ medium, scales = "free_x", space = "free_x") + ggtitle("Printability")
  save_plot(p_all, "printability_all_media", 3 * fig_width + key_w, fig_height + 3)
  print(p_all)
  
  list2env(list(pr_raw = pr_raw, pr_meshes = mesh, pr_points = pr_points, pr_status = pr_status),
           envir = globalenv())
  invisible()
}
run_printability()

# =============================================================================
# 6c. F) RHEOLOGY  (viscosity vs shear rate, thixotropy / recovery, gelation)
# =============================================================================
# Folder layout (found recursively, folder names do not matter, only file names):
#   <rheo_dir>/Viscosity vs shear/<ink folder>/TS_Viscosity_<ink>_<medium>_<rep>.csv
#   <rheo_dir>/Thixotropy/<ink folder>/TS_Thixotropy_<ink>_<medium>_<rep>.csv
#   <rheo_dir>/Gelation/<ink folder>/TS_Gelation_<ink>_<medium>_<rep>.csv     (phase 1: 1 min isothermal at 40 °C)
#   <rheo_dir>/Gelation/<ink folder>/TS_Gelation_<ink>_<medium>_<rep>_1.csv   (phase 2: 40 -> 15 °C in 5 min)
#     Phase 2 time restarts at 0; it is shifted by the last time of phase 1 so the curve is continuous (1 + 5 = 6 min).
# ('_' inside the ink name stands for '.', e.g. 30k-c13_3MC0Gel_PM3_1 -> 30k-c13.3MC0Gel (PM3), replicate 1)
# Files are the Anton Paar RheoCompass text export (UTF-16, tab separated, one block per interval).
# All rheology curves are MEAN ± SD of the replicates:
#   viscosity / thixotropy: replicates averaged point by point (same measuring points in every replicate);
#   gelation: every replicate is interpolated onto a common time grid (0 to the end of the shortest
#             replicate, step 'gel_grid_step') and then averaged, so replicates with a different number
#             of points or slightly different time stamps are averaged correctly.
# No t-test/ANOVA is computed: with n = 2 a test has no meaningful power; curves are compared descriptively.
rheo_dir       <- '/Users/tomas_serra/Library/CloudStorage/GoogleDrive-tomasfserra1@gmail.com/My Drive/Tese Drive/Current Thesis Files/Ink Characterizations/Rheology'
rheo_in_groups <- c("viscosity", "thixotropy", "gelation")   # techniques added to the comparison groups
# --- gelation ---
gel_expected_min <- c(1, 5)  # duration (min) of phase 1 (isothermal, 6 x 10 s) and phase 2 (cooling 40 -> 15 °C, 30 x 10 s) -> total 6 min
gel_grid_step  <- 0.02     # min; common time grid for averaging the replicates
gel_time_unit  <- "auto"   # only used if a file has no unit row: "s", "min" or "auto"
# (the unit is normally READ from the unit row of each file, e.g. [s] or [min])
#  Replicates left out of the mean curves WITHOUT deleting the files. Write the phase-1 file name
#  without ".csv"; both phase files of that replicate are excluded. Example:
#   gel_exclude <- c("TS_Gelation_30k-c13_3MC0Gel_PM3_3", "TS_Gelation_90k-c40MC0Gel_PM3_1")
gel_exclude    <- c()
thix_n_end     <- 3        # recovery = mean viscosity of the last N points of interval 3 / last N points of interval 1
rheo_line_size <- 0.5
rheo_point_size <- 1.6

read_rheo_file <- function(f) {
  # Reads an Anton Paar export in any of the formats seen so far: UTF-16 (with BOM) or UTF-8, separated by
  # tabs, semicolons or commas, also files re-saved by Excel/Numbers (leading tab missing, columns glued together).
  x <- readBin(f, "raw", file.info(f)$size)
  txt <- NA_character_
  if (length(x) >= 2 && ((x[1] == as.raw(0xff) && x[2] == as.raw(0xfe)) || (x[1] == as.raw(0xfe) && x[2] == as.raw(0xff))))
    txt <- iconv(list(x), from = "UTF-16", to = "UTF-8")
  if (is.na(txt)) txt <- tryCatch(iconv(rawToChar(x), from = "UTF-8", to = "UTF-8"), error = function(e) NA_character_)
  if (is.na(txt)) txt <- tryCatch(iconv(rawToChar(x), from = "latin1", to = "UTF-8"), error = function(e) NA_character_)
  if (is.na(txt)) return(NULL)
  txt <- sub("^\ufeff", "", txt)
  ln <- strsplit(txt, "\r?\n")[[1]]
  is_data <- function(l) grepl("^[\t;, ]*[0-9]+[\t;, ]+[-0-9.,]+", l)
  iv <- 0; hdr <- NULL; tunit <- NA_character_; rows <- list()
  for (l in ln) {
    if (grepl("^Interval and data points", l)) { iv <- iv + 1; next }
    if (grepl("^Interval data:", l)) {
      hdr <- strsplit(sub("^Interval data:[\t;, ]*", "", l), "[\t;,]+")[[1]]; hdr <- hdr[nzchar(hdr)]; next
    }
    if (is.null(hdr)) next
    if (is.na(tunit) && grepl("\\[", l) && !is_data(l)) {          # unit row, e.g. [s] [s] [°C] [Pa]
      u <- tolower(regmatches(l, gregexpr("\\[[^]]*\\]", l))[[1]]); u <- u[u %in% c("[s]", "[min]", "[h]")]
      if (length(u)) { tunit <- gsub("[][]", "", u[1]); next }
    }
    if (is_data(l)) {
      sep <- if (grepl("\t", l)) "\t" else if (lengths(regmatches(l, gregexpr(";", l))) >= lengths(regmatches(l, gregexpr(",", l)))) ";" else ","
      v <- strsplit(l, sep, fixed = TRUE)[[1]]
      if (length(v) && !nzchar(v[1])) v <- v[-1]
      strict <- length(v) >= length(hdr) - 1 && length(v) <= length(hdr) &&
        grepl("^[0-9]+$", v[1]) && !any(grepl("[0-9] +[-0-9]", v))
      if (!strict) v <- strsplit(trimws(gsub(if (sep == ",") "[\t;, ]+" else "[\t; ]+", " ", l)), " ", fixed = TRUE)[[1]]   # glued columns -> split on blanks
      length(v) <- length(hdr)
      rows[[length(rows) + 1]] <- c(as.character(iv), v)
    }
  }
  if (!length(rows)) return(NULL)
  d <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(d) <- c("interval", hdr)
  for (k in setdiff(names(d), "Status")) d[[k]] <- suppressWarnings(as.numeric(gsub(",", ".", d[[k]])))
  attr(d, "time_unit") <- tunit
  d
}

rh_col <- function(d, pat) {
  i <- grep(pat, names(d), ignore.case = TRUE)[1]
  if (is.na(i)) { message("Column '", pat, "' not found in rheology file -> NA"); return(rep(NA_real_, nrow(d))) }
  d[[i]]
}
# unit of the Time column ("s", "min" or "h"), read from the file's unit row
rh_time_unit <- function(d) {
  tu <- attr(d, "time_unit"); tm <- rh_col(d, "^Time")
  if (!is.null(tu) && !is.na(tu) && tu %in% c("s", "min", "h")) return(tu)
  if (gel_time_unit != "auto") return(gel_time_unit)
  if (max(tm, na.rm = TRUE) > 5) "s" else "min"
}

collect_rheo <- function(test_pat) {
  files <- list.files(rheo_dir, pattern = paste0("^TS_", test_pat, "_.*\\.csv$"), recursive = TRUE,
                      full.names = TRUE, ignore.case = TRUE)
  out <- list()
  for (f in files) {
    b <- sub("\\.csv$", "", basename(f), ignore.case = TRUE)
    b <- sub(paste0("^TS_", test_pat, "_"), "", b, ignore.case = TRUE)
    rp <- as.integer(sub(".*_([0-9]+)$", "\\1", b)); b <- sub("_[0-9]+$", "", b)
    medium <- sub(".*_", "", b); ink <- canon_ink(gsub("_", ".", sub("_[^_]*$", "", b)))
    lab <- paste0(ink, " (", medium, ")")
    if (!lab %in% all_labels) { message("Rheology file not matched to an ink (skipped): ", basename(f)); next }
    d <- read_rheo_file(f)
    if (is.null(d)) { message("No data rows found in: ", basename(f)); next }
    d$label <- lab; d$medium <- medium; d$rep <- rp
    out[[length(out) + 1]] <- d
  }
  if (!length(out)) return(NULL)
  cols <- Reduce(intersect, lapply(out, names))
  do.call(rbind, lapply(out, function(d) d[, cols]))
}

# gelation: two files per replicate (<rep>.csv = phase 1, <rep>_1.csv = phase 2), joined into one curve
collect_gelation <- function() {
  files <- list.files(rheo_dir, pattern = "^TS_Gelation_.*\\.csv$", recursive = TRUE,
                      full.names = TRUE, ignore.case = TRUE)
  segs <- list()
  for (f in files) {
    b   <- sub("^TS_Gelation_", "", sub("\\.csv$", "", basename(f), ignore.case = TRUE), ignore.case = TRUE)
    tok <- strsplit(b, "_", fixed = TRUE)[[1]]
    n_end <- 0                                             # number of trailing numeric tokens (1: phase 1, 2: phase 2)
    while (n_end < length(tok) && grepl("^[0-9]+$", tok[length(tok) - n_end])) n_end <- n_end + 1
    k <- length(tok) - n_end                               # position of the medium
    if (!(n_end %in% 1:2) || k < 2) { message("Gelation file name not understood (skipped): ", basename(f)); next }
    medium <- tok[k]
    ink    <- canon_ink(paste(tok[1:(k - 1)], collapse = "."))
    lab    <- paste0(ink, " (", medium, ")")
    if (!lab %in% all_labels) { message("Gelation file not matched to an ink (skipped): ", basename(f)); next }
    d <- read_rheo_file(f)
    if (is.null(d)) { message("No data rows found in: ", basename(f)); next }
    segs[[length(segs) + 1]] <- data.frame(
      label = lab, medium = medium, rep = as.integer(tok[k + 1]),
      rep_id = paste0("TS_Gelation_", paste(tok[1:(k + 1)], collapse = "_")),
      phase = if (n_end == 2) as.integer(tok[k + 2]) + 1L else 1L, file = basename(f),
      t_unit = rh_time_unit(d),
      row = seq_len(nrow(d)), time = rh_col(d, "^Time"), temp = rh_col(d, "Temperature"),
      Gp = rh_col(d, "Storage"), Gpp = rh_col(d, "Loss Modulus"), stringsAsFactors = FALSE)
  }
  if (!length(segs)) return(NULL)
  s <- do.call(rbind, segs)
  out <- list()
  for (r in split(s, s$rep_id)) {
    r <- r[order(r$phase, r$row), ]
    ph <- sort(unique(r$phase))
    if (!1L %in% ph) message("Gelation: phase-1 file missing for ", r$rep_id[1], " (curve starts at phase 2)")
    if (!2L %in% ph) message("Gelation: phase-2 file (_1) missing for ", r$rep_id[1], " (curve stops after phase 1)")
    tm <- ifelse(r$t_unit == "s", r$time / 60, ifelse(r$t_unit == "h", r$time * 60, r$time))   # -> minutes
    r$t_min <- NA_real_; off <- 0
    for (q in ph) {                                        # each phase continues where the previous one ended
      w <- r$phase == q
      r$t_min[w] <- tm[w] + off
      off <- max(r$t_min[w], na.rm = TRUE)
    }
    r$has_phase1 <- 1L %in% ph
    r$excluded <- r$rep_id[1] %in% gel_exclude
    out[[length(out) + 1]] <- r
  }
  s <- do.call(rbind, out)
  bad <- setdiff(gel_exclude, s$rep_id)
  if (length(bad)) message("gel_exclude entries that match no gelation replicate: ", paste(bad, collapse = ", "))
  s$label_f <- factor(s$label, levels = all_labels)
  s <- s[order(s$label_f, s$rep, s$phase, s$row), ]; s$label_f <- NULL
  s
}
# end of phase 1 (min), median over replicates -> dashed line in the plots
gel_p1_end <- function(d) {
  x <- tapply(d$t_min[d$phase == 1], d$rep_id[d$phase == 1], max, na.rm = TRUE)
  if (length(x)) median(x) else NA_real_
}
# gelation mean ± SD: every replicate interpolated onto a common time grid, then averaged
gel_average <- function(d) {
  out <- list()
  for (lab in unique(as.character(d$label))) {
    reps <- split(d[d$label == lab, ], d$rep_id[d$label == lab])
    tmax <- min(sapply(reps, function(r) max(r$t_min, na.rm = TRUE)))
    grid <- seq(0, tmax, by = gel_grid_step)
    itp <- function(r, v) {
      k <- is.finite(r$t_min) & is.finite(r[[v]])
      if (sum(k) < 2) return(rep(NA_real_, length(grid)))
      approx(r$t_min[k], r[[v]][k], xout = grid, ties = mean, rule = 2)$y
    }
    mk <- function(v) {
      m <- do.call(cbind, lapply(reps, itp, v = v))
      list(mean = rowMeans(m, na.rm = TRUE),
           sd = if (ncol(m) > 1) apply(m, 1, sd, na.rm = TRUE) else rep(NA_real_, nrow(m)))
    }
    gp <- mk("Gp"); gpp <- mk("Gpp"); tp <- mk("temp")
    out[[length(out) + 1]] <- data.frame(
      label = lab, point = seq_along(grid), time_mean = grid, temp_mean = tp$mean,
      Gp_mean = gp$mean, Gp_sd = gp$sd, Gpp_mean = gpp$mean, Gpp_sd = gpp$sd,
      n = length(reps), medium = sub(".*\\((.*)\\)$", "\\1", lab), stringsAsFactors = FALSE)
  }
  o <- do.call(rbind, out); o$label <- factor(o$label, levels = all_labels)
  o[order(o$label, o$point), ]
}

agg_ms <- function(df, by, vars) {
  m <- aggregate(df[vars], df[by], mean, na.rm = TRUE)
  s <- aggregate(df[vars], df[by], function(x) if (sum(!is.na(x)) > 1) sd(x, na.rm = TRUE) else NA_real_)
  n <- aggregate(df[vars[1]], df[by], function(x) sum(!is.na(x)))
  names(m)[-seq_along(by)] <- paste0(vars, "_mean"); names(s)[-seq_along(by)] <- paste0(vars, "_sd"); names(n)[-seq_along(by)] <- "n"
  out <- Reduce(function(a, b) merge(a, b, by = by), list(m, s, n))
  out$medium <- sub(".*\\((.*)\\)$", "\\1", out$label)
  out$label <- factor(out$label, levels = all_labels)
  out[order(out$label, out[[by[length(by)]]]), ]
}
rh_ncap <- function(df) {
  nn <- tapply(df$n, as.character(df$label), max); r <- range(nn)
  paste0("Data shown as mean ± SD (n = ", if (r[1] == r[2]) r[1] else paste0(r[1], "–", r[2]), ")")
}
rh_scale <- function(strip) scale_colour_manual(values = ink_colors, name = "Ink", labels = strip_lab(strip))
rh_theme <- function() list(base_theme() + theme(legend.position = "right"), key_guides())

plot_visc <- function(df, strip_medium = TRUE, cap_w = 12) {
  df$lo <- pmax(df$visc_mean - df$visc_sd, df$visc_mean / 100); df$hi <- df$visc_mean + df$visc_sd
  tt <- unique(data.frame(label = as.character(df$label), T = round(df$temp_mean), stringsAsFactors = FALSE))
  tt <- tapply(strip_lab(strip_medium)(tt$label), tt$T, paste, collapse = ", ")
  tnote <- paste0("Temperature: ", paste0(names(tt), " °C (", unname(tt), ")", collapse = "; "), ".")
  ggplot(df, aes(rate_mean, visc_mean, colour = label)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, alpha = 0.3, linewidth = 0.3, na.rm = TRUE) +
    geom_line(linewidth = rheo_line_size) + geom_point(size = rheo_point_size) +
    scale_x_log10() + scale_y_log10() + annotation_logticks(sides = "bl", size = 0.25) + rh_scale(strip_medium) +
    labs(x = expression("Shear rate (s"^-1*")"), y = "Viscosity (Pa·s)",
         caption = cap_out(c(paste0(rh_ncap(df), ". Shear-rate ramp 0.1–1000 1/s."), tnote), cap_w)) +
    rh_theme()
}
plot_thix <- function(df, strip_medium = TRUE, cap_w = 12) {
  # as before (line + points + SD bars); the only change: points with no value are dropped and every ink is drawn in time order as ONE
  # continuous line, so the curve has no gap (the data themselves are not changed)
  df <- df[is.finite(df$time_mean) & is.finite(df$visc_mean) & df$visc_mean > 0, ]
  df <- df[order(df$label, df$time_mean), ]
  r1 <- max(df$time_mean[df$interval == 1]); r3 <- min(df$time_mean[df$interval == 3])
  p <- ggplot(df, aes(time_mean, visc_mean, colour = label, group = label))
  if (is.finite(r1) && is.finite(r3))
    p <- p + annotate("rect", xmin = r1, xmax = r3, ymin = 0, ymax = Inf, fill = "grey80", alpha = 0.6) +
    geom_text(data = data.frame(x = c(min(df$time_mean) + (r1 - min(df$time_mean)) / 2, (r1 + r3) / 2, r3 + (max(df$time_mean) - r3) / 2),
                                lab = c("pre-printing\n(0.1 s\u207B\u00B9)", "nozzle\n(700 s\u207B\u00B9)", "post-printing\n(0.1 s\u207B\u00B9)")),
              aes(x = x, y = Inf, label = lab), inherit.aes = FALSE, vjust = 1.2, size = 3, colour = "grey30")
  p + geom_errorbar(aes(ymin = pmax(visc_mean - visc_sd, visc_mean / 100), ymax = visc_mean + visc_sd),
                    width = 0, alpha = 0.3, linewidth = 0.3, na.rm = TRUE) +
    geom_line(linewidth = rheo_line_size) + geom_point(size = rheo_point_size) +
    scale_y_log10(expand = expansion(mult = c(0.05, 0.25))) + annotation_logticks(sides = "l", size = 0.25) + rh_scale(strip_medium) +
    labs(x = "Time (s)", y = "Viscosity (Pa·s)",
         caption = cap_out(paste0(rh_ncap(df), ". Shear rate: 0.1 s\u207B\u00B9 (60 s), 700 s\u207B\u00B9 (10 s, nozzle, grey), 0.1 s\u207B\u00B9 (60 s)."), cap_w)) + rh_theme()
}
gel_phase_marks <- function(t_range) {                    # t_range in SECONDS
  p1 <- if (exists("rh_gel_p1")) rh_gel_p1 * 60 else NA_real_
  if (!is.finite(p1)) return(list())
  list(geom_vline(xintercept = p1, linetype = "dashed", colour = "grey50", linewidth = 0.4),
       geom_text(data = data.frame(x = c((t_range[1] + p1) / 2, p1 + (t_range[2] - p1) / 2),
                                   txt = c("40 °C", "40 → 15 °C")),
                 aes(x = x, y = Inf, label = txt), inherit.aes = FALSE, vjust = 1.3, size = 3, colour = "grey30"))
}
gel_caption <- "Gelation: 1 min at 40 °C, then cooling to 15 °C. Solid: G′; dotted: G″; grey: temperature."
plot_gel <- function(df, strip_medium = TRUE, cap_w = 12) {
  g1 <- data.frame(label = df$label, medium = df$medium, t = df$time_mean, m = df$Gp_mean,  type = "G′ (storage)")
  g2 <- data.frame(label = df$label, medium = df$medium, t = df$time_mean, m = df$Gpp_mean, type = "G″ (loss)")
  d <- rbind(g1, g2); d <- d[is.finite(d$m) & d$m > 0, ]
  d$t <- d$t * 60                                          # minutes -> seconds
  d$grp <- paste(d$label, d$type)
  # LOGARITHMIC modulus axis, so that the diluted inks (G' and G'' of about 0.1 Pa) and the stiff inks (hundreds to thousands of Pa) are all visible;
  # temperature (10-40 °C) on the right axis, drawn on the same log scale: 10 °C = bottom, 40 °C = top
  ylo <- 10^floor(log10(min(d$m))); yhi <- max(d$m) * 1.35; t_lo <- 10; t_hi <- 40
  lr     <- log10(yhi) - log10(ylo)
  to_y   <- function(T) 10^(log10(ylo) + (T - t_lo) / (t_hi - t_lo) * lr)
  from_y <- function(y) t_lo + (log10(y) - log10(ylo)) / lr * (t_hi - t_lo)
  dt <- df[is.finite(df$temp_mean), c("label", "medium", "time_mean", "temp_mean")]
  dt$t <- dt$time_mean * 60; dt$y <- to_y(pmin(pmax(dt$temp_mean, t_lo), t_hi))
  ggplot(d, aes(t, m, colour = label, linetype = type, group = grp)) +
    geom_line(data = dt, aes(t, y, group = label), inherit.aes = FALSE, colour = "grey60", linewidth = 0.6) +
    gel_phase_marks(range(d$t, na.rm = TRUE)) +
    geom_line(linewidth = 0.8) +
    scale_x_continuous(breaks = function(l) seq(0, ceiling(l[2] / 60) * 60, 60), expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_log10(limits = c(ylo, yhi), expand = expansion(0),
                  sec.axis = sec_axis(~ from_y(.), name = "Temperature (°C)", breaks = seq(t_lo, t_hi, 10))) +
    annotation_logticks(sides = "l", size = 0.25) +
    rh_scale(strip_medium) +
    scale_linetype_manual(values = c("G′ (storage)" = "solid", "G″ (loss)" = "22"), name = NULL) +
    guides(linetype = guide_legend(override.aes = list(linewidth = 0.8))) +
    labs(x = "Time (s)", y = "G′, G″ (Pa)",
         caption = cap_out(paste0("Mean of the replicates (n = ", max(df$n), "). ", gel_caption), cap_w)) + rh_theme()
}
plot_recovery <- function(df, strip_medium = TRUE, cap_w = 12, reps = NULL, show_values = TRUE) {
  vt <- df$rec_mean + ifelse(is.finite(df$rec_sd), df$rec_sd, 0)
  if (!is.null(reps) && nrow(reps) > 0) {
    rmax <- tapply(reps$rec, as.character(reps$label), max, na.rm = TRUE)
    vt <- pmax(vt, unname(rmax[as.character(df$label)]), na.rm = TRUE)
  }
  df$.vy <- vt + 0.06 * max(vt, na.rm = TRUE)
  p <- ggplot(df, aes(label, rec_mean, fill = label)) +
    geom_hline(yintercept = 100, linetype = "dotted", colour = "grey30") +
    geom_col(colour = "black", linewidth = 0.3, width = 0.7) +
    geom_errorbar(aes(ymin = rec_mean - rec_sd, ymax = rec_mean + rec_sd), width = 0.2, na.rm = TRUE, alpha = sd_alpha)
  if (!is.null(reps)) p <- p + geom_point(data = reps, aes(label, rec), inherit.aes = FALSE, shape = 21, fill = "white", size = 1.8)
  if (isTRUE(show_values))
    p <- p + geom_text(data = df, aes(x = label, y = .vy, label = sprintf("%.0f", rec_mean)),
                       inherit.aes = FALSE, angle = 90, hjust = 0, vjust = 0.5, size = 2.3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.30)))
  p + scale_fill_manual(values = ink_colors, name = NULL, labels = disp,
                        guide = if (.key_mode == "none") "none" else guide_legend(keywidth = unit(0.45, "cm"), keyheight = unit(0.45, "cm"))) + scale_x_discrete() +
    labs(x = NULL, y = "Viscosity recovery (%)",
         caption = cap_out(paste0(rh_ncap(df), ". Recovery = viscosity after / before the nozzle interval × 100; dotted line: 100 %. No test (n = 2)."), cap_w)) +
    base_theme() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "right")
}

run_rheology <- function() {
  if (!dir.exists(rheo_dir)) { message("Rheology folder not found (is the drive mounted?): ", rheo_dir); return(invisible()) }
  rh_raw <<- list(); rh_avg <<- list(); rh_rec <<- NULL; rh_rec_rep <<- NULL
  # ---- read + average ----
  v <- collect_rheo("Viscosity")
  if (!is.null(v)) {
    d <- data.frame(label = v$label, rep = v$rep, point = rh_col(v, "Point"), rate = rh_col(v, "Shear Rate"),
                    visc = rh_col(v, "Viscosity") / 1000, temp = rh_col(v, "Temperature"))
    rh_raw$viscosity <<- d; rh_avg$viscosity <<- agg_ms(d, c("label", "point"), c("rate", "visc", "temp"))
  } else message("No viscosity files found.")
  t <- collect_rheo("Thixotropy")
  if (!is.null(t)) {
    d <- data.frame(label = t$label, rep = t$rep, interval = as.numeric(t$interval), point = rh_col(t, "Point"),
                    time = rh_col(t, "^Time"), visc = rh_col(t, "Viscosity") / 1000, stress = rh_col(t, "Shear Stress"))
    d$tp <- d$interval * 1000 + d$point
    rh_raw$thixotropy <<- d
    a <- agg_ms(d, c("label", "tp"), c("time", "visc", "stress")); a$interval <- floor(a$tp / 1000); rh_avg$thixotropy <<- a
    rec <- do.call(rbind, lapply(split(d, list(d$label, d$rep), drop = TRUE), function(x) {
      e <- function(k) { y <- x$visc[x$interval == k]; if (!length(y)) NA_real_ else mean(tail(y, thix_n_end)) }
      data.frame(label = x$label[1], rep = x$rep[1], visc_pre_Pas = e(1), visc_nozzle_Pas = e(2), visc_post_Pas = e(3),
                 rec = e(3) / e(1) * 100)
    }))
    rec$label  <- factor(rec$label, levels = all_labels)
    rec$medium <- factor(sub(".*\\((.*)\\)$", "\\1", as.character(rec$label)), levels = all_media)
    rh_rec_rep <<- rec
    rs <- agg_ms(transform(rec, label = as.character(label)), c("label"), c("visc_pre_Pas", "visc_nozzle_Pas", "visc_post_Pas", "rec"))
    rh_rec <<- rs
    write.csv(rec, out_path("rheology_thixotropy_recovery_per_replicate.csv"), row.names = FALSE)
    write.csv(rs,  out_path("rheology_thixotropy_recovery_summary.csv"), row.names = FALSE)
  } else message("No thixotropy files found.")
  gl <- collect_gelation()
  rh_gel_all <<- gl; rh_gel_p1 <<- NA_real_
  if (!is.null(gl)) {
    write.csv(gl, out_path("rheology_gelation_raw_data.csv"), row.names = FALSE)
    # ---- duration check: every replicate should last 1 + 5 = 6 min ----
    mx <- function(v) if (length(v) && any(is.finite(v))) max(v[is.finite(v)]) else NA_real_
    dur <- do.call(rbind, lapply(split(gl, gl$rep_id), function(x) data.frame(
      label = x$label[1], rep_id = x$rep_id[1], time_unit = paste(unique(paste0("phase", x$phase, ":", x$t_unit)), collapse = " "),
      phase1_end_min = mx(x$t_min[x$phase == 1]), total_min = mx(x$t_min),
      n_points_phase1 = sum(x$phase == 1), n_points_phase2 = sum(x$phase == 2), excluded = x$excluded[1])))
    write.csv(dur, out_path("rheology_gelation_durations.csv"), row.names = FALSE)
    message("Gelation duration per replicate (expected phase 1 = ", gel_expected_min[1], " min, total = ", sum(gel_expected_min), " min):")
    print(dur[, c("rep_id", "time_unit", "phase1_end_min", "total_min", "n_points_phase1", "n_points_phase2")], digits = 3, row.names = FALSE)
    off <- dur[is.na(dur$total_min) | abs(dur$total_min - sum(gel_expected_min)) > 0.25 |
                 is.na(dur$phase1_end_min) | abs(dur$phase1_end_min - gel_expected_min[1]) > 0.1, ]
    if (nrow(off)) message("!! These gelation replicates do NOT last ", sum(gel_expected_min), " min (", paste(gel_expected_min, collapse = " + "), ") -> check the files / time unit: ",
                           paste(off$rep_id, collapse = ", "))
    d <- gl[!gl$excluded & gl$has_phase1, ]
    if (nrow(d)) {
      rh_gel_p1 <<- gel_p1_end(d)
      rh_raw$gelation <<- d; rh_avg$gelation <<- gel_average(d)
      fin <- do.call(rbind, lapply(split(d, d$rep_id), function(x) {
        k <- tail(x[order(x$t_min), ], 5)
        data.frame(label = x$label[1], rep_id = x$rep_id[1], Gp_end_Pa = mean(k$Gp), Gpp_end_Pa = mean(k$Gpp),
                   tan_delta_end = mean(k$Gpp / k$Gp))
      }))
      write.csv(fin, out_path("rheology_gelation_end_values_per_replicate.csv"), row.names = FALSE)
      write.csv(agg_ms(fin, "label", c("Gp_end_Pa", "Gpp_end_Pa", "tan_delta_end")),
                out_path("rheology_gelation_end_values_summary.csv"), row.names = FALSE)
    } else message("Gelation: no replicate left after exclusions.")
  } else message("No gelation files found.")
  for (k in names(rh_avg)) write.csv(rh_avg[[k]], out_path(paste0("rheology_", k, "_mean_curves.csv")), row.names = FALSE)
  
  # ---- plots: one per medium + all media ----
  fns <- list(viscosity = plot_visc, thixotropy = plot_thix, gelation = plot_gel)
  for (k in names(rh_avg)) {
    df <- rh_avg[[k]]; w <- fig_width + 6
    for (m in all_media) {
      dm <- droplevels(df[df$medium == m, ]); if (!nrow(dm)) next
      p <- fns[[k]](dm, TRUE, w) + ggtitle(tm(tech_lab[[k]], m))
      save_plot(p, paste0("rheology_", k, "_", safe_name(m)), w, fig_height + 2); print(p)
    }
    w3 <- 3 * fig_width
    p <- fns[[k]](df, TRUE, w3) + facet_wrap(~ medium, nrow = 1) + ggtitle(tech_lab[[k]])
    save_plot(p, paste0("rheology_", k, "_all_media"), w3, fig_height + 3); print(p)
  }
  if (!is.null(rh_rec)) {
    rr <- rh_rec; rr$medium <- sub(".*\\((.*)\\)$", "\\1", as.character(rr$label)); rr$label <- factor(rr$label, levels = all_labels)
    rr$label <- droplevels(rr$label)
    for (m in all_media) {
      dm <- droplevels(rr[rr$medium == m, ]); if (!nrow(dm)) next
      p <- plot_recovery(dm, TRUE, fig_width + 4, rh_rec_rep[rh_rec_rep$label %in% dm$label, ]) + ggtitle(tm("Viscosity recovery", m))
      save_plot(p, paste0("rheology_recovery_", safe_name(m)), fig_width + 4 + key_w, fig_height + 1.5); print(p)
    }
    p <- plot_recovery(rr, TRUE, 3 * fig_width, rh_rec_rep) + facet_grid(. ~ medium, scales = "free_x", space = "free_x") + ggtitle("Viscosity recovery")
    save_plot(p, "rheology_recovery_all_media", 3 * fig_width + key_w, fig_height + 1.5); print(p)
  }
  invisible()
}
run_rheology()

# comparison-group plots for the rheology techniques (called from run_comparisons)
do_rheo_group <- function(g, tech) {
  if (tech %in% cmp_skip && tech != "thixotropy") return(invisible())
  df <- if (exists("rh_avg")) rh_avg[[tech]] else NULL
  if (is.null(df)) { message("[", g$id, " / ", tech, "] no rheology data"); return(invisible()) }
  labs <- g$labels[g$labels %in% as.character(df$label)]
  miss <- setdiff(g$labels, labs)
  if (length(miss)) message("[", g$id, " / ", tech, "] no data for: ", paste(miss, collapse = ", "))
  if (length(labs) < 2) { message("[", g$id, " / ", tech, "] fewer than 2 inks with data -> skipped"); return(invisible()) }
  d <- df[as.character(df$label) %in% labs, ]; d$label <- factor(as.character(d$label), levels = labs)
  w <- fig_width + 1
  ttl <- function(w) wrap_parts(g$title, w, chars_per_cm = 4)
  fn <- switch(tech, viscosity = plot_visc, thixotropy = plot_thix, gelation = plot_gel)
  if (!(tech == "thixotropy" && "thixotropy" %in% cmp_skip)) {
    p <- fn(d, FALSE, w) + ggtitle(tt(tech_lab[[tech]], g$title, w))
    save_plot(p, paste0("cmp_", g$id, "_rheology_", tech), w, fig_height + 1.5); print(p)
  }
  if (tech == "thixotropy" && !is.null(rh_rec)) {
    r <- rh_rec[as.character(rh_rec$label) %in% labs, ]; r$label <- factor(as.character(r$label), levels = labs)
    rp <- rh_rec_rep[as.character(rh_rec_rep$label) %in% labs, ]; rp$label <- factor(as.character(rp$label), levels = labs)
    w2 <- 3.4 + 1.45 * length(labs)
    p <- plot_recovery(r, FALSE, w2, rp) + ggtitle(tt("Viscosity recovery", g$title, w2))
    save_plot(p, paste0("cmp_", g$id, "_rheology_recovery"), w2, fig_height + 1.5); print(p)
  }
  invisible()
}

# =============================================================================
# 7. E) COMPARISON GROUPS (Results & Discussion)
# =============================================================================
#  For every group in 'comparison_groups' and every test listed there:
#   compression -> Young's modulus (mean ± SD, t-test vs reference) + mean stress-strain curves
#   water       -> water content (mean ± SD, t-test vs reference)
#   swelling    -> swelling degree vs time (mean ± SD, no test)
#   printability-> Pr dot plot (adequate range 0.9-1.1, dotted line Pr = 1); inks with no mesh shown as 'not printable'
#   rheology    -> viscosity, thixotropy (+ recovery), gelation (mean ± SD, no test)
#  The reference (first label, or ref = ...) is tested against every other ink of the group.
#  Plots are saved in Comparisons/<group folder>/.

run_comparisons <- function() {
  .key_mode <<- "none"; on.exit(.key_mode <<- "side")      # comparison plots carry NO key and NO ink names; the keys are saved as images
  ptabs <- list()
  save_key(all_labels, "cmp_key_all_inks", folder = file.path("Comparisons", "Keys"))      # the ONE key of all inks
  take <- function(df, labs) {                     # subset in group order
    df <- df[as.character(df$label) %in% labs, ]
    df$label <- factor(as.character(df$label), levels = labs)
    df
  }
  for (g in comparison_groups) {
    ref <- if (!is.null(g$ref)) g$ref else g$labels[1]
    for (tech in c(g$tests, if (exists("rheo_in_groups")) rheo_in_groups)) {
      if (tech %in% c("viscosity", "thixotropy", "gelation")) { do_rheo_group(g, tech); next }
      if (tech == "printability") {
        if (!exists("pr_points")) { message("[", g$id, " / printability] no printability data"); next }
        labs <- g$labels[g$labels %in% as.character(pr_points$label)]
        nop  <- g$labels[g$labels %in% as.character(pr_status$label[pr_status$n_meshes == 0])]
        nod  <- setdiff(g$labels, c(labs, nop))
        if (length(nod)) message("[", g$id, " / printability] no folder/data for: ", paste(nod, collapse = ", "))
        if (length(labs) + length(nop) < 1) { message("[", g$id, " / printability] nothing to plot"); next }
        lev  <- g$labels[g$labels %in% c(labs, nop)]            # group order, not printable inks included
        npd  <- pr_status[as.character(pr_status$label) %in% nop, ]
        w_p  <- 3.4 + 1.45 * length(lev)
        p <- make_pr_plot(pr_points[as.character(pr_points$label) %in% labs, ], lev, npd,
                          strip_medium = FALSE, cap_w = w_p) +
          ggtitle(tt("Printability", g$title, w_p))
        save_plot(p, paste0("cmp_", g$id, "_printability"), w_p, fig_height + 1.5)
        print(p)
        next
      }
      avail <- switch(tech,
                      compression = if (exists("E_sum")) as.character(E_sum$label) else character(0),
                      water       = if (exists("wc")) as.character(wc$label) else character(0),
                      swelling    = as.character(res$label), character(0))
      labs <- g$labels[g$labels %in% avail]
      miss <- setdiff(g$labels, labs)
      if (length(miss)) message("[", g$id, " / ", tech, "] no data for: ", paste(miss, collapse = ", "))
      if (length(labs) < 2) { message("[", g$id, " / ", tech, "] fewer than 2 inks with data -> skipped"); next }
      others <- setdiff(labs, ref)
      pairs  <- if (ref %in% labs) lapply(others, function(o) c(ref, o)) else list()
      n_l    <- length(labs)
      w_bar  <- 3.4 + 1.45 * n_l
      ttl    <- function(w) wrap_parts(g$title, w, chars_per_cm = 4)
      ref_note <- paste0("Each ink compared with the reference: ", ref, ".")
      
      bar_group <- function(df, mean_col, sd_col, rep_df, rep_col, test_fun, ylab, name, zoom = NULL, breaks = waiver(), room = 0, rel = TRUE, gap, step, tick) {
        df$panel <- "all"
        sm <- data.frame(panel = "all", label = df$label, top = df[[mean_col]] + df[[sd_col]])
        bb <- build_brackets(pairs, sm, test_fun, "panel", gap = gap, step = step, tick = tick, rel = rel)
        if (!is.null(bb$tab)) { bb$tab$group <- g$id; bb$tab$technique <- ylab; ptabs[[length(ptabs) + 1]] <<- bb$tab }
        p <- make_bar_plot(df, mean_col, sd_col, ylab = ylab,
                           caption = cap_out(c("Data shown as mean ± SD (n = 3).", ref_note, stat_legend(),
                                               ns_caption(bb$tab, strip = FALSE)), w_bar),
                           strip_medium = FALSE, rep_df = rep_df, rep_col = rep_col, br_df = bb$br,
                           ylim_zoom = zoom, ybreaks = breaks, room = room) + ggtitle(tt(sub(" \\(.*\\)$", "", ylab), g$title, w_bar))
        save_plot(p, paste0("cmp_", g$id, "_", name), w_bar, fig_height + 1.5)
        print(p)
      }
      
      if (tech == "compression") {
        dE <- take(E_sum, labs); rE <- take(comp_fits, labs)
        bar_group(dE, "E_mean", "E_sd", rE, "E_kPa", E_test, "Young's modulus (kPa)", "youngs_modulus",
                  gap = 0.05, step = 0.10, tick = 0.02)
        dc <- take(comp_avg, labs)
        if (!("stress_strain" %in% cmp_skip)) {
          w_c <- fig_width + 1
          p <- make_curve_plot(dc, strip_medium = FALSE, cap_w = w_c) + ggtitle(tt("Compression curves", g$title, w_c))
          save_plot(p, paste0("cmp_", g$id, "_stress_strain"), w_c, fig_height + 1)
          print(p)
        }
      }
      if (tech == "water") {
        dW <- take(wc, labs)
        rW <- if (exists("wc_rep") && !is.null(wc_rep)) take(wc_rep, labs) else NULL
        bar_group(dW, "water_content_pct", "water_content_sd", rW, "value", wc_test, "Water content (%)", "water_content",
                  zoom = wc_y_min, breaks = seq(wc_y_min, 100, 2), room = wc_label_room, rel = FALSE,
                  gap = wc_bracket_gap, step = wc_bracket_step, tick = wc_bracket_tick)
      }
      if (tech == "swelling") {
        dS <- take(res, labs)
        w_s <- fig_width + 1
        p <- make_swelling_plot(dS, strip_medium = FALSE, cap_w = w_s) + ggtitle(tt("Swelling degree", g$title, w_s))
        save_plot(p, paste0("cmp_", g$id, "_swelling"), w_s, fig_height + 1.5)
        print(p)
      }
    }
  }
  if (length(ptabs)) {
    allp <- do.call(rbind, lapply(ptabs, function(t) t[, c("group", "technique", "a", "b", "t", "df", "p_value", "significance")]))
    write.csv(allp, out_path("comparison_groups_ttests.csv"), row.names = FALSE)
    print(allp, digits = 3, row.names = FALSE)
  }
  invisible()
}
if (isTRUE(run_comparison_plots)) run_comparisons()


# =============================================================================
# 7b. pH OF THE INKS  (merged from pH_plots.R, so that one script makes every figure)
# =============================================================================
#  Input: ph_csv, blocks of rows -> a medium name, then one row per ink with its pH.
#  Output: the per-medium plot, the all-media plot and the object used as panel (b)
#  of the printability panel. One reading per ink, so there are no error bars.
ph_plot_obj <- NULL
run_ph <- function() {
  if (!nzchar(ph_csv) || !file.exists(ph_csv)) { message("pH csv not found: ", ph_csv); return(invisible()) }
  l <- readLines(ph_csv, warn = FALSE, encoding = "UTF-8")
  l <- iconv(l, "UTF-8", "UTF-8", sub = ""); l <- sub("^\ufeff", "", l)
  raw <- read.csv(text = paste(l, collapse = "\n"), header = FALSE, stringsAsFactors = FALSE,
                  fill = TRUE, col.names = paste0("V", 1:8))
  medium <- NA; rows <- list(); missing <- character()
  for (i in seq_len(nrow(raw))) {
    a <- trimws(raw$V1[i]); b <- trimws(raw$V2[i])
    if (!nzchar(a)) next
    if (a %in% all_media) { medium <- a; next }
    v <- suppressWarnings(as.numeric(b))
    if (is.na(v)) { missing <- c(missing, paste0(a, " (", medium, ")")); next }
    rows[[length(rows) + 1]] <- data.frame(medium = medium, ink = canon_ink(a), pH = v)
  }
  if (!length(rows)) { message("No pH readings found in ", ph_csv); return(invisible()) }
  d <- do.call(rbind, rows)
  d$label <- paste0(d$ink, " (", d$medium, ")")
  unknown <- setdiff(d$label, all_labels)
  if (length(unknown)) message("pH: ink not in 'inks' (no colour): ", paste(unknown, collapse = ", "))
  d <- d[d$label %in% all_labels, ]
  d$medium <- factor(d$medium, levels = all_media)
  d$label  <- factor(d$label,  levels = all_labels)
  d <- d[order(d$label), ]
  ph <<- d
  write.csv(d[, c("medium", "ink", "pH")], out_path("pH_table.csv"), row.names = FALSE)
  if (length(missing)) message("No pH reading for: ", paste(missing, collapse = ", "))
  
  mk <- function(df, show_values = TRUE) {
    p <- ggplot(df, aes(label, pH, fill = label)) +
      geom_hline(yintercept = 7, linetype = "dotted", colour = "grey40") +
      geom_col(colour = "black", linewidth = 0.4, width = 0.7)
    if (show_values) p <- p + geom_text(aes(label = sprintf("%.2f", pH)), vjust = -0.4, size = 2.5)
    p + scale_fill_manual(values = ink_colors, name = NULL, labels = disp,
                          guide = if (.key_mode == "none") "none" else
                            guide_legend(keywidth = unit(0.45, "cm"), keyheight = unit(0.45, "cm"))) +
      scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2), expand = expansion(mult = c(0, 0.02))) +
      labs(x = NULL, y = "pH",
           caption = cap_out("One reading per ink (n = 1); dotted line: neutral pH.", panel_width_cm)) +
      base_theme() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                           legend.position = "right")
  }
  for (m in levels(droplevels(d$medium))) {
    dm <- droplevels(subset(d, medium == m)); if (!nrow(dm)) next
    p <- mk(dm) + ggtitle(tm("pH", m))
    save_plot(p, paste0("pH_", safe_name(m)), max(fig_width * 0.6, 3 + 1.8 * nrow(dm)) + key_w, fig_height + 1.5)
  }
  p_all <- mk(d) + facet_grid(. ~ medium, scales = "free_x", space = "free_x") + ggtitle("pH of the inks")
  save_plot(p_all, "pH_all_media", 3 * fig_width + key_w, fig_height + 1.5)
  ph_plot_obj <<- mk(d, show_values = FALSE) +
    facet_grid(. ~ medium, scales = "free_x", space = "free_x")
  invisible()
}
run_ph()







# =============================================================================
# 8. FIGURE PANELS FOR CHAPTER 4  (one image per characterisation technique)
# =============================================================================
#  Each panel has two rows:
#    row 1 : (a) the curves or dots of every ink, wide, and (b) the quantity the
#            text argues from, narrower  (widths panel_widths, default 2 : 1)
#    row 2 : (c) the ink key, a horizontal strip of key_ncol columns
#  Saved at the width they are USED at (panel_width_cm), so LaTeX does not shrink them.
#  Two versions of every panel are written to <out_dir>/Panels/:
#    <Technique>.png           -> (a) grouped by printing medium
#    <Technique>_all_inks.png  -> (a) all inks on one set of axes
#  Panels (b) carry ONE x and ONE y axis: the media are separated by dotted lines.

if (!requireNamespace("patchwork", quietly = TRUE))
  stop("The panels need the 'patchwork' package: install.packages('patchwork')")

#  the dots of the curve plots are too heavy at panel size: shrink every point layer
shrink_points <- function(p, size = panel_point_size) {
  if (is.null(size)) return(p)
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomPoint")) {
      cur <- p$layers[[i]]$aes_params$size
      p$layers[[i]]$aes_params$size <- if (is.null(cur)) size else min(cur, size)
    }
  }
  p
}

panel_strip <- function(p, right_pad = 4, keep_legend = FALSE, blank_strips = FALSE,
                        small_points = TRUE) {
  if (isTRUE(small_points)) p <- shrink_points(p)
  p + theme(legend.position = if (keep_legend) p$theme$legend.position else "none",
            plot.title = element_blank(), plot.caption = element_blank(),
            plot.subtitle = element_blank(),
            text        = element_text(size = panel_base_size),
            axis.text   = element_text(size = panel_base_size - 1),
            axis.title  = element_text(size = panel_base_size),
            strip.text  = if (blank_strips) element_blank() else
              element_text(size = panel_base_size, face = "bold"),
            strip.clip  = "off",
            plot.tag    = element_text(size = panel_base_size + 2, face = "bold"),
            plot.tag.position = c(0, 1),
            panel.spacing.x = unit(9, "pt"),
            plot.margin = margin(8, right_pad, 2, 2))
}

#  ink key as a horizontal strip under the two plots: the names carry their own
#  printing-medium prefix, so the entries only need the colour and the name
panel_key <- function(labs, ncol = key_ncol, text_size = key_text_size) {
  labs <- unique(labs[labs %in% names(ink_colors)])
  n <- length(labs); if (!n) return(NULL)
  nr <- ceiling(n / ncol)
  d  <- data.frame(lab = disp(labs), col = unname(ink_colors[labs]),
                   col_i = ((seq_len(n) - 1) %% ncol) + 1,
                   row_i = ((seq_len(n) - 1) %/% ncol) + 1,
                   stringsAsFactors = FALSE)
  ggplot(d) +
    geom_tile(aes(x = col_i - 0.44, y = -row_i, fill = col), width = 0.07, height = 0.55,
              colour = "black", linewidth = 0.3) +
    geom_text(aes(x = col_i - 0.39, y = -row_i, label = lab), hjust = 0, size = text_size) +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0.45, ncol + 0.55), ylim = c(-nr - 0.6, -0.4), expand = FALSE) +
    theme_void() + theme(plot.margin = margin(2, 2, 2, 2))
}

panel_save <- function(a, b, name, key_labs, widths = panel_widths, ncol_key = key_ncol,
                       height_cm = NULL, keep_legend = integer(0), blank_strips = 2L,
                       small_points = 1L) {
  top <- list(a, b); keep <- !vapply(top, is.null, logical(1))
  if (!any(keep)) { message("[panel] ", name, ": no data -> not saved"); return(invisible(NULL)) }
  top <- top[keep]
  w   <- widths[seq_len(length(top))]
  kl  <- seq_along(top) %in% keep_legend
  bs  <- seq_along(top) %in% blank_strips
  sp  <- seq_along(top) %in% small_points      # only the curve plots get the small dots
  top <- Map(function(p, t, k, b, s) panel_strip(p, keep_legend = k, blank_strips = b,
                                                 small_points = s) + labs(tag = t),
             top, paste0("(", letters[seq_along(top)], ")"), kl, bs, sp)
  row1 <- patchwork::wrap_plots(top, nrow = 1, widths = w)
  kp   <- panel_key(key_labs, ncol = ncol_key)
  if (is.null(kp)) {
    pl <- row1
    h  <- if (is.null(height_cm)) panel_height_cm else height_cm
  } else {
    nkey <- ceiling(length(unique(key_labs[key_labs %in% names(ink_colors)])) / ncol_key)
    kp <- kp + labs(tag = paste0("(", letters[length(top) + 1], ")")) +
      theme(plot.tag = element_text(size = panel_base_size + 2, face = "bold"),
            plot.tag.position = c(0, 1), plot.margin = margin(4, 2, 2, 2))
    key_cm <- 0.52 * nkey + 0.25
    h  <- (if (is.null(height_cm)) panel_height_cm else height_cm) + key_cm
    pl <- patchwork::wrap_plots(list(row1, kp), ncol = 1,
                                heights = c(h - key_cm, key_cm))
  }
  d <- file.path(out_dir, panel_dir); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(d, paste0(name, ".png")), pl, width = panel_width_cm,
         height = h, units = "cm", dpi = panel_dpi, bg = "white")
  message("[panel] saved: ", file.path(panel_dir, paste0(name, ".png")))
  invisible(pl)
}

#  separators between printing media for a plot with ONE pair of axes
medium_marks <- function(labs_in_order) {
  m <- sub(".*\\((.*)\\)$", "\\1", as.character(labs_in_order))
  r <- rle(m); ends <- cumsum(r$lengths)
  list(if (length(ends) > 1)
    geom_vline(xintercept = head(ends, -1) + 0.5, linetype = "dotted",
               colour = "grey55", linewidth = 0.4))
}

#  translucent bands for the stiffness windows reported in the literature
#  (Chapter 2): 1-3 kPa in 3D and 12-15 kPa in 2D
stiffness_bands <- function(p, ymax = NULL, show_values = TRUE) {
  #  the two windows are drawn ON TOP of the bars, translucent, so that they stay
  #  visible across the whole panel; their limits are repeated on the right axis
  p <- p +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = E_band_3D[1], ymax = E_band_3D[2],
             fill = E_band_3D_fill, alpha = 0.45, colour = NA) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = E_band_2D[1], ymax = E_band_2D[2],
             fill = E_band_2D_fill, alpha = 0.45, colour = NA) +
    annotate("segment", x = -Inf, xend = Inf,
             y = c(E_band_3D, E_band_2D), yend = c(E_band_3D, E_band_2D),
             linetype = "dashed", colour = "grey25", linewidth = 0.3)
  p <- p + scale_y_continuous(
    expand = expansion(mult = c(0, if (isTRUE(show_values)) 0.5 else 0.12)),
    sec.axis = sec_axis(~ ., breaks = c(mean(E_band_3D), mean(E_band_2D)),
                        labels = c(paste0("3D\n", E_band_3D[1], "-", E_band_3D[2]),
                                   paste0("2D\n", E_band_2D[1], "-", E_band_2D[2])),
                        name = NULL))
  suppressMessages(p)
}

#  row of each ink whose x is closest to 'target'
nearest_by <- function(df, key, target, value_cols) {
  df <- df[is.finite(df[[key]]), ]
  do.call(rbind, lapply(split(df, as.character(df$label)), function(x) {
    i <- which.min(abs(x[[key]] - target)); if (!length(i)) return(NULL)
    data.frame(label = as.character(x$label[i]), medium = as.character(x$medium[i]),
               x[i, value_cols, drop = FALSE], row.names = NULL)
  }))
}

make_panels_all <- function() {
  eqx <- facet_grid(. ~ medium, scales = "free_x", space = "free_x")   # equal bar width
  
  # ---------------- (1) FLOW CURVE ----------------------------------------
  try({
    if (!is.null(rh_avg$viscosity)) {
      labs_v <- all_labels[all_labels %in% as.character(rh_avg$viscosity$label)]
      a_f <- plot_visc(rh_avg$viscosity, TRUE, panel_width_cm) + facet_grid(. ~ medium) +
        scale_x_log10(breaks = c(0.1, 10, 1000), labels = c("0.1", "10", "1000"),
                      expand = expansion(mult = 0.10))
      a_1 <- plot_visc(rh_avg$viscosity, FALSE, panel_width_cm) +
        scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000), labels = c("0.1", "1", "10", "100", "1000"))
      v <- nearest_by(rh_avg$viscosity, "rate_mean", visc_ref_rate, c("visc_mean", "visc_sd"))
      v$label <- factor(v$label, levels = labs_v); v <- v[order(v$label), ]
      b <- ggplot(v, aes(label, visc_mean, colour = label)) +
        geom_linerange(aes(ymin = pmax(visc_mean - visc_sd, visc_mean / 50), ymax = visc_mean + visc_sd),
                       linewidth = 0.6, alpha = sd_alpha, na.rm = TRUE) +
        geom_point(size = 2.4) +
        geom_text(aes(y = (visc_mean + ifelse(is.finite(visc_sd), visc_sd, 0)) * 1.25,
                      label = sprintf("%.2f", visc_mean)),
                  colour = "black", angle = 90, hjust = 0, vjust = 0.5, size = 2.3) +
        scale_colour_manual(values = ink_colors, guide = "none") +
        scale_y_log10(expand = expansion(mult = c(0.05, 0.35))) + annotation_logticks(sides = "l", size = 0.25) +
        medium_marks(levels(v$label)) +
        labs(x = NULL, y = paste0("Viscosity at ", visc_ref_rate, " s\u207B\u00B9 (Pa\u00B7s)")) +
        base_theme() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      panel_save(a_f, b, "Flow_curve", labs_v)
      panel_save(a_1, b, "Flow_curve_all_inks", labs_v)
    }
  })
  
  # ---------------- (2) THIXOTROPY ----------------------------------------
  try({
    if (!is.null(rh_avg$thixotropy)) {
      labs_t <- all_labels[all_labels %in% as.character(rh_avg$thixotropy$label)]
      a_f <- plot_thix(rh_avg$thixotropy, TRUE, panel_width_cm) + facet_grid(. ~ medium)
      a_1 <- plot_thix(rh_avg$thixotropy, FALSE, panel_width_cm)
      b <- NULL
      if (!is.null(rh_rec)) {
        rr <- rh_rec
        rr$medium <- factor(sub(".*\\((.*)\\)$", "\\1", as.character(rr$label)), levels = all_media)
        rr$label  <- factor(as.character(rr$label), levels = all_labels)
        rr <- droplevels(rr[order(rr$label), ])
        rp <- rh_rec_rep
        if (!is.null(rp)) {
          rp$medium <- factor(sub(".*\\((.*)\\)$", "\\1", as.character(rp$label)), levels = all_media)
          rp$label  <- factor(as.character(rp$label), levels = levels(rr$label))
          rp <- rp[!is.na(rp$label), ]
        }
        b <- plot_recovery(rr, TRUE, panel_width_cm, rp) + eqx
      }
      panel_save(a_f, b, "Thixotropy", labs_t)
      panel_save(a_1, b, "Thixotropy_all_inks", labs_t)
    }
  })
  
  # ---------------- (3) GELATION ------------------------------------------
  try({
    if (!is.null(rh_avg$gelation)) {
      labs_g <- all_labels[all_labels %in% as.character(rh_avg$gelation$label)]
      a_f <- plot_gel(rh_avg$gelation, TRUE, panel_width_cm) + facet_grid(. ~ medium)
      a_1 <- plot_gel(rh_avg$gelation, FALSE, panel_width_cm)
      gg <- rh_avg$gelation[is.finite(rh_avg$gelation$temp_mean) & is.finite(rh_avg$gelation$Gp_mean), ]
      b <- NULL
      if (nrow(gg)) {
        d <- do.call(rbind, lapply(gel_ref_temps, function(T) {
          x <- nearest_by(gg, "temp_mean", T, c("Gp_mean", "Gp_sd"))
          if (is.null(x)) return(NULL); x$T <- paste0(T, " \u00B0C"); x
        }))
        d$label <- factor(d$label, levels = labs_g)
        d$T     <- factor(d$T, levels = paste0(gel_ref_temps, " \u00B0C"))
        d <- d[order(d$label), ]
        b <- ggplot(d, aes(label, Gp_mean)) +
          geom_line(aes(group = label), colour = "grey55", linewidth = 0.5) +
          geom_point(aes(fill = label, shape = T), colour = "black", size = 2.4, stroke = 0.4) +
          scale_fill_manual(values = ink_colors, guide = "none") +
          scale_shape_manual(values = c(21, 24), name = NULL,
                             guide = guide_legend(override.aes = list(fill = "grey70", size = 2.4))) +
          scale_y_log10() + annotation_logticks(sides = "l", size = 0.25) +
          medium_marks(levels(d$label)) +
          labs(x = NULL, y = "G\u2032 (Pa)") + base_theme() +
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                legend.position = "top", legend.direction = "horizontal",
                legend.justification = c(0, 1), legend.box.margin = margin(0, 0, -6, 0),
                legend.background = element_blank(), legend.key = element_blank(),
                legend.margin = margin(0, 0, 0, 0), legend.title = element_blank(),
                legend.text = element_text(size = panel_base_size - 1.5))
      }
      panel_save(a_f, b, "Gelation", labs_g, keep_legend = 2)
      panel_save(a_1, b, "Gelation_all_inks", labs_g, keep_legend = 2)
    }
  })
  
  # ---------------- (4) COMPRESSION ---------------------------------------
  try({
    if (exists("comp_avg") && exists("E_sum")) {
      labs_c <- all_labels[all_labels %in% as.character(E_sum$label)]
      a_f <- make_curve_plot(comp_avg, TRUE, panel_width_cm) + facet_grid(. ~ medium)
      a_1 <- make_curve_plot(comp_avg, FALSE, panel_width_cm)
      b <- make_bar_plot(E_sum, "E_mean", "E_sd", ylab = "Young's modulus (kPa)", caption = NULL,
                         rep_df = if (exists("comp_fits")) comp_fits else NULL, rep_col = "E_kPa",
                         br_df = if (exists("bE") && !is.null(bE$br)) bE$br else NULL,
                         show_values = TRUE, value_fmt = "%.1f") + eqx
      b <- stiffness_bands(b, ymax = max(E_sum$E_mean + E_sum$E_sd, na.rm = TRUE) * 1.35)
      panel_save(a_f, b, "Compression", labs_c)
      panel_save(a_1, b, "Compression_all_inks", labs_c)
    }
  })
  
  # ---------------- (5) DEGRADATION (mass loss) ---------------------------
  #  (a) every ink over the 10 days; (b) the same curves without the diluted inks,
  #  which lose five times more mass and compress the axis of the others
  try({
    if (exists("res")) {
      labs_s <- all_labels[all_labels %in% as.character(res$label)]
      a_f <- make_swelling_plot(res, TRUE, panel_width_cm) + facet_grid(. ~ medium)
      a_1 <- make_swelling_plot(res, FALSE, panel_width_cm)
      #  (b) mass loss at two fixed times, the early plateau and the end of the assay
      ml_days <- c(3, max(res$day, na.rm = TRUE))
      ml <- do.call(rbind, lapply(ml_days, function(dd) {
        x <- res[res$day == dd, c("label", "medium", "swelling_pct", "swelling_sd")]
        if (!nrow(x)) return(NULL)
        x$day <- factor(paste0("day ", dd), levels = paste0("day ", ml_days)); x
      }))
      b <- NULL
      if (!is.null(ml) && nrow(ml)) {
        if (isTRUE(report_mass_loss)) ml$swelling_pct <- -ml$swelling_pct
        ml$label <- factor(as.character(ml$label), levels = labs_s)
        b <- ggplot(ml, aes(x = label, y = swelling_pct)) +
          geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
          geom_errorbar(aes(ymin = swelling_pct - swelling_sd, ymax = swelling_pct + swelling_sd),
                        width = 0.25, linewidth = 0.4, colour = "grey45") +
          geom_point(aes(fill = label, shape = day), colour = "black", size = 2.4, stroke = 0.4) +
          scale_fill_manual(values = ink_colors, guide = "none") +
          scale_shape_manual(values = c(22, 23), name = NULL,
                             guide = guide_legend(override.aes = list(fill = "grey70", size = 2.4))) +
          medium_marks(levels(droplevels(ml$label))) +
          labs(x = NULL, y = if (isTRUE(report_mass_loss)) "Mass loss (%)" else "Swelling degree (%)") +
          base_theme() +
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                legend.position = "top", legend.direction = "horizontal",
                legend.justification = c(0, 1), legend.box.margin = margin(0, 0, -6, 0),
                legend.background = element_blank(), legend.key = element_blank(),
                legend.margin = margin(0, 0, 0, 0), legend.title = element_blank(),
                legend.text = element_text(size = panel_base_size - 1.5))
      }
      panel_save(a_f, b, "Degradation", labs_s, keep_legend = 2)
      panel_save(a_1, b, "Degradation_all_inks", labs_s, keep_legend = 2)
    }
  })
  
  # ---------------- (6) WATER CONTENT + PRINTABILITY ----------------------
  try({
    if (exists("wc") && exists("pr_points")) {
      labs_w <- all_labels[all_labels %in% as.character(wc$label)]
      a <- make_bar_plot(wc, "water_content_pct", "water_content_sd", ylab = "Water content (%)",
                         caption = NULL, rep_df = if (exists("wc_rep")) wc_rep else NULL, rep_col = "value",
                         br_df = if (exists("bw")) bw$br else NULL,
                         ylim_zoom = wc_y_min, ybreaks = seq(wc_y_min, 100, 2), room = wc_label_room,
                         show_values = TRUE, value_fmt = "%.1f") + eqx
      pr_lev <- all_labels[all_labels %in% unique(c(as.character(pr_points$label),
                                                    if (exists("pr_status")) as.character(pr_status$label[pr_status$n_meshes == 0])))]
      b <- make_pr_plot(pr_points, pr_lev,
                        if (exists("pr_status")) pr_status[pr_status$n_meshes == 0, ] else NULL,
                        TRUE, panel_width_cm) + eqx
      panel_save(a, b, "Water_content_and_printability", all_labels, blank_strips = c(1L, 2L))
    }
  })
  
  # ---------------- (7) pH (appendix figure, on its own) ------------------
  try({
    if (exists("ph_plot_obj") && !is.null(ph_plot_obj)) {
      p <- ph_plot_obj
      panel_save(p, NULL, "pH_appendix", all_labels, widths = 1, blank_strips = integer(0))
    }
  })
  invisible()
}
if (isTRUE(make_panels)) make_panels_all()

message("Done. Files saved in: ", normalizePath(out_dir))

