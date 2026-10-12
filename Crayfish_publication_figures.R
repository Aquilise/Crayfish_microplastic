# CRAYFISH MANUSCRIPT FIGURES
# Input: final Excel workbook, GLM_dataset sheet.
# Run this entire script in RStudio. A file selection window opens if no path is set.
# First install missing packages once:
# install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "patchwork", "ragg", "svglite"))
# Outputs: editable PDF/SVG and 600-dpi PNG, plus plotted-data CSV files.
# This script does not run statistical models or add unverified significance marks.
# Verified against the workbook structure; not executed here because R is unavailable.

required <- c("readxl", "dplyr", "tidyr", "ggplot2", "patchwork", "ragg", "svglite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install these packages first: ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({library(readxl); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)})

# USER SETTINGS ---------------------------------------------------------------
excel_file <- ""       # Example: "C:/Users/Lenovo/Documents/Das et al.,RVS(5).xlsx"
output_directory <- "" # Blank = Crayfish_figures beside your Excel workbook.
font_family <- "sans"  # Or "Arial" on Windows.
digestive_label <- "GIT" # Define the pooled digestive sample in the manuscript.
make_raman_figure <- FALSE # Enable only after reconciling Raman IDs/count summaries.
export_tiff <- FALSE    # Optional 600-dpi LZW TIFF; large files.
png_dpi <- 600

if (!nzchar(excel_file)) excel_file <- file.choose()
if (!file.exists(excel_file)) stop("Excel file not found: ", excel_file)
if (!nzchar(output_directory)) output_directory <- file.path(dirname(excel_file), "Crayfish_figures")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
if (!"GLM_dataset" %in% excel_sheets(excel_file)) stop("Missing GLM_dataset sheet.")

# IMPORT AND VALIDATE ---------------------------------------------------------
dat <- read_excel(excel_file, sheet = "GLM_dataset", skip = 1, .name_repair = "unique")
names(dat) <- trimws(names(dat))
needed <- c("ID", "Site", "Individual", "Organ", "Tissue_g", "MP_count", "Presence",
            "Fibers", "Fragments", "Films", "Foam", "Red", "Green", "Black", "Blue",
            "Yellow", "White", "Grey", "Small", "Medium", "Large")
if (length(setdiff(needed, names(dat)))) stop("Missing columns: ", paste(setdiff(needed, names(dat)), collapse = ", "))
dat <- dat %>% filter(!is.na(ID)) %>% select(all_of(needed))
numeric_columns <- setdiff(needed, c("Site", "Individual", "Organ"))
dat <- dat %>% mutate(across(all_of(numeric_columns), as.numeric),
                      Site = trimws(as.character(Site)),
                      Individual = trimws(as.character(Individual)),
                      Organ = trimws(as.character(Organ)))
dat$Organ[dat$Organ %in% c("Gills", "GILLS", "gills")] <- "Gill"
if (anyNA(dat)) stop("Missing or nonnumeric required observations: inspect source records; do not replace missing values with zeros automatically.")
if (any(dat$Tissue_g <= 0)) stop("Organ mass must be positive and expressed in grams.")
if (any(dat$MP_count < 0 | dat$MP_count %% 1 != 0)) stop("MP_count must contain nonnegative integer counts.")
if (any(!dat$Site %in% c("S1", "S2", "S3")) || any(!dat$Organ %in% c("Gill", "GIT"))) stop("Unexpected site or organ label.")
if (any(dat$Presence != as.integer(dat$MP_count > 0))) stop("Presence disagrees with MP_count.")
if (anyDuplicated(dat[c("Individual", "Organ")])) stop("Duplicate animal-organ records.")
if (any(table(dat$Individual) != 2L)) stop("Each crayfish must have two organ records.")
prefix_site <- c(L = "S1", K = "S2", H = "S3")
if (any(is.na(prefix_site[substr(dat$Individual, 1, 1)])) ||
    any(dat$Site != unname(prefix_site[substr(dat$Individual, 1, 1)]))) stop("Site/sample-prefix mapping is inconsistent.")
for (cols in list(c("Fibers", "Fragments", "Films", "Foam"),
                  c("Red", "Green", "Black", "Blue", "Yellow", "White", "Grey"),
                  c("Small", "Medium", "Large"))) {
  if (any(as.matrix(dat[cols]) < 0) || any(as.matrix(dat[cols]) %% 1 != 0)) stop("Invalid particle characteristic counts.")
  if (any(rowSums(dat[cols]) != dat$MP_count)) stop("Particle characteristic totals do not agree with MP_count.")
}
dat <- dat %>% mutate(Site = factor(Site, levels = c("S1", "S2", "S3")),
                      Organ = factor(Organ, levels = c("Gill", "GIT")),
                      MP_per_g = MP_count / Tissue_g)
site_n <- dat %>% distinct(Site, Individual) %>% count(Site, name = "Animals")
site_labels <- setNames(paste0(site_n$Site, "\n", "n = ", site_n$Animals), as.character(site_n$Site))
organ_labels <- c(Gill = "Gills", GIT = digestive_label)
organ_palette <- c(Gill = "#0072B2", GIT = "#D55E00")
site_palette <- c(S1 = "#0072B2", S2 = "#009E73", S3 = "#CC79A7")

# STYLE AND EXPORT ------------------------------------------------------------
theme_paper <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(axis.title = element_text(colour = "#222222"),
          axis.text = element_text(colour = "#333333"),
          axis.line = element_line(linewidth = 0.4, colour = "#555555"),
          axis.ticks = element_line(linewidth = 0.35),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 11, margin = margin(b = 7)),
          legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = 9), legend.key.size = grid::unit(4, "mm"),
          plot.title = element_text(face = "bold", size = 12, margin = margin(b = 8)),
          plot.caption = element_text(size = 8, colour = "#555555", hjust = 0, margin = margin(t = 10)),
          plot.tag = element_text(face = "bold", size = 15),
          plot.margin = margin(10, 12, 10, 10))
}
save_figure <- function(plot, filename, width, height) {
  path <- file.path(output_directory, filename)
  ggsave(paste0(path, ".png"), plot, width = width, height = height, units = "in",
         dpi = png_dpi, device = ragg::agg_png, bg = "white")
  if (capabilities("cairo")) {
    ggsave(paste0(path, ".pdf"), plot, width = width, height = height,
           device = grDevices::cairo_pdf, family = font_family, bg = "white")
  } else {
    ggsave(paste0(path, ".pdf"), plot, width = width, height = height,
           device = grDevices::pdf, bg = "white", useDingbats = FALSE)
  }
  ggsave(paste0(path, ".svg"), plot, width = width, height = height,
         device = svglite::svglite, bg = "white")
  if (export_tiff) ggsave(paste0(path, ".tiff"), plot, width = width, height = height,
                          dpi = png_dpi, device = ragg::agg_tiff, compression = "lzw", bg = "white")
}
wilson <- function(x, n, z = 1.959964) {
  p <- x/n; denominator <- 1 + z^2/n
  centre <- (p + z^2/(2*n))/denominator
  half <- z * sqrt(p*(1-p)/n + z^2/(4*n^2))/denominator
  data.frame(lower = 100*(centre-half), upper = 100*(centre+half))
}

# FIGURE 2: OCCURRENCE ---------------------------------------------------------
individual <- dat %>% group_by(Site, Individual) %>%
  summarise(positive = any(MP_count > 0), .groups = "drop")
overall <- individual %>% group_by(Site) %>%
  summarise(n = n(), positive = sum(positive), .groups = "drop") %>%
  mutate(percent = 100*positive/n, label = paste0(positive, "/", n))
overall <- bind_cols(overall, wilson(overall$positive, overall$n))
organ_occ <- dat %>% group_by(Site, Organ) %>%
  summarise(n = n(), positive = sum(Presence), .groups = "drop") %>%
  mutate(percent = 100*positive/n, label = paste0(positive, "/", n))
organ_occ <- bind_cols(organ_occ, wilson(organ_occ$positive, organ_occ$n))
p_overall <- ggplot(overall, aes(Site, percent)) +
  geom_col(aes(fill = Site), width = .52, alpha = .9) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .12, linewidth = .5) +
  geom_text(aes(y = upper + 4, label = label), size = 3.2, family = font_family) +
  scale_fill_manual(values = site_palette, guide = "none") +
  scale_x_discrete(labels = site_labels) +
  scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 110), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Occurrence (%)", title = "Either compartment positive") + theme_paper()
dodge <- position_dodge(width = .68)
p_organs <- ggplot(organ_occ, aes(Site, percent, fill = Organ, group = Organ)) +
  geom_col(position = dodge, width = .6, alpha = .9) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = dodge, width = .14, linewidth = .5) +
  geom_text(aes(y = upper + 4, label = label), position = dodge, size = 3.1, family = font_family) +
  scale_fill_manual(values = organ_palette, labels = organ_labels) +
  scale_x_discrete(labels = site_labels) +
  scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 110), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Occurrence (%)", title = "Occurrence by compartment") + theme_paper()
figure2 <- p_overall + p_organs + plot_layout(widths = c(1, 1.25)) +
  plot_annotation(tag_levels = "A", caption = "Labels: positive/examined animals. Error bars: 95% Wilson intervals.\nIntervals are descriptive; paired-organ comparisons require the statistical model.",
                  theme = theme(plot.caption = element_text(size = 8, hjust = 0)))
save_figure(figure2, "Figure_2_occurrence", 7.2, 4.4)
write.csv(overall, file.path(output_directory, "Occurrence_individuals.csv"), row.names = FALSE)
write.csv(organ_occ, file.path(output_directory, "Occurrence_organs.csv"), row.names = FALSE)

# NEW ABUNDANCE FIGURE: PAIRED OBSERVATIONS ------------------------------------
set.seed(20261012)
animal_shift <- dat %>% distinct(Individual) %>% mutate(shift = runif(n(), -.09, .09))
paired <- dat %>% left_join(animal_shift, by = "Individual") %>%
  mutate(x = as.numeric(Organ) + shift)
segments <- paired %>% select(Site, Individual, Organ, x, MP_per_g) %>%
  pivot_wider(names_from = Organ, values_from = c(x, MP_per_g))
abundance_plot <- ggplot() +
  geom_segment(data = segments, aes(x = x_Gill, xend = x_GIT, y = MP_per_g_Gill, yend = MP_per_g_GIT),
               colour = "#A0A0A0", alpha = .45, linewidth = .35) +
  geom_boxplot(data = dat, aes(as.numeric(Organ), MP_per_g, group = Organ, fill = Organ),
               width = .48, alpha = .25, linewidth = .55, outlier.shape = NA) +
  geom_point(data = paired, aes(x, MP_per_g, fill = Organ), shape = 21,
             size = 2.4, stroke = .35, colour = "white", alpha = .9) +
  facet_wrap(~ Site, nrow = 1, labeller = as_labeller(site_labels)) +
  scale_fill_manual(values = organ_palette, labels = organ_labels, guide = "none") +
  scale_x_continuous(breaks = c(1, 2), labels = unname(organ_labels), limits = c(.6, 2.4)) +
  scale_y_continuous(expand = expansion(mult = c(.025, .08))) +
  labs(x = NULL, y = "Suspected particles/g wet organ mass",
       caption = "Each point is one organ sample; connecting lines identify samples from the same crayfish.\nBoxes show medians and interquartile ranges; whiskers extend to values within 1.5 times the IQR. Zero counts are included.") +
  theme_paper()
save_figure(abundance_plot, "Figure_abundance_paired", 7.2, 4.5)
abundance_summary <- dat %>% group_by(Site, Organ) %>%
  summarise(n = n(), particles = sum(MP_count), mean = mean(MP_per_g),
            sample_SD = sd(MP_per_g), median = median(MP_per_g), .groups = "drop")
write.csv(abundance_summary, file.path(output_directory, "Abundance_summary.csv"), row.names = FALSE)

# FIGURE 4: COMPOSITION --------------------------------------------------------
composition_data <- function(columns) {
  dat %>% select(Site, all_of(columns)) %>%
    pivot_longer(-Site, names_to = "Category", values_to = "Count") %>%
    group_by(Site, Category) %>% summarise(Count = sum(Count), .groups = "drop") %>%
    group_by(Site) %>% mutate(total = sum(Count), percent = 100*Count/total) %>% ungroup() %>%
    mutate(Category = factor(Category, levels = columns))
}
shape_data <- composition_data(c("Fragments", "Foam", "Fibers", "Films"))
color_data <- composition_data(c("Red", "Green", "Black", "Blue", "Yellow", "White", "Grey"))
size_data <- composition_data(c("Small", "Medium", "Large"))
shape_palette <- c(Fragments = "#0072B2", Foam = "#E69F00", Fibers = "#009E73", Films = "#CC79A7")
color_palette <- c(Red = "#C94C4C", Green = "#379B70", Black = "#333333", Blue = "#347DB1",
                   Yellow = "#E8BC45", White = "#FFFFFF", Grey = "#A0A0A0")
size_palette <- c(Small = "#173F5F", Medium = "#3E8EAB", Large = "#A8DADC")
particle_totals <- dat %>% group_by(Site) %>% summarise(total = sum(MP_count), .groups = "drop")
particle_labels <- setNames(paste0(particle_totals$Site, "\n", particle_totals$total, " particles"), as.character(particle_totals$Site))
composition_plot <- function(data, palette, title, labels = names(palette), legend_rows = 2) {
  ggplot(data, aes(Site, percent, fill = Category)) +
    geom_col(width = .65, colour = "#777777", linewidth = .25, position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = palette, labels = labels, drop = FALSE) +
    scale_x_discrete(labels = particle_labels) +
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 100), expand = expansion(mult = c(0, 0))) +
    guides(fill = guide_legend(nrow = legend_rows, byrow = TRUE)) +
    labs(x = NULL, y = "Composition (%)", title = title) + theme_paper() +
    theme(legend.text = element_text(size = 8), legend.key.size = grid::unit(3, "mm"))
}
p_shape <- composition_plot(shape_data, shape_palette, "Morphotype", legend_rows = 2)
p_color <- composition_plot(color_data, color_palette, "Color", legend_rows = 3)
p_size <- composition_plot(size_data, size_palette, "Size class", labels = c("50–1500 µm", "1500–3000 µm", "3000–5000 µm"), legend_rows = 3)
figure4 <- p_shape + p_color + p_size + plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A", caption = "Percentages describe stereomicroscopically identified suspected particles, not only Raman-confirmed particles.",
                  theme = theme(plot.caption = element_text(size = 8, hjust = 0)))
save_figure(figure4, "Figure_4_particle_composition", 7.2, 4.8)
write.csv(bind_rows(mutate(shape_data, Attribute = "Morphotype"),
                    mutate(color_data, Attribute = "Color"), mutate(size_data, Attribute = "Size")),
          file.path(output_directory, "Particle_composition.csv"), row.names = FALSE)

# OPTIONAL FIGURE 5A: RAMAN IDENTIFICATION RECORDS -----------------------------
if (make_raman_figure) {
  poly <- read_excel(excel_file, sheet = "Polymers_stations", skip = 1)
  names(poly) <- trimws(names(poly))
  poly <- poly %>% filter(!is.na(Sample), !is.na(Polymer)) %>%
    transmute(Sample = as.character(Sample), Polymer = trimws(as.character(Polymer)))
  if (anyDuplicated(poly$Sample)) stop("Repeated Raman particle IDs remain. Reconcile the identifiers before exporting Figure 5A.")
  polymer_names <- c("Polyethylene" = "PE", "Polyethylene terephthalate" = "PET",
                     "Polystyrene" = "PS", "Low-density polyethylene" = "LDPE",
                     "Polypropylene" = "PP", "Polycarbonate" = "PC")
  poly$Polymer <- unname(polymer_names[poly$Polymer])
  if (anyNA(poly$Polymer)) stop("Unrecognized Raman polymer name.")
  poly_summary <- poly %>% count(Polymer, name = "Records") %>%
    mutate(percent = 100*Records/sum(Records), label = sprintf("%d (%.1f%%)", Records, percent),
           Polymer = reorder(Polymer, Records))
  p_poly <- ggplot(poly_summary, aes(Records, Polymer)) +
    geom_col(fill = "#0072B2", width = .65) +
    geom_text(aes(label = label), hjust = -.12, size = 3.4, family = font_family) +
    scale_x_continuous(expand = expansion(mult = c(0, .32))) +
    labs(x = "Raman identification records", y = NULL,
         caption = "Composition of the analysed Raman subset; not an estimate for all suspected particles.") + theme_paper()
  save_figure(p_poly, "Figure_5A_Raman_composition", 4.8, 3.5)
  write.csv(poly_summary, file.path(output_directory, "Raman_composition.csv"), row.names = FALSE)
}

# CAPTION DRAFTS AND PROVENANCE ------------------------------------------------
writeLines(c(
  "Figure 2. Occurrence of visually suspected microplastic particles in spiny-cheek crayfish. (A) Individuals with particles in either sampled compartment. (B) Gill and digestive-sample occurrence. Labels give positive/examined animals; error bars show descriptive 95% Wilson confidence intervals.",
  "Figure abundance. Tissue-normalized suspected particle concentrations in paired gill and digestive samples. Points represent individual organ observations and lines connect samples from the same crayfish. Boxes show medians and interquartile ranges. Whiskers extend to observations within 1.5 times the interquartile range; all observations, including zeros, are shown.",
  "Figure 4. Composition of stereomicroscopically identified suspected particles by (A) morphotype, (B) color, and (C) size class. Denominators are the total suspected particle counts at each site. Site mapping: S1 Lipno, S2 Kořensko, S3 Hněvkovice.",
  "The optional Raman chart is disabled until Raman identifiers and count summaries are reconciled. No figure contains inferential significance marks."
), file.path(output_directory, "Figure_captions.txt"))
capture.output(sessionInfo(), file = file.path(output_directory, "R_session_info.txt"))
message("Finished. Figures and plotted-data tables are saved in: ", normalizePath(output_directory, winslash = "/"))
message("Animals: ", n_distinct(dat$Individual), "; organ samples: ", nrow(dat), "; suspected particles: ", sum(dat$MP_count))
print(figure2)
print(abundance_plot)
print(figure4)
