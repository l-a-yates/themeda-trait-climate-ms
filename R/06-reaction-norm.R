# 06-reaction-norm.R -- Reaction norms (Fig 2, upper): per-site mean trait
# values at HT and LT (grey lines, raw data with no model adjustment) and the
# modelled population mean with 95% credible interval (red). Means are taken
# on the log scale and back-transformed.

source("R/00-functions.R")
source("R/01-data-prep.R")

theme_set(theme_classic(base_size = 14) +
          theme(axis.text   = element_text(size = rel(0.9)),
                legend.text = element_text(size = rel(0.95))))

grey_col <- "grey60"
red_col  <- RColorBrewer::brewer.pal(3, "Set1")[1]

dir.create("out/plots", showWarnings = FALSE, recursive = TRUE)

# per-site means
site_trait_means <- map_dfr(traits$trait, function(trait) {
  themeda %>%
    filter(!is.na(.data[[paste0("log_", trait)]])) %>%
    group_by(site, treat) %>%
    summarise(value = exp(mean(.data[[paste0("log_", trait)]])), .groups = "drop") %>%
    mutate(trait = trait)
})

# modelled population means: exp(b_Intercept -/+ 0.5 * b_treatLTvsHT)
pop_means <- map_dfr(traits$trait, function(trait) {
  b <- as_tibble(read_fit(trait))
  tibble(trait = trait,
         HT = exp(b$b_Intercept - 0.5 * b$b_treatLTvsHT),
         LT = exp(b$b_Intercept + 0.5 * b$b_treatLTvsHT))
}) %>%
  pivot_longer(c(HT, LT), names_to = "treat", values_to = "value") %>%
  group_by(trait, treat) %>%
  mean_qi(value, .width = 0.95) %>%
  ungroup()

rn_panel <- function(trait) {
  g <- filter(site_trait_means, trait == .env$trait)
  r <- filter(pop_means, trait == .env$trait)
  units <- traits$units[traits$trait == trait]
  ggplot(mapping = aes(x = treat)) +
    geom_line(aes(y = value, group = site), data = g,
              colour = grey_col, alpha = 0.6, linewidth = 0.4) +
    geom_point(aes(y = value), data = g, colour = grey_col, alpha = 0.5, size = 0.8) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper, group = 1), data = r,
                fill = red_col, alpha = 0.25) +
    geom_line(aes(y = value, group = 1), data = r, colour = red_col, linewidth = 1.1) +
    geom_point(aes(y = value), data = r, colour = red_col, size = 2) +
    labs(x = NULL,
         y = if (units == "") NULL else parse(text = paste0("(", units, ")")),
         title = parse(text = traits$short_label[traits$trait == trait])) +
    theme(plot.title   = element_text(hjust = 0.5, size = rel(0.9)),
          axis.title.y = element_text(size = rel(0.8)))
}

# one row per trait group, centred on the widest row
max_cols <- max(table(traits$group))
rows <- map(levels(traits$group), function(g) {
  panels <- map(traits$trait[traits$group == g], rn_panel)
  pad <- (max_cols - length(panels)) / 2
  if (pad == 0) return(ggarrange(plotlist = panels, nrow = 1))
  ggarrange(plotlist = c(list(NULL), panels, list(NULL)), nrow = 1,
            widths = c(pad, rep(1, length(panels)), pad))
})
body <- ggarrange(rows[[1]], NULL, rows[[2]], NULL, rows[[3]],
                  ncol = 1, heights = c(1, 0.1, 1, 0.1, 1))

legend <- get_legend(
  ggplot(tibble(treat = c("HT", "LT", "HT", "LT"), y = c(0, 1, 0, 1),
                series = rep(c("Site-level", "Population-level"), each = 2)),
         aes(treat, y, colour = series, group = series)) +
    geom_line(linewidth = 1.1) +
    scale_colour_manual(values = c("Site-level" = grey_col, "Population-level" = red_col),
                        breaks = c("Site-level", "Population-level")) +
    labs(colour = NULL) +
    theme(legend.position = "top"))

p_final <- ggarrange(legend, body, ncol = 1, heights = c(0.5, 12)) %>%
  annotate_figure(bottom = text_grob("Growth-temperature treatment", size = 12))

save_fig("out/plots/reaction_norm", p_final, width = 15, height = 10,
         limitsize = FALSE)
