# 03-plots.R -- Climate-slope, partial-effect and treatment-effect figures.
#   fig_<group>             Figs 3-5: slope intervals above partial effects
#   fig_delta_<group>       Figs S2-S4: slope intervals above LT - HT differences
#   slope_summary_<scale>   Fig S1: slopes for every trait x climate predictor
#   treat_effect_forest     Fig 2 (lower): treatment effect per trait
#   combined_log_<trait>_delta  all three rows for a single trait

source("R/00-functions.R")
source("R/01-data-prep.R")

theme_set(theme_classic(base_size = 14) +
          theme(axis.text   = element_text(size = rel(0.9)),
                legend.text = element_text(size = rel(0.95))))
brewer_set1 <- RColorBrewer::brewer.pal(3, "Set1")
set.seed(2024)   # partial residuals and ribbons use random posterior draws

# settings
ribbon_prob  <- 0.95        # credible ribbon width
delta_points <- "partial"   # per-site LT - HT points: "partial" residuals or "observed" values
slope_summary_scale <- "propdiff"
  # Fig S1 scale. "propdiff": HT, LT, avg as exp(b) - 1, delta as
  # exp(b_LT) - exp(b_HT); "ratio": all four as exp(b) - 1; "log": raw slopes.

dir.create("out/plots", showWarnings = FALSE, recursive = TRUE)

hidden <- element_text(color = "transparent")

# Slope-interval, partial-effect and LT - HT panels for one trait, one of each
# per climate predictor, with a shared y range within each row.
trait_panels <- function(fit, post, trait) {
  y_lab <- traits$label[traits$trait == trait]
  int <- partial <- delta <- list()
  for (pred in preds) {
    first  <- pred == preds[1]
    y_text <- if (first) element_text() else hidden

    int[[pred]] <- plot_intervals(post, pred, brewer_set1, transform = \(x) exp(x) - 1) +
      labs(y = if (first) parse(text = "Prop.~change") else "") +
      theme(axis.text.y = y_text)

    partial[[pred]] <- plot_predictor(fit, pred, pred_scales[[trait]], themeda,
                                      prob = ribbon_prob, alpha = 0.2) +
      labs(y = if (first) parse(text = y_lab) else "") +
      theme(legend.position = "none", axis.text.y = y_text)

    delta[[pred]] <- plot_delta(fit, pred, pred_scales[[trait]], themeda,
                                points = delta_points, prob = ribbon_prob) +
      labs(y = if (first) parse(text = paste0("Delta~", y_lab)) else "") +
      theme(legend.position = "none", axis.text.y = y_text)
  }
  list(int     = apply_common_ylim(int),
       partial = apply_common_ylim(partial),
       delta   = apply_common_ylim(delta))
}

# One trait group as a figure: per trait a row of climate columns, each with
# the slope intervals above the `row` panel ("partial" or "delta").
group_figure <- function(panels, row) {
  cells <- list()
  for (trait in names(panels)) {
    bottom <- trait == last(names(panels))
    for (pred in preds) {
      x_lab <- if (bottom) parse(text = climate$label[climate$pred == pred]) else ""
      cells <- c(cells, list(ggarrange(panels[[trait]]$int[[pred]],
                                       panels[[trait]][[row]][[pred]] + labs(x = x_lab),
                                       ncol = 1, heights = c(1, 2), align = "v")))
    }
  }
  legend <- get_legend(panels[[1]]$int[[1]] +
                         guides(color = guide_legend(nrow = 1)) +
                         theme(legend.position = "top"))
  n <- length(panels)
  body <- ggarrange(plotlist = cells, nrow = n, ncol = length(preds),
                    labels = "AUTO", label.y = 1.03, legend = "none")
  ggarrange(legend, body, ncol = 1, heights = c(0.7, n * 5))
}

# --- per-trait and grouped figures -------------------------------------------
slope_draws <- list()
treat_draws <- list()

for (g in levels(traits$group)) {
  panels <- list()
  for (trait in traits$trait[traits$group == g]) {
    cat("Plotting", trait, "\n")
    fit  <- read_fit(trait)
    post <- extract_posteriors(fit)
    b    <- as_tibble(fit)[["b_treatLTvsHT"]]
    slope_draws[[trait]] <- post %>% mutate(.draw = row_number(), trait = trait)
    treat_draws[[trait]] <- tibble(.draw = seq_along(b), trait = trait,
                                   quantity = "pct_change", value = exp(b) - 1)

    p <- trait_panels(fit, post, trait)
    panels[[trait]] <- p

    columns <- map(preds, ~ ggarrange(p$int[[.x]],
                                      p$partial[[.x]] + labs(x = "") + theme(axis.text.x = hidden),
                                      p$delta[[.x]],
                                      ncol = 1, heights = c(1, 2, 2), align = "v"))
    save_fig(paste0("out/plots/combined_log_", trait, "_delta"),
             ggarrange(plotlist = columns, nrow = 1), width = 16, height = 9)
  }

  stem   <- str_replace_all(tolower(g), " ", "_")
  height <- 5 * length(panels) + 0.7
  save_fig(paste0("out/plots/fig_", stem), group_figure(panels, "partial"),
           width = 16, height = height)
  save_fig(paste0("out/plots/fig_delta_", stem), group_figure(panels, "delta"),
           width = 16, height = height)
}

# --- Fig 2 (lower): treatment effect -----------------------------------------
p_treat <- plot_treat_effect(bind_rows(treat_draws), traits, "pct_change",
                             label = "Treatment effect", traits_on = "x") +
  scale_y_continuous(labels = scales::label_percent())
save_fig("out/plots/treat_effect_forest", p_treat, width = 8, height = 5)

# --- Fig S1: slope summary ---------------------------------------------------
types <- c("HT", "LT", "avg", "delta")

slope_summary <- bind_rows(slope_draws) %>%
  pivot_longer(-c(.draw, trait), names_to = c("pred", ".value"), names_sep = "_") %>%
  mutate(delta = switch(slope_summary_scale,
                        propdiff = exp(LT) - exp(HT),
                        ratio    = exp(delta) - 1,
                        log      = delta),
         across(c(HT, LT, avg), ~ if (slope_summary_scale == "log") .x else exp(.x) - 1)) %>%
  pivot_longer(all_of(types), names_to = "type", values_to = "value") %>%
  mutate(pred  = factor(pred, levels = preds),
         type  = factor(type, levels = types),
         trait = factor(trait, levels = traits$trait))

p_slope_summary <- ggplot(slope_summary, aes(x = type, y = value, color = type)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_pointinterval(.width = c(0.5, 0.95)) +
  geom_text(aes(y = upper, label = "*"),
            data = credible_stars(slope_summary, trait, pred, type),
            vjust = 0.3, size = 6, color = "black", show.legend = FALSE) +
  scale_color_manual(values = slope_colours(brewer_set1)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  coord_cartesian(clip = "off") +
  facet_grid(trait ~ pred, scales = "free_y",
             labeller = labeller(
               trait = as_labeller(set_names(traits$short_label, traits$trait), label_parsed),
               pred  = as_labeller(set_names(climate$strip, climate$pred), label_parsed))) +
  labs(x = NULL, color = NULL,
       y = parse(text = if (slope_summary_scale == "log") "beta~(log~scale)"
                        else "Prop.~change~per~SD")) +
  theme(legend.position = "none",
        axis.text.x   = element_blank(),
        axis.ticks.x  = element_blank(),
        axis.line.x   = element_blank(),
        strip.text.y  = element_text(angle = 0),
        panel.spacing = grid::unit(0.4, "lines"))

save_fig(paste0("out/plots/slope_summary_", slope_summary_scale), p_slope_summary,
         width = 11, height = 1.5 * nrow(traits) + 1, limitsize = FALSE)
