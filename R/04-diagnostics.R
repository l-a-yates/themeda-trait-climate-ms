# 04-diagnostics.R -- Convergence diagnostics (R-hat, bulk and tail ESS,
# divergent transitions) and posterior-predictive checks for every fit.
# Outputs in out/plots/checks/:
#   diagnostics_per_param_all.csv   one row per trait x parameter
#   diagnostics_summary.csv         one row per trait, with pass/fail
#   diagnostics_rhat, diagnostics_ess   across-trait summaries
#   ppcheck_{dens,intervals,pit,skew}_<trait>   checks grouped by treatment

source("R/00-functions.R")
source("R/01-data-prep.R")
library(posterior)

theme_set(theme_classic())

# settings
rhat_thresh <- 1.01
ess_thresh  <- 400    # for both bulk and tail ESS
ppc_ndraws  <- 500

dir.create("out/plots/checks", showWarnings = FALSE, recursive = TRUE)

skewness <- function(y) mean((y - mean(y))^3) / sd(y)^3

passes <- function(max_rhat, min_ess_bulk, min_ess_tail, n_divergent) {
  max_rhat < rhat_thresh & min_ess_bulk > ess_thresh & min_ess_tail > ess_thresh &
    (is.na(n_divergent) | n_divergent == 0)
}

# --- convergence -------------------------------------------------------------
per_param <- map_dfr(traits$trait, function(trait) {
  fit  <- read_fit(trait)
  conv <- summarise_draws(as_draws_array(fit$fit), default_convergence_measures()) %>%
    filter(!variable %in% c("lp__", "lprior"))
  tibble(trait, parameter = conv$variable, rhat = conv$rhat,
         ess_bulk = conv$ess_bulk, ess_tail = conv$ess_tail,
         n_divergent_total = sum(nuts_params(fit, pars = "divergent__")$Value))
})

summary_tbl <- per_param %>%
  group_by(trait) %>%
  summarise(max_rhat     = max(rhat,     na.rm = TRUE),
            min_ess_bulk = min(ess_bulk, na.rm = TRUE),
            min_ess_tail = min(ess_tail, na.rm = TRUE),
            n_divergent  = first(n_divergent_total),
            .groups      = "drop") %>%
  mutate(pass = passes(max_rhat, min_ess_bulk, min_ess_tail, n_divergent)) %>%
  arrange(match(trait, traits$trait))

write_csv(per_param,   "out/plots/checks/diagnostics_per_param_all.csv")
write_csv(summary_tbl, "out/plots/checks/diagnostics_summary.csv")
print(summary_tbl)

# R-hat and ESS by trait, excluding per-site random effects (they stay in the CSV)
plot_dat <- per_param %>%
  filter(!str_detect(parameter, "^r_")) %>%
  mutate(trait = factor(trait, levels = traits$trait))

diag_plot <- function(d, x, threshold) {
  ggplot(d, aes(x = {{ x }}, y = trait)) +
    geom_vline(xintercept = threshold, linetype = "dashed", colour = "red") +
    geom_jitter(width = 0, height = 0.15, alpha = 0.55, size = 1.4) +
    scale_y_discrete(limits = rev) +
    labs(y = NULL, subtitle = paste("Threshold (dashed):", threshold,
                                    "-- per-site random effects excluded")) +
    theme(panel.grid.major.y = element_line(colour = "grey92"))
}

p_rhat <- diag_plot(plot_dat, rhat, rhat_thresh) +
  labs(x = expression(hat(R)), title = "Convergence: R-hat by trait")
save_fig("out/plots/checks/diagnostics_rhat", p_rhat, width = 7, height = 5)

ess_long <- plot_dat %>%
  pivot_longer(c(ess_bulk, ess_tail), names_to = "metric", values_to = "ess") %>%
  mutate(metric = recode(metric, ess_bulk = "ESS-bulk", ess_tail = "ESS-tail"))
p_ess <- diag_plot(ess_long, ess, ess_thresh) +
  scale_x_log10() +
  facet_wrap(~ metric, nrow = 1) +
  labs(x = "Effective sample size (log scale)", title = "Convergence: ESS by trait")
save_fig("out/plots/checks/diagnostics_ess", p_ess, width = 10, height = 5)

# --- posterior-predictive checks ---------------------------------------------
for (trait in traits$trait) {
  cat("PP checks for", trait, "\n")
  fit <- read_fit(trait)
  ppc <- function(type, ...) pp_check(fit, type = type, group = "treat", ndraws = ppc_ndraws, ...)
  out <- function(name) paste0("out/plots/checks/ppcheck_", name, "_", trait)
  save_fig(out("dens"),      ppc("dens_overlay_grouped"), width = 6, height = 4)
  save_fig(out("intervals"), ppc("intervals_grouped"),    width = 6, height = 4)
  save_fig(out("pit"),       ppc("pit_ecdf_grouped"),     width = 6, height = 4)
  save_fig(out("skew"),      ppc("stat_grouped", stat = skewness), width = 6, height = 4)
}

fails <- summary_tbl$trait[!summary_tbl$pass]
cat(if (length(fails) == 0) "\nAll traits passed the diagnostic thresholds.\n"
    else paste0("\nTraits with diagnostic issues: ", paste(fails, collapse = ", "), "\n"))
