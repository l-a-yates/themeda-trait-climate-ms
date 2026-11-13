# 07-resid.R -- Residual SD (sigma) from the heteroscedastic fits, for the
# supplement. sigma is the residual SD on the log-trait scale, so it reads as
# proportional (CV-like) variability on the trait scale.
#   sigma_HT, sigma_LT  population-level sigma per treatment
#   sigma_effect        LT vs HT difference in log sigma, log(sigma_LT / sigma_HT)
#   site_sd             among-site SD of log sigma
#   sigma_dw            sqrt of the residual variance used in the variance
#                       decomposition (cell-size weighted over sites x treatments)
# Outputs in out/plots/: resid_sigma_by_treat, resid_sigma_effect_<scale>,
# resid_site_spread_<scale>, resid_site_sigma_caterpillar,
# resid_sigma_summary.csv, resid_site_sigma.csv

source("R/00-functions.R")
source("R/01-data-prep.R")

theme_set(theme_classic(base_size = 14) +
          theme(axis.text   = element_text(size = rel(0.9)),
                legend.text = element_text(size = rel(0.95))))
brewer_set1 <- RColorBrewer::brewer.pal(3, "Set1")

# settings
resid_scale <- "fold"
  # scale for sigma_effect and site_sd. "fold": exp(x), a ratio of SDs;
  # "pct": (exp(x) - 1) * 100, % change in SD; "log": as modelled.
ci_widths <- c(0.66, 0.95)

dir.create("out/plots", showWarnings = FALSE, recursive = TRUE)

on_scale <- function(x) switch(resid_scale, log = x, fold = exp(x), pct = (exp(x) - 1) * 100)

# --- sigma draws per trait ---------------------------------------------------
sigma_list <- list()
site_list  <- list()

for (trait in traits$trait) {
  cat("Residual sigma for", trait, "\n")
  fit   <- read_fit(trait)
  draws <- as_draws_matrix(fit)
  b0    <- as.numeric(draws[, "b_sigma_Intercept"])
  bLT   <- as.numeric(draws[, "b_sigma_treatLTvsHT"])

  sigma_list[[trait]] <- tibble(
    .draw        = seq_along(b0),
    trait        = trait,
    sigma_HT     = exp(b0 - 0.5 * bLT),
    sigma_LT     = exp(b0 + 0.5 * bLT),
    sigma_effect = bLT,
    site_sd      = as.numeric(draws[, "sd_site__sigma_Intercept"]),
    sigma_dw     = sqrt(suppressMessages(variance_partition(fit))$V_resid))

  # per-site sigma, averaged over treatments (site effect in log sigma)
  r_cols <- grep("^r_site__sigma\\[", colnames(draws), value = TRUE)
  site_list[[trait]] <- map_dfr(r_cols, function(col) {
    r <- as.numeric(draws[, col])
    tibble(.draw = seq_along(r), trait = trait,
           site = sub("^r_site__sigma\\[(.+),Intercept\\]$", "\\1", col),
           sigma_avg = exp(b0 + r), r_site = r)
  })
}

# long: .draw, trait, quantity, raw (as modelled), value (display scale)
sigma_draws <- bind_rows(sigma_list) %>%
  pivot_longer(c(sigma_HT, sigma_LT, sigma_effect, site_sd, sigma_dw),
               names_to = "quantity", values_to = "raw") %>%
  mutate(value = if_else(quantity %in% c("sigma_effect", "site_sd"), on_scale(raw), raw))

site_sigma <- bind_rows(site_list)

# --- Fig 1: residual SD by treatment -----------------------------------------
by_treat <- sigma_draws %>%
  filter(quantity %in% c("sigma_HT", "sigma_LT")) %>%
  left_join(select(traits, trait, group), by = "trait") %>%
  mutate(treat = str_remove(quantity, "sigma_"),
         trait = factor(trait, levels = rev(traits$trait)))

p_by_treat <- ggplot(by_treat, aes(x = value, y = trait, colour = treat)) +
  stat_pointinterval(.width = ci_widths, position = position_dodge(width = 0.5)) +
  scale_y_discrete(labels = function(x) parse(text = set_names(traits$short_label, traits$trait)[x])) +
  scale_colour_manual(name = NULL,
                      values = c(HT = brewer_set1[1], LT = brewer_set1[2]),
                      labels = c(HT = "High temp (HT)", LT = "Low temp (LT)")) +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  labs(x = "Residual SD (log-response scale)", y = NULL) +
  theme(legend.position  = "top",
        strip.background = element_blank(),
        panel.spacing    = grid::unit(0.7, "lines"))
save_fig("out/plots/resid_sigma_by_treat", p_by_treat, width = 8, height = 6)

# --- Fig 2: treatment effect on residual SD ----------------------------------
p_effect <- plot_treat_effect(sigma_draws, traits, "sigma_effect",
                              label = switch(resid_scale,
                                             fold = "Residual SD ratio (LT / HT)",
                                             pct  = "Change in residual SD (%)",
                                             log  = "log(SD ratio), LT vs HT"),
                              traits_on = "y", .width = ci_widths,
                              vline_at = if (resid_scale == "fold") 1 else 0)
save_fig(paste0("out/plots/resid_sigma_effect_", resid_scale), p_effect, width = 8, height = 6)

# --- Fig 3: among-site variation in residual SD ------------------------------
# on the fold scale 1 = no among-site variation; other scales have no reference
p_spread <- plot_treat_effect(sigma_draws, traits, "site_sd",
                              label = switch(resid_scale,
                                             fold = "Among-site fold in residual SD",
                                             pct  = "Among-site spread in residual SD (%)",
                                             log  = "Among-site SD of log-sigma"),
                              traits_on = "y", .width = ci_widths,
                              vline_at = if (resid_scale == "fold") 1 else NA)
save_fig(paste0("out/plots/resid_site_spread_", resid_scale), p_spread, width = 8, height = 6)

# --- Fig 4: per-site residual SD ---------------------------------------------
# sites in one order across all traits: mean over traits of the median site effect
site_order <- site_sigma %>%
  group_by(trait, site) %>%
  summarise(m = median(r_site), .groups = "drop") %>%
  group_by(site) %>%
  summarise(m = mean(m), .groups = "drop") %>%
  arrange(m) %>%
  pull(site)

p_caterpillar <- site_sigma %>%
  mutate(site  = factor(site, levels = site_order),
         trait = factor(trait, levels = traits$trait)) %>%
  ggplot(aes(x = sigma_avg, y = site)) +
  stat_pointinterval(.width = ci_widths) +
  facet_wrap(~ trait, scales = "free_x", ncol = 3,
             labeller = as_labeller(set_names(traits$short_label, traits$trait), label_parsed)) +
  labs(x = "Residual SD (log-response scale)", y = "Site (origin)")
save_fig("out/plots/resid_site_sigma_caterpillar", p_caterpillar,
         width = 14, height = 3.2 * ceiling(nrow(traits) / 3), limitsize = FALSE)

# --- tables ------------------------------------------------------------------
sigma_summary <- bind_rows(
  sigma_draws %>% group_by(trait, quantity) %>%
    mean_qi(value, .width = ci_widths) %>% mutate(.point_type = "mean"),
  sigma_draws %>% group_by(trait, quantity) %>%
    median_qi(value, .width = ci_widths) %>% mutate(.point_type = "median")
) %>%
  ungroup() %>%
  mutate(scale = resid_scale) %>%
  select(trait, quantity, scale, .point_type, .width, value, .lower, .upper) %>%
  arrange(quantity, .point_type, .width, match(trait, traits$trait))
write_csv(sigma_summary, "out/plots/resid_sigma_summary.csv")

site_summary <- site_sigma %>%
  group_by(trait, site) %>%
  mean_qi(sigma_avg, .width = ci_widths) %>%
  left_join(site_sigma %>% group_by(trait, site) %>%
              summarise(r_site_med = median(r_site), .groups = "drop"),
            by = c("trait", "site")) %>%
  ungroup() %>%
  arrange(match(trait, traits$trait), site) %>%
  select(trait, site, .width, sigma_avg, .lower, .upper, r_site_med)
write_csv(site_summary, "out/plots/resid_site_sigma.csv")

# --- console summary: sigma per treatment and the LT / HT ratio (95% CrI) ----
sigma_draws %>%
  filter(quantity %in% c("sigma_HT", "sigma_LT", "sigma_effect")) %>%
  group_by(trait, quantity) %>%
  mean_qi(raw, .width = 0.95) %>%
  ungroup() %>%
  mutate(across(c(raw, .lower, .upper), ~ if_else(quantity == "sigma_effect", exp(.x), .x)),
         est = sprintf("%.3f [%.3f, %.3f]", raw, .lower, .upper)) %>%
  select(trait, quantity, est) %>%
  pivot_wider(names_from = quantity, values_from = est) %>%
  arrange(match(trait, traits$trait)) %>%
  rename(`sigma_LT / sigma_HT` = sigma_effect) %>%
  print(n = Inf)
