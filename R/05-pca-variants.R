# 05-pca-variants.R -- Variance decomposition with the four climate predictors
# replaced by the leading k = 2, 3, 4 principal components of 20 bioclimatic
# variables, compared with the climate-predictor fits (Fig 6).
# PC fits are cached in out/fits/pca/; figures and the summary table go to
# out/plots/pca/.

source("R/00-functions.R")
source("R/01-data-prep.R")

theme_set(theme_classic())

K_SET <- c(2, 3, 4)

dir.create("out/fits/pca",  showWarnings = FALSE, recursive = TRUE)
dir.create("out/plots/pca", showWarnings = FALSE, recursive = TRUE)

# --- PCA of site climate -----------------------------------------------------
# Unrotated PCA of Moisture_Index and the 19 bioclim variables, one row per
# site. Each PC is signed so a reference variable loads positively.
env      <- read_csv("data/site_bioclim.csv", show_col_types = FALSE)
env_vars <- c("Moisture_Index", grep("_bio[0-9]+$", names(env), value = TRUE))
site_env <- env %>%
  group_by(Location) %>%
  summarise(across(all_of(env_vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

pca      <- prcomp(site_env[env_vars], scale. = TRUE)
sign_ref <- c(PC1 = "MAT_bio1", PC2 = "MAP_bio12",
              PC3 = "TempRange_bio7", PC4 = "PSeasonality_bio15")
pc_scores <- tibble(location = site_env$Location)
for (pc in names(sign_ref)) {
  pc_scores[[pc]] <- pca$x[, pc] * sign(pca$rotation[sign_ref[[pc]], pc])
}

var_explained <- pca$sdev^2 / sum(pca$sdev^2)
cat("Variance explained:", sprintf("PC%d %.1f%%", 1:4, 100 * var_explained[1:4]), "\n")

themeda_pca <- left_join(themeda, pc_scores, by = "location")
stopifnot(!anyNA(themeda_pca$PC1))

# --- fit the PC models (compiled once per k; cached fits are skipped) --------
pca_fit_file <- function(k, trait) sprintf("out/fits/pca/m_pca%d_log_%s.rds", k, trait)

for (k in K_SET) {
  pcs <- paste0("PC", seq_len(k))
  fit <- NULL
  for (i in seq_len(nrow(traits))) {
    trait <- traits$trait[i]
    if (file.exists(pca_fit_file(k, trait))) next
    cat(sprintf("Fitting pca%d %s\n", k, trait))

    d <- themeda_pca %>%
      select(site, treat, all_of(pcs), y = all_of(paste0("log_", trait))) %>%
      filter(!is.na(y)) %>%
      mutate(across(all_of(pcs), ~ as.numeric(scale(.x))),
             treat = factor(treat, levels = c("HT", "LT")))
    contrasts(d$treat) <- treat_contrast

    fit <- if (is.null(fit)) {
      brm(model_formula(pcs), data = d, prior = prior_fixed, iter = 4000,
          control = list(adapt_delta = 0.95), save_pars = save_pars(all = TRUE),
          cores = 1, seed = traits$seed[i])
    } else {
      update(fit, newdata = d, recompile = FALSE, iter = 4000,
             control = list(adapt_delta = 0.95), save_pars = save_pars(all = TRUE),
             cores = 1, seed = traits$seed[i])
    }
    saveRDS(fit, pca_fit_file(k, trait))
  }
}

# --- variance decomposition for every configuration --------------------------
configs <- c("climate", paste0("pca", K_SET))

vp_all <- map_dfr(traits$trait, function(trait) {
  bind_rows(
    variance_partition(read_fit(trait)) %>% mutate(config = "climate"),
    map_dfr(K_SET, ~ variance_partition(readRDS(pca_fit_file(.x, trait)),
                                        preds = paste0("PC", seq_len(.x))) %>%
              mutate(config = paste0("pca", .x)))
  ) %>% mutate(trait = trait)
})

comp_levels <- c("resid", "site", "GxE", "G", "E")
vp_colours  <- c(G = "#66C2A5", E = "#FC8D62", GxE = "#8DA0CB",
                 site = "#E78AC3", resid = "#DFF1C5")

vp_summary <- vp_all %>%
  group_by(config, trait) %>%
  summarise(across(all_of(paste0("p_", comp_levels)), mean), .groups = "drop") %>%
  arrange(match(config, configs), match(trait, traits$trait))
write_csv(vp_summary, "out/plots/pca/varpart_pca_summary.csv")

# traits reversed so the first trait group sits at the top of the bars
bar_long <- vp_summary %>%
  pivot_longer(starts_with("p_"), names_to = "component", values_to = "proportion") %>%
  mutate(component = factor(str_remove(component, "^p_"), levels = comp_levels),
         config    = factor(config, levels = configs),
         trait     = factor(trait, levels = rev(traits$trait)))

make_bars <- function(df) {
  ggplot(df, aes(x = proportion, y = trait, fill = component)) +
    geom_col(position = "stack", width = 0.7) +
    facet_wrap(~ config, nrow = 1) +
    scale_fill_manual(values = vp_colours) +
    guides(fill = guide_legend(reverse = TRUE)) +
    labs(x = "Proportion of variance", y = "Trait", fill = "Component") +
    theme(legend.position = "bottom")
}

# --- figures -----------------------------------------------------------------
p_bars <- make_bars(bar_long)
save_fig("out/plots/pca/var_decomp_pca_bars", p_bars, width = 11, height = 5)

p_bars_nr <- make_bars(drop_resid_renorm(bar_long, groups = c("config", "trait"))) +
  labs(x = "Proportion of variance (excl. residual)")
save_fig("out/plots/pca/var_decomp_pca_bars_no_resid", p_bars_nr, width = 11, height = 5)

# each PC configuration against the climate configuration
p_agree <- bar_long %>%
  pivot_wider(names_from = config, values_from = proportion) %>%
  pivot_longer(all_of(paste0("pca", K_SET)), names_to = "pca_config", values_to = "pca_prop") %>%
  ggplot(aes(x = climate, y = pca_prop, colour = component)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2, alpha = 0.85) +
  facet_wrap(~ pca_config, nrow = 1) +
  scale_colour_manual(values = vp_colours) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Proportion (climate basis)", y = "Proportion (PC basis)", colour = "Component") +
  theme(legend.position = "bottom")

p_compare <- ggarrange(p_bars, p_agree, ncol = 1, labels = "AUTO",
                       heights = c(1.25, 1), common.legend = TRUE, legend = "bottom")
save_fig("out/plots/pca/var_decomp_pca_compare", p_compare, width = 13, height = 11)
