# 00-functions.R -- Packages, model specification and shared helpers for the
# Themeda trait x climate analysis.

library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(brms)
library(tidybayes)
library(ggdist)
library(modelr)

# --- model specification -----------------------------------------------------

# Treatment sum contrast: HT = -0.5, LT = +0.5. Climate coefficients b_<X> are
# then treatment-averaged slopes and b_treatLTvsHT:<X> are LT - HT differences.
treat_contrast <- matrix(c(-0.5, 0.5), nrow = 2,
                         dimnames = list(c("HT", "LT"), "LTvsHT"))

# y ~ treat*X1 + treat*X2 + ... + (1 | site), sigma ~ treat + (1 | site)
model_formula <- function(preds) {
  rhs <- paste0("treat*", preds, collapse = " + ")
  bf(as.formula(paste("y ~", rhs, "+ (1 | site)")),
     sigma ~ treat + (1 | site))
}

# Priors common to all traits
prior_fixed <- c(
  set_prior("normal(0, 1)",         class = "b"),
  set_prior("student_t(3, 0, 5)",   class = "Intercept"),
  set_prior("student_t(3, 0, 2.5)", class = "sd"),
  set_prior("student_t(3, 0, 2.5)", class = "Intercept", dpar = "sigma"),
  set_prior("student_t(3, 0, 2.5)", class = "sd",        dpar = "sigma")
)

fit_file <- function(trait) paste0("out/fits/m_log_", trait, ".rds")

read_fit <- function(trait) {
  f <- fit_file(trait)
  if (!file.exists(f)) stop("Missing ", f, ": run R/02-models.R first")
  readRDS(f)
}

# --- posterior summaries -----------------------------------------------------

# Climate slopes per draw: <X>_avg (treatment-averaged), <X>_delta (LT - HT),
# and the per-treatment slopes <X>_HT and <X>_LT.
extract_posteriors <- function(fit, preds = c("MT", "MP", "AR", "PR")) {
  draws <- as_tibble(fit)
  map_dfc(preds, function(p) {
    avg   <- draws[[paste0("b_", p)]]
    delta <- draws[[paste0("b_treatLTvsHT:", p)]]
    stopifnot(!is.null(avg), !is.null(delta))
    tibble("{p}_avg"   := avg,
           "{p}_delta" := delta,
           "{p}_HT"    := avg - 0.5 * delta,
           "{p}_LT"    := avg + 0.5 * delta)
  })
}

# Variance decomposition per draw into treatment (E), climate (G),
# treatment x climate (GxE), site and residual components. Fixed-effect
# components are the variance of their linear predictor over the observed
# design; site is the site-intercept variance; residual is the mean modelled
# sigma^2 over the observed site x treatment cells. Returns V_* and
# proportions p_*.
variance_partition <- function(fit, preds = c("MT", "MP", "AR", "PR")) {
  d <- fit$data
  d$treat <- factor(d$treat, levels = c("HT", "LT"))
  contrasts(d$treat) <- treat_contrast
  X <- model.matrix(as.formula(paste("~", paste0("treat*", preds, collapse = " + "))), d)
  draws <- as_draws_matrix(fit)

  fe_var <- function(cols) {
    lp <- X[, cols, drop = FALSE] %*% t(as.matrix(draws[, paste0("b_", cols), drop = FALSE]))
    colMeans(lp^2) - colMeans(lp)^2
  }

  cells   <- count(d, site, treat)
  r_sigma <- as.matrix(draws[, sprintf("r_site__sigma[%s,Intercept]", as.character(cells$site))])
  log_sigma <- as.numeric(draws[, "b_sigma_Intercept"]) +
    outer(as.numeric(draws[, "b_sigma_treatLTvsHT"]), ifelse(cells$treat == "LT", 0.5, -0.5)) +
    r_sigma

  V <- tibble(
    V_E     = fe_var("treatLTvsHT"),
    V_G     = fe_var(preds),
    V_GxE   = fe_var(paste0("treatLTvsHT:", preds)),
    V_site  = as.numeric(draws[, "sd_site__Intercept"])^2,
    V_resid = as.numeric(exp(2 * log_sigma) %*% (cells$n / sum(cells$n)))
  )

  # the three fixed-effect blocks should add up to the full fixed-effect variance
  all_fe <- c("treatLTvsHT", preds, paste0("treatLTvsHT:", preds))
  ratio  <- fe_var(all_fe) / (V$V_E + V$V_G + V$V_GxE)
  message(sprintf("Additivity check: mean ratio = %.4f", mean(ratio)))

  V_total <- rowSums(V)
  bind_cols(tibble(.draw = seq_len(nrow(V))), V,
            rename_with(V / V_total, ~ sub("^V_", "p_", .x)))
}

# Drop the residual component from a long variance-proportion table and
# re-normalise the rest to sum to 1 within `groups`.
drop_resid_renorm <- function(df, groups,
                              keep_levels = c("site", "GxE", "G", "E")) {
  df %>%
    filter(component != "resid") %>%
    group_by(across(all_of(groups))) %>%
    mutate(proportion = proportion / sum(proportion)) %>%
    ungroup() %>%
    mutate(component = factor(as.character(component), levels = keep_levels))
}

# Intervals whose 95% credible interval excludes zero, with the upper end of
# the interval as the asterisk position.
credible_stars <- function(draws, ...) {
  draws %>%
    group_by(...) %>%
    summarise(lower = quantile(value, 0.025),
              upper = quantile(value, 0.975), .groups = "drop") %>%
    filter(lower > 0 | upper < 0)
}

# --- partial-effect and treatment-difference plots ---------------------------

# Per-site mean of a raw climate predictor
site_means <- function(data, pred) {
  data %>%
    group_by(site) %>%
    summarise(.raw_pred = mean(.data[[pred]]), .groups = "drop")
}

# Population-level predicted means (response scale) across the observed range
# of `pred`, other predictors at their means; `pred` back on the raw scale.
predicted_draws <- function(fit, pred, sc, raw_site, ndraws) {
  x_range <- (range(raw_site$.raw_pred) - sc["centre"]) / sc["scale"]
  fit$data %>%
    data_grid(x  = seq(x_range[1], x_range[2], length.out = 300),
              MT = mean(MT), MP = mean(MP), AR = mean(AR), PR = mean(PR),
              treat = c("HT", "LT")) %>%
    mutate("{pred}" := x, x = NULL) %>%
    add_epred_draws(fit, re_formula = NA, ndraws = ndraws) %>%
    mutate(.epred = exp(.epred),
           "{pred}" := .data[[pred]] * sc["scale"] + sc["centre"])
}

# Partial residuals on the response scale (`.value`): each observation's
# deviation from its fitted value, added to the population-level prediction at
# its site's `pred` value with the other predictors at their means.
partial_residuals <- function(fit, pred, sc, raw_site, ndraws) {
  fit$data %>%
    as_tibble() %>%
    left_join(raw_site, by = "site") %>%
    mutate(x = (.raw_pred - sc["centre"]) / sc["scale"],
           MT = mean(MT), MP = mean(MP), AR = mean(AR), PR = mean(PR)) %>%
    mutate("{pred}" := x, x = NULL) %>%
    {mutate(., .epred = colMeans(posterior_epred(fit, ndraws = ndraws)))} %>%
    {mutate(., .epred_marginal = colMeans(posterior_epred(fit, newdata = ., re_formula = NA,
                                                          ndraws = ndraws)))} %>%
    mutate(.value = exp(.epred_marginal) + exp(y) - exp(.epred),
           "{pred}" := .raw_pred)
}

# Partial effect of one climate predictor: modelled HT and LT means with
# credible ribbons, and partial residuals as points.
plot_predictor <- function(fit, pred, scales, data, prob = 0.95,
                           ndraws = 1000, alpha = 0.4) {
  sc       <- scales[[pred]]
  raw_site <- site_means(data, pred)
  points   <- partial_residuals(fit, pred, sc, raw_site, ndraws)

  predicted_draws(fit, pred, sc, raw_site, ndraws) %>%
    ggplot(aes(x = .data[[pred]])) +
    stat_ribbon(aes(y = .epred, fill = treat), .width = prob, alpha = alpha) +
    stat_summary(aes(y = .epred, col = treat), fun = mean, geom = "line") +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    geom_point(aes(y = .value, col = treat), data = points) +
    labs(y = "y")
}

# Treatment difference (LT - HT) along one climate predictor: modelled
# difference with credible ribbon, and per-site differences (+- SE) as points,
# from partial residuals (`points = "partial"`) or observed values ("observed").
plot_delta <- function(fit, pred, scales, data, points = "partial",
                       prob = 0.95, ndraws = 1000, alpha = 0.4) {
  sc       <- scales[[pred]]
  raw_site <- site_means(data, pred)

  ribbon <- predicted_draws(fit, pred, sc, raw_site, ndraws) %>%
    ungroup() %>%
    select(.draw, all_of(pred), treat, .epred) %>%
    pivot_wider(names_from = treat, values_from = .epred) %>%
    mutate(.delta = LT - HT)

  obs <- if (points == "partial") {
    partial_residuals(fit, pred, sc, raw_site, ndraws)
  } else {
    fit$data %>%
      as_tibble() %>%
      left_join(raw_site, by = "site") %>%
      mutate(.value = exp(y), "{pred}" := .raw_pred)
  }

  site_delta <- obs %>%
    group_by(site) %>%
    filter(any(treat == "HT") & any(treat == "LT")) %>%
    summarise(
      delta = mean(.value[treat == "LT"]) - mean(.value[treat == "HT"]),
      se    = sqrt(var(.value[treat == "LT"]) / sum(treat == "LT") +
                   var(.value[treat == "HT"]) / sum(treat == "HT")),
      "{pred}" := mean(.data[[pred]]),
      .groups = "drop")

  ggplot(ribbon, aes(x = .data[[pred]])) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_linerange(aes(ymin = delta - se, ymax = delta + se),
                   data = filter(site_delta, !is.na(se)),
                   color = "black", alpha = 0.5) +
    geom_point(aes(y = delta), data = site_delta, color = "black") +
    stat_ribbon(aes(y = .delta), .width = prob, alpha = alpha, fill = "grey70") +
    stat_summary(aes(y = .delta), fun = mean, geom = "line", color = "black") +
    labs(y = "LT - HT")
}

# --- slope and treatment-effect plots ----------------------------------------

slope_colours <- function(brewer_set1) {
  c(HT = brewer_set1[1], LT = brewer_set1[2], avg = brewer_set1[3], delta = "grey30")
}

# Slope intervals for one climate predictor: HT, LT, average and delta, with an
# asterisk above each interval whose 95% CI excludes zero.
plot_intervals <- function(post, pred, brewer_set1, transform = identity) {
  types <- c("HT", "LT", "avg", "delta")
  d <- post %>%
    select(all_of(paste0(pred, "_", types))) %>%
    pivot_longer(everything(), names_to = "type", values_to = "value") %>%
    mutate(value = transform(value),
           type  = factor(str_remove(type, paste0(pred, "_")), levels = types))

  ggplot(d, aes(x = type, y = value, color = type)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    stat_pointinterval(.width = c(0.5, 0.95), point_size = 3,
                       interval_size_range = c(0.8, 1.5)) +
    geom_text(aes(y = upper, label = "*"), data = credible_stars(d, type),
              vjust = 0.3, size = 8, color = "black", show.legend = FALSE) +
    scale_color_manual(name = NULL,
                       values = slope_colours(brewer_set1),
                       breaks = types,
                       labels = expression("High Temperature (HT)",
                                           "Low Temperature (LT)",
                                           "Average (1/2(HT + LT))",
                                           Delta*"-slope")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none",
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x  = element_blank())
}

# Give a list of plots a shared y range
apply_common_ylim <- function(plot_list) {
  y_range <- map(plot_list, function(p) {
    layers <- ggplot_build(p)$data
    range(unlist(map(layers, ~ unlist(.x[intersect(c("y", "ymin", "ymax"), names(.x))]))),
          na.rm = TRUE)
  })
  shared <- range(unlist(y_range))
  map(plot_list, ~ .x + coord_cartesian(ylim = shared, clip = "off"))
}

# Posterior of one per-trait quantity as half-eye intervals, facetted by trait
# group. `draws` is long with columns .draw, trait, quantity, value.
# `traits_on` puts traits on the x axis (groups side by side) or y axis
# (groups stacked).
plot_treat_effect <- function(draws, traits, quantity, label = NULL,
                              traits_on = c("x", "y"),
                              palette = c("Gas exchange"             = "#66C2A5",
                                          "Leaf economic"            = "#FC8D62",
                                          "Morphology and phenology" = "#8DA0CB"),
                              vline_at = 0, .width = c(0.5, 0.95),
                              scale = 0.7, slab_alpha = 0.6) {
  traits_on    <- match.arg(traits_on)
  label_vec    <- set_names(traits$short_label, traits$trait)
  trait_levels <- if (traits_on == "y") rev(traits$trait) else traits$trait

  d <- draws %>%
    filter(.data$quantity == .env$quantity, trait %in% traits$trait) %>%
    left_join(select(traits, trait, group), by = "trait") %>%
    mutate(trait = factor(trait, levels = trait_levels))

  trait_labels <- function(x) parse(text = label_vec[x])
  ref_line <- if (is.na(vline_at)) NULL
              else if (traits_on == "y") geom_vline(xintercept = vline_at, linetype = "dashed")
              else geom_hline(yintercept = vline_at, linetype = "dashed")

  layout <- if (traits_on == "y") {
    list(aes(x = value, y = trait, fill = group, colour = group),
         scale_y_discrete(labels = trait_labels),
         facet_wrap(~ group, ncol = 1, scales = "free_y"),
         labs(x = label, y = NULL),
         theme())
  } else {
    list(aes(x = trait, y = value, fill = group, colour = group),
         scale_x_discrete(labels = trait_labels),
         facet_wrap(~ group, nrow = 1, scales = "free_x"),
         labs(x = NULL, y = label),
         theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  }

  ggplot(d) +
    layout[[1]] +
    ref_line +
    stat_halfeye(point_interval = mean_qi, .width = .width, scale = scale,
                 slab_alpha = slab_alpha, normalize = "groups") +
    layout[2:5] +
    scale_fill_manual(values = palette, drop = FALSE) +
    scale_colour_manual(values = palette, drop = FALSE) +
    theme(legend.position  = "none",
          strip.background = element_blank(),
          panel.spacing    = grid::unit(0.7, "lines"))
}

# --- output ------------------------------------------------------------------

# Save a ggplot as both .pdf and .png, on a white background; any extension on
# `filename` is ignored.
save_fig <- function(filename, plot = last_plot(),
                     width = NA, height = NA, dpi = 300, bg = "white", ...) {
  stem <- sub("\\.(pdf|png)$", "", filename)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, dpi = dpi, bg = bg, ...)
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = dpi, bg = bg, ...)
  invisible(c(pdf = paste0(stem, ".pdf"), png = paste0(stem, ".png")))
}
