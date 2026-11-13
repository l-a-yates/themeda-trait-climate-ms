# 02-models.R -- Fit the main model to each trait and save the fits to
# out/fits/m_log_<trait>.rds. Stan code is compiled for the first trait and
# reused for the rest.

source("R/00-functions.R")
source("R/01-data-prep.R")

dir.create("out/fits", showWarnings = FALSE, recursive = TRUE)

# log(trait) ~ treat*MT + treat*MP + treat*AR + treat*PR + (1 | site),
# sigma ~ treat + (1 | site)
formula <- model_formula(preds)

fit <- NULL
for (i in seq_len(nrow(traits))) {
  trait   <- traits$trait[i]
  control <- list(adapt_delta = if (trait == "Fl") 0.99 else 0.95)
  cat("Fitting", trait, "\n")

  fit <- if (is.null(fit)) {
    brm(formula, data = trait_data[[trait]], prior = prior_fixed,
        iter = 4000, control = control, save_pars = save_pars(all = TRUE),
        cores = 1, seed = traits$seed[i])
  } else {
    update(fit, newdata = trait_data[[trait]], recompile = FALSE,
           iter = 4000, control = control, save_pars = save_pars(all = TRUE),
           cores = 1, seed = traits$seed[i])
  }
  saveRDS(fit, fit_file(trait))
}
