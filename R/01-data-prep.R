# 01-data-prep.R -- Trait and climate metadata, the cleaned data set
# (`themeda`), and the per-trait model data used by every later script.

# Traits in figure order. `column` is the raw data column; `seed` is the MCMC
# seed for that trait's fit; labels and units are plotmath strings.
traits <- tribble(
  ~trait, ~column,                ~group,                     ~short_label,     ~units,                         ~seed,
  "A",    "A",                    "Gas exchange",             "A[sat]",         "mu*mol~CO[2]~m^{-2}~s^{-1}",   2025,
  "gs",   "gs",                   "Gas exchange",             "g[s]",           "mol~H[2]*O~m^{-2}~s^{-1}",     2026,
  "WUE",  "A_gs",                 "Gas exchange",             "WUE",            "mu*mol~mol^{-1}",              2035,
  "PN",   "PNUE",                 "Gas exchange",             "PNUE",           "mu*mol~CO[2]~g^{-1}~N~s^{-1}", 2036,
  "LM",   "LMA",                  "Leaf economic",            "LMA",            "g~m^{-2}",                     2034,
  "LDMC", "LDMC",                 "Leaf economic",            "LDMC",           "g~g^{-1}",                     2028,
  "Na",   "N_area",               "Leaf economic",            "N[area]",        "g~m^{-2}",                     2029,
  "N",    "Nitrogen",             "Leaf economic",            "N[mass]",        "'%'",                          2037,
  "Rm",   "R_mass",               "Leaf economic",            "R[mass]",        "mu*mol~CO[2]~g^{-1}~s^{-1}",   2027,
  "LW",   "Leaf_width",           "Morphology and phenology", "Leaf~width",     "mm",                           2032,
  "Ll",   "Leaf_length",          "Morphology and phenology", "Leaf~length",    "cm",                           2031,
  "Ti",   "Tillers",              "Morphology and phenology", "Tiller~number",  "",                             2030,
  "Fl",   "Flowering_time_weeks", "Morphology and phenology", "Flowering~time", "weeks",                        2033
) %>%
  mutate(group = factor(group, levels = unique(group)),
         label = if_else(units == "", short_label,
                         paste0(short_label, "~(", units, ")")))

# Climate-of-origin predictors: mean annual temperature (BIO1), mean annual
# precipitation (BIO12), annual temperature range (BIO7), precipitation
# seasonality (BIO15).
climate <- tribble(
  ~pred, ~column,              ~label,                               ~strip,
  "MT",  "MAT_bio1",           "Mean~annual~temperature~('°C')",     "MAT",
  "MP",  "MAP_bio12",          "Mean~annual~precipitation~(mm)",     "MAP",
  "AR",  "TempRange_bio7",     "Annual~temperature~range~('°C')",    "AR",
  "PR",  "PSeasonality_bio15", "Precipitation~seasonality",          "PR"
)
preds <- set_names(climate$pred)

themeda <- read_csv("data/themeda_traits.csv", show_col_types = FALSE) %>%
  # LDMC outlier: ID 55831-b is ~0.037 g/g (next smallest ~0.25), a data error
  mutate(LDMC = replace(LDMC, which(ID == "55831-b"), NA_real_)) %>%
  select(location = Location,
         treat = Growth_Temperature,
         all_of(set_names(traits$column, traits$trait)),
         all_of(set_names(climate$column, climate$pred))) %>%
  mutate(across(all_of(traits$trait), log, .names = "log_{.col}"),
         treat = factor(treat),
         site = factor(abbreviate(location)))

# Per-trait model data: log trait as `y`, climate predictors z-scored over the
# rows observed for that trait, then PR residualised on MT and MP, and AR on
# MT, MP and PR. Treatment is sum-coded at +-0.5 (HT = -0.5, LT = +0.5).
# `pred_scales` keeps each predictor's raw mean and SD for back-transforming.
trait_data  <- list()
pred_scales <- list()
for (trait in traits$trait) {
  d <- themeda %>%
    select(site, treat, all_of(preds), y = all_of(paste0("log_", trait))) %>%
    filter(!is.na(y))

  pred_scales[[trait]] <- map(d[preds], ~ c(centre = mean(.x), scale = sd(.x)))

  d <- d %>%
    mutate(across(all_of(preds), ~ as.numeric(scale(.x))),
           PR = resid(lm(PR ~ MT + MP)),
           AR = resid(lm(AR ~ MT + MP + PR)),
           treat = factor(treat, levels = c("HT", "LT")))
  contrasts(d$treat) <- treat_contrast
  trait_data[[trait]] <- d
}
