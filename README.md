# Trait–climate relations in *Themeda triandra*

Analysis code for Jacob et al., *Trait–climate relations in Themeda triandra, a widely distributed C4 grass and crop wild relative*.

Fifteen accessions of *T. triandra* were grown at two growth temperatures
(HT 30 °C, LT 20 °C). For each of 13 log-transformed traits, a Bayesian mixed
model (brms) relates the trait to growth temperature, four climate-of-origin
predictors and their interactions, with a site random intercept and a
heteroscedastic residual SD:

```r
bf(log(trait) ~ treat*MT + treat*MP + treat*AR + treat*PR + (1 | site),
   sigma ~ treat + (1 | site))
```

## Data

The trait data will be available in the Supplementary Material of the paper.
To run the code, place the two input files in `data/`:

- `themeda_traits.csv`: trait measurements and site climate
- `site_bioclim.csv`: bioclimatic variables per site (for the PCA comparison)

## Running

R packages are pinned with renv: run `renv::restore()` once. Then, from the
repo root:

| Script | Output (in `out/`) | ms figure |
|---|---|---|
| `R/02-models.R` | model fits, one per trait | |
| `R/03-plots.R` | `fig_<group>`, `fig_delta_<group>`, `slope_summary_propdiff`, `treat_effect_forest` | Figs 2 (lower), 3–5, S1–S4 |
| `R/04-diagnostics.R` | convergence and posterior-predictive checks | |
| `R/05-pca-variants.R` | PC-basis fits, `var_decomp_pca_bars` | Fig 6 |
| `R/06-reaction-norm.R` | `reaction_norm` | Fig 2 (upper) |
| `R/07-resid.R` | residual SD summaries | Supplement |

Run `02-models.R` first; the other scripts load its outputs and can run in any
order. `R/00-functions.R` contains the packages, model specification and shared
helpers. `R/01-data-prep.R` contains the trait and climate metadata and data
preparation.

## Licence

Code is released under the MIT licence (see `LICENSE`).
