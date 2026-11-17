# Subevent inversion results: 2024-12-05-mww70-off-coast-of-northern

- Selected subevents: **4** | best misfit **1.4992**
- Total moment: **Mw 7.04** (catalog Mw 7.0)

## L-curve (misfit vs subevents)
```
n  best_misfit  reduction_%  delta_aic  chains  chains_within_1%  island_bests
1       2.3857            -      125.8      32                32  2.3853,2.3857,2.3875,2.3879
2       1.8898         20.8       52.6      32                25  1.8862,1.8905,1.8923,1.8935
3       1.6189         14.3       10.7      32                17  1.6261,1.6298,1.6390,1.6590
4       1.4992          7.4        0.0      32                24  1.5029,1.5982,1.6097,1.6191
5       1.4271          4.8        0.3      32                 2  1.4952,1.5866,1.6736,1.7213
misfit components of each best model: penalized = data_misfit x clvd_factor x stress_factor (evaluate_misfit_components.sh)
n    penalized  data_misfit  clvd_factor  stress_factor  worst_axis_deg  p_misfit  sh_misfit  regional_misfit
N=1     2.3857       2.3795       1.0026         1.0000             0.1     2.774      1.522            4.211
N=2     1.8898       1.8084       1.0450         1.0000            29.5     1.834      1.566            2.878
N=3     1.6189       1.5775       1.0262         1.0000            29.9     1.492      1.444            2.664
N=4     1.4992       1.4348       1.0449         1.0000            29.8     1.251      1.392            2.670
N=5     1.4271       1.3398       1.0651         1.0000            29.7     1.130      1.333            3.179
selection: AIC = ln(misfit) / sigma^2 + 2 k N, sigma = 0.05, k = 10 parameters per subevent, so a subevent costs 5% of misfit (ln(misfit) has to fall by 0.05, the misfit by 4.9%, to tie); delta_aic = AIC above the lowest; every N with delta_aic below 2 is tied, and the smallest of them is selected
NOTE: a tie - N=4 (delta_aic 0.0), N=5 (delta_aic 0.3) are equally supported (delta_aic below 2); N=4 is selected as the smallest: report both
lowest AIC: N=4
SUGGESTED number of subevents (smallest N with delta_aic below 2): 4
```

## Best model

| # | centroid (s) | X east (km) | Y north (km) | duration (s) | depth (km) | Mw |
|---|---|---|---|---|---|---|
| E1 | 5.22 | 0.00 | 0.00 | 3.74 | 3.0 | 6.58 |
| E2 | 8.45 | 6.47 | 2.45 | 8.05 | 3.0 | 6.87 |
| E3 | 12.77 | 19.36 | 4.19 | 3.48 | 3.1 | 6.50 |
| E4 | 15.37 | 47.66 | -13.49 | 5.46 | 16.0 | 6.73 |

## Figures

- `subevents.pdf/png` — map (lon/lat, 95% error bars), depth section, STFs
- `fits_P.pdf`, `fits_Pvel.pdf`, `fits_SH.pdf`, `fits_rayl.pdf` — waveform fits
- `misfit_evolution.pdf` — misfit vs sample number, all chains
- `histoplot_py.pdf` — posterior histograms (gated ensemble)
- `lcurve.pdf` — model selection: best misfit per number of subevents (star: the selected number, the smallest of those tied with the lowest AIC; dashed line: a misfit 5% lower per added subevent, the cost of a subevent; hollow points: numbers that cannot be selected, their last subevent lies on the bound of the centroid range; small open circles: island bests) and, to its right, the factors of the misfit

## Convergence of the chains
```
chains: 32; best misfit 1.4992
chains whose best is within 1% of the global best: 24
chains whose best is within 2% of the global best: 32
chains whose best is within 5% of the global best: 32
split R-hat of the gated chains (1 = converged, > 1.1 = chains disagree):
  E1: centroid time 2.07, duration 1.86, depth 1.23
  E2: centroid time 1.60, x east 1.35, y north 1.29, duration 1.40, depth 1.13
  E3: centroid time 1.92, x east 1.87, y north 1.88, duration 1.27, depth 1.04
  E4: centroid time 1.82, x east 1.37, y north 1.50, duration 1.87, depth 1.55
largest split R-hat: 2.07
```

